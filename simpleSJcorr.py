#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 Matthew Stone <mrstone3@wisc.edu>
# Distributed under terms of the MIT license.

"""
Simple splice junction correction
"""

import argparse
import pysam
import pandas as pd
import pyranges as pr


def get_splice_junctions(read, as_pyrange=False):
    """
    Get positions of internal splice junctions from CIGAR string.

    Skips (N) are considered to represent introns. True for GMAP alignments,
    may not be otherwise.

    Arguments
    ---------
    read : pysam.AlignedSegment
        Sequence to parse splice junctions from
    as_pyrange : bool, optional
        If true, return a PyRanges table of splice junctions instead of a
        list of coordinates

    Returns
    -------
    SJs : <tuple of tuple of int> OR <pr.PyRanges>
        1) List of (start, end) pairs of splice junctions OR
        2) PyRanges object with (chromosome, start, end) columns
    """

    splice_sites = []
    curr_pos = read.pos

    for op, oplen in read.cigartuples:
        # Track current position wrt reference by
        # adding matching bases or deleted bases
        if op == 0 or op == 2:
            curr_pos += oplen

        # Assume skipped bases indicate intron
        elif op == 3:
            splice_begin = curr_pos + 1
            curr_pos += oplen
            splice_end = curr_pos
            splice_sites.append((splice_begin, splice_end))

    # Return splice junctions in specified format
    if as_pyrange:
        df = pd.DataFrame(splice_sites, columns=['Start', 'End'])
        df['Chromosome'] = read.reference_name
        return pr.PyRanges(df[['Chromosome', 'Start', 'End']])
    else:
        return tuple(splice_sites)


def correct_splice_junctions(read, dists, window=5):
    """
    Corrects splice junction coordinates by updating CIGAR string
    
    Arguments
    ---------
    read : pysam.AlignedSegment
        Read to correct
    dists : list of tuples of int
        Distances to nearest start/end of a reference splice junction
    window : int, optional
        Maximum distance to correct

    Returns
    -------
    read : pysam.AlignedSegment
        Read with corrected CIGAR string
    """

    # TODO: think about refactoring to not copy/paste the SJ parsing logic
    cigarops = []
    curr_SJ = 0
    curr_pos = read.pos

    # turn the list into an iterator so we can look forward when correcting ends
    cigartuples = iter(read.cigartuples)

    for op, oplen in cigartuples:
        # Track current position wrt reference by
        # adding matching bases or deleted bases
        if op == 0 or op == 2:
            curr_pos += oplen

        # Assume skipped bases indicate intron
        elif op == 3:
            splice_begin = curr_pos + 1  # BAM is 0-start, TxClean bed is 1
            curr_pos += oplen
            splice_end = curr_pos

            # TODO: unlikely for now, but add logic to catch if total correction
            # length exceeds intron length
            start_dist, end_dist = dists[curr_SJ]
            if 1 <= abs(start_dist) <= window:
                # Case 1) Ref pos > tx pos
                # Pad end of last exon with deletion (D) ops, and shorten
                # intron skip accordingly
                if start_dist > 0:
                    # If last op was a del, just add to its length
                    if cigarops[-1][0] == 2:
                        cigarops[-1][1] += start_dist

                    # Otherwise add a deletion
                    else:
                        cigarops.append([2, start_dist])

                    # then truncate the intron skip
                    oplen -= start_dist
                    
                # Case 2) Ref pos < tx pos
                # Replace last match/delete operations with insertions, and
                # extend intron skip accordingly
                else:
                    # track remaining length to be corrected, and amount of
                    # insertion currently added
                    correction_len = abs(start_dist)
                    ins_len = 0
                    while correction_len > 0:
                        last_op, last_oplen = cigarops.pop()
                        # if we hit an insertion, just extend it
                        if last_op == 1:
                            last_oplen += correction_len + ins_len
                            correction_len = 0
                        # otherwise try to trim match/del 
                        elif last_op in [0, 2]:
                            # if the correction fits in the last operation,
                            # trim the last op's length and replace with insertion
                            if last_oplen >= correction_len:
                                last_oplen -= correction_len
                                cigarops.append([last_op, last_oplen])
                                cigarops.append([1, correction_len + ins_len])
                                correction_len = 0
                            # otherwise drop the op and move to the next
                            else:
                                correction_len -= last_oplen
                                ins_len += last_oplen
                        else:
                            raise Exception("Correcting into previous intron")

                    # Extend intron skip
                    oplen += abs(start_dist)

            if 1 <= abs(end_dist) <= window:
                # 1) ref pos > tx pos => extend intron
                if end_dist > 0:
                    correction_len = end_dist
                    ins_len = 0

                    # extend the intron skip first
                    cigarops.append([op, oplen + end_dist])

                    # then truncate the next operations
                    while correction_len > 0:
                        next_op, next_oplen = next(cigartuples)
                        if next_op == 1:
                            next_oplen += correction_len + ins_len
                            correction_len = 0
                            
                        elif next_op in [0, 2]:
                            if next_oplen >= correction_len:
                                next_oplen -= correction_len
                                cigarops.append([1, correction_len + ins_len])
                                correction_len = 0
                            # otherwise drop the op and move to the next
                            else:
                                correction_len -= next_oplen
                                ins_len += next_oplen

                    op, oplen = next_op, next_oplen

                # 2) ref pos < tx pos => truncate intron
                else:
                    # First, truncate intron skip and add to cigar op list
                    oplen -= abs(end_dist)
                    cigarops.append([op, oplen])

                    # Next, get next op. If it's an deletion, extend it,
                    # otherwise add a new deletion before it
                    next_op, next_oplen = next(cigartuples)
                    if next_op == 2:
                        next_oplen += abs(end_dist)
                    else:
                        cigarops.append([2, abs(end_dist)])

                    # and set up the next op to be added back
                    op, oplen = next_op, next_oplen

            curr_SJ += 1

        # note: adding as list, not tuple, so we can update dels in place above
        cigarops.append([op, oplen])

    OP_MAP = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H'}
    cigarstring = ''
    for op, oplen in cigarops:
        cigarstring += '{0}{1}'.format(oplen, OP_MAP[op])

    read.cigarstring = cigarstring


def simpleSJcorr(read, ref_SJs, window=5, as_unit=False):
    """
    Correct splice junctions in a read relative to a reference set of junctions

    Arguments
    ---------
    read : pysam.AlignedSegment
        Aligned sequence to correct
    SJs : pr.PyRanges
        Table of reference splice junctions (Chromosome, Start, End)
    window : int
        Maximum distance for correction
    as_unit : bool
        Don't correct start and end independently (not yet implemented)
    """

    read_SJs = get_splice_junctions(read, as_pyrange=True)
    read_starts = pr.PyRanges(chromosomes=read_SJs.Chromosome, starts=read_SJs.Start, ends=read_SJs.Start + 1)
    read_ends = pr.PyRanges(chromosomes=read_SJs.Chromosome, starts=read_SJs.End, ends=read_SJs.End + 1)

    ref_starts = pr.PyRanges(chromosomes=ref_SJs.Chromosome, starts=ref_SJs.Start, ends=ref_SJs.Start + 1)
    ref_ends = pr.PyRanges(chromosomes=ref_SJs.Chromosome, starts=ref_SJs.End, ends=ref_SJs.End + 1)

    nearest_start = read_starts.nearest(ref_starts, suffix='_ref')
    start_dist = nearest_start.Start_ref - nearest_start.Start
    nearest_end = read_ends.nearest(ref_ends, suffix='_ref')
    end_dist = nearest_end.Start_ref - nearest_end.Start

    dists = list(zip(start_dist, end_dist))

    correct_splice_junctions(read, dists, window)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bam',
                        help='CCS/transcripts to correct (BAM format).')
    parser.add_argument('splice_junctions',
                        help='Reference splice junctions to correct against. '
                        'Expected bed-like (chrom, start, end).')
    parser.add_argument('fout',
                        help='Corrected sequences (BAM format).')
    parser.add_argument('-w', '--window', type=int, default=5,
                        help='Permissible window for correction [5 bp].')
    args = parser.parse_args()

    # Open sequence BAM
    bam = pysam.AlignmentFile(args.bam)
    fout = pysam.AlignmentFile(args.fout, 'wb', template=bam)

    # Load splice junctions
    df = pd.read_table(args.splice_junctions,
                       names=['Chromosome', 'Start', 'End'],
                       usecols=range(3))
    ref_SJs = pr.PyRanges(df)

    for read in bam:
        simpleSJcorr(read, ref_SJs)
        fout.write(read)


if __name__ == '__main__':
    main()
