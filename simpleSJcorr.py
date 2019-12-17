#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 Matthew Stone <mrstone3@wisc.edu>
# Distributed under terms of the MIT license.

"""
Simple splice junction correction
"""

import argparse
import logging
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


# TODO: refactor below into class or helper functions
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

    # turn the list into an iterator so we can look forward when correcting ends
    cigartuples = iter(read.cigartuples)

    n_intron = 0
    for op, oplen in cigartuples:
        # If we reach an intron, attempt correction
        if op == 3:
            n_intron += 1

            if n_intron >= 7:
                continue

            # TODO: unlikely for now, but add logic to catch if total correction
            # length exceeds intron length

            # Get distances from start and end to nearest reference start/end
            start_dist, end_dist = dists[curr_SJ]

            # Correct start
            if 1 <= abs(start_dist) <= window:
                # Case 1) Ref pos > tx pos
                # Pad end of last exon with deletion (D) ops, and shorten
                # intron skip accordingly
                if start_dist > 0:
                    logging.info("Truncating start of intron {0} by {1} bp".format(n_intron, abs(start_dist)))
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
                    logging.info("Extending start of intron {0} by {1} bp".format(n_intron, abs(start_dist)))
                    # track remaining length to be corrected, and amount of
                    # insertion currently added
                    correction_len = abs(start_dist)
                    ins_len = 0

                    while correction_len > 0:
                        last_op, last_oplen = cigarops.pop()
                        # if we hit an insertion, just extend it
                        if last_op == 1:
                            last_oplen += ins_len
                            correction_len = 0
                        # otherwise try to trim match/del 
                        elif last_op == 0:
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
                        elif last_op == 2:
                            # if the correction fits in the last operation,
                            # trim the last op's length
                            if last_oplen >= correction_len:
                                last_oplen -= correction_len
                                cigarops.append([last_op, last_oplen])
                                cigarops.append([1, ins_len])
                                correction_len = 0
                            # otherwise drop the op and move to the next
                            else:
                                correction_len -= last_oplen
                        else:
                            msg = "Skipping read {0}: correcting into previous intron"
                            logging.warning(msg.format(read.qname))
                            return

                    # Extend intron skip
                    oplen += abs(start_dist)
                    

            if 1 <= abs(end_dist) <= window:
                # 1) ref pos > tx pos => extend intron
                if end_dist > 0:
                    logging.info("Extending end of intron {0} by {1} bp".format(n_intron, abs(end_dist)))
                    correction_len = end_dist
                    ins_len = 0

                    # extend the intron skip first
                    cigarops.append([op, oplen + end_dist])

                    # then truncate the next operations
                    while correction_len > 0:
                        next_op, next_oplen = next(cigartuples)
                        if next_op == 1:
                            next_oplen += ins_len
                            correction_len = 0
                            
                        elif next_op == 0:
                            if next_oplen >= correction_len:
                                next_oplen -= correction_len
                                cigarops.append([1, correction_len + ins_len])
                                correction_len = 0
                            # otherwise drop the op and move to the next
                            else:
                                correction_len -= next_oplen
                                ins_len += next_oplen
                        elif next_op == 2:
                            if next_oplen >= correction_len:
                                next_oplen -= correction_len
                                cigarops.append([1, ins_len])
                                correction_len = 0
                            # otherwise drop the op and move to the next
                            else:
                                correction_len -= next_oplen
                        else:
                            msg = "Skipping read {0}: correcting into next intron"
                            logging.warning(msg.format(read.qname))
                            return

                    op, oplen = next_op, next_oplen

                # 2) ref pos < tx pos => truncate intron
                else:
                    logging.info("Truncating end of intron")
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

    update_cigarops(read, cigarops)


def update_cigarops(read, cigarops):
    """
    Updates a read's cigar string with a specified list of cigar operations

    It's easier to work with the list of cigar tuples but pysam doesn't permit
    updating the read.cigartuples attribute directly - it's parsed from the 
    cigar string each time. So, we build our list of desired operations, then
    create the corresponding CIGAR string and update the read.

    Updates the read in place.
    
    Arguments
    ---------
    read : pysam.AlignedSegment
        Read to correct
    cigarops : list of tuples of int
        List of cigartuples, i.e. [(op, oplen)]
    """

    # Construct the CIGAR string corresponding to our list of operations
    OP_MAP = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H',
              6: 'P', 7: '=', 8: 'X'}
    cigarstring = ''
    for op, oplen in cigarops:
        cigarstring += '{0}{1}'.format(oplen, OP_MAP[op])

    # Cache the old cigar string and query length before updating
    old_cigar, old_qlen = read.cigarstring, read.infer_query_length()
    old_tuples = read.cigartuples
    read.cigarstring = cigarstring

    # And reset back to old alignment if it doesn't match
    if read.infer_query_length() != old_qlen:
        msg = 'Skipping read {0}: inconsistent query length after CIGAR update'
        logging.warning(msg.format(read.qname))
        read.cigarstring = old_cigar


def merge_nearby_introns(read, window=5):
    """
    Merge introns that have insufficient aligned sequence between them.

    Updates the read in place.
    
    Arguments
    ---------
    read : pysam.AlignedSegment
        Read to correct
    window : int, optional
        SJ correction window; introns merged if they are within 2*window
    """

    cigarops = []
    prev_intron = None
    curr_ops = []

    for op, oplen in read.cigartuples:
        if op == 3:
            if prev_intron is None:
                # Mark first intron and add all preceding ops to our list
                prev_intron = (op, oplen)
                cigarops += curr_ops
                curr_ops = []
            else:
                # Count query-consuming operations
                dist = sum([x[1] for x in curr_ops if x[0] in [0, 1]])
                if dist <= 2 * window:
                    logging.info("MERGING")

                    # Merge previous and current intron
                    prev_intron = (3, prev_intron[1] + oplen)

                    # Replace intermediate operations with an insertion
                    # TODO: figure out correct operation to include here
                    curr_ops = [(1, dist)]
                else:
                    # Add last intron to our list of tracked operations
                    cigarops.append(prev_intron)
                    # Then add the operations between it and the current intron
                    cigarops += curr_ops

                    # Reset our trackers
                    curr_ops = []
                    prev_intron = (op, oplen)
        else:
            # Track all operations between introns
            curr_ops.append((op, oplen))

    cigarops.append(prev_intron)
    cigarops += curr_ops

    update_cigarops(read, cigarops)


def simpleSJcorr(read, ref_SJs, window=5, merge=False, as_unit=False):
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
    merge : bool
        Merge splice junctions within 2*window
    as_unit : bool
        Don't correct start and end independently (not yet implemented)
    """

    if merge:
        merge_nearby_introns(read, window)

    read_SJs = get_splice_junctions(read, as_pyrange=True)
    # skip reads with no introns
    if len(read_SJs) == 0:
        return

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
    parser.add_argument('-m', '--merge', action='store_true', default=False,
                        help='Merge introns that are closer than 2*window.')
    parser.add_argument('--log', help="Log file")
    args = parser.parse_args()

    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        level=logging.WARNING,
        datefmt='%Y-%m-%d %H:%M:%S')
    if args.log is not None:
        logging.getLogger().addHandler(logging.FileHandler(args.log))

    # Open sequence BAM
    bam = pysam.AlignmentFile(args.bam)
    fout = pysam.AlignmentFile(args.fout, 'wb', template=bam)

    # Load splice junctions
    df = pd.read_table(args.splice_junctions,
                       names=['Chromosome', 'Start', 'End'],
                       usecols=range(3))
    ref_SJs = pr.PyRanges(df)

    for read in bam:
        simpleSJcorr(read, ref_SJs, window=args.window, merge=args.merge)
        fout.write(read)


if __name__ == '__main__':
    main()
