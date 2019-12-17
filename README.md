# simpleSJcorr

Simple splice junction correction.

## Summary

This program corrects the boundaries of all splice junctions in a long-read
RNA-sequencing read with respect to a set of reference splice junctions. If the
start or end coordinate of an observed splice junction is within a permissible
window of a reference junction (by default 5 bp), the program edits the read
alignment so the boundary matches the reference.

Specifically, the read sequence remains unaltered. Correction involves editing
the CIGAR string of the read so the inferred lengths of the query sequence and
its alignment to the reference remain the same, while the inferred positions of
the splice junction boundaries are shifted to match the corresponding reference
positions. In short, if the intron needs to be truncated, the corresponding
skip (N) operation is shortened and the exon alignment is extended with
deletions. If the intron needs to be extended, the skip operation is lengthened
and aligned bases are replaced with insertions. 

## Algorithm details

1. Splice junction/intron boundaries are parsed from the read's CIGAR string.
   Per the definition of CIGAR operations in the [SAM
   spec](https://samtools.github.io/hts-specs/SAMv1.pdf), matches and deletions
   count towards progress along the reference sequence, and skips (N) are
   considered introns.
2. Optionally (using `-m` flag), merge adjacent introns if the intervening
   sequence is shorter than twice the correction window. This prevents cases
   where extending one intron may overlap with another.  
   Specifically, any intervening query sequence is converted to an insertion
   and placed before the start of the merged intron, and any corresponding
   reference sequence is converted to a skip to maintain the start and end
   positions of the merged intron.
3. Start and end positions of each splice junction are matched to the nearest
   start or end, respectively, of a reference splice junction. If the distance
   to the reference boundary is within the user-specified threshold (default 5 
   bp), the alignment is corrected to match the reference position.
4. Correction of an intron takes two forms - truncation and extension.  
   Truncation is the simpler operation; the length of the intron is shortened
   by shortening the length of the corresponding skip operation, and the length of
   the reference alignment is maintained by adding a deletion of corresponding
   length to the start/end of the intron.  
   Extension consumes any matches or deletions beyond the intron boundary up to
   the correction distance, replacing matches with insertions to conserve query
   sequence length, and then extends the intron (skip) by the correction 
   distance.

## Usage

```
simpleSJcorr.py -h
usage: simpleSJcorr.py [-h] [-w WINDOW] [-m] [--log LOG]
                       bam splice_junctions fout

Simple splice junction correction

positional arguments:
  bam                   CCS/transcripts to correct (BAM format).
  splice_junctions      Reference splice junctions to correct against.
                        Expected bed-like (chrom, start, end).
  fout                  Corrected sequences (BAM format).

optional arguments:
  -h, --help            show this help message and exit
  -w WINDOW, --window WINDOW
                        Permissible window for correction [5 bp].
  -m, --merge           Merge introns that are closer than 2*window.
  --log LOG             Log file
```
