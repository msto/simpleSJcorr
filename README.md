# simpleSJcorr

Simple splice junction correction.

If the start or end coordinate of a read's splice junction is within a
permissible window of a reference junction (by default 5 bp), the script
"corrects" the read alignment so the coordinate matches the reference.
Correction involves updating the CIGAR string of the alignment; the sequence
remains unaltered. If the intron needs to be extended, the SKIP (N) operation
is extended and aligned bases (M or D) are replaced with insertions (I). If the
intron needs to be truncated, the SKIP operation is shortened and the exon
alignment is extended with deletions.

