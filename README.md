Localalign
===

This package implements a Smith-Waterman style local alignment algorithm. You
can align a query sequence to a reference. The scoring functions can be based
on a matrix, or simple identity. Weights can be adjusted for match/mismatch
and gaps, with gap extention penalties. Additionally, the gap penalty can be
subject to a decay to prioritize long gaps.

The input files are FASTA format sequences, or strings of sequences.
