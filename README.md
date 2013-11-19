swalign
===

This package implements a Smith-Waterman style local alignment algorithm. You
can align a query sequence to a reference. The scoring functions can be based
on a matrix, or simple identity. Weights can be adjusted for match/mismatch
and gaps, with gap extention penalties. Additionally, the gap penalty can be
subject to a decay to prioritize long gaps.

The input files are FASTA format sequences, or strings of sequences.

Here is some skeleton code to get you started:

	import swalign
	# choose your own values here… 2 and -1 are common.
	match = 2
	mismatch = -1
	scoring = swalign.NucleotideScoringMatrix(match, mismatch)

	sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...
	alignment = sw.align('ACACACTA','AGCACACA')
	alignment.dump()

For other uses, see the script in bin/swalign.