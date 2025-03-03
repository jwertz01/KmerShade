# KmerShade
Color genes/ sequences over locus.

This pipeline takes in two fastas, one with sample sequences over a locus, and the other with sub-sequences (e.g. genes) within the locus.
It divides the sample sequences into k-mers and colors each k-mer according to which gene it maps most closely to.
