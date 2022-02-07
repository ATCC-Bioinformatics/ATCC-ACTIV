# ATCC-ACTIV
Codebase for scripts written for the collaboration with ACTIV-TRACE working group

This repository is intended to house ATCC-written code designed to work with ACTIV TRACE on the SARS-CoV-2 activity.

**SPDI_parsimony.py**

Designed to alter different variant calling reporting formats into a consistent and context-removed format to better find agreement between different groups.
It is based on NIH's SPDI normalization results and expects an input format that is consistent with that produced by the ACTIV-TRACE group.

It was principally written to remove identical leading/trailing bases between Ref and Alt alleles, handle Indel and SNP calls in the same variant record, as well as force multiple adjacent SNPs to individual rows per SNP call. 
It can readily be used as a template for general VCF format normalization.
