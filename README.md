LDSC (LD SCore) v 0.00001
======================

Copyright (c) 2014 Brendan Bulik-Sullivan & Hilary Finucane

(Warning: under very active development)


What is LD Score?
--------------

LD Score is a measure of the amount of linkage disequilibrium (LD) around a SNP. 
LD Score allows one to properly take LD into account when attempting to make 
inferences about the genetic architecture of complex disease using GWAS summary 
statistics.

What can I do with LDSC?
---------------------

1. Estimate LD Score (and other, more exotic sums-of-moments of the genotype distribution).
2. Quantify the inflation in GWAS test statistics from confounding bias.
3. Estimate heritability from GWAS summary statistics.
4. Estimate partitioned heritability from GWAS summary statistcs.
5. Estimate genetic covariance and correlation from GWAS summary statistics.


Installation
------------

0. Download python 2.7
1. Download dependencies: numpy, scipy, pandas, bitarrary, progressbar, argparse
	(e.g., via pip install)
2. Download dependencies for unittests: unittest, nose, nose_parameterized
3. git clone https://github.com/bulik/ldsc.git
4. cd to repository root directory, run tests by typing nosetests