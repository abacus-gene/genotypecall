# genotypecall

This program calculates the genotype-calling error rate for both
homozygotes and heterozygotes as a function of read depth (n) and
base-calling error rate (epsilon).  Genotype is inferred using ML
based on reads for one sample (one individual), and the true
single-read error rate epsilon is assumed to be known and given.
Error rates are assumed to be the same among reads, independent of the
true nucleotides.
   
The model is one of simple binomial sampling.  See equations (4.12),
(4.13) and (4.14) in the mathematical notes of Li (2011), for a more
complex model.  Note that a homozygote has genotype 11: with 1 to be
the true base, and 0 the error.  Data is k, the number of 1's on p.20.

To compile the program, using something like the following

      cl -Ox genotypecall.c tools.c 
      gcc -o genotypecall -O3 genotypecall.c tools.c -lm 
      icc -o genotypecall -O3 genotypecall.c tools.c -lm 
      genotypecall


Li, H., 2011 A statistical framework for SNP calling, mutation
discovery, association mapping and population genetical parameter
estimation from sequencing data. Bioinformatics 27: 2987-2993.

