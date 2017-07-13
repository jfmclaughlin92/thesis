# thesis

A repository of scripts used in JF McLaughlin's MS thesis. Versions are as used in the thesis.

find_chrom.py Using BLASTn output (as an .xml), find loci mapped to a specific locus (user-specified) and create an output file excluding them. Initially written to remove Z-linked loci from UCE datasets. Requirements: BioPython, phyluce

biallelicSNPs.py Find non-biallelic SNPs and create a new vcf file without them. Necessary for preparing files for analyses limited to biallelic data.

randomSNP.py When thinning to 1 SNP per locus, select a random SNP (rather than just the first SNP listed).

bootstrap_dadi.py resamples individuals (with replacement) to create bootstrapped datasets and runs each a user-specified number of times with DaDi.

ngapi_dadi.py resamples individuals without replacement within a given subsampled dataset, and similar to bootstrap_dadi.py, runs each a specified number of times in DaDi. Used for investigating impact of low sample size on estimates of demographic parameters.
