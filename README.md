# ABC-GWAS
This repository contains scripts to run analysis for the modules in the ABC-GWAS web resource. The folders containing the relevant files should be modified by the user accordingly. There are three scripts in this folder:
1. eQTL_analysis.R: This script is used to run eQTL analysis for the TCGA data. Imputed genotypes, gene expresssion and copy number segmentation data are the necessary input datasets.
2. motif_analysis.R: Given a set of SNPs, this script computes permutation test p-values for the differential binding affinity of a collection of transcription factors.
3. compute_corr_gene_TF.R: This script requires the results from eQTL and motif analysis from the first two scripts. Here, the aim is to compute the correlation coefficients between the eQTL genes and candidate transcription factors for each SNP.
