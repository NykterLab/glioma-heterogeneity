#!/usr/bin/env Rscript
# Script for running Battenberg
# Give arguments in order: normal bam path, normal sample name, tumor bam path, tumor sample name,
# breakpoints file path, working directory path, minimun ploidy, maximum ploidy, number of threads, is male

library(foreach)
library(doParallel)
library(parallel)
library(GenomicRanges)
library(VariantAnnotation)
library(Battenberg)
library(ASCAT)
library(copynumber)

args <- commandArgs(trailingOnly = TRUE)

normal_bam <- normalizePath(args[1])
normal_sample <- args[2]
tumor_bam <- normalizePath(args[3])
tumor_sample <- args[4]
svfilepath <- normalizePath(args[5])
work_dir <- normalizePath(args[6])
min_ploidy <- args[7]
max_ploidy <- args[8]
nthreads <- args[9]
ismale <- args[10]

system( paste0("mkdir -p ", work_dir), intern=T)
setwd( work_dir)

imputeinfofile = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/imputation/impute_info_mod.txt"
g1000prefix = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"

g1000prefix_ac = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
gccorrectprefix = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/GC_correction_hg38/1000G_GC_chr"
problemloci = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/probloci/probloci.txt"

beaglejar = "/data/tools/beagle5/beagle.12Jul19.0df.jar"
beagleref = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/beagle5/chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz"
beagleplink = "/data/tools/battenberg_hg38/Battenberg_hg38_reference/beagle5/plink.chrCHROMNAME.GRCh38.map"

battenberg(analysis="paired", # or cell_line
		tumourname=tumor_sample,
		normalname=normal_sample,
		tumour_data_file=tumor_bam,
		normal_data_file=normal_bam,
		ismale=ismale,
		nthreads=nthreads,
		imputeinfofile=imputeinfofile,
		g1000prefix=g1000prefix,
		g1000allelesprefix=g1000prefix_ac,	   
		gccorrectprefix=gccorrectprefix, 
		problemloci=problemloci,
		min_ploidy=min_ploidy,
		max_ploidy=max_ploidy,
		allelecounter_exe="alleleCounter",
		skip_preprocessing=FALSE,
		skip_allele_counting=FALSE,
		skip_phasing=FALSE,
		usebeagle=TRUE,
		beaglejar=beaglejar,
		beagleref=beagleref,
		beagleplink=beagleplink,		   
		GENOMEBUILD="hg38",		   
		prior_breakpoints_file=svfilepath)

