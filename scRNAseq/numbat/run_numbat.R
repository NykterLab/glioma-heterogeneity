#!/usr/bin/env Rscript
# Script for processing WGS copy number segments and running Numbat
# Author: Serafiina Jaatinen

library(glue)
library(stringr)
library(data.table)
library(dplyr)
library(vcfR)
library(Matrix)
library(numbat)

setwd("/data/project/scRNAseq/numbat")

load("ref_annot_seurat.Robj")
load("ref_raw_rna_seurat.Robj")
ref_agg = aggregate_counts(ref_raw_rna, ref_annot)

load("GB2_raw_rna_seurat.Robj")

L1_ac = read.table("L1_allele_counts.tsv", sep="\t", header=T, stringsAsFactors=F)
L2_ac = read.table("L2_allele_counts.tsv", sep="\t", header=T, stringsAsFactors=F)
L3_ac = read.table("L3_allele_counts.tsv", sep="\t", header=T, stringsAsFactors=F)

# barcode ids
L1_ac$cell = paste0("L1_", L1_ac$cell)
L2_ac$cell = paste0("L2_", L2_ac$cell)
L3_ac$cell = paste0("L3_", L3_ac$cell)

ac = bind_rows(L1_ac, L2_ac, L3_ac)


# Numbat github issue #153 for using copy numbe segments from multiple samples
# Take union of events across samples, overlapping events: 
L1_sc = read.table("cD-GB2-L1_subclones.txt", sep="\t", header=T, stringsAsFactors=F)
L2_sc = read.table("cD-GB2-L2_subclones.txt", sep="\t", header=T, stringsAsFactors=F)
L3_sc = read.table("cD-GB2-L3_subclones.txt", sep="\t", header=T, stringsAsFactors=F)
# Ploidy ~4 in GB2, possible CNV states:
# neutral 2+2, 3+1
# del 2+1
# amp > 4
# loh X+0
# bamp 3+3, 4+4
# bdel 1+1

L1_sc = L1_sc[, c("chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]
L2_sc = L2_sc[, c("chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]
L3_sc = L3_sc[, c("chr", "startpos", "endpos", "nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]

format_seg = function(sc) {
        seg = data.frame(CHROM=integer(), seg=character(), seg_start=integer(), seg_end=integer(), cnv_state=character())
        for (i in 1:nrow(sc)) {
                # first subclone
                add_seg = data.frame(CHROM=sc[i, "chr"], seg=paste0("seg_", sc[i, "chr"], "_", i),
                                                        seg_start=sc[i, "startpos"], seg_end=sc[i, "endpos"], cnv_state=NA)

                if (sc[i, "nMaj1_A"] + sc[i, "nMin1_A"]==4) {
                        add_seg$cnv_state[1] = "neu"
                } else if (sc[i, "nMaj1_A"] + sc[i, "nMin1_A"]>4) {
                        add_seg$cnv_state[1] = "amp"
                        if (sc[i, "nMaj1_A"]==sc[i, "nMin1_A"]) {
                                add_seg$cnv_state[1] = "bamp"
                        }
                } else if (sc[i, "nMaj1_A"] + sc[i, "nMin1_A"]<4) {
                        add_seg$cnv_state[1] = "del"
                        if (sc[i, "nMaj1_A"] == sc[i, "nMin1_A"]) {
                                add_seg$cnv_state[1] = "bdel"
                        }
                }
                if (sc[i, "nMin1_A"] == 0 | sc[i, "nMaj1_A"]==0) {
                        add_seg$cnv_state[1] = "loh"
                }
                # second subclone
                if (!is.na(sc[i, "frac2_A"])) {
                        add_seg = rbind(add_seg, add_seg)
                        add_seg[2, "seg"] = paste0(add_seg[2, "seg"], "_sc")

                        if (sc[i, "nMaj2_A"] + sc[i, "nMin2_A"]==4) {
                                add_seg$cnv_state[2] = "neu"
                        } else if (sc[i, "nMaj2_A"] + sc[i, "nMin2_A"]>4) {
                                add_seg$cnv_state[2] = "amp"
                                if (sc[i, "nMaj2_A"]==sc[i, "nMin2_A"]) {
                                        add_seg$cnv_state[2] = "bamp"
                                }
                        } else if (sc[i, "nMaj2_A"] + sc[i, "nMin2_A"]<4) {
                                add_seg$cnv_state[2] = "del"
                                if (sc[i, "nMaj2_A"] == sc[i, "nMin2_A"]) {
                                        add_seg$cnv_state[2] = "bdel"
                                }
                        }
                        if (sc[i, "nMin2_A"] == 0 | sc[i, "nMaj2_A"]==0) {
                                add_seg$cnv_state[2] = "loh"
                        }
                        if (add_seg[1, "cnv_state"]==add_seg[2, "cnv_state"]) {
                                add_seg = add_seg[1,]
                        }
                }
                seg = rbind(seg, add_seg)
        }       
        seg = seg[!seg$CHROM=="X", ]
        seg$CHROM = as.integer(seg$CHROM)
        return(seg)
}

L1_seg = format_seg(L1_sc)
L2_seg = format_seg(L2_sc)
L3_seg = format_seg(L3_sc)

L1_seg$seg = paste0("L1_", L1_seg$seg)
L2_seg$seg = paste0("L2_", L2_seg$seg)
L3_seg$seg = paste0("L3_", L3_seg$seg)

seg = bind_rows(L1_seg, L2_seg, L3_seg)
seg = seg[!duplicated(seg[, c("CHROM", "seg_start", "seg_end", "cnv_state")]),]

write.table(seg, "GB2_consensus_segs.txt", row.names=F, quote=F, sep="\t")

# Run Numbat
out = run_numbat(GB2_raw_rna,
                ref_agg,
                ac,
                segs_consensus_fix=seg,
                genome="hg38",
                #t = 1e-5, # lower t for complex copy numbers and subclones, higher t for controlling FP 
                ncores=4,
                skip_nj=T,
                max_iter = 10,
                check_convergence=T,
                tau=0.6, # increase tau to get fewer clones in general
                min_genes=10,
                plot=T,
                out_dir="multisample_results",
                diploid_chroms=c(5, 11), # based on WGS
                call_clonal_loh=T)

