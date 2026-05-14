# Script for public data processing
# Author: Seraiina Jaatinen

library(Seurat) # v5

ref_data_path = "/reference_data/"

# Richards 2021
richards_raw = read.csv(paste0(ref_data_path, "Richards_NatureCancer_GBM_scRNAseq_counts.csv"))
rownames(richards_raw) = richards_raw[,1]
richards_raw = richards_raw[,-1]
richards = CreateSeuratObject(richards_raw, project="Richards", min.cells=3, min.features=200)

# Neftel 2019 GSE131928
# neftel2 not included
neftel_raw1 = ReadMtx(mtx=paste0(ref_data_path, "IDHwtGBM.processed.10X.counts.mtx"), cells=paste0(ref_data_path, "cells1.proc.tsv"), features=paste0(ref_data_path, "genes1.tsv"))
neftel1 = CreateSeuratObject(neftel_raw1, project="Neftel_1", min.cells=3, min.features=200)

# Pombo Antunes 2021 GSE163120
# human recurrent GBM
antunes_raw1 = read.csv(paste0(ref_data_path, "GSM4972210_Human.GBM.R1_2_3_4_4nc.filtered.gene.bc.matrix.csv"))
# human primary GBM
antunes_raw2 = read.csv(paste0(ref_data_path, "GSM4972211_Human.GBM.ND1_2_3_4_5_6_7.filtered.gene.bc.matrix.csv"))
rownames(antunes_raw1) = antunes_raw1[,1]
antunes_raw1 = antunes_raw1[,-1]
rownames(antunes_raw2) = antunes_raw2[,1]
antunes_raw2 = antunes_raw2[,-1]
antunes1 = CreateSeuratObject(antunes_raw1, project="Pombo_Antunes_R", min.cells=3, min.features=200)
antunes2 = CreateSeuratObject(antunes_raw2, project="Pombo_Antunes_ND", min.cells=3, min.features=200)
# add samplenames
antunes1_annot = read.csv(paste0(ref_data_path, "GSM4972210_annot.Human.GBM.R1_2_3_4_4nc.csv"))
antunes2_annot = read.csv(paste0(ref_data_path, "GSM4972211_annot.Human.GBM.ND1_2_3_4_5_6_7.csv"))
antunes1_annot$cell = sub("-", ".", antunes1_annot$cell)
all(antunes1_annot$cell == colnames(antunes1))
colnames(antunes1_annot) = c("x", "y", "ident", "cluster", "cell", "orig.ident")
antunes1 = AddMetaData(antunes1, antunes1_annot$orig.ident, col.name="orig.ident")
#
antunes2_annot$cell = sub("-", ".", antunes2_annot$cell)
all(antunes2_annot$cell == colnames(antunes2))
colnames(antunes2_annot) = c("x", "y", "ident", "cluster", "cell", "orig.ident")
antunes2 = AddMetaData(antunes2, antunes2_annot$orig.ident, col.name="orig.ident")

# Slyper 2019 GSE140819
# seurat object produced in Seurat v4 where hdf5r library worked
# slyper_raw = Read10X_h5(paste0(ref_data_path, "GSM4186981_MGH125_fresh_channel1_raw_feature_bc_matrix.h5"))
# slyper = CreateSeuratObject(slyper_raw, project="Slyper", min.cells=3, min.features=200)
# save(slyper, file=paste0(paste0(ref_data_path, "Slyper_raw_seurat_obj.Robj")))
load(paste0(ref_data_path, "Slyper_raw_seurat_obj.Robj"))
slyper = UpdateSeuratObject(slyper)


# Filter
richards_ft = richards
richards_ft$percent.mt = PercentageFeatureSet(richards_ft, pattern="^MT-")
richards_ft = subset(richards_ft, subset=percent.mt < 5)

neftel_ft = neftel1
neftel_ft$percent.mt = PercentageFeatureSet(neftel_ft, pattern="^MT-")
neftel_ft = subset(neftel_ft, subset=percent.mt < 5)

antunes1_ft = antunes1
antunes1_ft$percent.mt = PercentageFeatureSet(antunes1_ft, pattern="^MT-")
antunes1_ft = subset(antunes1_ft, subset=percent.mt < 5)
antunes1_ft = SetIdent(antunes1_ft, value=antunes1_ft@meta.data$orig.ident)

antunes2_ft = antunes2
antunes2_ft$percent.mt = PercentageFeatureSet(antunes2_ft, pattern="^MT-")
antunes2_ft = subset(antunes2_ft, subset=percent.mt < 5)
antunes2_ft = SetIdent(antunes2_ft, value=antunes2_ft@meta.data$orig.ident)

slyper_ft = slyper
slyper_ft$percent.mt = PercentageFeatureSet(slyper_ft, pattern="^MT-")
slyper_ft = subset(slyper_ft, subset=percent.mt < 5)

richards_ft$cell_barcode_orig = rownames(richards_ft@meta.data)
antunes1_ft$cell_barcode_orig = rownames(antunes1_ft@meta.data)
antunes2_ft$cell_barcode_orig = rownames(antunes2_ft@meta.data)
neftel_ft$cell_barcode_orig = rownames(neftel_ft@meta.data)
slyper_ft$cell_barcode_orig = rownames(slyper_ft@meta.data)

richards_ft$author = "Richards"
antunes1_ft$author = "Antunes"
antunes2_ft$author = "Antunes"
neftel_ft$author = "Neftel"
slyper_ft$author = "Slyper"

# add preliminary annotations analyzed in Tiihonen et al. 2026
load(paste0(ref_data_path, "seurat_metadata_GBM.RData"))
richards_ft$preliminary_cell_types = seurat_metadata_GBM[richards_ft$cell_barcode_annot, "preliminary_cell_types"]
antunes1_ft$preliminary_cell_types = seurat_metadata_GBM[antunes1_ft$cell_barcode_annot, "preliminary_cell_types"]
antunes2_ft$preliminary_cell_types = seurat_metadata_GBM[antunes2_ft$cell_barcode_annot, "preliminary_cell_types"]
neftel_ft$preliminary_cell_types = seurat_metadata_GBM[neftel_ft$cell_barcode_annot, "preliminary_cell_types"]
slyper_ft$preliminary_cell_types = seurat_metadata_GBM[slyper_ft$cell_barcode_annot, "preliminary_cell_types"]

# separate layers per samples
richards_ft_list = SplitObject(richards_ft, split.by="ident")
neftel_ft_list = SplitObject(neftel_ft, split.by="ident")
antunes1_ft_list = SplitObject(antunes1_ft, split.by="ident")
antunes2_ft_list = SplitObject(antunes2_ft, split.by="ident")
slyper_ft_list = SplitObject(slyper_ft, split.by="ident")

# save reference lists for numbat
save(richards_ft_list, file=paste0(ref_data_path, "richards_ft_list.Robj"))
save(neftel_ft_list, file=paste0(ref_data_path, "neftel_ft_list.Robj"))
save(antunes1_ft_list, file=paste0(ref_data_path, "antunes1_ft_list.Robj"))
save(antunes2_ft_list, file=paste0(ref_data_path, "antunes2_ft_list.Robj"))
save(slyper_ft_list, file=paste0(ref_data_path, "slyper_ft_list.Robj"))

