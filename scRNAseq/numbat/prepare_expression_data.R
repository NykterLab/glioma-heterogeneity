#!/usr/bin/env Rscript
# Script for processing Seurat v5 expression data for Numbat copy number analysis
# Author: Serafiina Jaatinen

library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(devtools)
library(SCP)
library(harmony)

sc_path = "/data/project/scRNAseq"
setwd(sc_path)

ref_data_path = "/reference_data/"

##################################################
# Preprocessed public data
load(paste0(ref_data_path, "richards_ft_list.Robj"))
load(paste0(ref_data_path, "neftel_ft_list.Robj"))
load(paste0(ref_data_path, "pombo1_ft_list.Robj"))
load(paste0(ref_data_path, "pombo2_ft_list.Robj"))
load(paste0(ref_data_path, "slyper_ft_list.Robj"))

richards_normal_list = lapply(richards_ft_list, function(x) subset(x, cells=WhichCells(x, expression=!preliminary_cell_types %in% c("cancer cells", "cancer cells (no CNV)", "proliferating cancer cells (no CNV)"))))

ref_merged = merge(slyper_ft_list[[1]], y=c(richards_ft_list, neftel_ft_list, pombo1_ft_list, pombo2_ft_list), 
					add.cell.ids=c(names(slyper_ft_list), names(richards_ft_list), names(neftel_ft_list), names(pombo1_ft_list), names(pombo2_ft_list)),
					merge.data=FALSE, project="Harmony_ref")

# Only normal cells for the reference
cells_to_keep = setdiff(Cells(ref_merged), WhichCells(ref_merged, expression = preliminary_cell_types %in% c("cancer cells", "cancer cells (no CNV)", "proliferating cancer cells (no CNV)")))
ref_normal = subset(ref_merged, cells=cells_to_keep)

ref_normal <- NormalizeData(ref_normal)
ref_normal <- FindVariableFeatures(ref_normal)
ref_normal <- ScaleData(ref_normal)
ref_normal <- RunPCA(ref_normal)

ref_normal_harmony = IntegrateLayers(object=ref_normal, method=HarmonyIntegration, orig.reduction="pca", new.reduction="harmony")

ref_normal_harmony <- FindNeighbors(ref_normal_harmony, reduction = "harmony", dims = 1:50)
ref_normal_harmony <- FindClusters(ref_normal_harmony, resolution = 0.5, cluster.name = "harmony_clusters")
ref_normal_harmony <- RunUMAP(ref_normal_harmony, dims = 1:50, reduction = "harmony", reduction.name = "umap.harmony")

sample_cols = iwanthue(48)
d = DimPlot(ref_normal_harmony, cols=sample_cols, reduction = "umap.harmony", group.by=c("orig.ident", "author", "preliminary_cell_types", "harmony_clusters"))
ggsave("seurat/ref_harmony_umap.pdf", d, width=17, height=10)

ref_normal_harmony = JoinLayers(ref_normal_harmony)

# Save raw counts and cluster annotations
ref_raw_rna = GetAssayData(ref_normal_harmony, assay="RNA", layer="counts")
ref_annot = data.frame(cell=rownames(ref_normal_harmony@meta.data), group=ref_normal_harmony@meta.data$seurat_clusters)
save(ref_raw_rna, file="numbat/ref_raw_rna_seurat.Robj")
save(ref_annot, file="numbat/ref_annot_seurat.Robj")

################################################
# In-house expression data

load("seurat/L1_norm_cellqc_filtered.Robj")
load("seurat/L2_norm_cellqc_filtered.Robj")
load("seurat/L3_norm_cellqc_filtered.Robj")

L1_norm = UpdateSeuratObject(L1_norm_cellqc_f)
L2_norm = UpdateSeuratObject(L2_norm_cellqc_f)
L3_norm = UpdateSeuratObject(L3_norm_cellqc_f)

L1_raw_rna = GetAssayData(L1_norm, assay="RNA", layer="counts")
L2_raw_rna = GetAssayData(L2_norm, assay="RNA", layer="counts")
L3_raw_rna = GetAssayData(L3_norm, assay="RNA", layer="counts")

save(L1_raw_rna, file="numbat/L1_raw_rna_seurat.Robj")
save(L2_raw_rna, file="numbat/L2_raw_rna_seurat.Robj")
save(L3_raw_rna, file="numbat/L3_raw_rna_seurat.Robj")

# Integrated, seurat v5 object
load("seurat/GB2_harmony.Robj")
raw_rna = GetAssayData(GB2_harmony, assay="RNA", layer="counts")
save(raw_rna, file="numbat/GB2_raw_rna_seurat.Robj")
