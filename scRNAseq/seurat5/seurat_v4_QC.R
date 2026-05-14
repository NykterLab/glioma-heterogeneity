# Script for Seurat QC and processing of scRNA-seq data (CellRanger output)
# Author: Serafiina Jaatinen

# Use Seurat v4.4.0 in QC because SCP library does not work with Seurat v5
library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(devtools)
library(SCP)

sc_path = "/data/project/scRNAseq/"
setwd(sc_path)

samplenames <- c("GB2_L1", "GB2_L2", "GB2_L3")

# Sample-specific QC parameters
db_rates <- c(0.0066, 0.066, 0.048)
UMI_thresholds <- c(3000, 1000, 1000)
gene_thresholds <- c(5000, 1500, 1500)
mito_thresholds <- c(20, 12, 8)

samples <- list()
norm_samples <- list()
cellqc_samples <- list()
filtered_samples <- list()

for (i in seq_along(samplenames)) {

  sample_name <- samplenames[i]
  sample_data <- Read10X(data.dir = paste0("cellranger/", samplenames[i], "/outs/filtered_feature_bc_matrix/"))
  samples[[i]] <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 3, min.features = 200)

  # Plot QC metrics
  samples[[i]]$percent.mt <- PercentageFeatureSet(samples[[i]], pattern = "^MT-")
  samples[[i]]$percent.ribo <- PercentageFeatureSet(samples[[i]], pattern = "^RP[SL]")

  v <- VlnPlot(samples[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
  ggsave(paste0("seurat/", sample_name, "_nFeature_nCount_percMt_QC.pdf"), plot = v, width = 10, height = 5)

  v_log <- VlnPlot(samples[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, log = TRUE)
  ggsave(paste0("seurat/", sample_name, "_nFeature_nCount_percMt_QC_log.pdf"), plot = v_log, width = 10, height = 5)

  # Normalize and preprocess
  norm_samples[[i]] <- NormalizeData(samples[[i]])
  norm_samples[[i]] <- FindVariableFeatures(norm_samples[[i]])
  norm_samples[[i]] <- ScaleData(norm_samples[[i]], features = rownames(norm_samples[[i]]))
  norm_samples[[i]] <- RunPCA(norm_samples[[i]], features = VariableFeatures(norm_samples[[i]]))
  norm_samples[[i]] <- RunUMAP(norm_samples[[i]], reduction = "pca", dims = 1:50, n.neighbors = 30)

  # CellQC diagnostics
  cellqc_samples[[i]] <- RunCellQC(norm_samples[[i]], return_filtered = FALSE, 
    qc_metrics = c("doublets", "umi", "gene", "mito", "ribo_mito_ratio", "outlier"),
    db_method = "scDblFinder",
    db_rate = db_rates[i],
    outlier_threshold = c("log10_nCount:lower:2.5", "log10_nCount:higher:5", "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5"),
    outlier_n = 1,
    UMI_threshold = UMI_thresholds[i],
    gene_threshold = gene_thresholds[i],
    mito_threshold = mito_thresholds[i])

  dp <- CellDimPlot(srt = cellqc_samples[[i]], group.by = c("CellQC", "db_qc", "umi_qc", "gene_qc", "mito_qc", "outlier_qc", "ribo_mito_ratio_qc"), reduction = "UMAP")
  ggsave(paste0("seurat/", sample_name, "_norm_cellQC_umap_all.pdf"), dp, width = 15, height = 15)

  csp <- CellStatPlot(srt = cellqc_samples[[i]], stat.by = c("db_qc", "umi_qc", "gene_qc", "mito_qc", "outlier_qc", "ribo_mito_ratio_qc"), plot_type = "upset", stat_level = "Fail")
  ggsave(paste0("seurat/", sample_name, "_norm_cellQC_stat.pdf"), csp, width = 7, height = 5)

  dp_feat <- FeatureDimPlot(srt = cellqc_samples[[i]], features = c("percent.mt", "percent.ribo", "nFeature_RNA", "nCount_RNA"), reduction = "UMAP")
  ggsave(paste0("seurat/", sample_name, "_norm_percmito_percribo_nfeat_ncount_umap.pdf"), dp_feat, width = 10, height = 10)

  # Filtering
  filtered_samples[[i]] <- RunCellQC(norm_samples[[i]],return_filtered = TRUE,
    qc_metrics = c("doublets", "gene", "mito", "ribo_mito_ratio"),
    db_method = "scDblFinder",
    db_rate = db_rates[i],
    gene_threshold = gene_thresholds[i],
    mito_threshold = mito_thresholds[i])

  # Save filtered object
  save(filtered_samples[[i]], file = paste0("seurat/", sample_name, "_norm_cellqc_filtered.Robj"))
}
