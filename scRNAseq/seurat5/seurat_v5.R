# Script for Seurat v5
# Author: Serafiina Jaatinen

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(devtools)
library(SCP)
library(hues)
library(scCustomize)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(purrr)
library(scales)
library(glue)
library(data.table)
library(ggtree)
library(tidygraph)
library(msigdbr)
library(openxlsx)
library(pheatmap)

sc_path = "/data/project/scRNAseq/"
setwd(sc_path)

load("seurat/GB2_L1_norm_cellqc_filtered.Robj")
load("seurat/GB2_L2_norm_cellqc_filtered.Robj")
load("seurat/GB2_L3_norm_cellqc_filtered.Robj")

script_path = "/data/scripts/scRNAseq/seurat5/"

source(paste0(script_path, "load_genesets.R"))


L1_f = UpdateSeuratObject(L1_norm_cellqc_f)
L2_f = UpdateSeuratObject(L2_norm_cellqc_f)
L3_f = UpdateSeuratObject(L3_norm_cellqc_f)

L1_f = Convert_Assay(seurat_object=L1_f, convert_to="V5")
L2_f = Convert_Assay(seurat_object=L2_f, convert_to="V5")
L3_f = Convert_Assay(seurat_object=L3_f, convert_to="V5")

L1_f$cell_barcode_orig = rownames(L1_f@meta.data)
L2_f$cell_barcode_orig = rownames(L2_f@meta.data)
L3_f$cell_barcode_orig = rownames(L3_f@meta.data)

L1_f$author = "Tampere"
L2_f$author = "Tampere"
L3_f$author = "Tampere"

L1_f$preliminary_cell_types = "cancer cells GB2"
L2_f$preliminary_cell_types = "cancer cells GB2"
L3_f$preliminary_cell_types = "cancer cells GB2"

sample_cols = iwanthue(48)

L1_norm = NormalizeData(L1_f)
L1_norm = FindVariableFeatures(L1_norm)
L1_norm = ScaleData(L1_norm, features=rownames(L1_norm))

L2_norm = NormalizeData(L2_f)
L2_norm = FindVariableFeatures(L2_norm)
L2_norm = ScaleData(L2_norm, features=rownames(L2_norm))

L3_norm = NormalizeData(L3_f)
L3_norm = FindVariableFeatures(L3_norm)
L3_norm = ScaleData(L3_norm, features=rownames(L3_norm))

L1_norm <- RunPCA(L1_norm)
L1_norm <- FindNeighbors(L1_norm, dims = 1:50, reduction = "pca")
L1_norm <- FindClusters(L1_norm, resolution = 0.8, cluster.name = "clusters")
L1_norm <- RunUMAP(L1_norm, dims = 1:50, reduction = "pca", reduction.name = "umap")

L2_norm <- RunPCA(L2_norm)
L2_norm <- FindNeighbors(L2_norm, dims = 1:50, reduction = "pca")
L2_norm <- FindClusters(L2_norm, resolution = 0.8, cluster.name = "clusters")
L2_norm <- RunUMAP(L2_norm, dims = 1:50, reduction = "pca", reduction.name = "umap")

L3_norm <- RunPCA(L3_norm)
L3_norm <- FindNeighbors(L3_norm, dims = 1:50, reduction = "pca")
L3_norm <- FindClusters(L3_norm, resolution = 0.8, cluster.name = "clusters")
L3_norm <- RunUMAP(L3_norm, dims = 1:50, reduction = "pca", reduction.name = "umap")

d = CellDimPlot(L1_norm, reduction="umap", group.by=c("clusters"))
ggsave("seurat/L1_norm_clusters_umap.pdf", d, width=15, height=10)

d = CellDimPlot(L2_norm, reduction="umap", group.by=c("clusters"))
ggsave("seurat/L2_norm_clusters_umap.pdf", d, width=15, height=10)

d = CellDimPlot(L3_norm, reduction="umap", group.by=c("clusters"))
ggsave("seurat/L3_norm_clusters_umap.pdf", d, width=15, height=10)

L1_norm = CellCycleScoring(L1_norm, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T)
L2_norm = CellCycleScoring(L2_norm, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T)
L3_norm = CellCycleScoring(L3_norm, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T)

save(L1_norm, file="seurat/L1_norm.Robj")
save(L2_norm, file="seurat/L2_norm.Robj")
save(L3_norm, file="seurat/L3_norm.Robj")


##########################
# Neftel subtypes

neftel_list = list(AC, OPC, NPC1, NPC2, MES1, MES2)
names(neftel_list) = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")

L1_norm = AddModuleScore(L1_norm, features=neftel_list, name=names(neftel_list), search=T)
L2_norm = AddModuleScore(L2_norm, features=neftel_list, name=names(neftel_list), search=T)
L3_norm = AddModuleScore(L3_norm, features=neftel_list, name=names(neftel_list), search=T)

colnames(L1_norm@meta.data)[16:21] = names(neftel_list)
colnames(L2_norm@meta.data)[16:21] = names(neftel_list)
colnames(L3_norm@meta.data)[16:21] = names(neftel_list)

L1_norm@meta.data$neftel_state = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")[max.col(L1_norm@meta.data[, c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")], ties.method="first")]
L2_norm@meta.data$neftel_state = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")[max.col(L2_norm@meta.data[, c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")], ties.method="first")]
L3_norm@meta.data$neftel_state = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")[max.col(L3_norm@meta.data[, c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")], ties.method="first")]


# Plot in Neftel space
source(paste0(script_path, "plot_neftel_space.R"))

L1_norm@meta.data$D = apply(L1_norm@meta.data, 1, calculate_D)
L2_norm@meta.data$D = apply(L2_norm@meta.data,1, calculate_D)
L3_norm@meta.data$D = apply(L3_norm@meta.data,1, calculate_D)

L1_norm@meta.data$x = apply(L1_norm@meta.data, 1, calculate_x)
L2_norm@meta.data$x = apply(L2_norm@meta.data, 1, calculate_x)
L3_norm@meta.data$x = apply(L3_norm@meta.data, 1, calculate_x)

p = plot_neftel_space(L1_norm@meta.data, "neftel_state") + scale_color_manual(values=neftel_pal)
ggsave("seurat/L1_norm_neftel_space.pdf", p,  width=6.5, height=4)
p = plot_neftel_space(L2_norm@meta.data, "neftel_state") + scale_color_manual(values=neftel_pal)
ggsave("seurat/L2_norm_neftel_space.pdf", p,  width=6.5, height=4)
p = plot_neftel_space(L3_norm@meta.data, "neftel_state") + scale_color_manual(values=neftel_pal)
ggsave("seurat/L3_norm_neftel_space.pdf", p,  width=6.5, height=4)


dp = CellDimPlot(srt=L1_norm, group.by="neftel_state", reduction="umap", palcolor=c("#A6CEE3","#1F78B4",   "#B2DF8A",   "#33A02C",   "#FDBF6F"))
ggsave("seurat/L1_norm_neftel_state_annot.pdf", dp, width=5, height=5)
dp = CellDimPlot(srt=L2_norm, group.by="neftel_state", reduction="umap", palcolor=c("#33A02C", "#A6CEE3"  , "#1F78B4" ,  "#B2DF8A", "#FDBF6F" ,  "#FF7F00"))
ggsave("seurat/L2_norm_neftel_state_annot.pdf", dp, width=5, height=5)
dp = CellDimPlot(srt=L3_norm, group.by="neftel_state", reduction="umap", palcolor=c("#FDBF6F", "#A6CEE3", "#B2DF8A", "#1F78B4","#33A02C" ,  "#FF7F00"))
ggsave("seurat/L3_norm_neftel_state_annot.pdf", dp, width=5, height=5)


# Plot genesets
features = c(sojka_list, mp_list[2:13], richards_list[c(1:6, 64:65, 68:79)], list(developmental_gsc_top=developmental_gsc_top, injury_gsc_top=injury_gsc_top, ncc=ncc))

L1_norm = AddModuleScore(L1_norm, features=features, name=names(features), search=T)
L2_norm = AddModuleScore(L2_norm, features=features, name=names(features), search=T)
L3_norm = AddModuleScore(L3_norm, features=features, name=names(features), search=T)

plot_dot_plot = function(seurat_obj, minval = NULL, maxval = NULL, module_names) {
  
  module_names <- unique(module_names)  # avoid duplicate factor levels
  
  # Get metadata + module scores
  df <- FetchData(seurat_obj, vars = c("clusters", module_names))
  # Long format
  
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = all_of(module_names),
      names_to = "module",
      values_to = "score"
    )
  
  # Compute mean score & % cells with expression > 0
  df_summary <- df_long %>%
    group_by(clusters, module) %>%
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      perc_cells = (sum(score > 0, na.rm = TRUE) / n()) * 100,
      .groups = "drop"
    )
  
  # Compute row/module ordering 
  module_order <- df_summary %>%
    group_by(module) %>%
    summarise(avg_score_all_clusters = mean(mean_score, na.rm = TRUE)) %>%
    arrange(avg_score_all_clusters) %>%       # low → high; reverse if desired
    pull(module)

  # Apply row order (reverse so highest is at top of plot)
  df_summary <- df_summary %>%
    mutate(module = factor(module, levels = rev(module_order)))
  
  # Auto-determine color scale limits if not supplied
  if (is.null(minval)) minval <- min(df_summary$mean_score, na.rm = TRUE)
  if (is.null(maxval)) maxval <- max(df_summary$mean_score, na.rm = TRUE)
  
  message(paste0("Color scale limits: ", round(minval, 3), " to ", round(maxval, 3)))
  
  # Plot
  ggplot(df_summary, aes(x = clusters, y = module)) +
    geom_point(aes(size = perc_cells, color = mean_score)) +
    scale_color_distiller(palette = "Spectral", limits = c(minval, maxval), oob = scales::squish) +
    scale_size_continuous(range = c(1, 10), name = "% expressing") +
    theme_linedraw(base_size = 14) +
    labs(x = "Cluster", y = "Module (Gene Set)", color = "Mean Module Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

p = plot_dot_plot(L1_norm, module_names=c(names(fts), "S.Score", "G2M.Score"), minval=-0.1, maxval=0.4)
ggsave("seurat/L1_norm_fts_dotplot.pdf", p, width=8, height=20)
p = plot_dot_plot(L2_norm, module_names=c(names(fts), "S.Score", "G2M.Score"), minval=-0.1, maxval=0.4)
ggsave("seurat/L2_norm_fts_dotplot.pdf", p, width=12, height=20)
p = plot_dot_plot(L3_norm, module_names=c(names(fts), "S.Score", "G2M.Score"), minval=-0.1, maxval=0.4)
ggsave("seurat/L3_norm_fts_dotplot.pdf", p, width=12, height=20)



# Differential expression and gene set enrichment
L1_norm = RunDEtest(L1_norm, group_by="clusters", assay="RNA", only.pos=F)
L1_norm = RunGSEA(L1_norm, group_by="clusters", TERM2GENE=hallmark_sets, DE_threshold="TRUE")
g = GSEAPlot(L1_norm, group_by="clusters", plot_type="comparison", db="custom", direction="both")
ggsave("seurat/L1_norm_gsea_dotplot_clusters.pdf", g, width=8, height=4)

L2_norm = RunDEtest(L2_norm, group_by="clusters", assay="RNA", only.pos=F)
L2_norm = RunGSEA(L2_norm, group_by="clusters", TERM2GENE=hallmark_sets, DE_threshold="TRUE")
g = GSEAPlot(L2_norm, group_by="clusters", plot_type="comparison", db="custom", direction="both")
ggsave("seurat/L2_norm_gsea_dotplot_clusters.pdf", g, width=10, height=6)

L3_norm = RunDEtest(L3_norm, group_by="clusters", assay="RNA", only.pos=F)
L3_norm = RunGSEA(L3_norm, group_by="clusters", TERM2GENE=hallmark_sets, DE_threshold="TRUE")
g = GSEAPlot(L3_norm, group_by="clusters", plot_type="comparison", db="custom", direction="both")
ggsave("seurat/L3_norm_gsea_dotplot_clusters.pdf", g, width=10, height=6)

