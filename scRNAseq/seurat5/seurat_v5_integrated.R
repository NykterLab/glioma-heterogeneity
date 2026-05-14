# Script for Seurat v5 analysis on integrated samples
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
library(rlang)

sc_path = "/data/project/scRNAseq/"
setwd(sc_path)

load("seurat/L1_norm.Robj")
load("seurat/L2_norm.Robj")
load("seurat/L3_norm.Robj")

script_path = "/data/scripts/scRNAseq/seurat5/"

source(paste0(script_path, "load_genesets.R"))

###########################
# Integration
GB2 = merge(L1_norm, y=c(L2_norm, L3_norm), 
                add.cell.ids=c("GB2_L1", "GB2_L2", "GB2_L3"),
                merge.data=FALSE, project="Harmony_GB2")

GB2 <- NormalizeData(GB2)
GB2 <- FindVariableFeatures(GB2)
GB2 <- ScaleData(GB2)

GB2 <- RunPCA(GB2)
GB2 <- FindNeighbors(GB2, dims = 1:50, reduction = "pca")
GB2 <- FindClusters(GB2, resolution = 0.5, cluster.name = "unintegrated_clusters")
GB2 <- RunUMAP(GB2, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")

d = DimPlot(GB2, cols=sample_cols, reduction = "umap.unintegrated", group.by=c("unintegrated_clusters", "orig.ident"))
ggsave("seurat/GB2_unintegrated_umap.pdf", d, width=10, height=5)

GB2_harmony = IntegrateLayers(object=GB2, method=HarmonyIntegration, orig.reduction="pca", new.reduction="harmony")

GB2_harmony <- FindNeighbors(GB2_harmony, reduction = "harmony", dims = 1:50)
GB2_harmony <- FindClusters(GB2_harmony, resolution = 0.8, cluster.name = "harmony_clusters")
GB2_harmony <- RunUMAP(GB2_harmony, dims = 1:50, reduction = "harmony", reduction.name = "umap.harmony")

d = DimPlot(GB2_harmony, cols=sample_cols, reduction = "umap.harmony", group.by=c("orig.ident", "harmony_clusters"))
ggsave("seurat/GB2_harmony_umap.pdf", d, width=10, height=5)

GB2_harmony = JoinLayers(GB2_harmony)

#########################
# QC
GB2_harmony = CellCycleScoring(GB2_harmony, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=T)

f = FeatureDimPlot(GB2_harmony, features=c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), reduction="umap.harmony", bg_cutoff=-Inf)
ggsave("seurat/GB2_harmony_qc_umap.pdf", f, width=15, height=10)

f = FeatureDimPlot(GB2_harmony, features=c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score"), reduction="pca", bg_cutoff=-Inf)
ggsave("seurat/GB2_harmony_qc_pca.pdf", f, width=15, height=10)

v = VlnPlot(GB2_harmony, features=c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
ggsave("seurat/GB2_harmony_clusters_qc_violin.pdf", v, width=15, height=10)

Idents(GB2_harmony) = "harmony_clusters"

#########################
# Neftel states
neftel_list = list(AC, OPC, NPC1, NPC2, MES1, MES2)
names(neftel_list) = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")

neftel_pal = c("MES1"="#A6CEE3", "MES2"="#1F78B4", "NPC1"="#B2DF8A", "NPC2"="#33A02C", "AC"="#FDBF6F", "OPC"="#FF7F00")

GB2_harmony = AddModuleScore(GB2_harmony, features=neftel_list, name=names(neftel_list), search=T)
colnames(GB2_harmony@meta.data)[16:21] = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")

GB2_harmony@meta.data$neftel_state = c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")[max.col(GB2_harmony@meta.data[, c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state")], ties.method="first")]

f = FeaturePlot(GB2_harmony, features=c("AC_state", "OPC_state", "NPC1_state", "NPC2_state", "MES1_state", "MES2_state"), reduction="umap.harmony", keep_scale="all")
ggsave("seurat/GB2_harmony_neftel_states.pdf", f, width=10, height=10)

# Plot in Neftel space
source(paste0(script_path, "plot_neftel_space.R"))

GB2_harmony@meta.data$D = apply(GB2_harmony@meta.data, 1, calculate_D)
GB2_harmony@meta.data$x = apply(GB2_harmony@meta.data, 1, calculate_x)

plot_neftel_space(GB2_harmony@meta.data, "neftel_state") + scale_color_manual(values=neftel_pal)

save(GB2_harmony, file="seurat/GB2_harmony.Robj")

##########################
# Plot gene sets
fts = c(sojka_list, mp_list[2:13], richards_list[c(1:6, 64:65, 68:79)], list(developmental_gsc_top=developmental_gsc_top, injury_gsc_top=injury_gsc_top, ncc=ncc))

GB2_harmony = AddModuleScore(GB2_harmony, features=fts, name=names(fts), search=T)

plot_dot_plot_facet_by_sample <- function(seurat_obj, module_names, color_palette = "Spectral", minval_shift = NULL, maxval_shift = NULL, cluster_col = "harmony_clusters_0.8") {

  # Extract metadata and module scores (including chosen cluster column)
  df <- FetchData(seurat_obj, vars = c("orig.ident", cluster_col, module_names))

  # Convert to long format
  df_long <- df %>%
    pivot_longer(cols = all_of(module_names), names_to = "module", values_to = "score")

  # Convert cluster_col to symbol for tidy evaluation
  cluster_sym <- sym(cluster_col)

  # Compute per-cluster, per-sample mean score and cell count
  df_summary <- df_long %>%
    group_by(orig.ident, !!cluster_sym, module) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(), .groups = "drop")

  # Compute total cell count per orig.ident (for percentage scaling)
  total_cells <- df_long %>%
    group_by(orig.ident) %>%
    summarise(total_cells = n(), .groups = "drop")

  df_summary <- df_summary %>%
    left_join(total_cells, by = "orig.ident") %>%
    mutate(perc_cells = (n_cells / total_cells) * 100)

  # --- Determine module order based on selected sample ---
  orig_levels <- sort(unique(df_summary$orig.ident))
  # Pick second sample (as in your code)
  last_ident <- orig_levels[2]

  module_order <- df_summary %>%
    filter(orig.ident == last_ident) %>%
    group_by(module) %>%
    summarise(mean_across_clusters = mean(mean_score, na.rm = TRUE)) %>%
    arrange(desc(mean_across_clusters)) %>%
    pull(module)

  # Apply module order (same for all panels)
  df_summary$module <- factor(df_summary$module, levels = module_order)

  # Rename cluster column for consistent plotting
  colnames(df_summary)[colnames(df_summary) == cluster_col] <- "cluster"

  # Order clusters numerically if possible
  df_summary$cluster <- factor(df_summary$cluster, levels = sort(unique(as.numeric(as.character(df_summary$cluster)))))

  # Auto-detect limits
  minval <- min(df_summary$mean_score, na.rm = TRUE)
  maxval <- max(df_summary$mean_score, na.rm = TRUE)
  if (!is.null(minval_shift)) minval <- minval_shift
  if (!is.null(maxval_shift)) maxval <- maxval_shift
  message(glue("Color scale limits: {round(minval, 3)} to {round(maxval, 3)}"))

  # Plot
  p <- ggplot(df_summary, aes(x = cluster, y = module)) +
    geom_point(aes(size = perc_cells, color = mean_score)) +
    scale_color_distiller(palette = color_palette, limits = c(minval, maxval), oob = scales::squish, direction = -1, name = "Mean Module Score") +
    scale_size_continuous(range = c(1, 10), guide = guide_legend(title = "% of cells in cluster")) +
    facet_grid(. ~ orig.ident, scales = "free_x", space = "free_x") +
    theme_linedraw(base_size = 13) +
    labs(x = "Cluster", y = "Module (Gene Set)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right",
      panel.spacing.x = unit(0.3, "lines"))

  return(p)
}

p <- plot_dot_plot_facet_by_sample(GB2_harmony, c("S.Score", "G2M.Score", names(fts)), minval_shift=-0.05, maxval_shift=0.2)
ggsave("seurat/GB2_integrated_geneset_dotplot_clusters.pdf", p, width = 17, height = 17)


# Look at gene sets in Numbat clones
nb = Numbat$new(out_dir="numbat/multisample_results", i=3)
clone_post = nb$clone_post
GB2_harmony$cna_clone = clone_post$clone_opt

p <- plot_dot_plot_facet_by_sample(GB2_harmony, c("S.Score", "G2M.Score", names(fts)), minval_shift=-0.1, maxval_shift=0.3, cluster_col="cna_clone")
ggsave("seurat/GB2_integrated_geneset_dotplot_cnaclones.pdf", p, width = 15, height = 17)


########################
# DE between clusters
# SCP approach 
GB2_harmony = RunDEtest(GB2_harmony, group_by="harmony_clusters_0.8", assay="RNA", only.pos=F)
GB2_harmony = RunGSEA(GB2_harmony, group_by="harmony_clusters_0.8", TERM2GENE=hallmark_sets, DE_threshold="TRUE") # "p_val_adj < 0.05 & avg_log2FC > 0.25")

g = GSEAPlot(GB2_harmony, group_by="harmony_clusters_0.8", plot_type="comparison", db="custom", direction="both")
ggsave("seurat/GB2_gsea_dotplot_scp_clusters0.8.pdf", g, width=12, height=4)

GB2_harmony = RunDEtest(GB2_harmony, group_by="cna_clone", assay="RNA", only.pos=F)
GB2_harmony = RunGSEA(GB2_harmony, group_by="cna_clone", TERM2GENE=hallmark_sets, DE_threshold="TRUE") # "p_val_adj < 0.05 & avg_log2FC > 0.25")

g = GSEAPlot(GB2_harmony, group_by="cna_clone", plot_type="comparison", db="custom", direction="both")
ggsave("seurat/GB2_gsea_dotplot_scp_cnaclones.pdf", g, width=12, height=4)


plot_SCP_GSEA_dotplot_facet_by_sample <- function(srt, gsea_tool = "GSEA_cluster_wilcox", color_palette = "Spectral", padj_cutoff = 0.05, zero_center = TRUE, cluster_col = "harmony_clusters_0.8") {

  # 1) extract gsea results
  res_list <- srt@tools[[gsea_tool]]$results
  if (is.null(res_list)) stop(paste0("No GSEA results found in srt@tools$", gsea_tool))

  res_df <- purrr::map_dfr(names(res_list), function(cl) {
    r <- res_list[[cl]]
    if (inherits(r, "gseaResult")) {
      df <- as.data.frame(r@result)
      df$cluster <- cl
      return(df)
    } else return(NULL)
  })
  if (nrow(res_df) == 0) stop("No valid GSEA result data found.")

  # normalize cluster names from GSEA ("10-custom" -> "10")
  res_df$cluster_clean <- gsub("-.*$", "", res_df$cluster)
  res_df$cluster_clean <- as.character(res_df$cluster_clean)

  # 2) prepare cluster map from metadata
  if (!cluster_col %in% colnames(srt@meta.data)) stop(paste("cluster_col", cluster_col, "not found"))

  # read metadata cluster vector as character (preserve order of appearance)
  meta_clusters <- as.character(srt@meta.data[[cluster_col]])
  meta_orig <- as.character(srt@meta.data$orig.ident)

  # compute counts per cluster x sample
  cluster_map <- tibble::tibble(orig.ident = meta_orig, cluster_clean = meta_clusters) %>%
    group_by(orig.ident, cluster_clean) %>%
    summarise(n_cells = n(), .groups = "drop")

  total_cells <- tibble::tibble(orig.ident = meta_orig) %>%
    group_by(orig.ident) %>%
    summarise(total_cells = n(), .groups = "drop")

  cluster_map <- left_join(cluster_map, total_cells, by = "orig.ident") %>%
    mutate(cluster_clean = as.character(cluster_clean),
           perc_cells = (n_cells / total_cells) * 100)

  # 3) join GSEA and cluster_map (many-to-many OK)
  res_df <- left_join(res_df, cluster_map, by = "cluster_clean", relationship = "many-to-many")

  # significance flag
  res_df$is_signif <- res_df$p.adjust < padj_cutoff

  # 4) determine desired cluster ordering from metadata
  # try numeric ordering if all cluster labels are numeric-like
  unique_meta_clusters <- unique(as.character(srt@meta.data[[cluster_col]]))

  # attempt numeric conversion
  numeric_clusters <- suppressWarnings(as.numeric(unique_meta_clusters))
  if (!any(is.na(numeric_clusters))) {
    # numeric ordering
    desired_levels <- as.character(sort(numeric_clusters))
    message("Using numeric ordering for clusters based on metadata.")
  } else if (is.factor(srt@meta.data[[cluster_col]])) {
    # if the metadata column is a factor, preserve its levels
    desired_levels <- as.character(levels(srt@meta.data[[cluster_col]]))
    message("Using factor levels from metadata for cluster ordering.")
  } else {
    # fallback: preserve the order of first appearance in metadata
    desired_levels <- unique_meta_clusters
    message("Using unique appearance order from metadata for cluster ordering.")
  }

  # restrict desired_levels to only those present in the GSEA data (otherwise ggplot will show empty x's)
  present_levels <- intersect(desired_levels, unique(res_df$cluster_clean))
  if (length(present_levels) == 0) {
    # if nothing intersects, fall back to ordering inferred from res_df
    present_levels <- sort(unique(res_df$cluster_clean))
    message("No overlap between metadata clusters and GSEA cluster names; using GSEA cluster order fallback.")
  }

  # apply ordering
  res_df$cluster_clean <- factor(res_df$cluster_clean, levels = present_levels)

  # 5) row ordering by second sample (or first if only one)
  orig_levels <- sort(unique(res_df$orig.ident))
  if (length(orig_levels) >= 2) {
    ref_sample <- orig_levels[2]
  } else {
    ref_sample <- orig_levels[1]
    message("Only one sample detected — using it for row ordering.")
  }

  row_order <- res_df %>%
    filter(orig.ident == ref_sample) %>%
    group_by(Description) %>%
    summarise(mean_NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_NES)) %>%
    pull(Description)

  # if some Description are missing in the reference sample, append them
  all_terms <- unique(res_df$Description)
  missing_terms <- setdiff(all_terms, row_order)
  row_order <- c(row_order, missing_terms)

  res_df$Description <- factor(as.character(res_df$Description), levels = row_order)

  # 6) color limits
  minval <- min(res_df$NES, na.rm = TRUE)
  maxval <- max(res_df$NES, na.rm = TRUE)
  if (zero_center) midval <- 0 else midval <- (minval + maxval) / 2


  # --- determine desired_levels robustly from metadata (numeric if possible) ---
  unique_meta_clusters <- unique(as.character(srt@meta.data[[cluster_col]]))

  numeric_try <- suppressWarnings(as.numeric(unique_meta_clusters))
  if (all(!is.na(numeric_try))) {
    desired_levels <- as.character(sort(numeric_try))
  } else if (is.factor(srt@meta.data[[cluster_col]])) {
    desired_levels <- as.character(levels(srt@meta.data[[cluster_col]]))
  } else {
    desired_levels <- unique_meta_clusters  # first-seen order
  }

  # keep only levels that are present in res_df (avoid empty levels)
  present_levels <- intersect(desired_levels, unique(as.character(res_df$cluster_clean)))

  # If present_levels ends up empty, fall back to res_df order
  if (length(present_levels) == 0) {
    present_levels <- unique(as.character(res_df$cluster_clean))
    message("No overlap between metadata clusters and GSEA results; using GSEA cluster order.")
  }

  # enforce factor levels in the data
  res_df$cluster_clean <- factor(as.character(res_df$cluster_clean), levels = present_levels)

  # 7) plot
  p <- ggplot(res_df, aes(x = cluster_clean, y = Description)) +
    # outline for significant
    geom_point(data = subset(res_df, is_signif), aes(size = perc_cells), shape = 21, stroke = 0.9, color = "black", fill = NA) +
    # main colored points
    geom_point(aes(size = perc_cells, color = NES), shape = 16) +
    facet_grid(. ~ orig.ident, scales = "free_x", space = "free_x") +
    scale_x_discrete(limits = present_levels) +   # <--- force order here
    scale_color_distiller(palette = color_palette, direction = -1, limits = c(minval, maxval), oob = scales::squish) +
    scale_size_continuous(range = c(2, 10), guide = guide_legend(title = "% of cells in cluster")) +
    theme_linedraw(base_size = 13) +
    labs(x = "Cluster", y = "Enriched Term", color = "NES (Enrichment Score)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "right", panel.spacing.x = unit(0.3, "lines"))

  return(p)
}

p <- plot_SCP_GSEA_dotplot_facet_by_sample(
  srt = GB2_harmony,
  gsea_tool = "GSEA_harmony_clusters_0.8_wilcox",
  padj_cutoff = 0.05,
  cluster_col="harmony_clusters_0.8"
)
ggsave("seurat/GB2_gsea_dotplot_clusters.pdf", p, width=19, height=17)

p <- plot_SCP_GSEA_dotplot_facet_by_sample(
  srt = GB2_harmony,
  gsea_tool = "GSEA_cna_clone_wilcox",
  padj_cutoff = 0.05,
  cluster_col="cna_clone"
)
ggsave("seurat/GB2_gsea_dotplot_cnaclones.pdf", p, width=19, height=17)
