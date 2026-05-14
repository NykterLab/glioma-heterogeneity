# Script for Numbat results postprocessing
# Author: Serafiina Jaatinen

library(ggplot2)
library(numbat)
library(dplyr)
library(glue)
library(data.table)
library(ggtree)
library(stringr)
library(tidygraph)
library(patchwork)
library(Seurat)
library(SCP)
library(RColorBrewer)

numbat_result_path = "data/project/scRNAseq/numbat/multisample_results/"
setwd(numbat_result_path)

nb = Numbat$new(out_dir=numbat_result_path, i=3)


annots = data.frame(cell=nb$clone_post$cell)
annots$sample = substr(annots$cell, 1, 8)

sample_pal = c("GB2_L1"="palevioletred1", "GB2_L2"="deepskyblue", "GB2_L3"="darkolivegreen2")

pdf(paste0(numbat_result_path, "phylo_heatmap_pmin0.7.pdf"), width=8, height=5)
nb$plot_phylo_heatmap(clone_bar=T, clone_stack=F, clone_line=F, show_phylo=F, p_min=0.7,
                      annot=annots, pal_annot=sample_pal, pal_clone=brewer.pal(12, "Paired"))
dev.off()

pdf(paste0(numbat_result_path, "clone_profiles.pdf"), width=15, height=10)
nb$plot_clone_profile()
dev.off()

pdf(paste0(numbat_result_path, "mut_history.pdf"), width=10, height=10)
nb$plot_mut_history()
dev.off()

pdf(paste0(numbat_result_path, "sc_tree.pdf"), width=10, height=10)
nb$plot_sc_tree()
dev.off()

joint_post = nb$joint_post # posterior prob of CNV states for each segment in each cell, both expression and allele based
clone_post = nb$clone_post # sc clone assignment and tumor vs normal classification posteriors
bulk_clones = nb$bulk_clones # clone-level pseudobulk profiles based on final cell lineage tree, rows are SNPs

joint_post$origin = substr(joint_post$cell, 1, 8)
joint_post$chr_start_end = paste0(joint_post$CHROM, "_", joint_post$seg_start, "_", joint_post$seg_end)

segs = joint_post %>% distinct(chr_start_end, .keep_all=T)
segs = segs[, c("CHROM", "seg_start", "seg_end", "chr_start_end", "seg_label")]

#############################
# Add clones to seurat object

load("../../seurat/integrated/GB2_harmony.Robj")

identical(clone_post$cell, rownames(GB2_harmony@meta.data))
GB2_harmony$cna_clone = clone_post$clone_opt

d = CellDimPlot(GB2_harmony,  reduction = "umap.harmony", group.by="cna_clone", palette=brewer.pal(12, "Paired"))
ggsave(paste0(numbat_result_path, "numbat_clones.pdf"), d, width=5, height=5)

#############################
# CNAs affecting mutation history and evolution
clones = clone_post[,2:3]
clones = unique(clones)
# clone 1 normal cells
mut_cl2 = unlist(strsplit(as.character(clones[clones$clone_opt==2, "GT_opt"]), ","))
mut_cl3 = unlist(strsplit(as.character(clones[clones$clone_opt==3, "GT_opt"]), ","))
mut_cl4 = unlist(strsplit(as.character(clones[clones$clone_opt==4, "GT_opt"]), ","))
mut_cl5 = unlist(strsplit(as.character(clones[clones$clone_opt==5, "GT_opt"]), ","))
mut_cl6 = unlist(strsplit(as.character(clones[clones$clone_opt==6, "GT_opt"]), ","))
mut_cl7 = unlist(strsplit(as.character(clones[clones$clone_opt==7, "GT_opt"]), ","))
mut_cl8 = unlist(strsplit(as.character(clones[clones$clone_opt==8, "GT_opt"]), ","))
mut_cl9 = unlist(strsplit(as.character(clones[clones$clone_opt==9, "GT_opt"]), ","))
mut_cl10 = unlist(strsplit(as.character(clones[clones$clone_opt==10, "GT_opt"]), ","))
mut_cl11 = unlist(strsplit(as.character(clones[clones$clone_opt==11, "GT_opt"]), ","))
mut_cl12 = unlist(strsplit(as.character(clones[clones$clone_opt==12, "GT_opt"]), ","))

# Clonal CNAs in mut_cl2
setdiff(mut_cl3, mut_cl2)
setdiff(mut_cl4, mut_cl3)
setdiff(mut_cl5, mut_cl3)
setdiff(mut_cl6, mut_cl2)
setdiff(mut_cl7, mut_cl2)
setdiff(mut_cl8, mut_cl7)
setdiff(mut_cl9, mut_cl7)
setdiff(mut_cl10, mut_cl7)
setdiff(mut_cl11, mut_cl7)
setdiff(mut_cl12, mut_cl2)

# clone3
GB2_harmony$L1_seg_6_221 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L1_seg_6_221"]
GB2_harmony$L3_seg_16_573 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_16_573"]
GB2_harmony$L2_seg_16_690 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L2_seg_16_690"]
# clone4
GB2_harmony$L1_seg_6_223 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L1_seg_6_223"]
# clone5
GB2_harmony$L3_seg_9_379 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_9_379"]
# clone6
GB2_harmony$L2_seg_19_765_amp = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L2_seg_19_765_amp"]
# clone7
GB2_harmony$L2_seg_1_13 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L2_seg_1_13"]
GB2_harmony$L3_seg_1_5 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_1_5"]
GB2_harmony$L3_seg_19_632 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_19_632"]
# clone8
GB2_harmony$L2_seg_19_765_bamp = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L2_seg_19_765_bamp"]
# clone9
GB2_harmony$L3_seg_11_462 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_11_462"]
# clone10
GB2_harmony$L3_seg_1_9 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_1_9"]
# clone11
GB2_harmony$L3_seg_5_180 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L3_seg_5_180"]
# clone12
GB2_harmony$L2_seg_5_213 = joint_post$p_cnv[match(rownames(GB2_harmony@meta.data), joint_post$cell) & joint_post$seg=="L2_seg_5_213"]


pdf(paste0(numbat_result_path, "mut_history_cnv_prob_order.pdf"), width=5, height=4.5)
FeaturePlot(GB2_harmony, features=c("L1_seg_6_221"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_16_573"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L2_seg_16_690"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L1_seg_6_223"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_9_379"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L2_seg_19_765_amp"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L2_seg_1_13"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_1_5"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_19_632"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L2_seg_19_765_bamp"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_11_462"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_1_9"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L3_seg_5_180"), reduction="umap.harmony", order=T) +  
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
FeaturePlot(GB2_harmony, features=c("L2_seg_5_213"), reduction="umap.harmony", order=T) + 
  scale_color_gradientn(colours = c("royalblue", "white", "red3"), values = scales::rescale(c(0, 0.9, 1)),limits = c(0, 1))
dev.off()

###########################
# Neftel space

source("/data/scripts/seurat5/plot_neftel_space.R")

clone_ids <- sort(unique(GB2_harmony$cna_clone))
GB2_harmony@meta.data[paste0("clone_", clone_ids)] <-
  sapply(clone_ids, function(cl) ifelse(GB2_harmony$cna_clone == cl, "Yes", "No"))


for (i in 2:length(clone_ids)) {
  pdf(paste0(numbat_result_path, "clones_separate_neftel_space_clone", clone_ids[i],".pdf"), width=5, height=4)
  p = plot_neftel_space(GB2_harmony@meta.data, paste0("clone_", clone_ids[i])) + 
    scale_color_manual(values = c("No" = "grey80", "Yes" = brewer.pal(12, "Paired")[i-1])) +
    scale_alpha_manual(values = c("No" = 0.2, "Yes" = 1), guide="none") +
    theme(legend.position = "none")
  print(p)
  dev.off()
}

###############################
# Plot clones in separate UMAPs
meta <- GB2_harmony@meta.data
umap_df = Embeddings(GB2_harmony, "umap.harmony") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
meta$UMAP_1 <- umap_df[,2]
meta$UMAP_2 <- umap_df[,3]

table(meta$cna_clone, meta$orig.ident)

pdf(paste0(numbat_result_path, "clones_separate_umap.pdf"), width=5, height=5)
for (i in 2:length(clone_ids)) {
  colname <- paste0("clone_", clone_ids[i])
  p <- ggplot(meta, aes(x = UMAP_1, y = UMAP_2,
                        color = meta[[colname]],
                        alpha = meta[[colname]])) +
    geom_point(size = 0.1) +
    scale_color_manual(values = c("No" = "grey80", "Yes" = brewer.pal(12, "Paired")[i-1])) +
    scale_alpha_manual(values = c("No" = 0.2, "Yes" = 1)) +
    theme_classic() +
    coord_fixed() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    ggtitle(colname)
  print(p)
}
dev.off()

##############################
# Stacked barplot for clone neftel states

neftel_pal = c("MES1"="#A6CEE3", "MES2"="#1F78B4", "NPC1"="#B2DF8A", "NPC2"="#33A02C", "AC"="#FDBF6F", "OPC"="#FF7F00")

plot_stacked_prop_cna <- function(seurat_obj, fill_col, facet_by_orig = FALSE) {
    meta <- seurat_obj@meta.data

    required_cols <- c("cna_clone", fill_col)
    if (facet_by_orig) {
        required_cols <- c(required_cols, "orig.ident")
    }

    missing_cols <- setdiff(required_cols, colnames(meta))
    if (length(missing_cols) > 0) {
        stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
    }
    # Build counting table
    if (facet_by_orig) {
        df <- meta %>%
            dplyr::count(orig.ident, cna_clone, .data[[fill_col]]) %>%
            dplyr::group_by(orig.ident, cna_clone) %>%
            dplyr::mutate(prop = n / sum(n)) %>%
            dplyr::ungroup()
    } else {
        df <- meta %>%
            dplyr::count(cna_clone, .data[[fill_col]]) %>%
            dplyr::group_by(cna_clone) %>%
            dplyr::mutate(prop = n / sum(n)) %>%
            dplyr::ungroup()
    }
    # Plot
    p <- ggplot(df, aes(x = cna_clone, y = prop, fill = .data[[fill_col]])) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(x = "CNA clone", y = "Proportion of cells", fill = fill_col) +
        scale_fill_manual(values=neftel_pal, drop=F) +
        theme_classic(base_size = 13) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))

    if (facet_by_orig) {
        p <- p + facet_wrap(~ orig.ident)
    }
    return(p)
}

GB2_harmony$cna_clone <- factor(GB2_harmony$cna_clone)

p = plot_stacked_prop_cna(seurat_obj = GB2_harmony, fill_col = "neftel_state", facet_by_orig = F)
ggsave(paste0(numbat_result_path, "stacked_clone_neftel.pdf"), width=7, height=5)

