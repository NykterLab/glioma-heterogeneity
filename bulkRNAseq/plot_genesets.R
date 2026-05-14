# Script for plotting bulk RNAseq

library(openxlsx)
library("org.Hs.eg.db")
library(pheatmap)
library(reshape2)
library(randomcoloR)
library(ggplot2)

setwd("data/project/bulkRNAseq/")

bulk = read.xlsx("normalized_expressions.xlsx")
rownames(bulk) = bulk[,1]
bulk = bulk[2:ncol(bulk)]
bulk = log2(bulk + 1)

source("../scRNAseq/load_genesets.R")

#########################
# Plot genesets
# look for gene aliases present in the dataset
check_bulk_geneset <- function(geneset, bulk) {

    missing_genes <- geneset[!geneset %in% rownames(bulk)]
    if (length(missing_genes) == 0) {
        message("All genes found in bulk rownames")
        return(geneset)
    }
    message("Missing genes: ", paste(missing_genes, collapse = ", "))

    gene_synonyms_sym <- tryCatch(
        AnnotationDbi::select(org.Hs.eg.db, keys = missing_genes, keytype = "SYMBOL", columns = c("SYMBOL", "ALIAS")),
        error = function(e) NULL)

    gene_synonyms_alias <- tryCatch(
        AnnotationDbi::select(org.Hs.eg.db, keys = missing_genes, keytype = "ALIAS", columns = c("SYMBOL", "ALIAS")),
        error = function(e) NULL)

    # Merge
    gene_synonyms <- dplyr::bind_rows(gene_synonyms_sym, gene_synonyms_alias)

    if (is.null(gene_synonyms) || nrow(gene_synonyms) == 0) {
        message("No valid SYMBOL or ALIAS found for missing genes.")
        geneset <- geneset[geneset %in% rownames(bulk)]
        return(geneset)
    }

    # Clean
    gene_synonyms <- gene_synonyms %>%
        dplyr::filter(!is.na(SYMBOL), !is.na(ALIAS)) %>%
        dplyr::distinct()

    # Build alias → bulk rowname match
    alias_map <- gene_synonyms %>%
        dplyr::mutate(
            replacement = dplyr::case_when(
                ALIAS %in% rownames(bulk) ~ ALIAS,
                SYMBOL %in% rownames(bulk) ~ SYMBOL,
                TRUE ~ NA_character_
            )
        ) %>%
        dplyr::filter(!is.na(replacement)) %>%
        dplyr::distinct(SYMBOL, replacement, .keep_all = TRUE)

    # Apply replacements
    for (i in seq_len(nrow(alias_map))) {
        geneset[geneset == alias_map$SYMBOL[i]] <- alias_map$replacement[i]
        geneset[geneset == alias_map$ALIAS[i]] <- alias_map$replacement[i]
    }

    # Report failures
    remaining_missing <- geneset[!geneset %in% rownames(bulk)]
    if (length(remaining_missing) > 0) {
        message("Still missing after alias matching: ", paste(remaining_missing, collapse = ", "))
    } else {
        message("All genes resolved into bulk rownames")
    }

    # Keep only genes present in bulk
    geneset <- geneset[geneset %in% rownames(bulk)]
    return(geneset)
}

calculate_score <- function(df) {
    # keep genes with variance > 0
    df <- df[apply(df, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
    cat("n_genes after filtering:", nrow(df), "\n")

    # z-scale per gene
    z_score <- t(scale(t(df)))
    cat("Any NaN in z_score?:", any(is.nan(z_score)), "\n")

    # activity score
    activity.z <- colSums(z_score / sqrt(nrow(df)), na.rm = TRUE)
    return(activity.z)
}

fts = c(sojka_list, mp_list[2:13], richards_list[c(1:6, 68:79)], 
    list(developmental_gsc_top=developmental_gsc_top, injury_gsc_top=injury_gsc_top, hypoxia=c("CA9", "VEGFA", "ADM", "PDK1"),
    cellular=ivy[ivy$ROI=="CT", "gene"], MVP=ivy[ivy$ROI=="CTmvp", "gene"], pseudopalisading=ivy[ivy$ROI=="CTpan", "gene"], leading_edge=ivy[ivy$ROI=="LE", "gene"]))

fts_checked = lapply(fts, function(gs) check_bulk_geneset(gs, bulk))

# tissue only
bulk_tissue = bulk[,c(1:12)]
scores_tissue <- lapply(fts, function(gs) {
    gs2 <- check_bulk_geneset(gs, bulk_tissue)
    if (length(gs2) == 0) {
        warning("Gene set is empty after alias correction")
        return(rep(NA, ncol(bulk_tissue)))
    }
    df <- bulk_tissue[gs2, , drop = FALSE]
    calculate_score(df)
})
scores_tissue <- do.call(rbind, scores_tissue)
rownames(scores_tissue) <- names(fts)

pdf("all_sets_bulk_tissue_rowscale.pdf", width=8, height=20)
pheatmap(scores_tissue[,c(5:8,1:4,9:12)], cluster_rows=F, cluster_cols=F, scale="row")
dev.off()

fts_supp = fts_checked
fts_main = fts_checked[c(1:19, 41:44)]

# test gene sets in one sample vs others
test_sample_vs_others <- function(df, groups) {
  # df = gene expression matrix (genes x samples)
  # groups = list of 3 groups of 4 column indices
  results <- matrix("", nrow = 1, ncol = ncol(df))
  colnames(results) <- colnames(df)
  
  for (g in groups) {
    for (i in g) {
      others <- setdiff(g, i)
      p <- tryCatch(
        t.test(df[, i], df[, others])$p.value,
        error = function(e) NA
      )
      if (!is.na(p) && p < 0.01) {
        results[1, i] <- "*"
      }
    }
  }
  return(results)
}

sig_mat <- lapply(fts_supp, function(gs) {
  gs2 <- check_bulk_geneset(gs, bulk)
  if (length(gs2) == 0) {
    return(rep("", ncol(bulk)))
  }
  
  df <- bulk[gs2, , drop = FALSE]
  df <- df[apply(df, 1, var) > 0, , drop = FALSE]
  z <- t(scale(t(df)))
  test_sample_vs_others(z, list(1:4,5:8,9:12))
})

sig_mat <- do.call(rbind, sig_mat)
rownames(sig_mat) <- names(fts_supp)

pdf("all_sets_bulk_tissue_rowscale_signif_supp.pdf", width=8, height=14)
pheatmap(scores_tissue[names(fts_supp), c(5:8,1:4,9:12)],cluster_rows = T,cluster_cols = FALSE,scale = "row",
        display_numbers = sig_mat[, c(5:8,1:4,9:12)],number_color = "black",fontsize_number = 14)
dev.off()

sig_mat <- lapply(fts_main, function(gs) {
  gs2 <- check_bulk_geneset(gs, bulk)
  if (length(gs2) == 0) {
    return(rep("", ncol(bulk)))
  }
  
  df <- bulk[gs2, , drop = FALSE]
  df <- df[apply(df, 1, var) > 0, , drop = FALSE]
  z <- t(scale(t(df)))
  test_sample_vs_others(z, list(1:4,5:8,9:12))
})

sig_mat <- do.call(rbind, sig_mat)
rownames(sig_mat) <- names(fts_main)

# cluster metaprograms and sojka sets separately
mat1 <- scores_tissue[names(fts_main)[1:7], ]
mat2 <- scores_tissue[names(fts_main)[8:19], ]
mat3 = scores_tissue[names(fts_main)[20:23], ]

ord1 <- hclust(dist(mat1))$order
ord2 <- hclust(dist(mat2))$order
ord3 = hclust(dist(mat3))$order

mat_new <- rbind(mat1[ord1, ], mat2[ord2, ], mat3[ord3, ])

sig_new <- rbind(sig_mat[1:7, ][ord1, ], sig_mat[8:19, ][ord2, ], sig_mat[20:23,][ord3, ])

pdf("all_sets_bulk_tissue_rowscale_signif_main.pdf", width=8, height=10)
pheatmap(mat_new[, c(5:8,1:4,9:12)], cluster_rows = FALSE, cluster_cols = FALSE, scale = "row",
  display_numbers = sig_new[, c(5:8,1:4,9:12)], number_color="black", fontsize_number=14)
dev.off()

############################
# Gene expression in genomically altered genes
alterations = read.table("altered_genes.txt", stringsAsFactors=F)[[1]]

pdf("bulk_rna/altered_genes.pdf", width=5, height=8)
pheatmap(bulk[alterations, c(5:8,1:4,9:15)], cluster_rows=F, cluster_cols=F, scale="row")
dev.off()

############################
# Deconvolution results visualization

deconv = read.table("deconvolution_celltype_proportions.txt", sep="\t", header=T, stringsAsFactors=F)

df_long <- melt(as.matrix(dc))
colnames(df_long) <- c("Category", "Box", "Proportion")
df_long$Category <- factor(df_long$Category)

n_cat <- length(levels(df_long$Category))
cols <- distinctColorPalette(n_cat)
cols_named <- setNames(cols, levels(df_long$Category))

p = ggplot(df_long, aes(x = Box, y = Proportion, fill = Category)) +
  #geom_bar(stat = "identity") +
  geom_col(position = "fill") +
  scale_fill_manual(values=cols_named) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic()

ggsave("deconvolution_boxplot.pdf", p, width=8, height=6)

