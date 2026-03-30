## ============================================================
## MAYA Pathway Activity Workflow
## Snakemake Compatible
## ============================================================

if (!requireNamespace("MAYA", quietly = TRUE)) {
    remotes::install_github("one-biosciences/maya")
}
library(MAYA)
library(Seurat)
library(msigdbr)
library(dplyr)
library(tibble)

## ============================================================
## Inputs / Outputs
## ============================================================

seurat_file   <- snakemake@input$seurat
clusters_file <- snakemake@input$subclone

heatmap1   <- snakemake@output$heatmap1
heatmap2   <- snakemake@output$heatmap2
umap_file  <- snakemake@output$umap
top_pdf    <- snakemake@output$top_genes
spec_pdf   <- snakemake@output$specificity

## ============================================================
## Load Seurat Object
## ============================================================

message("Loading Seurat object")

seurat_obj <- tryCatch(
  readRDS(seurat_file),
  error = function(e) stop("Failed to load Seurat object: ", e$message)
)

count_mat <- GetAssayData(seurat_obj,
                          assay = "RNA",
                          layer = "counts")

message("Matrix dimension: ",
        nrow(count_mat), " genes x ",
        ncol(count_mat), " cells")

## ============================================================
## Load Subclone Metadata
## ============================================================

message("Loading metadata")

meta <- read.table(clusters_file) |>
  rename(
    cells_id = V1,
    subclone = V2
  ) |>
  mutate(barcode = sub(".*_", "", cells_id)) |>
  filter(barcode %in% colnames(seurat_obj)) |>
  column_to_rownames("barcode")

## ============================================================
## Pathway Activity Analysis
## ============================================================

message("Running pathway analysis")

activity_summary <- MAYA_pathway_analysis(
  expr_mat = count_mat,
  modules_list = "hallmark",
  is_logcpm = FALSE
)

activity_mat_scale <- scale_0_1(activity_summary$activity_matrix)

modules <- names(activity_summary$PCA_obj)

## ============================================================
## Heatmap 1 (cluster rows + cols)
## ============================================================

message("Generating heatmap 1")

png(heatmap1,
    width = 1200,
    height = 1000,
    res = 150)

plot_heatmap_activity_mat(
  activity_mat = activity_mat_scale,
  meta = meta,
  annot_name = "subclone",
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  clustering_distance = "euclidean",
  clustering_method = "ward.D2"
)

dev.off()

## ============================================================
## Heatmap 2 (rows only)
## ============================================================

message("Generating heatmap 2")

png(heatmap2,
    width = 1200,
    height = 1000,
    res = 150)

plot_heatmap_activity_mat(
  activity_mat = activity_mat_scale,
  meta = meta,
  annot_name = "subclone",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance = "euclidean",
  clustering_method = "ward.D2"
)

dev.off()

## ============================================================
## UMAP
## ============================================================

message("Generating UMAP")

png(umap_file,
    width = 1000,
    height = 800,
    res = 150)

plot_umap_annot(
  umap = activity_summary$umap,
  labels = as.factor(meta$subclone),
  title = "Clone annotation - HALLMARK"
)

dev.off()

## ============================================================
## LogCPM normalization
## ============================================================

message("Running logCPM")

options(future.globals.maxSize = 8 * 1024^3)

logcpm <- logcpmNormalization(count_mat)

## ============================================================
## Top contributing genes (PDF)
## ============================================================

message("Plotting top contributing genes")

pdf(top_pdf,
    width = 10,
    height = 8)

for (mod in modules) {
  try({
    plot_heatmap_pathway_top_contrib_genes(
      expr_mat = logcpm,
      PCA_object = activity_summary$PCA_obj,
      module = mod,
      n = 20,
      meta = meta,
      annot_name = "subclone"
    )
  }, silent = TRUE)
}

dev.off()

## ============================================================
## Pathway specificity (PDF)
## ============================================================

message("Plotting pathway specificity")

pdf(spec_pdf,
    width = 8,
    height = 6)

for (mod in modules) {
  try({
    plot_pathway_specificity(
      PCA_object = activity_summary$PCA_obj,
      module = mod,
      meta = meta,
      annot_name = "subclone"
    )
  }, silent = TRUE)
}

dev.off()

message("Workflow completed successfully")
