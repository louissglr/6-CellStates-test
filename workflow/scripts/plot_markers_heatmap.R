library(CoGAPS)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)
library(tidyr)
library(tibble)

# Snakemake inputs/outputs
cogaps_file  <- snakemake@input[["cogaps_rds"]]
markers_file <- snakemake@input[["markers_excel"]]
output_pdf   <- snakemake@output[["pdf"]]

sample_name <- snakemake@wildcards[["sample"]]
n_patterns  <- snakemake@wildcards[["npatterns"]]

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# =========================
# Lecture des données
# =========================
cogaps <- readRDS(cogaps_file)

# Lecture du fichier markers (avec colonne Rank)
markers_all <- read.xlsx(markers_file, sheet = 1) |> as_tibble()
scores_all  <- read.xlsx(markers_file, sheet = 2) |> as_tibble()

# =========================
# Feature loadings (A-matrix)
# =========================
A_matrix <- as_tibble(cogaps@featureLoadings, rownames = "Gene") |>
  pivot_longer(
    cols = starts_with("Pattern"),
    names_to = "Pattern",
    values_to = "amplitude_features"
  )

# =========================
# Top 20 markers par pattern
# =========================
top_markers <- markers_all |>
  pivot_longer(
    cols = starts_with("Pattern"),
    names_to = "Pattern",
    values_to = "Gene"
  ) |>
  filter(!is.na(Gene)) |>
  filter(Rank <= 10)

# Jointure avec les feature loadings
top_markers <- top_markers |>
  left_join(A_matrix, by = c("Gene", "Pattern"))

# =========================
# Construction matrices pour Heatmap featureLoadings
# =========================
feature_loadings <- as.data.frame(cogaps@featureLoadings)

top_genes <- markers_all |>
  arrange(Rank) |>
  select(-Rank) |>
  lapply(function(x) x[10:0])

patterns <- colnames(feature_loadings)
genes_ordered <- unlist(top_genes, use.names = FALSE)

mat <- feature_loadings[genes_ordered, patterns]
mat <- as.matrix(mat)
mat_log10 <- log10(mat + 1e-6)

# =========================
# Construction matrices pour Heatmap Scores
# =========================
# Les gènes sont identiques à ceux utilisés pour featureLoadings
scores_matrix <- scores_all |>
  column_to_rownames("Gene") |>  # Mettre la colonne Gene en rownames
  select(all_of(patterns))       # S'assurer que l'ordre des colonnes est correct

mat_scores <- scores_matrix[genes_ordered, ]
mat_scores <- as.matrix(mat_scores)
mat_scores_log10 <- log10(mat_scores + 1e-6)

# =========================
# Annotations Heatmap
# =========================
pattern_colors <- structure(rainbow(length(patterns)), names = patterns)

ha <- HeatmapAnnotation(
  Pattern = patterns,
  col = list(Pattern = pattern_colors),
  annotation_name_side = "left"
)

# =========================
# Export PDF
# =========================
pdf(output_pdf, width = 12, height = 14)

# Heatmaps featureLoadings
Heatmap(
  mat,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Top 10 marker genes per pattern",
  row_title = "Genes",
  top_annotation = ha,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  name = "featureLoadings",
  row_gap = unit(1, "mm")
)

Heatmap(
  mat_log10,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Top 10 marker genes per pattern (log10)",
  row_title = "Genes",
  top_annotation = ha,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  name = "log10(featureLoadings)",
  row_gap = unit(1, "mm")
)

score_col_fun <- colorRamp2(
  breaks = c(min(mat_scores, na.rm = TRUE), max(mat_scores, na.rm = TRUE)),
  colors = c("red", "blue")  # rouge = bas, vert = haut
)


# Heatmaps Scores
Heatmap(
  mat_scores,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Top 10 marker genes per pattern (PatternScores)",
  row_title = "Genes",
  top_annotation = ha,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  name = "PatternScores",
  row_gap = unit(1, "mm"),
  col = score_col_fun
)

dev.off()
message("Heatmaps saved: ", output_pdf)
