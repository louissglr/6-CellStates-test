library(ggplot2)
library(dplyr)
library(ComplexUpset)

# ============================================================
# ================== UTILITAIRE PDF VIDE ====================
# ============================================================

create_empty_pdf <- function(filename, title, message) {
  pdf(filename, width = 10, height = 7)
  plot.new()
  text(0.5, 0.6, title, cex = 1.4, font = 2)
  text(0.5, 0.4, message, cex = 1)
  dev.off()
}

# ============================================================
# ====================== VARIABLES ===========================
# ============================================================

sample_name <- snakemake@wildcards[["sample"]]
n_patterns  <- as.integer(snakemake@wildcards[["npatterns"]])

input_csv_enrichment  <- snakemake@input[["csv_enrich"]]
input_csv_overrep     <- snakemake@input[["csv_overr"]]

output_pdf_enrichment <- snakemake@output[["pdf_inter_enrich"]]
output_pdf_overrep    <- snakemake@output[["pdf_inter_overr"]]

output_csv_enrichment <- snakemake@output[["csv_inter_enrich"]]
output_csv_overrep    <- snakemake@output[["csv_inter_overr"]]

# ============================================================
# =============== 1) HALLMARKS ENRICHMENT ====================
# ============================================================

hallmarks_df  <- read.csv(input_csv_enrichment, stringsAsFactors = FALSE)
hallmarks_sig <- hallmarks_df[hallmarks_df$padj < 0.05, ]

if (nrow(hallmarks_sig) == 0) {

  message("⚠️ Enrichment: aucun hallmark significatif")

  create_empty_pdf(
    output_pdf_enrichment,
    paste0("UpSet plot - ", sample_name,
           " (", n_patterns, " patterns) - Enrichment"),
    "Aucun hallmark significatif.\nUpSet plot non applicable."
  )

  write.csv(
    data.frame(combination = character(0), hallmarks = character(0)),
    output_csv_enrichment,
    row.names = FALSE,
    quote = TRUE
  )

} else {

  hallmarks_binary <- as.data.frame.matrix(
    table(hallmarks_sig$gene.set, hallmarks_sig$Pattern) > 0
  )

  hallmarks_binary <- cbind(
    gene.set = rownames(hallmarks_binary),
    hallmarks_binary
  )
  rownames(hallmarks_binary) <- NULL

  pattern_cols <- colnames(hallmarks_binary)[-1]

  # ---------------- UpSet plot ----------------
  if (length(pattern_cols) >= 2) {

    upset_plot_enrichment <- ComplexUpset::upset(
      hallmarks_binary,
      intersect = pattern_cols,
      name = "Hallmarks"
    ) +
      ggtitle(paste0(
        "UpSet plot - ", sample_name,
        " (", n_patterns, " patterns) - Enrichment"
      )) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )

    ggsave(
      filename = output_pdf_enrichment,
      plot = upset_plot_enrichment,
      width = 10,
      height = 7
    )

  } else {

    message("⚠️ Enrichment: < 2 patterns significatifs, UpSet non généré")

    create_empty_pdf(
      output_pdf_enrichment,
      paste0("UpSet plot - ", sample_name,
             " (", n_patterns, " patterns) - Enrichment"),
      "Moins de deux patterns significatifs.\nUpSet plot non applicable."
    )
  }

  # ---------------- CSV intersections ----------------
  combinations <- apply(
    hallmarks_binary[, pattern_cols, drop = FALSE],
    1,
    function(x) paste(pattern_cols[which(x)], collapse = "|")
  )

  valid_idx <- combinations != ""

  result_df_enrichment <- if (any(valid_idx)) {
    gene_sets <- hallmarks_binary$gene.set[valid_idx]
    split_list <- split(gene_sets, combinations[valid_idx])

    data.frame(
      combination = names(split_list),
      hallmarks   = sapply(split_list, paste, collapse = ", "),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(combination = character(0), hallmarks = character(0))
  }

  write.csv(
    result_df_enrichment,
    output_csv_enrichment,
    row.names = FALSE,
    quote = TRUE
  )
}

# ============================================================
# ========== 2) HALLMARK OVERREPRESENTATION ==================
# ============================================================

overrep_df  <- read.csv(input_csv_overrep, stringsAsFactors = FALSE)
overrep_sig <- overrep_df[overrep_df$padj < 0.05, ]

if (nrow(overrep_sig) == 0) {

  message("⚠️ Overrepresentation: aucun hallmark significatif")

  create_empty_pdf(
    output_pdf_overrep,
    paste0("UpSet plot - ", sample_name,
           " (", n_patterns, " patterns) - Overrepresentation"),
    "Aucun hallmark significatif.\nUpSet plot non applicable."
  )

  write.csv(
    data.frame(combination = character(0), hallmarks = character(0)),
    output_csv_overrep,
    row.names = FALSE,
    quote = TRUE
  )

} else {

  overrep_binary <- as.data.frame.matrix(
    table(overrep_sig$gene.set, overrep_sig$Pattern) > 0
  )

  overrep_binary <- cbind(
    gene.set = rownames(overrep_binary),
    overrep_binary
  )
  rownames(overrep_binary) <- NULL

  pattern_cols_overrep <- colnames(overrep_binary)[-1]

  # ---------------- UpSet plot ----------------
  if (length(pattern_cols_overrep) >= 2) {

    upset_plot_overrep <- ComplexUpset::upset(
      overrep_binary,
      intersect = pattern_cols_overrep,
      name = "Hallmarks"
    ) +
      ggtitle(paste0(
        "UpSet plot - ", sample_name,
        " (", n_patterns, " patterns) - Overrepresentation"
      )) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )

    ggsave(
      filename = output_pdf_overrep,
      plot = upset_plot_overrep,
      width = 10,
      height = 7
    )

  } else {

    message("⚠️ Overrepresentation: < 2 patterns significatifs, UpSet non généré")

    create_empty_pdf(
      output_pdf_overrep,
      paste0("UpSet plot - ", sample_name,
             " (", n_patterns, " patterns) - Overrepresentation"),
      "Moins de deux patterns significatifs.\nUpSet plot non applicable."
    )
  }

  # ---------------- CSV intersections ----------------
  combinations_overrep <- apply(
    overrep_binary[, pattern_cols_overrep, drop = FALSE],
    1,
    function(x) paste(pattern_cols_overrep[which(x)], collapse = "|")
  )

  valid_idx_overrep <- combinations_overrep != ""

  result_df_overrep <- if (any(valid_idx_overrep)) {
    gene_sets <- overrep_binary$gene.set[valid_idx_overrep]
    split_list <- split(gene_sets, combinations_overrep[valid_idx_overrep])

    data.frame(
      combination = names(split_list),
      hallmarks   = sapply(split_list, paste, collapse = ", "),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(combination = character(0), hallmarks = character(0))
  }

  write.csv(
    result_df_overrep,
    output_csv_overrep,
    row.names = FALSE,
    quote = TRUE
  )
}

message("✅ Upset plot Enrichment et Overrepresentation terminés")
