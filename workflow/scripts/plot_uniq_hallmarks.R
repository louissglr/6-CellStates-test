suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(viridis)
  library(tidytext)
  library(ggplot2)
})


input_csv_enrichment <- snakemake@input[["csv_enrich"]]
input_csv_overr      <- snakemake@input[["csv_overr"]]
output_png_enrich    <- snakemake@output[["png_enrich"]]
output_png_overr     <- snakemake@output[["png_overr"]]

sample_name <- snakemake@wildcards[["sample"]]
n_patterns  <- as.integer(snakemake@wildcards[["npatterns"]])

padj_cutoff <- 0.05
threshold   <- -10 * log10(padj_cutoff)


safe_read <- function(path){
  if(!file.exists(path)){
    stop(paste("Input file does not exist:", path))
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if(nrow(df) == 0){
    warning(paste("Empty dataframe:", path))
  }
  return(df)
}


generate_combined_plot <- function(df, title_prefix){

  required_cols <- c("padj", "gene.set", "Pattern")
  if(!all(required_cols %in% colnames(df))){
    stop("Missing required columns in dataframe.")
  }

  df <- df %>% filter(!is.na(padj)) %>% filter(padj < padj_cutoff)
  if(nrow(df) == 0){
    warning(paste0("No significant hallmarks for ", title_prefix))
    return(NULL)
  }

  df <- df %>%
    mutate(
      neg_log_padj = -10 * log10(padj),
      gene.set_short = str_remove(gene.set, "^HALLMARK_"),
      Pattern_num = suppressWarnings(as.integer(str_extract(Pattern, "\\d+")))
    ) %>%
    filter(!is.na(Pattern_num)) %>%
    mutate(
      Pattern = factor(
        Pattern,
        levels = paste0("Pattern_", sort(unique(Pattern_num)))
      )
    ) %>%
    select(-Pattern_num) %>%
    group_by(gene.set) %>%
    mutate(
      uniqueness = ifelse(n_distinct(Pattern) == 1, "Unique", "Non-unique")
    ) %>%
    ungroup()

  p_combined <- ggplot(
    df,
    aes(
      x = reorder_within(gene.set_short, neg_log_padj, Pattern),
      y = neg_log_padj,
      fill = uniqueness
    )
  ) +
    geom_col() +
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed", linewidth = 0.6) +
    scale_x_reordered() +
    coord_flip() +
    facet_grid(Pattern ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Unique" = "red3", "Non-unique" = "grey70")) +
    labs(
      title = paste0("Significant Hallmarks per Pattern (Unique vs Non-Unique) - ",
                     title_prefix, " ", sample_name, " (", n_patterns, " patterns)"),
      x = "",
      y = expression(-10*log[10](padj)),
      fill = "Hallmark type"
    ) +
    theme_classic() +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold"),
      axis.text.y = element_text(size = 5),
      legend.position = "right"
    )

  return(p_combined)
}

# Read CSV
df_enrich <- safe_read(input_csv_enrichment)
df_overr  <- safe_read(input_csv_overr)

# Generate plots
p_enrich <- generate_combined_plot(df_enrich, "Enrichment")
p_overr  <- generate_combined_plot(df_overr,  "Overrepresentation")


png(output_png_enrich, width = 10, height = 8, units = "in", res = 300)
if(!is.null(p_enrich)) print(p_enrich)
dev.off()

png(output_png_overr, width = 10, height = 8, units = "in", res = 300)
if(!is.null(p_overr)) print(p_overr)
dev.off()