suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
})

contributions_file <- snakemake@input[["contributions_file"]]
clusters_file <- snakemake@input[["subclone"]]
kw_output_csv <- snakemake@output[["kw_csv"]]
pairwise_output_csv <- snakemake@output[["pairwise_csv"]]

# Metadata 
sample_name <- snakemake@wildcards[["sample"]]
n_patterns <- as.integer(snakemake@wildcards[["npatterns"]])

contrib <- read.table(
  contributions_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) |>
  mutate(barcode = sub("^.*#", "", barcode))

# ---- Read clusters ----
clusters <- read.table(
  clusters_file,
  header = FALSE,
  sep = "\t",
  col.names = c("barcode", "subclone"),
  stringsAsFactors = FALSE
) |>
  mutate(
    barcode = sub("^.*?_", "", barcode),
    subclone = factor(subclone)
  )

# ---- Merge + long format ----
contrib_long <- contrib |>
  left_join(clusters, by = "barcode") |>
  pivot_longer(
    cols = starts_with("Pattern_"),
    names_to = "Pattern",
    values_to = "Weight"
  ) |>
  drop_na(subclone)

# ---- Kruskal-Wallis + correction BH ----
kw_results <- contrib_long |>
  group_by(Pattern) |>
  kruskal_test(Weight ~ subclone) |>
  adjust_pvalue(method = "BH") |>
  add_significance("p.adj") |>
  mutate(
    sample = sample_name,
    n_patterns = n_patterns
  )

write.csv(kw_results, kw_output_csv, row.names = FALSE)
message("Kruskal-Wallis results written to: ", kw_output_csv)

# Pairwise Wilcoxon post-hoc (per Pattern) 
pairwise_results <- contrib_long %>%
  group_by(Pattern) %>%
  summarise(
    pw = list(pairwise.wilcox.test(
      x = Weight,
      g = subclone,
      p.adjust.method = "BH"
    )),
    .groups = "drop"
  ) %>%
  mutate(
    pw_df = map(pw, \(x) {
      pval_mat <- x$p.value
      df <- as.data.frame(as.table(pval_mat))
      colnames(df) <- c("group1", "group2", "p.adj")
      df %>%
        drop_na(p.adj) %>%
        mutate(
          significance = case_when(
            p.adj < 0.001 ~ "***",
            p.adj < 0.01  ~ "**",
            p.adj < 0.05  ~ "*",
            TRUE          ~ "ns"
          ),
          sample = sample_name,
          n_patterns = n_patterns
        )
    })
  ) %>%
  select(Pattern, pw_df) %>%
  unnest(pw_df)

write.csv(pairwise_results, pairwise_output_csv, row.names = FALSE)
message("Pairwise Wilcoxon results written to: ", pairwise_output_csv)
