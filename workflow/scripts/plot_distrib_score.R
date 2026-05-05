suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tibble)
  library(forcats)
  library(ggplot2)
  library(ggridges)
})

contrib_file <- snakemake@input[["contrib"]] 
subclones_files <- snakemake@input[["subclones"]]
output_png   <- snakemake@output[["png"]]
sample_name <- snakemake@wildcards[["sample"]]

contrib <- read.table(contrib_file, header = TRUE)

subclones <- read.table(subclones_files) |>
  rename(
    cells_id = V1,
    subclone = V2
  ) |>
  mutate(barcode = sub(".*_", "", cells_id)) |>
  filter(barcode %in% contrib$barcode) |>
  column_to_rownames("barcode")

df <- contrib |>
  inner_join(
    subclones |> rownames_to_column("barcode"),
    by = "barcode")


df_long <- df |>
  pivot_longer(
    cols = starts_with("Pattern_"),
    names_to = "pattern",
    values_to = "score"
  ) |>
  mutate(
    pattern = factor(pattern, levels = paste0("Pattern_", sort(as.numeric(sub("Pattern_", "", unique(pattern))))))
  )

p <- ggplot(df_long, aes(x = score, y = factor(subclone), fill = factor(subclone))) +
  geom_density_ridges(alpha = 0.6, scale = 1.2) +
  facet_wrap(~ pattern) +   
  labs(
    x = "Weight Score",
    y = "Subclone",
    fill = "Subclone"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line()
  )

png(output_png, width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()