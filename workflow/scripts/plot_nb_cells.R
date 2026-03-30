suppressPackageStartupMessages({
library(tidyverse)
library(dplyr)
library(tibble)
library(forcats)
})

contrib_file <- snakemake@input[["contrib"]] 
subclones_files <- snakemake@input[["subclones"]]
output_pdf  <- snakemake@output[["pdf"]]
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
    by = "barcode"
  )

df_prop <- df |>
  group_by(subclone, dominant_pattern) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(subclone) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

## Ordre pattern
df_prop <- df_prop |>
  mutate(
    dominant_pattern = fct_reorder(
      dominant_pattern,
      as.numeric(gsub("Pattern_", "", dominant_pattern))
    )
  )

## Ordre subclone
df_prop <- df_prop |>
  mutate(
    subclone = factor(
      subclone,
      levels = sort(unique(subclone))
    )
  )



df_counts <- df |>
  count(subclone, name = "n_cells") |>
  mutate(
    subclone = factor(
      subclone,
      levels = levels(df_prop$subclone)
    )
  )

## --------------------------------------------------
## Plot
## --------------------------------------------------

p <- ggplot(df_prop, aes(x = subclone, y = prop, fill = dominant_pattern)) +
  geom_col(width = 0.8) +
  coord_flip() +
  
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1),
    expand = c(0, 0)
  ) +
  
  scale_x_discrete(
    breaks = levels(df_prop$subclone)
  ) +
  
  geom_text(
    data = df_counts,
    aes(
      x = subclone,
      y = 1.02,
      label = paste0("n = ", n_cells)
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.5
  ) +
  
  labs(
    x = "Subclone",
    y = "Proportion of cells",
    fill = "Dominant pattern"
  ) +
  
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )


pdf(output_pdf, width = 10, height = 8)

print(p)
dev.off()
