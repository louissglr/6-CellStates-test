suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tibble)
  library(forcats)
})

contrib_file <- snakemake@input[["contrib"]] 
subclones_files <- snakemake@input[["subclones"]]

output_png  <- snakemake@output[["png1"]]       # Plot subclone -> pattern
output_png2 <- snakemake@output[["png2"]]     # Plot pattern -> subclone

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

df_prop <- df_prop |>
  mutate(
    dominant_pattern = fct_reorder(
      dominant_pattern,
      as.numeric(gsub("Pattern_", "", dominant_pattern))
    ),
    subclone = factor(
      subclone,
      levels = sort(unique(subclone))
    )
  )

df_counts <- df |>
  count(subclone, name = "n_cells") |>
  mutate(
    subclone = factor(subclone, levels = levels(df_prop$subclone))
  )

p1 <- ggplot(df_prop, aes(x = subclone, y = prop, fill = dominant_pattern)) +
  geom_col(width = 0.8) +
  coord_flip() +
  
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1),
    expand = c(0, 0)
  ) +
  
  geom_text(
    data = df_counts,
    aes(x = subclone, y = 1.02, label = paste0("n = ", n_cells)),
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

## Save plot 1
png(output_png, width = 10, height = 8, units = "in", res = 300)
print(p1)
dev.off()


df_prop_pattern <- df |>
  group_by(dominant_pattern, subclone) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(dominant_pattern) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

## Ordres cohérents
df_prop_pattern <- df_prop_pattern |>
  mutate(
    dominant_pattern = fct_reorder(
      dominant_pattern,
      as.numeric(gsub("Pattern_", "", dominant_pattern))
    ),
    subclone = factor(
      subclone,
      levels = levels(df_prop$subclone)
    )
  )

## Annotation n cellules / pattern
df_counts_pattern <- df |>
  count(dominant_pattern, name = "n_cells") |>
  mutate(
    dominant_pattern = fct_reorder(
      dominant_pattern,
      as.numeric(gsub("Pattern_", "", dominant_pattern))
    )
  )

p2 <- ggplot(df_prop_pattern, aes(x = dominant_pattern, y = prop, fill = subclone)) +
  geom_col(width = 0.8) +
  coord_flip() +
  
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1),
    expand = c(0, 0)
  ) +
  
  geom_text(
    data = df_counts_pattern,
    aes(x = dominant_pattern, y = 1.02, label = paste0("n = ", n_cells)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.5
  ) +
  
  labs(
    x = "Dominant pattern",
    y = "Proportion of cells",
    fill = "Subclone"
  ) +
  
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

## Save plot 2
png(output_png2, width = 10, height = 8, units = "in", res = 300)
print(p2)
dev.off()