suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

# Inputs / outputs
contributions_file <- snakemake@input[["contributions_file"]]
clusters_file <- snakemake@input[["subclone"]]

output_png1 <- snakemake@output[["png1"]]
output_png2 <- snakemake@output[["png2"]]

# -------------------
# Lecture des fichiers
# -------------------
contrib <- read.table(
  contributions_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

clusters <- read.table(
  clusters_file,
  header = FALSE,
  sep = "\t",
  col.names = c("barcode", "subclone"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    barcode = sub("^.*?_", "", barcode),
    subclone = factor(subclone)
  )

# -------------------
# Merge
# -------------------
contrib2 <- contrib %>%
  left_join(clusters, by = "barcode")

# -------------------
# Format long
# -------------------
contrib_long <- contrib2 %>%
  pivot_longer(
    cols = starts_with("Pattern_"),
    names_to = "Pattern",
    values_to = "Weight"
  ) %>%
  mutate(
    subclone = factor(subclone),
    Pattern = factor(Pattern)
  )

# -------------------
# Boxplot 1 : comparaison des patterns
# -------------------
bxp <- ggboxplot(
  contrib_long,
  x = "Pattern",
  y = "Weight",
  color = "Pattern",
  palette = "npg",
  outlier.size = 0.5,
  legend = "none"
) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format"
  ) +
  geom_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "{p.adj.signif}",
    hide.ns = TRUE
  )

# -------------------
# Boxplot 2 : subclones par pattern
# -------------------
bxp2 <- ggboxplot(
  contrib_long,
  x = "subclone",
  y = "Weight",
  color = "subclone",
  palette = "npg",
  facet.by = "Pattern",
  outlier.size = 0.5,
  legend = "none"
) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y.npc = "top"
  ) +
  geom_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "{p.adj.signif}",
    hide.ns = TRUE
  )

# -------------------
# Export PNG
# -------------------
ggsave(
  filename = output_png1,
  plot = bxp,
  width = 14,
  height = 10,
  dpi = 150
)

ggsave(
  filename = output_png2,
  plot = bxp2,
  width = 14,
  height = 10,
  dpi = 150
)
