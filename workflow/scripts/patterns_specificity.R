suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggpubr)
})

contributions_file <- snakemake@input[["contributions_file"]]
clusters_file <- snakemake@input[["subclone"]]
output_pdf <- snakemake@output[["pdf"]]
sample_name <- snakemake@wildcards[["sample"]]
n_patterns <- as.integer(snakemake@wildcards[["npatterns"]])

# Lecture des fichiers
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

# Merge des contributions et clusters
contrib2 <- contrib %>%
  left_join(clusters, by = "barcode")

# Conversion en format long + correction de l'ordre des patterns
contrib_long <- contrib2 %>%
  pivot_longer(
    cols = starts_with("Pattern_"),
    names_to = "Pattern",
    values_to = "Weight"
  ) %>%
  mutate(
    subclone = factor(subclone),
    Pattern = factor(
      Pattern,
      levels = paste0("Pattern_", 1:length(unique(Pattern)))
    )
  )

# -------------------
# Version 1 : subclone vs Weight
# -------------------
v1 <- ggplot(contrib_long, aes(x = subclone, y = Weight, fill = Pattern)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Subclone",
    y = "Weight",
    title = paste(
      "Pattern distribution per subclone\nSample:",
      sample_name,
      "- n_patterns:",
      n_patterns
    )
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------
# Version 2 : Pattern vs Weight
# -------------------
v2 <- ggplot(contrib_long, aes(x = Pattern, y = Weight, fill = subclone)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Pattern",
    y = "Weight",
    title = paste(
      "Subclone distribution per pattern\nSample:",
      sample_name,
      "- n_patterns:",
      n_patterns
    )
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# -------------------
# Version 3 : distribution des poids par pattern
# -------------------
v3 <- ggplot(contrib_long, aes(x = Pattern, y = Weight, fill = Pattern)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = "Pattern",
    y = "Weight",
    title = paste(
      "Weight Distribution per Pattern\nSample:",
      sample_name,
      "- n_patterns:",
      n_patterns
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# -------------------
# Plots ggpubr avec tests statistiques
# -------------------
contrib_long$Pattern <- as.factor(contrib_long$Pattern)

# Comparaison des patterns
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

# Comparaison des subclones par pattern
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
# Sauvegarde des plots dans l'ordre souhaité
# -------------------
pdf(output_pdf, width = 14, height = 10)

print(v3)   # anciennement v3 -> en premier
print(v2)
print(v1)   # anciennement v1 -> en troisième
print(bxp)
print(bxp2)

dev.off()

message("End: ", output_pdf)
