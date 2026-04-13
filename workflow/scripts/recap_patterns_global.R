library(CoGAPS)
library(dplyr)
library(stringr)


cogaps_files <- snakemake@input[["cogaps_files"]]
output_csv <- snakemake@output[["csv"]]

# Loop over all files and extract info

recap_list <- lapply(cogaps_files, function(file) {
  
  fname <- basename(file)
  parts <- str_match(fname, "(.*)\\.npatterns-(\\d+)\\.cogaps-object\\.Rds")
  sample_name <- parts[2]
  expected_patterns <- as.integer(parts[3])

  cogaps <- readRDS(file)
  obtained_patterns <- ncol(cogaps@sampleFactors)
  
  mean_chisq <- getMeanChiSq(cogaps)

  data.frame(
    sample = sample_name,
    expected_patterns = expected_patterns,
    obtained_patterns = obtained_patterns,
    meanChisqValue = mean_chisq
  )
})

# Combine all runs
recap_table <- bind_rows(recap_list)
write.csv(recap_table, output_csv, row.names = FALSE)
