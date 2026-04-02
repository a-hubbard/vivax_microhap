# Create heatmaps showing read counts for each sample/target

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--read_counts", help = "CSV file containing read counts"), 
  make_option("--shared_functions", help = "Path to shared functions file"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
if (interactive()) {
  arg$read_counts <- "../../results/AmpSeq/serialdil2/rl_read_counts.csv"
  arg$shared_functions <- "../shared_functions.R"
  arg$out_base <- "../../results/figs/suppfig1"
}

# Access shared functions
source(arg$shared_functions)

# Read in data ---------------------------------------------------------
rep_read_counts <- read_csv(
  arg$read_counts, 
  col_types = cols(
    .default = col_character(), 
    parasitemia = col_integer(), 
    n_read = col_double()
  ), 
  progress = FALSE
)

# Summarize for each sample --------------------------------------------
sample_read_counts <- rep_read_counts %>%
  # Keep protocol "v2" over "v1" when both present
  filter(parasitemia != 100 | protocol == "v2") %>%
  group_by(sample_id, specimen_id, parasitemia, locus) %>%
  summarize(n_read = mean(n_read), .groups = "drop")

# Plot and save --------------------------------------------------------
fig <- sample_read_counts %>%
  select(-sample_id) %>%
  mutate(specimen_id = str_remove(specimen_id, "pv")) %>%
  mutate(
    sample_id = str_c("Specimen: ", specimen_id, ", Parasitemia: ", parasitemia)
  ) %>%
  arrange(parasitemia, desc(specimen_id)) %>%
  mutate(
    sample_id = factor(sample_id, levels = unique(sample_id), ordered = TRUE)
  ) %>%
  read_count_heatmap()
w <- 15
h <- 7
ggsave(
  str_c(arg$out_base, ".pdf"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".png"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".eps"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
