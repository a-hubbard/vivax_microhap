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
}

# Access shared functions
source(arg$shared_functions)

# Read in data ---------------------------------------------------------
read_counts <- read_csv(
  arg$read_counts, 
  col_types = cols(
    .default = col_character(), 
    n_read = col_double()
  ), 
  progress = FALSE
)

# Plot and save --------------------------------------------------------
fig <- read_count_heatmap(read_counts)
w <- 15
h <- 11
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
