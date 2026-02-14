# Create scatterplot of parasite density versus read count

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(dplyr)
library(ggplot2)
library(magrittr)
library(optparse)                                                                
library(readr)
library(stringr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--rl_pct_good_amp_parasitemia", 
    help = 
      str_c(
        "CSV file containing percent of loci successfully amplified and ", 
        "parasitemia values"
      )
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg <- list(
    s_read_counts_parasitemia = 
      "../../results/AmpSeq/serialdil2/rl_pct_good_amp_parasitemia.csv"
  )
}

# Set default ggplot2 theme
theme_set(theme_bw())
# Access shared functions
source("../shared_functions.R")

# Read in read counts and parasitemia values ---------------------------
rl_pct_good_amp_parasitemia <- read_csv(
    arg$s_read_counts_parasitemia, 
    col_types = cols(
      .default = col_character(), 
      n_read = col_double(), 
      parasitemia = col_double()
    ), 
    progress = FALSE
  )

# Plot parasitemia versus sample read count ----------------------------
fig <- rl_pct_good_amp_parasitemia %>%
  pct_good_amp_vs_parasitemia()
w <- 5
h <- 4.5
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
