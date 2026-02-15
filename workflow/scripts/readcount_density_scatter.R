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
  make_option("--shared_functions", help = "Path to shared functions file"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg <- list(
    rl_pct_good_amp_parasitemia = 
      "../../results/AmpSeq/serialdil2/rl_pct_good_amp_parasitemia.csv", 
    shared_functions = "../shared_functions.R"
  )
}

# Set default ggplot2 theme
theme_set(theme_bw())
# Access shared functions
source(arg$shared_functions)

# Read in read counts and parasitemia values ---------------------------
rl_pct_good_amp_parasitemia <- read_csv(
    arg$rl_pct_good_amp_parasitemia, 
    col_types = cols(
      .default = col_character(), 
      parasitemia = col_integer(), 
      pct_loci_good_amp = col_double()
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
