# Create plots of cardinality versus effective cardinality

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
    "--marker_metrics", 
    help = 
      str_c(
        "TSV file contained cardinality and effective cardinality for ", 
        "multiple panels and populations"
      )
  ), 
  make_option("--shared_functions", help = "Path to shared functions file"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg <- list(
    marker_metrics = "../../results/paneljudge_marker_metrics.tsv", 
    shared_functions = "../shared_functions.R", 
    out_base = "../../results/figs/fig3"
  )
}

# Set default ggplot2 theme
theme_set(theme_bw())
# Access shared functions
source(arg$shared_functions)

# Read in marker metrics -----------------------------------------------
marker_metrics <- read_tsv(
    arg$marker_metrics, 
    col_types = cols(
      .default = col_character(), 
      cardinality = col_integer(), 
      eff_cardinality = col_double()
    ), 
    progress = FALSE
  )

# Plot parasitemia versus sample read count ----------------------------
fig <- marker_metrics %>%
  card_vs_effcard()
w <- 8
h <- 7.5
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
