# Plot the confidence interval width of rhat estimates

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
    "--rel_sim_res", 
    help = "TSV file with results of paneljudge relatedness simulations"
  ), 
  make_option("--shared_functions", help = "Path to shared functions file"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg <- list(
    rel_sim_res = "../../results/rel_sim_res_all.tsv", 
    shared_functions = "../shared_functions.R", 
    out_base = "../../results/figs/suppfig4"
  )
}

# Set default ggplot2 theme
theme_set(theme_bw())
# Access shared functions
source(arg$shared_functions)

# Read in relatedness simulation results -------------------------------
rel_sim_res <- read_tsv(
    arg$rel_sim_res, 
    col_types = cols(
      .default = col_double(), 
      panel = col_character(), 
      pop = col_character(), 
      pair_num = col_integer(), 
      pop_lbl = col_character(), 
      panel_lbl = col_character()
    ), 
    progress = FALSE
  )

# Plot rhat confidence interval widths ---------------------------------
fig <- rel_sim_res %>%
  rhat_ci_width_boxplot()
w <- 10
h <- 6
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
