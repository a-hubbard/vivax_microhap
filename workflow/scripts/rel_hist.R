# Create histograms showing relatedness values

# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(patchwork)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--rel_csv", help = "CSV file containing relatedness values"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(rel_csv = "../../results/relatedness/sample_rel.csv")

rel_hist <- function(rel_res) {
  ggplot(data = rel_res, mapping = aes(x = estimate)) +
    geom_histogram(bins = 50) +
    labs(x = "Relatedness", y = "No. of Sample Pairs") +
    theme_bw()
}

# Read data ------------------------------------------------------------
rel_res <- read_csv(
  arg$rel_csv, 
  col_types = cols(
    .default = col_double(), 
    sample_a = col_character(), 
    sample_b = col_character()
  ), 
  progress = FALSE
)

# Plot and save --------------------------------------------------------
full_hist <- rel_hist(rel_res)
nozero_hist <- filter(rel_res, estimate != 0) %>%
  rel_hist()
# Note the outer parentheses are necessary to make the tagging work
fig <- (full_hist | nozero_hist) +
  plot_annotation(tag_levels = "A")
w <- 6
h <- 3
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
