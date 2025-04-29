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
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
if (interactive()) {
  arg$read_counts <- "../../results//AmpSeq/uci1223/sl_read_counts.csv"
}

read_count_heatmap_bysample <- function(reads) {
  reads %>%
    ggplot(mapping = aes(x = locus, y = sample_id, fill = n_read)) +
    geom_tile() +
    scale_fill_fermenter(
      palette = "YlGnBu", 
      name = "Reads", 
      breaks = c(10, 100, 1000, 2000)
    ) +
    labs(x = "Locus", y = "Sample") +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    )
}

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
fig <- read_count_heatmap_bysample(read_counts)
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
