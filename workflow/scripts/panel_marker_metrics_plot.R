# Bar plots showing marker metrics from panel design

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(patchwork)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--window_stats", 
    help = "TSV file containing stats for each window from the panel design"
  ), 
  make_option(
    "--targets", 
    help = "CSV file containing targets selected for final panel"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   window_stats = "../../resources/filtered_windows_tab.txt", 
#   targets = "../../results/trgs2filter/good_amp.csv"
# )

metric_barplot <- function(stats, metric, label_y_axis = TRUE) {
  theme_tweaks <- NULL
  if (! label_y_axis) {
    theme_tweaks <- theme(
      axis.text.y = element_blank(), 
      axis.title.y = element_blank()
    )
  }
  stat_lbl_key <- c(
    "mean_TD" = "Mean Tajima's D", 
    "mean_ND" = "Mean Nucleotide Diversity", 
    "mean_FST" = "Mean FST"
  )
  stats %>%
    ggplot(mapping = aes(x = .data[[metric]], y = target)) +
    geom_col() +
    labs(x = stat_lbl_key[metric], y = "Target") +
    theme_bw() +
    theme_tweaks
}

# Read in data ---------------------------------------------------------
window_stats <- read_tsv(
    arg$window_stats, 
    col_types = cols(
      .default = col_character(), 
      mean_TD = col_double(), 
      count_TD = col_integer(), 
      mean_ND = col_double(), 
      count_ND = col_integer(), 
      mean_FST = col_double(), 
      count_FST = col_integer()
    ), 
    progress = FALSE
  ) %>%
  # This column is empty
  select(-`...10`)
targets <- read_csv(
  arg$targets, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)

# Filter window stats to selected targets ------------------------------
window_stats <- window_stats %>%
  # Reconstruct target names in format used in the selected targets file
  mutate(version = "v1") %>%
  unite(target, chromosome, version, wind.start, wind.end)
target_stats <- targets %>%
  left_join(window_stats, by = "target") %>%
  # Remove targets that were added after sliding window analysis (e.g., 
  # purported drug resistance markers) and therefore don't have stats
  filter(! is.na(mean_TD))

# Create and save barplots
fig <- (
    metric_barplot(target_stats, "mean_ND") | 
    metric_barplot(target_stats, "mean_FST", label_y_axis = FALSE)
  ) +
  plot_annotation(tag_levels = "A")
w <- 8
h <- 9
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
