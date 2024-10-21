# Box plots showing read counts by locus for the different protocols

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
    "--readcounts_metadata", 
    help = "CSV file containing read counts and metadata for May 2022 runs"
  ), 
  make_option(
    "--selected_loci", 
    help = "CSV file identifying loci selected for final panel"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   readcounts_metadata = "../../results/may2022_readcounts_metadata.csv", 
#   selected_loci = "../../results/loci2filter/good_amp.csv"
# )

protocol_boxplot <- function(read_counts, dna_source) {
  if (dna_source == "DBS") {
    mapping <- aes(x = Treatment, y = mean_read_count, color = Extraction)
    color_scale <- scale_color_discrete(name = "Extraction\nMethod")
  } else if (dna_source == "Whole blood") {
    mapping <- aes(x = Treatment, y = mean_read_count)
    color_scale <- NULL
  }
  read_counts %>%
    filter(Source == dna_source) %>%
    ggplot(mapping = mapping) +
    geom_boxplot() +
    scale_y_continuous(trans = "log10") +
    color_scale +
    labs(x = "Enrichment Method", y = "Mean Locus Read Count") +
    theme_bw()
}

# Read in data ---------------------------------------------------------
reads_by_locus <- read_csv(
    arg$readcounts_metadata, 
    col_types = cols(
      .default = col_character(), 
      n_read = col_integer(), 
      CT = col_double(), 
      CT_Pv = col_double()
    ), 
    progress = FALSE
  ) %>%
  # Filter out runs with both SWGA and target pre-amplification
  filter(Treatment != "SWGA & Targ. Pre-amp.")
selected_loci <- read_csv(
  arg$selected_loci, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)

# Filter loci and compute mean read count across samples ---------------
mean_read_counts <- selected_loci %>%
  # Use left_join to filter to selected loci
  left_join(reads_by_locus, by = "locus") %>%
  # Make all missing data explicit
  complete(
    SID, 
    Source, 
    Treatment, 
    Extraction, 
    locus, 
    fill = list(n_read = 0)
  ) %>%
  group_by(locus, Treatment, Extraction, Source) %>%
  summarize(mean_read_count = mean(n_read), .groups = "drop") %>%
  # Add 1 to all reads to enable use of log10 transformation
  mutate(mean_read_count = mean_read_count + 1)

# Make plots and save --------------------------------------------------
dbs_plot <- protocol_boxplot(mean_read_counts, "DBS")
wb_plot <- protocol_boxplot(mean_read_counts, "Whole blood")
fig <- (dbs_plot | wb_plot) +
  plot_annotation(tag_levels = "A")
w <- 9
h <- 5
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
