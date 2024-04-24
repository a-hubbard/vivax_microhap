# Bar plots showing percent OT reads for the different protocols

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
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))

protocol_barplot <- function(pct_ot_reads, dna_source) {
  if (dna_source == "DBS") {
    facet_cmd <- facet_grid(rows = vars(Extraction), cols = vars(Treatment))
  } else if (dna_source == "Whole blood") {
    facet_cmd <- facet_wrap(vars(Treatment))
  }
  pct_ot_reads %>%
    filter(Source == dna_source) %>%
    ggplot(mapping = aes(x = SID, y = pct_read_ot)) +
    geom_col() +
    facet_cmd +
    ylim(0, 100) +
    labs(x = "Sample", y = "On-Target Reads (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

# Read in data ---------------------------------------------------------
pct_ot_reads <- read_csv(
    arg$readcounts_metadata, 
    col_types = cols(
      .default = col_character(), 
      n_read_total = col_integer(), 
      n_read_ot = col_integer(), 
      pct_read_ot = col_double(), 
      CT = col_double(), 
      CT_Pv = col_double()
    ), 
    progress = FALSE
  ) %>%
  # Filter out runs with both sWGA and target pre-amplification
  filter(Treatment != "sWGA & Targ. Pre-amp.")

# Make plots and save --------------------------------------------------
dbs_plot <- protocol_barplot(pct_ot_reads, "DBS")
wb_plot <- protocol_barplot(pct_ot_reads, "Whole blood")
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
