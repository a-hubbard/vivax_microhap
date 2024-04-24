# Combine read counts and metadata from May 2022 run into one CSV

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--read_counts", 
    help = "TSV file containing read counts separated by target"
  ), 
  make_option(
    "--total_read_counts", 
    help = "TSV file containing total read counts from each well"
  ), 
  make_option(
    "--metadata", 
    help = "CSV file containing metadata for the May 2022 experiments"
  ), 
  make_option("--out", help = "CSV to contain output")
)
arg <- parse_args(OptionParser(option_list = opts))

# Read in data ---------------------------------------------------------
# Read counts for each target
read_counts <- read_tsv(
  arg$read_counts, 
  col_names = c("target", "protocol_well", "n_read"), 
  col_types = cols(.default = col_character(), n_read = col_integer()), 
  progress = FALSE
)
# Total read counts from each well
total_read_counts <- read_tsv(
  arg$total_read_counts, 
  col_names = c("protocol_well", "n_read_total"), 
  col_types = cols(.default = col_character(), n_read_total = col_integer()), 
  progress = FALSE
)
metadata <- read_csv(
  arg$metadata, 
  col_types = cols(
    .default = col_character(), 
    CT = col_double(), 
    CT_Pv = col_double()
  ), 
  progress = FALSE
)

# Compute percent on-target (OT) reads ---------------------------------
ot_reads <- read_counts %>%
  group_by(protocol_well) %>%
  summarize(n_read_ot = sum(n_read), .groups = "drop")
pct_ot_reads <- total_read_counts %>%
  left_join(ot_reads, by = "protocol_well") %>%
  mutate(pct_read_ot = (n_read_ot / n_read_total) * 100)

# Join metadata and write ----------------------------------------------
pct_ot_reads <- pct_ot_reads %>%
  mutate(is_preamp = str_detect(protocol_well, "preAmped")) %>%
  mutate(well = str_sub(protocol_well, start = -3L)) %>%
  select(-protocol_well) %>%
  separate_wider_position(well, c(well_letter = 1, well_num = 2)) %>%
  mutate(well_num = as.integer(well_num)) %>%
  mutate(well_num = as.character(well_num)) %>%
  unite(well, well_letter, well_num, sep = "") %>%
  left_join(metadata, by = c("well" = "Well")) %>%
  mutate(
    Treatment = if_else(
      Enrichment == "sWGA", 
      if_else(is_preamp, "sWGA & Targ. Pre-amp.", "sWGA"), 
      if_else(is_preamp, "Targ. Pre-amp.", "None")
    )
  ) %>%
  select(-is_preamp, -Enrichment) %>%
  relocate(SID, Source, Extraction, Treatment, .before = 1) %>%
  write_csv(arg$out)
