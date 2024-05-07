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
metadata <- read_csv(
  arg$metadata, 
  col_types = cols(
    .default = col_character(), 
    CT = col_double(), 
    CT_Pv = col_double()
  ), 
  progress = FALSE
)

# Join metadata and write ----------------------------------------------
read_counts <- read_counts %>%
  mutate(is_preamp = str_detect(protocol_well, "preAmped")) %>%
  mutate(well = str_sub(protocol_well, start = -3L)) %>%
  separate_wider_position(well, c(well_letter = 1, well_num = 2)) %>%
  mutate(well_num = as.integer(well_num)) %>%
  mutate(well_num = as.character(well_num)) %>%
  unite(well, well_letter, well_num, sep = "") %>%
  left_join(metadata, by = c("well" = "Well")) %>%
  mutate(
    Treatment = if_else(
      Enrichment == "sWGA", 
      if_else(is_preamp, "SWGA & Targ. Pre-amp.", "SWGA"), 
      if_else(is_preamp, "Targ. Pre-amp.", "None")
    )
  ) %>%
  select(-is_preamp, -Enrichment) %>%
  relocate(SID, Source, Extraction, Treatment, .before = 1) %>%
  write_csv(arg$out)
