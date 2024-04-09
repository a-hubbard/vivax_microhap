# Convert aligned sequences to DNAbin objects and subset by population

# Load required libraries
library(ape)
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--alignments", 
    help = "RDS file containing tibble with sequence alignments"
  ), 
  make_option("--mh_csv", help = "CSV file containing microhaplotype data"), 
  make_option(
    "--out", 
    help = "RDS file to contain tibble with population-level DNAbin objects"
  )
)
arg <- parse_args(OptionParser(option_list = opts))

subset_by_pop <- function(dnabin, pop_id, sample_pop_key) {
  samples <- sample_pop_key %>%
    filter(Population == pop_id) %$%
    sample_id
  dnabin[samples, ]
}

# Read in data ---------------------------------------------------------
mh_align <- read_rds(arg$alignments)
mh_meta <- read_csv(
  arg$mh_csv, 
  col_types = cols(
    .default = col_character(), 
    Lat = col_double(), 
    Long = col_double(), 
    Year = col_integer(), 
    `% callable` = col_double(), 
    Fws = col_double(), 
    F_MISS = col_double()
  ), 
  progress = FALSE
)

# Make sample/population key -------------------------------------------
sample_pop_key <- mh_meta %>%
  select(sample_id, Population) %>%
  distinct() %>%
  group_by(Population) %>%
  mutate(n_samp = n()) %>%
  ungroup() %>%
  filter(n_samp >= 15) %>%
  select(-n_samp)

# Convert to DNAbin, subset by population, and write -------------------
mh_align %>%
  mutate(pop_id = list(unique(sample_pop_key$Population))) %>%
  unnest(pop_id) %>%
  mutate(seqs_dnabin = map(trg_seqs_aligned, ape::as.DNAbin)) %>%
  select(-trg_seqs_aligned) %>%
  mutate(
    seqs_dnabin = map2(seqs_dnabin, pop_id, subset_by_pop, sample_pop_key)
  ) %>%
  write_rds(arg$out)
