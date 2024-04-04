# Reformat primer list and filter out targets with known issues

# Load required libraries
library(seqinr)
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(optparse)
library(readxl)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--sd_id_file", help = "SeekDeep ID file"), 
  make_option("--primer_info", help = "Excel file containing primer info"), 
  make_option("--primers_tsv", help = "TSV file to contain primers"), 
  make_option(
    "--fwd_primers", 
    help = "Path of FASTA file to contain forward primers"
  ), 
  make_option(
    "--rev_primers", 
    help = "Path of FASTA file to contain reverse primers"
  )
)
arg <- parse_args(OptionParser(option_list = opts))

# Read in data ---------------------------------------------------------
primers_all <- read_tsv(
  arg$sd_id_file, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)
# This Excel also has the primer sequences, so it would be possible to 
# get all of the necessary information from this file. However, the 
# adapter sequences have already been removed from the SeekDeep ID 
# file, making it quicker to just take the primer sequences from there.
targets_include_exclude <- read_excel(
  arg$primer_info, 
  sheet = "Pv_Include_Exclude"
)

# Filter target list ---------------------------------------------------
primers_filtered <- primers_all %>%
  filter(! (target %in% targets_include_exclude$Exclude)) %>%
  # We decided to use the OLD samples, which were not sequenced for the 
  # SNP markers, so these are removed
  filter(! str_detect(target, "SNP"))

# Write to disk --------------------------------------------------------
# TSV file
primers_filtered %>%
  write_tsv(arg$primers_tsv)
# FASTA files
seqinr::write.fasta(
  as.list(primers_filtered$forward), 
  primers_filtered$target, 
  arg$fwd_primers
)
seqinr::write.fasta(
  as.list(primers_filtered$reverse), 
  primers_filtered$target, 
  arg$rev_primers
)
