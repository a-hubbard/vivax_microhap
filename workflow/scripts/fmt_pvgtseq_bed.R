# Reformat PvGTSeq target spreadsheet into BED

# Load required libraries
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(dplyr)
library(magrittr)
library(optparse)
library(readr)
library(readxl)
library(stringr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--target_info", help = "Excel file containing target info"), 
  make_option(
    "--chrom_key", 
    help = str_c(
      "TSV file containing key for matching chromosome names to designations ", 
      "used in reference genome"
    )
  ), 
  make_option("--panel_bed", help = "Path for BED file describing panel")
)
arg <- parse_args(OptionParser(option_list = opts))
if (interactive()) {
  arg <- list(
    target_info = "../../resources/S4_Table.xlsx"
  )
}

# Read in data ---------------------------------------------------------
target_info <- read_excel(arg$target_info, sheet = "S4_Table") %>%
  select(chromosome, start, end, amplicon_name) %>%
  relocate(amplicon_name, .after = end) %>%
  # Remove mitochondria target
  filter(chromosome != "PvP01_MIT_v1")

# Read and join chromosome name key ------------------------------------
chrom_key <- read_tsv(
    arg$chrom_key, 
    col_names = c("chrom_name", "chrom_assembly_name"), 
    col_types = cols(.default = col_character()), 
    progress = FALSE
  )
target_info <- target_info %>%
  left_join(chrom_key, by = c("chromosome" = "chrom_name")) %>%
  select(-chromosome) %>%
  relocate(chrom_assembly_name, .before = 1)

# Add BED required columns and write to disk ---------------------------
target_info %>%
  # The positive strand is recorded for all targets, because the input 
  # spreadsheet does not include strand information and the paneljudge 
  # analyses will be unaffected by strand designation
  mutate(score = ".", strand = "+") %>%
  write_tsv(arg$panel_bed, col_names = FALSE)
