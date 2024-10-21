# Get haplotypes from aligned sequences, saving haplotye IDs to tidy CSV

# Load required libraries
library(ape)
library(pegas)
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
    "--haplotypes", 
    help = "RDS file containing tibble with haplotype objects"
  ), 
  make_option(
    "--mh_csv_w_hap_ids", 
    help = "CSV file to contain microhaplotype data, with haplotype IDs added"
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   alignments = "../../results/microhap/MalariaGEN/mg_microhap_aligned.rds"
# )

alignment2hap <- function(msa_alignment) {
  quiet_hap <- quietly(pegas::haplotype)
  msa_alignment %>%
    ape::as.DNAbin() %>%
    quiet_hap()
}

match_samples2haps <- function(alignments, haplotypes) {
  # Get sample names from msa alignment object
  sample_ids <- alignments %>%
    Biostrings::unmasked() %>%
    names()
  # Get list of sequence indices corresponding to each haplotype
  indices <- attr(haplotypes, "index")
  # Get list of sample IDs corresponding to each haplotype
  sample_ids_list <- vector("list", length(indices))
  for (i in seq_along(indices)) {
    sample_ids_list[[i]] <- sample_ids[indices[[i]]]
  }
  # Create tibble that matches haplotype and sample IDs
  tibble(
      hap_id = dimnames(haplotypes)[[1]], 
      sample_id = sample_ids_list
    ) %>%
    unnest(sample_id)
}

# Read in aligned sequences and convert to haplotype objects -----------
mh_aligned <- read_rds(arg$alignments)
mh_hap <- mh_aligned %>%
  mutate(hap_out = map(locus_seqs_aligned, alignment2hap)) %>%
  mutate(hap = map(hap_out, pluck, "result")) %>%
  mutate(hap_warn = map(hap_out, pluck, "warnings")) %>%
  select(-hap_out) 

# Print warnings obtained during conversion ----------------------------
# The expected set of warnings is:
# 1) "some sequences of different lengths were assigned to the same 
# haplotype;" 
# 2) "some sequences were not assigned to the same haplotype because of 
# ambiguities;" and 
# 3) "no segregating site detected with these options."
# 1 is the expected outcome of using trailingGapsAsN = TRUE, 2 arises 
# because sometimes a sequence with ambiguities cannot be definitively 
# assigned to a single non-ambiguous haplotype (i.e., it is compatible 
# with two or more), and 3 arises from the underlying data and is not 
# an issue with the analysis.
print("The warnings obtained during the haplotype conversion are:")
mh_hap %>%
  select(locus, hap_warn) %>%
  unnest(hap_warn) %$%
  unique(hap_warn)

# Save haplotype tibble to disk ----------------------------------------
mh_hap %>%
  select(-locus_seqs_aligned, -hap_warn) %>%
  write_rds(arg$haplotypes)

# Match sample and haplotype IDs for all loci --------------------------
sample_hap_key <- mh_hap %>%
  select(-chrom, -hap_warn) %>%
  mutate(sample_hap_key = map2(locus_seqs_aligned, hap, match_samples2haps)) %>%
  select(-locus_seqs_aligned, -hap) %>%
  unnest(sample_hap_key)

# Read microhaplotype CSV, join haplotype IDs, and write back to disk --
mh_data <- read_csv(
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
  ) %>%
  left_join(sample_hap_key, by = c("locus", "sample_id")) %>%
  write_csv(arg$mh_csv_w_hap_ids)
