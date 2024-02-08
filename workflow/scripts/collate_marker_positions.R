# Create GFF file specifying marker positions

# Load required libraries ----------------------------------------------
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(dplyr)
library(magrittr)
library(optparse)
library(readr)
library(stringr)

# Parse arguments ------------------------------------------------------
# opts <- list(
#   make_option("--sample_names", help = "TSV file containing sample IDs"), 
#   make_option("--out", help = "Path of text file to contain output")
# )
# arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
arg <- list(
  sam_info = "../../resources/saminfo.txt", 
  snp_info = "../../resources/Pviv42_New-Amplicons.txt"
)

sam_info <- read_csv(
    arg$sam_info, 
    col_names = c("locus", "chrom", "start_pos", "end_pos", "length?")
  ) %>%
  select(-`length?`)

snp_info <- read_tsv(arg$snp_info) %>%
  select(LocusID, CHROM, FWD_PRIMER_POS, REV_PRIMER_POS) %>%
  rename(
    locus = LocusID, 
    chrom = CHROM, 
    start_pos = FWD_PRIMER_POS, 
    end_pos = REV_PRIMER_POS
  ) %>%
  mutate(locus = str_c("PvSNP_", locus))

seekdeep_id <- read_tsv("~/Alfred-Liz-PV-microhaplotype-files/vera-2023-5-19/supp_files/id_file.txt")
n_distinct(seekdeep_id$target)

vera_readcount <- read_tsv(
  "../../resources/vera_readcount.txt", 
  col_names = c("locus", "sample_id", "readcount")
)
n_distinct(vera_readcount$locus)
n_distinct(vera_readcount$sample_id)

marker_info <- bind_rows(sam_info, snp_info)
unique(marker_info$locus)

setdiff(seekdeep_id$target, vera_readcount$locus)
setdiff(vera_readcount$locus, seekdeep_id$target)

setdiff(vera_readcount$locus, marker_info$locus)
setdiff(marker_info$locus, vera_readcount$locus)
