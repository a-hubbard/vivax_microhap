# Convert microhaplotype data from CSV to genind format

# Load required libraries
library(hubpopgen)
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(dplyr)
library(magrittr)
library(optparse)
library(readr)
library(stringr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--mh_csv", 
    help = "CSV containing microhaplotype data in tidy format"
  ), 
  make_option(
    "--mh_gd", 
    help = "RDS file to contain genind of microhaplotype data"
  )
)
arg <- parse_args(OptionParser(option_list = opts))

# Read in microhaplotype data ------------------------------------------
mh_data <- read_csv(
    arg$mh_csv, 
    col_types = cols(
      .default = col_character(), 
      n_read = col_integer(), 
      n_samp_wdata = col_integer(), 
      n_trg_wdata = col_integer(), 
      ParasiteDensity = col_double(), 
      CT_Pv = col_double(), 
      CT_Pf18s = col_double(), 
      CT_PfvarATS = col_double(), 
      Age = col_integer(), 
      n_hap = col_integer(), 
      moi = col_integer()
    ), 
    progress = FALSE
  ) %>%
  select(target, sample_id, ASV, pop)

# Convert to genind object and save ------------------------------------
mh_data %>%
  filter(! is.na(pop)) %>%
  # Required by adegenet
  mutate(target = str_replace_all(target, "\\.", "_")) %>%
  rename(sample_ID = sample_id, locus = target, allele = ASV) %>%
  hubpopgen::tib2genind() %>%
  write_rds(arg$mh_gd)
