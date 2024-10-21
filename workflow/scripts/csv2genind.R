# Convert microhaplotype data from CSV to genind format

# Load required libraries
library(adegenet)
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(tidyverse)

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
# Arguments used for development
# arg <- list(
#   mh_csv = "../../results/microhap/MalariaGEN/mg_microhap_w_hap_ids.csv"
# )

# Read in microhaplotype data ------------------------------------------
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
  select(locus, sample_id, hap_id, Population, Country, Site)

# Convert to genind object and save ------------------------------------
mh_data_wide <- mh_data %>%
  # Required by adegenet
  mutate(locus = str_replace_all(locus, "\\.", "_")) %>%
  mutate(Population = as.factor(Population)) %>%
  pivot_wider(values_from = "hap_id", names_from = "locus") %>%
  column_to_rownames("sample_id")
mh_data_wide %>%
  select(-Population, -Country, -Site) %>%
  # Note that while it is not mentioned in the documentation, it 
  # would appear based on empirical testing that this does not alter 
  # the order of the samples
  adegenet::df2genind(
    ploidy = 1, 
    strata = select(mh_data_wide, Population, Country, Site)
  ) %>%
  write_rds(arg$mh_gd)
