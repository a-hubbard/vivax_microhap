# Convert microhaplotype data from CSV to genind format

# Load required libraries
library(hubpopgen)
library(pegas)
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
  select(target, sample_id, hap_id, geog) %>%
  # Remove samples without defined geographies, which means I decided 
  # the sample size was insufficient for population-level analysis
  filter(! is.na(geog))

# Convert to genind object and save ------------------------------------
mh_data %>%
  # Required by adegenet
  mutate(target = str_replace_all(target, "\\.", "_")) %>%
  rename(
    sample_ID = sample_id, 
    locus = target, 
    allele = hap_id, 
    pop = geog
  ) %>%
  hubpopgen::tib2genind() %>%
  write_rds(arg$mh_gd)
