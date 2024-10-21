# Filter microhaplotype data for selection

# Load required libraries ----------------------------------------------
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--mh_in", 
    help = "CSV file containing unfiltered microhaplotype data"
  ), 
  make_option(
    "--loci2filter", 
    help = "CSV file containing list of loci to filter out"
  ), 
  make_option(
    "--mh_out", 
    help = "CSV file to contain filtered microhaplotype data"
  )
)
arg <- parse_args(OptionParser(option_list = opts))

# Read in data ---------------------------------------------------------
mh_data <- read_csv(
  arg$mh_in, 
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
loci2filter <- read_csv(
    arg$loci2filter, 
    col_types = cols(.default = col_character()), 
    progress = FALSE
  ) %$%
  locus

# Filter data and write to disk ----------------------------------------
mh_data %>%
  filter(! (locus %in% loci2filter)) %>%
  write_csv(arg$mh_out)
