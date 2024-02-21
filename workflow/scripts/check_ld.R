# Check for linkage disequilibrium over all loci and save results

# Load required libraries ----------------------------------------------
library(poppr)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(readr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--data_gd", 
    help = "RDS file containing genetic data in genind format"
  ), 
  make_option("--seed", type = "integer", help = "Random number seed"), 
  make_option("--out", help = "Path of RDS file to contain output")
)
arg <- parse_args(OptionParser(option_list = opts))

set.seed(arg$seed)

# Read data, check for LD, and write results ---------------------------
read_rds(arg$data_gd) %>%
  poppr::clonecorrect(strata = NA) %>%
  poppr::ia(sample = 999, plot = FALSE) %>%
  write_rds(arg$out)
