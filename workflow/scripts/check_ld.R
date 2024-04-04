# Check for linkage disequilibrium and save results

# Load required libraries ----------------------------------------------
library(adegenet)
library(poppr)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(tidyverse)

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
data_gd <- read_rds(arg$data_gd)
adegenet::setPop(data_gd) <- ~Population
data_gd %>%
  poppr::clonecorrect(strata = ~Population) %>%
  poppr::pair.ia(sample = 1000, plot = FALSE) %>%
  write_rds(arg$out)
