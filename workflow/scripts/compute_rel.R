# Estimate IBD-based relatedness with Dcifer

# Load required libraries ----------------------------------------------
library(dcifer)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--gen_data", help = "CSV containing genetic data"), 
  make_option("--seed", type = "integer", help = "Random number seed"), 
  make_option(
    "--rel_mat", 
    help = "Path of RDS file to contain relatedness matrix"
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   gen_data = 
#     "../../results/microhap/MalariaGEN/mg_microhap_filtered4popgen.csv", 
#   seed = 1
# )

set.seed(arg$seed)

# Read data ------------------------------------------------------------
gen_data <- dcifer::readDat(arg$gen_data, "sample_id", "locus", "hap_id")

# Compute relatedness and save -----------------------------------------
moi <- rep_len(1, length(gen_data))
a_freq <- dcifer::calcAfreq(gen_data, moi)
dcifer::ibdDat(gen_data, moi, a_freq, confint = TRUE) %>%
  write_rds(arg$rel_mat)
