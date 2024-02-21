# Identify targets in LD with other targets

# Load required libraries ----------------------------------------------
library(hubpopgen)
library(poppr)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(dplyr)
library(magrittr)
library(optparse)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--gen_data", 
    help = "CSV file containing genetic data in tidy format"
  ), 
  make_option(
    "--pair_ld_thres", 
    type = "double", 
    help = "rbarD threshold to identify pairs of targets that are in LD"
  ), 
  make_option(
    "--signif_thres", 
    type = "double", 
    help = "p-value threshold to consider overall LD as significant"
  ), 
  make_option("--seed", type = "integer", help = "Random number seed"), 
  make_option(
    "--high_ld_trgs", 
    help = "Path of CSV file to contain list of targets"
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   gen_data = "../../results/microhap/microhap.csv", 
#   pair_ld_thres = 0.6, 
#   signif_thres = 0.01, 
#   seed = 1
# )

set.seed(arg$seed)

# Read genetic data ----------------------------------------------------
mh_data <- read_csv(
  arg$gen_data, 
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
)

# Save chromosome info for later ---------------------------------------
chrom_info <- mh_data %>%
  select(target, chrom) %>%
  distinct()

# Prepare microhaplotype data for conversion to genind -----------------
mh_data <- mh_data %>%
  select(target, sample_id, ASV, pop) %>%
  filter(! is.na(pop)) %>%
  # Required by adegenet
  mutate(target = str_replace_all(target, "\\.", "_")) %>%
  rename(sample_ID = sample_id, locus = target, allele = ASV)

# Estimate initial LD values -------------------------------------------
# Compute pairwise LD
pairwise_ld_res <- mh_data %>%
  hubpopgen::tib2genind() %>%
  poppr::clonecorrect(strata = NA) %>%
  poppr::pair.ia(plot = FALSE)
# Find target pairs that are in high LD
high_ld_trg_pairs <- as_tibble(pairwise_ld_res, rownames = "target_pair") %>%
  separate(target_pair, c("target_a", "target_b"), sep = ":") %>%
  # Join chromosome information
  left_join(chrom_info, by = c("target_a" = "target")) %>%
  rename(chrom_a = chrom) %>%
  left_join(chrom_info, by = c("target_b" = "target")) %>%
  rename(chrom_b = chrom) %>%
  # Regardless of what the test shows, we know that targets on 
  # different chromosomes cannot be physically linked
  filter(chrom_a == chrom_b) %>%
  select(-chrom_a, -chrom_b) %>%
  mutate(rbarD = as.double(rbarD)) %>%
  filter(rbarD >= arg$pair_ld_thres)
# Compute overall LD for first check
overall_ld_res <- mh_data %>%
  hubpopgen::tib2genind() %>%
  poppr::clonecorrect(strata = NA) %>%
  poppr::ia(sample = 999, plot = FALSE)

# Iteratively remove targets until LD becomes non-significant ----------
trgs2filter <- character()
while (overall_ld_res[["p.rD"]] < arg$signif_thres) {

  # Filter out any targets already identified for removal
  high_ld_trg_pairs_filtered <- high_ld_trg_pairs %>%
    filter(! (target_a %in% trgs2filter) & ! (target_b %in% trgs2filter))
  # Compute occurrence of each target in high LD pairs by counting the 
  # appearances of each target in each "side" of the pairs and then 
  # summing the two numbers. This approach is valid because all pairs 
  # should be unique.
  high_ld_trg_a <- high_ld_trg_pairs_filtered %>%
    select(-target_b) %>%
    rename(target = target_a) %>%
    group_by(target) %>%
    summarize(trg_prev = n(), .groups = "drop")
  high_ld_trg_b <- high_ld_trg_pairs_filtered %>%
    select(-target_a) %>%
    rename(target = target_b) %>%
    group_by(target) %>%
    summarize(trg_prev = n(), .groups = "drop")
  high_ld_trgs <- bind_rows(high_ld_trg_a, high_ld_trg_b) %>%
    group_by(target) %>%
    summarize(trg_prev = sum(trg_prev), .groups = "drop") %>%
    arrange(desc(trg_prev))
  # End loop if there are no more high LD targets
  if (identical(nrow(high_ld_trgs), 0L)) {
    warning(
      "All targets in high LD pairs (identified as rbarD >= ", 
      arg$pair_ld_thres, 
      ") have been removed and overall LD is still significant. Consider ", 
      "lowering high_ld_pair_thres."
    )
    break
  }
  # Add target in most high LD pairs to trgs2filter
  trgs2filter <- append(trgs2filter, high_ld_trgs$target[[1]])
  print(
    str_c(
      high_ld_trgs$target[[1]], 
      " is in ", 
      high_ld_trgs$trg_prev[[1]], 
      " high LD pair(s) and will be removed next."
    )
  )

  # Filter data and recompute overall LD
  overall_ld_res <- mh_data %>%
    filter(! locus %in% trgs2filter) %>%
    hubpopgen::tib2genind() %>%
    poppr::clonecorrect(strata = NA) %>%
    poppr::ia(sample = 999, plot = FALSE)
  print(
    str_c(
      "Overall LD has been recalculated. rbarD is ", 
      signif(overall_ld_res[["rbarD"]]), 
      " with a p-value of ", 
      overall_ld_res[["p.rD"]], 
      "."
    )
  )

}

# Write high LD targets to disk ----------------------------------------
tibble(target = trgs2filter) %>%
  write_csv(arg$high_ld_trgs)
