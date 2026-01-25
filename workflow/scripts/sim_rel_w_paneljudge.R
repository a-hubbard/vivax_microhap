# Simulate relatedness estimates for the specified panel of targets

# Load required libraries ----------------------------------------------
library(paneljudge)
# These libraries will be referenced without the library name and so 
# should be loaded second
library(dplyr)
library(magrittr)
library(optparse)                                                                
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--bed", help = "BED file describing panel"), 
  make_option(
    "--af", 
    help = "TSV file containing allele frequencies to use in simulation"
  ), 
  make_option(
    "--shared_functions", 
    help = "R file containing shared functions"
  ), 
  make_option(
    "--n_pair", 
    type = "integer", 
    default = 30, 
    help = "Number of pairs to simulate relatedness for"
  ), 
  make_option(
    "--seed", 
    type = "integer", 
    default = 1, 
    help = "Seed value for random number generation"
  ), 
  make_option("--out", help = "TSV to contain output")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$bed <- "../../results/primer_mapping/captured_seq_coords_good_amp.bed"
  arg$af <- "../../results/af/af_PvGAP_AF.tsv"
  arg$shared_functions <- "../shared_functions.R"
  arg$n_pair <- 1
  arg$seed <- 1
}

source(arg$shared_functions)
set.seed(arg$seed)

# Simulated relatedness estimates with paneljudge
sim_rel_w_panel <- function(af_long, distances, n_pair = 5) {

  # Convert allele frequencies to matrix
  af_wide <- af_dcifer2wide(af_long)
  # Join to distances, to ensure order matches
  distances_af_wide <- distances %>%
    left_join(af_wide, by = "target_id")
  af_matrix <- distances_af_wide %>%
    select(-target_id, -distance) %>%
    as.matrix()
  rownames(af_matrix) <- distances_af_wide$target_id

  ### Copied from paneljudge vignette with light modifications ###
  # Define data-generating r values
  rs <- c("0.01"=0.01, "0.25"=0.25, "0.50"=0.50, "0.75"=0.75, "0.99"=0.99)
  k <- 5 # Data-generating switch rate parameter value
  n <- n_pair # Number of pairs to per simulate per r in rs
  fs <- af_matrix # example marker allele frequencies
  ds <- distances_af_wide$distance # distances between marker mid-points
  ### End of copied code ###

  # Simulate relatedness estimates
  # The warn_fs argument of these functions doesn't seem to work
  quiet_simulate_Ys <- quietly(paneljudge::simulate_Ys)
  quiet_compute_r_and_k_CIs <- quietly(paneljudge::compute_r_and_k_CIs)
  # Also, this function returns a lot of these warnings: "Some 
  # per-marker allele exceed per-marker non-zero allele frequencies. 
  # Data are permissible due to non-zero epsilon." Based on my digging, 
  # I think isn't a concern, so I suppress the warnings.
  quiet_estimate_r_and_k <- quietly(paneljudge::estimate_r_and_k)
  expand.grid(datagen_r = rs, pair_num = 1:n) %>%
    as_tibble() %>%
    mutate(
      Ys = map(
        datagen_r, 
        ~ quiet_simulate_Ys(fs, ds, k, .x, warn_fs = FALSE)[[1]]
      )
    ) %>%
    mutate(
      krhat = map(
        Ys, 
        ~ quiet_estimate_r_and_k(fs, ds, .x, warn_fs = FALSE)[[1]]
      )
    ) %>%
    mutate(
      ci = map(
        krhat, 
        ~ quiet_compute_r_and_k_CIs(
          fs, 
          ds, 
          khat = .x['khat'], 
          rhat = .x['rhat'], 
          warn_fs = FALSE)[[1]]
      )
    ) %>%
    mutate(rhat = map_dbl(krhat, pluck, "rhat")) %>%
    mutate(ci_2_5 = map_dbl(ci, ~ .x["rhat", "2.5%"])) %>%
    mutate(ci_97_5 = map_dbl(ci, ~ .x["rhat", "97.5%"])) %>%
    select(-Ys, -krhat, -ci)

}

# Read in BED file describing panel markers ----------------------------
marker_info <- read_tsv(
    arg$bed, 
    col_types = cols(
      .default = col_character(), 
      start_pos = col_integer(), 
      end_pos = col_integer()
    ), 
    col_names = c(
      "chrom", 
      "start_pos", 
      "end_pos", 
      "target_id", 
      "score", 
      "strand"
    ), 
    progress = FALSE
  ) %>%
  select(-score, -strand) %>%
  # Filter to diversity loci
  filter(str_detect(target_id, "PvP01"))

# Read in allele frequencies -------------------------------------------
af <- read_tsv(
    arg$af, 
    col_types = cols(.default = col_character(), freq = col_double()), 
    progress = FALSE
  )

# Compute between marker distances -------------------------------------
marker_info <- marker_info %>%
  # For this calculation, it doesn't actually matter whether this gets 
  # the chromosomes in the "official" order
  arrange(chrom, start_pos) %>%
  # Compute coordinate in the middle of each target
  mutate(midpoint = (start_pos + end_pos) / 2) %>%
  mutate(
    # If this target is on the same chromosome as the next one, compute 
    # distance as the difference in midpoints. Otherwise, set to Inf.
    distance = if_else(
      chrom == lead(chrom) & ! is.na(lead(chrom)), 
      lead(midpoint) - midpoint, 
      Inf
    )
  )

# Simulate relatedness estimates and write to disk ---------------------
distances <- marker_info %>%
  select(target_id, distance)
sim_rel_w_panel(af, distances, n_pair = arg$n_pair) %>%
  write_tsv(arg$out)
