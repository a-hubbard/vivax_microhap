# Create boxplots of concordance among replicates

# Load required libraries ----------------------------------------------
library(DescTools)
# These libraries will be referenced without the library name and so 
# should be loaded second
library(dplyr)
library(ggplot2)
library(magrittr)
library(optparse)                                                                
library(purrr)
library(readr)
library(stringr)
library(tidyr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--allele_table", help = "TSV file containing microhaplotypes"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
if (interactive()) {
  arg$allele_table <- 
    "../../results/microhap/serialdil2/serialdil2_microhap.tsv"
  arg$out_base <- "../../results/figs/suppfig3"
}

# Compute the concordance of each pair of replicates
compute_rep_ccc <- function(rep_a_data, rep_b_data, incl_missing = FALSE) {
  rep_a_data <- rep_a_data %>%
    rename(rep_a_wsaf = wsaf)
  rep_b_data <- rep_b_data %>%
    rename(rep_b_wsaf = wsaf)
  if (incl_missing) {
    full_join(rep_a_data, rep_b_data, by = c("locus", "cigar")) %>%
      replace_na(list(rep_a_wsaf = 0, rep_b_wsaf = 0)) %$%
      DescTools::CCC(rep_a_wsaf, rep_b_wsaf) %$%
      rho.c %$%
      est
  } else {
    inner_join(rep_a_data, rep_b_data, by = c("locus", "cigar")) %$%
      DescTools::CCC(rep_a_wsaf, rep_b_wsaf) %$%
      rho.c %$%
      est
  }
}

# Read in data ---------------------------------------------------------
allele_table <- read_tsv(
  arg$allele_table, 
  col_types = cols(
    .default = col_character(), 
    parasitemia = col_integer(), 
    n_read = col_integer(), 
    wsaf = col_double()
  ), 
  progress = FALSE
)

# Add arbitrary replicate numbers --------------------------------------
rep_nums <- allele_table %>%
  distinct(specimen_id, parasitemia, rep_id) %>%
  group_by(specimen_id, parasitemia) %>%
  mutate(rep_num = row_number()) %>%
  ungroup()
allele_table <- allele_table %>%
  left_join(rep_nums, by = c("specimen_id", "parasitemia", "rep_id"))

# Build data frame with all combinations of replicates -----------------
rep_data <- allele_table %>%
  select(-rep_id, -n_read) %>%
  nest(allele_data = c(locus, cigar, wsaf))
rep_combos <- expand.grid(rep_a = 1:3, rep_b = 1:3) %>%
  as_tibble() %>%
  filter(rep_a != rep_b)
rep_data_combos <- rep_data %>%
  rename(rep_a = rep_num) %>%
  rename(rep_a_data = allele_data) %>%
  left_join(rep_combos, by = "rep_a", relationship = "many-to-many") %>%
  left_join(
    rep_data, 
    by = c("specimen_id", "parasitemia", "rep_b" = "rep_num")
  ) %>%
  rename(rep_b_data = allele_data)

# Compute CCC with and without FPs and FNs -----------------------------
rep_ccc_nomiss <- rep_data_combos %>%
  mutate(ccc = map2_dbl(rep_a_data, rep_b_data, compute_rep_ccc)) %>%
  mutate(missing_cond = "No FPs and FNs")
rep_ccc_inclmiss <- rep_data_combos %>%
  mutate(
    ccc = map2_dbl(
      rep_a_data, 
      rep_b_data, 
      compute_rep_ccc, 
      incl_missing = TRUE
    )
  ) %>%
  mutate(missing_cond = "Incl. FPs and FNs")
rep_ccc <- bind_rows(rep_ccc_nomiss, rep_ccc_inclmiss)

# Plot and save --------------------------------------------------------
fig <- rep_ccc %>%
  arrange(parasitemia) %>%
  mutate(
    parasitemia = factor(
      parasitemia, 
      levels = unique(parasitemia), 
      ordered = TRUE
    )
  ) %>%
  ggplot(mapping = aes(x = parasitemia, y = ccc)) +
  geom_boxplot() +
  facet_wrap(vars(missing_cond)) +
  labs(
    x = expression("Parasite Density (parasites" ~ "/" ~ mu ~ "L)"), 
    y = "Concordance Correlation Coefficient"
  )
w <- 8
h <- 5.5
ggsave(
  str_c(arg$out_base, ".pdf"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".png"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
