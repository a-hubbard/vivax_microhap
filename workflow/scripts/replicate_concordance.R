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

# Create combinations and compute CCC with no missing data -------------
rep_data_long <- allele_table %>%
  select(specimen_id, parasitemia, rep_num, locus, cigar, wsaf)
rep_ccc_nomiss <- rep_data_long %>%
  # Join with itself to add all combinations of replicate pairs
  inner_join(
    rep_data_long,
    by = c("specimen_id", "parasitemia", "locus", "cigar"),
    suffix = c("_a", "_b"),
    relationship = "many-to-many"
  ) %>%
  # Only keep unique pairs (avoid duplicates and self-comparisons)
  filter(rep_num_a < rep_num_b) %>%
  group_by(parasitemia) %>%
  summarise(
    ccc = list(DescTools::CCC(wsaf_a, wsaf_b)$rho.c),
    .groups = "drop"
  ) %>%
  unnest(ccc) %>%
  mutate(missing_cond = "No FPs and FNs")

# Create combinations and compute CCC, including missing data ----------
rep_data_long_inclmiss <- rep_data_long %>%
  complete(
    nesting(specimen_id, parasitemia, locus, cigar), 
    rep_num, 
    fill = list(wsaf = 0)
  )
rep_ccc_inclmiss <- rep_data_long_inclmiss %>%
  # Join with itself to add all combinations of replicate pairs
  inner_join(
    rep_data_long_inclmiss,
    by = c("specimen_id", "parasitemia", "locus", "cigar"),
    suffix = c("_a", "_b"),
    relationship = "many-to-many"
  ) %>%
  # Only keep unique pairs (avoid duplicates and self-comparisons)
  filter(rep_num_a < rep_num_b) %>%
  group_by(parasitemia) %>%
  summarise(
    ccc = list(DescTools::CCC(wsaf_a, wsaf_b)$rho.c),
    .groups = "drop"
  ) %>%
  unnest(ccc) %>%
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
  ggplot(mapping = aes(x = parasitemia, y = est)) +
  geom_point(size = 3) +
  geom_errorbar(
    mapping = aes(ymin = lwr.ci, ymax = upr.ci), 
    width = 0, 
    linewidth = 1
  ) +
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
