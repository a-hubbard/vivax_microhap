# Scatter plots showing percent OT reads versus parasite density

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--readcounts_metadata", 
    help = "CSV file containing read counts and metadata for May 2022 runs"
  ), 
  make_option(
    "--total_read_counts", 
    help = "TSV file containing total read counts from each well"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))

# Read in data ---------------------------------------------------------
# Read counts for each locus, with metadata
read_counts <- read_csv(
    arg$readcounts_metadata, 
    col_types = cols(
      .default = col_character(), 
      n_read = col_integer(), 
      CT = col_double(), 
      CT_Pv = col_double()
    ), 
    progress = FALSE
  ) %>%
  # Filter out runs with both SWGA and target pre-amplification
  filter(Treatment != "SWGA & Targ. Pre-amp.") %>%
  # Merge source and extraction columns
  mutate(
    source_extraction = if_else(
      Source == "Whole blood", 
      "Whole blood", 
      if_else(Extraction == "Chelex", "DBS+Chelex", "DBS+Kit")
    )
  ) %>%
  select(-Source, -Extraction)
# Total read counts from each well
total_read_counts <- read_tsv(
  arg$total_read_counts, 
  col_names = c("protocol_well", "n_read_total"), 
  col_types = cols(.default = col_character(), n_read_total = col_integer()), 
  progress = FALSE
)

# Compute percent on-target (OT) reads ---------------------------------
reads_sample_data <- read_counts %>%
  group_by(protocol_well, SID, source_extraction, Treatment, CT_Pv) %>%
  summarize(n_read_ot = sum(n_read), .groups = "drop") %>%
  left_join(total_read_counts, by = "protocol_well") %>%
  mutate(pct_read_ot = (n_read_ot / n_read_total) * 100)

# Compute parasite density ---------------------------------------------
reads_sample_data <- reads_sample_data %>%
  mutate(
    para_dens = if_else(
      source_extraction == "Whole blood", 
      (10^((CT_Pv - 41.663)/-3.289))/5, 
      (10^((CT_Pv - 41.663)/-3.289))
    )
  )

# Make plots and save --------------------------------------------------
fig <- reads_sample_data %>%
  ggplot(mapping = aes(x = para_dens, y = pct_read_ot)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(source_extraction), cols = vars(Treatment)) +
  labs(x = "Parasite Density", y = "On-Target Reads (%)") +
  theme_bw()
w <- 8
h <- 5
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
