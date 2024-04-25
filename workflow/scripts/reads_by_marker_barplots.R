# Box plots showing read counts by marker for the different protocols

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
    "--selected_trgs", 
    help = "CSV file identifying targets selected for final panel"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   readcounts_metadata = "../../results/may2022_readcounts_metadata.csv", 
#   selected_trgs = "../../results/trgs2filter/good_amp.csv"
# )

# Read in data ---------------------------------------------------------
reads_by_marker <- read_csv(
    arg$readcounts_metadata, 
    col_types = cols(
      .default = col_character(), 
      n_read = col_integer(), 
      CT = col_double(), 
      CT_Pv = col_double()
    ), 
    progress = FALSE
  ) %>%
  # Filter to protocol conditions of interest
  filter(
    Source == "DBS", 
    Extraction == "Kit", 
    Treatment %in% c("sWGA", "Targ. Pre-amp.")
  ) %>%
  select(-Source, -Extraction)
selected_trgs <- read_csv(
  arg$selected_trgs, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)

# Filter markers and compute mean read count across samples ------------
mean_read_counts <- selected_trgs %>%
  # Use left_join to filter to selected markers
  left_join(reads_by_marker, by = "target") %>%
  # Make all missing data explicit
  complete(SID, Treatment, target, fill = list(n_read = 0)) %>%
  group_by(target, Treatment) %>%
  summarize(mean_read_count = mean(n_read), .groups = "drop")

# Make plots and save --------------------------------------------------
fig <- mean_read_counts %>%
  ggplot(mapping = aes(x = target, y = mean_read_count)) +
  geom_col() +
  facet_wrap(vars(Treatment), ncol = 1) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Target", y = "Mean Read Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
w <- 10
h <- 6
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
