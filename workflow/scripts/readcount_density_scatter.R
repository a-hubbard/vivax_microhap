# Create scatterplot of parasite density versus read count

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(dplyr)
library(ggplot2)
library(magrittr)
library(optparse)                                                                
library(readr)
library(stringr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--s_read_counts_parasitemia", 
    help = "CSV file containing read counts and parasitemia values"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg <- list(
    s_read_counts_parasitemia = 
      "../../results/AmpSeq/uci1223/s_read_counts_parasitemia.csv"
  )
}

# Set default ggplot2 theme
theme_set(theme_bw())

# Read in read counts and parasitemia values ---------------------------
sample_mean_total_read_counts <- read_csv(
    arg$s_read_counts_parasitemia, 
    col_types = cols(
      .default = col_character(), 
      ct = col_double(), 
      n_read = col_double(), 
      parasitemia = col_double()
    ), 
    progress = FALSE
  )

# Plot parasitemia versus sample read count ----------------------------
fig <- sample_mean_total_read_counts %>%
  ggplot(mapping = aes(x = parasitemia, y = n_read)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = expression("Parasite Density (parasites" ~ "/" ~ mu ~ "L)"), 
    y = "Mean Reads per Replicate"
  )
w <- 5
h <- 4.5
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
