# Create heatmaps showing relatedness values

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(patchwork)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--rel", help = "CSV file containing relatedness values"), 
  make_option(
    "--scale", 
    help = 
      "String specifying whether to make heatmaps at country or site level"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))

rel_heatmap <- function(rel_tib, metric) {
  if (identical(metric, "mean_r")) {
    fill_scale <- scale_fill_distiller(
      palette = "YlGnBu", 
      name = "Mean r"
    )
  } else if (identical(metric, "frac_high_r")) {
    fill_scale <- scale_fill_distiller(
      palette = "YlGnBu", 
      name = "Frac. High r"
    )
  }
  ggplot(
      data = rel_tib, 
      mapping = aes(x = group_a, y = group_b, fill = .data[[metric]])
    ) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    fill_scale +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.25, "in")
    )
}

# Read data ------------------------------------------------------------
if (identical(arg$scale, "site")) {
  rel <- read_csv(
      arg$rel, 
      col_types = cols(
        .default = col_character(), 
        mean_r = col_double(), 
        frac_signif = col_double(), 
        frac_high_r = col_double()
      )
    ) %>%
    rename(group_a = site_a, group_b = site_b)
} else if (identical(arg$scale, "country")) {
  rel <- read_csv(
      arg$rel, 
      col_types = cols(
        .default = col_character(), 
        mean_r = col_double()
      )
    ) %>%
    rename(group_a = country_a, group_b = country_b)
}

# Convert group columns to a factor to order plot labels ---------------
if (identical(arg$scale, "site")) {
  group_order <- c(
    "Oddar Meanchey", 
    "Ho Chi Min", 
    "Binh Phuoc", 
    "Dak O", 
    "Krong Pa"
  )
} else if (identical(arg$scale, "country")) {
  group_order <- c(
    "Colombia", 
    "Ethiopia", 
    "Afghanistan", 
    "Cambodia", 
    "Vietnam", 
    "Indonesia"
  )
}
rel <- rel %>%
  mutate(
    group_a = factor(group_a, levels = group_order), 
    group_b = factor(group_b, levels = group_order)
  )

# Plot and save --------------------------------------------------------
if (identical(arg$scale, "site")) {
  # Note the outer parentheses are necessary to make the tagging work
  fig <- (
      rel_heatmap(rel, "mean_r") | rel_heatmap(rel, "frac_high_r")
    ) +
    plot_annotation(tag_levels = "A")
  w <- 7.5
  h <- 4.5
} else if (identical(arg$scale, "country")) {
  fig <- rel_heatmap(rel, "mean_r")
  w <- 5
  h <- 5
}
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
