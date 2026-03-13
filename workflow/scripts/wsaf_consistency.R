# Create rainbow plots showing the consistency of WSAFs

# Load required libraries ----------------------------------------------
library(RColorBrewer)
# These libraries will be referenced without the library name and so 
# should be loaded second
library(dplyr)
library(ggplot2)
library(magrittr)
library(optparse)                                                                
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
  arg$out_base <- "../../results/figs/suppfig2"
}

# Create stacked bar plots showing WSAFs for each sample and locus
# This functions is a lightly-modified version of code shared with me 
# by Andres Aranda-Diaz, and the idea for this type of plot was 
# originally developed by Nicholas Hathaway
plot_sample_locus_wsafs <- function(alleles) {

  # Get unique locus values
  loci <- unique(alleles$locus)
  # Initialize an empty list to store color mappings
  color_list <- list()
  # Iterate over each locus and generate color mappings
  for (locus_oi in loci) {
    locus_alleles <- alleles %>% 
      filter(locus == locus_oi) %$%
      unique(cigar)
    
    # Generate rainbow colors for this locus
    locus_colors <- RColorBrewer::brewer.pal(length(locus_alleles), "Dark2")
    names(locus_colors) <- locus_alleles
    
    # Append to the main list
    color_list <- c(color_list, locus_colors)
  }

  # Combine all named vectors into one named vector
  combined_colors <- unlist(color_list)

  # Define y-axis positioning variables
  alleles <- alleles %>%
    arrange(desc(specimen_id), parasitemia, desc(rep_num)) %>%
    mutate(
      row_lbl = factor(row_lbl, levels = unique(row_lbl), ordered = TRUE)
    ) %>%
    # This effectively enumerates the row_lbls, starting at one
    mutate(row_lbl_index = as.numeric(row_lbl))%>%
    group_by(row_lbl_index, locus) %>%
    mutate(
      y_min = row_lbl_index - 1 + lag(cumsum(wsaf), default = 0),
      y_max = row_lbl_index - 1 + cumsum(wsaf)
    ) %>%
    ungroup() %>% 
    mutate(locus = factor(locus),
           locus_num = as.numeric(locus)
    )
  y_axis_lbls <- alleles %>%
    distinct(row_lbl, row_lbl_index)
  x_axis_lbls <- alleles %>%
    distinct(locus, locus_num)

  # Plot with y-axis labeled by row_lbl
  alleles %>%
    ggplot() +
    geom_rect(
      aes(
        ymin = y_min,
        ymax = y_max,
        xmin = locus_num - 0.4,
        xmax = locus_num + 0.4,
        fill = cigar
      )
    ) +
    scale_fill_manual(values = combined_colors) +
    geom_rect(
      aes(
        ymin = y_min,
        ymax = y_max,
        xmin = locus_num - 0.4,
        xmax = locus_num + 0.4
      ),
      fill = NA,
      color = "black",
      linewidth = 0.1
    ) +
    scale_y_continuous(
      breaks = y_axis_lbls$row_lbl_index - 0.5,
      labels = y_axis_lbls$row_lbl, 
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      breaks = x_axis_lbls$locus_num, 
      labels = x_axis_lbls$locus, 
      expand = c(0, 0) 
    ) +
    labs(x = "Locus", y = "Sample") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 5, angle = 90, hjust = 1), 
      axis.text.y = element_text(size = 5, hjust = 0.5), 
      axis.ticks.y = element_blank(), 
      panel.grid = element_blank(), 
      panel.background = element_blank(), 
      plot.background = element_blank()
    )
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

# Plot and save --------------------------------------------------------
fig <- allele_table %>%
  mutate(
    row_lbl = str_c(
      "Sample: ", 
      specimen_id, 
      ", Parasite Copy Number: ", 
      parasitemia, 
      ", Rep: ", 
      rep_num
    )
  ) %>%
  plot_sample_locus_wsafs()
w <- 10
h <- 15
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
