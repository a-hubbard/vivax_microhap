# Create heatmap of read counts
read_count_heatmap <- function(reads, 
                               fill_scale = scale_fill_fermenter(
                                 palette = "YlGnBu", 
                                 name = "Reads", 
                                 breaks = c(10, 100, 1000, 2000)
                               ), 
                               theme_tweaks = NULL, 
                               y_lbl = "Sample", 
                               title_lbl = NULL) {
  ggplot(
      data = reads, 
      mapping = aes(x = locus, y = sample_id, fill = n_read)
    ) +
    geom_tile() +
    fill_scale +
    labs(x = "Locus", y = y_lbl, title = title_lbl) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    ) +
    theme_tweaks
}

# Convert allele frequency from Dcifer format to wide format
af_dcifer2wide <- function(af_long) {
  af_wide <- af_long %>%
    group_by(target_id) %>%
    mutate(allele_id = str_c("Allele_", row_number())) %>%
    ungroup() %>%
    select(-seq) %>%
    pivot_wider(names_from = "allele_id", values_from = "freq") %>%
    mutate(across(-target_id, ~replace_na(.x, 0)))

  af_wide
}

# Create a scatterplot of percent good amplification versus parasitemia
pct_good_amp_vs_parasitemia <- function(data_tib) {
  data_tib %>%
    ggplot(mapping = aes(x = parasitemia, y = pct_loci_good_amp)) +
    geom_hline(yintercept = 75, color = "red", linewidth = 0.5) +
    geom_point() +
    geom_smooth(method = "lm", formula = "y ~ log10(x)") +
    scale_x_continuous(trans = "log10") +
    labs(
      x = expression("Parasite Density (parasites" ~ "/" ~ mu ~ "L)"), 
      y = expression("Percentage of Loci with" >= "10 Reads")
    )
}
