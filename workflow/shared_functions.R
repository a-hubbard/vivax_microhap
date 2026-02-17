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

# Visualize cardinality versus effective cardinality
card_vs_effcard <- function(marker_metrics) {
  scale_max <- max(c(marker_metrics$cardinality, marker_metrics$eff_cardinality))
  scale_lim <- c(0, scale_max)
  marker_metrics %>%
    ggplot(mapping = aes(x = cardinality, y = eff_cardinality)) +
    geom_point(shape = 95, size = 3) +
    geom_abline(slope = 1, intercept = 0) +
    facet_grid(rows = vars(panel), cols = vars(pop_lbl)) +
    scale_x_continuous(limits = scale_lim) +
    scale_y_continuous(limits = scale_lim) +
    labs(x = "Cardinality", y = "Effective Cardinality")
}

# Compute and plot 95% CI widths of rhat by data-generating r
rhat_ci_width_boxplot <- function(rel_sim_res) {
  rel_sim_res %>%
    mutate(ci_width = ci_97_5 - ci_2_5) %>%
    filter(datagen_r %in% c(0.01, 0.15, 0.25, 0.5, 0.99)) %>%
    mutate(datagen_r = as.character(datagen_r)) %>%
    ggplot(mapping = aes(y = ci_width, x = datagen_r, color = panel)) +
    geom_boxplot() +
    facet_wrap(vars(pop_lbl)) +
    labs(
      x = expression("Data-generating" ~ italic(r)), 
      y = expression("95% CI Width of" ~ italic(hat(r))), 
      color = "Panel"
    ) +
    theme(legend.position = "bottom")
}

# Compute and plot RMSEs for each data-generating r
plot_rhat_rmse <- function(rel_sim_res) {
  rel_sim_res %>%
    group_by(panel, pop_lbl, datagen_r) %>%
    summarize(rmse = sqrt(mean((rhat - datagen_r)^2)), .groups = "drop") %>%
    ggplot(mapping = aes(x = datagen_r, y = rmse, color = panel)) +
    geom_line() +
    geom_point() +
    facet_wrap(vars(pop_lbl)) +
    labs(
      x = expression("Data-generating" ~ italic(r)), 
      y = expression("RMSE of" ~ italic(hat(r))), 
      color = "Panel"
    ) +
    theme(legend.position = "bottom")
}
