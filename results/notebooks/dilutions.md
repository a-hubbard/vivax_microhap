Effect of Dilution on Sequencing Yield
================
Alfred Hubbard

``` r
hap_setdiff <- function(hap_tib_a, hap_tib_b) {
  # This check is necessary because it is possible reads were obtained 
  # for the diluted case and not the undiluted case
  if (is.null(hap_tib_b)) {
    hap_tib_a
  } else {
    setdiff(hap_tib_a, hap_tib_b)
  }
}

# Read microhaplotype data ---------------------------------------------
mh_data <- read_csv(
    inp$mh_data,  
    col_types = cols(
      .default = col_character(), 
      n_read = col_integer(), 
      ParasiteDensity = col_double(), 
      CT_Pv = col_double(), 
      CT_Pf18s = col_double(), 
      CT_PfvarATS = col_double()
    ), 
    progress = FALSE
  ) %>%
  mutate(
    Dilution = factor(
      Dilution, 
      c("Neat", "1:10", "0.111111111", "0.736111111")
    )
  )

# Check the consistency of genotypes between dilutions -----------------
mh_data %>%
  # This is a check for contamination, so extremely low read count 
  # haplotypes are removed
  filter(n_read >= 10) %>%
  select(UCI_SID, Dilution, target, ASV) %>%
  nest(haplotypes = ASV) %>%
  pivot_wider(names_from = Dilution, values_from = haplotypes) %>%
  rename(undiluted_haps = Neat) %>%
  pivot_longer(
    c(`1:10`, `0.111111111`, `0.736111111`), 
    names_to = "dilution", 
    values_to = "diluted_haps"
  ) %>%
  # Filter out rows without any data in the dilution
  filter(map_lgl(diluted_haps, ~ ! is.null(.x))) %>%
  mutate(hap_diff = map2(diluted_haps, undiluted_haps, hap_setdiff)) %>%
  # Compute fraction of haplotypes that are in the dilution but not in 
  # the undiluted data
  mutate(n_hap_diluted = map_int(diluted_haps, nrow)) %>%
  mutate(n_hap_notneat = map_int(hap_diff, nrow)) %>%
  mutate(notneat_frac = n_hap_notneat / n_hap_diluted) %>%
  mutate(
    dilution = factor(dilution, c("1:10", "0.111111111", "0.736111111"))
  ) %>%
  group_by(dilution) %>%
  summarize(avg_notneat_frac = mean(notneat_frac), .groups = "drop")
```

    ## # A tibble: 3 × 2
    ##   dilution    avg_notneat_frac
    ##   <fct>                  <dbl>
    ## 1 1:10                   0.312
    ## 2 0.111111111            0.424
    ## 3 0.736111111            0.556

``` r
# Plot read counts, facetted by dilution -------------------------------
mh_data %>%
  group_by(UCI_SID, target, Dilution) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  # Add explicit missing data
  complete(UCI_SID, target, Dilution, fill = list(n_read = 0)) %>%
  ggplot(mapping = aes(x = target, y = UCI_SID, fill = n_read)) +
  geom_tile() +
  facet_wrap(vars(Dilution), ncol = 1) +
  scale_fill_fermenter(
    palette = "YlGnBu", 
    name = "Reads", 
    breaks = c(10, 100, 1000, 2000)
  ) +
  labs(x = "Marker", y = "Sample") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.position = "bottom", 
    legend.key.width = unit(0.35, "in")
  )
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/dilutions_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />
