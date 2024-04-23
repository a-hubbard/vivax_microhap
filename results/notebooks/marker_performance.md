Performance of Marker Panel in May 2022 Run
================
Alfred Hubbard

For reference, this table gives some of the metadata associated with the
samples used in this sequencing run:

``` r
read_count_heatmap <- function(reads, 
                               fill_scale = scale_fill_fermenter(
                                 palette = "YlGnBu", 
                                 name = "Reads", 
                                 breaks = c(10, 100, 1000, 2000)
                               ), 
                               theme_tweaks = NULL) {
  ggplot(
      data = reads, 
      mapping = aes(x = target, y = sample_id, fill = n_read)
    ) +
    geom_tile() +
    facet_wrap(vars(treatment_source), ncol = 1) +
    fill_scale +
    labs(x = "Marker", y = "Sample") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    ) +
    theme_tweaks
}

# Read in data ---------------------------------------------------------
mh_combined <- read_csv(
    inp$mh_combined, 
    col_types = cols(
      .default = col_integer(), 
      Primer = col_character(), 
      SeqRun = col_character(), 
      SID = col_character(), 
      Source = col_character(), 
      Dilution = col_character(), 
      Treatment = col_character(), 
      Primer_set = col_character(), 
      CT_Pv = col_double(), 
      CT_PfvarATS = col_double(), 
      CT_Pf18s = col_double(), 
      ParasiteDensity = col_double(), 
      Mixed = col_character(), 
      Median_reads_primer = col_double(), 
      Prop_value = col_double(), 
      Sequence = col_character(), 
      mh_allele_snps = col_character(), 
      mh_loci = col_character(), 
      TajimasD = col_double(), 
      Pnorm = col_double(), 
      Pbet = col_double()
    ), 
    progress = FALSE
  ) %>%
  select(-`...1`) %>%
  filter(SeqRun == "may2022", Primer_set == "OLD", Dilution == "Neat")

# Read, join, and print metadata ---------------------------------------
gbpcd_meta <- read_csv(
  inp$gbpcd_meta, 
  col_types = cols(.default = col_character(), Age = col_integer()), 
  progress = FALSE
)
mh_combined %>%
  left_join(gbpcd_meta, by = c("SID" = "UCI_SID")) %>%
  select(SID, DBS_ID, Age, Sex) %>%
  distinct()
```

    ## # A tibble: 10 × 4
    ##    SID           DBS_ID              Age Sex  
    ##    <chr>         <chr>             <int> <chr>
    ##  1 GBPCD_P04_D04 ABOB-202105-45067    10 F    
    ##  2 GBPCD_P06_A11 VI89-202106-901       4 M    
    ##  3 GBPCD_P07_H02 UKUN-202009-3918     27 M    
    ##  4 GBPCD_P08_D11 ABOB-202010-07273    27 F    
    ##  5 GBPCD_P08_F06 UKUN-202010-5083     19 M    
    ##  6 GBPCD_P13_F07 ABOB-202102-33880    12 F    
    ##  7 AH-46         <NA>                 NA <NA> 
    ##  8 AH-47         <NA>                 NA <NA> 
    ##  9 AH-48         <NA>                 NA <NA> 
    ## 10 AH-51         <NA>                 NA <NA>

``` r
# Calculate total reads ------------------------------------------------
total_reads <- mh_combined %>%
  rename(target = Primer, sample_id = SID) %>%
  unite(treatment_source, Treatment, Source) %>%
  group_by(target, sample_id, treatment_source) %>%
  summarize(n_read = sum(Reads), .groups = "drop") %>%
  complete(target, sample_id, treatment_source, fill = list(n_read = 0))

# Plot read counts for each sample/target for the DBS treatments -------
whole_blood_samples <- mh_combined %>%
  filter(str_detect(Source, "Whole blood")) %$%
  unique(SID)
total_reads %>%
  filter(
    ! sample_id %in% whole_blood_samples, 
    ! str_detect(treatment_source, "Whole blood")
  ) %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/marker_performance_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

These heatmaps show the total reads obtained for each sample/target
combination for all of the DBS treatments.

``` r
# Plot read counts for each sample/target for whole blood trials -------
total_reads %>%
  filter(
    sample_id %in% whole_blood_samples, 
    str_detect(treatment_source, "Whole blood")
  ) %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/marker_performance_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

These heatmaps show the same thing, except for the library prep methods
that utilized whole blood.

The point of this analysis is to identify any markers that consistently
don’t amplify and should not be included in the MalariaGEN analyses.
These plots indicate three markers that consistently yield poor results.
A marker set without these problematic markers is created for use in
subsequent analyses.

``` r
# Filter targets to analysis set and save ------------------------------
total_reads %>%
  select(target) %>%
  distinct() %>%
  filter(
    ! target %in% c(
      "PvP01_07_v1_500001_500200", 
      "PvP01_14_v1_1887501_1887700", 
      "PvP01_14_v1_2759601_2759800"
    )
  ) %>%
  write_csv(out$filtered_trgs)
```
