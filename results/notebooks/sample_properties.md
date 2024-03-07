Properties of Pv4 Samples After Variant Filtering
================
Alfred Hubbard

# Proportion of Missingness by Sample

``` r
# Read in and join sample information ----------------------------------
metadata <- read_tsv(
  inp$sample_metadata, 
  col_types = cols(
    .default = col_character(), 
    Lat = col_double(), 
    Long = col_double(), 
    Year = col_integer(), 
    `% callable` = col_double(), 
    `QC pass` = col_logical(), 
    `Is returning traveller` = col_logical()
  ), 
  progress = FALSE
)
coi <- read_tsv(
  inp$sample_coi, 
  col_types = cols(
    .default = col_character(), 
    Fws = col_double()
  ), 
  progress = FALSE
)
prop_miss <- read_tsv(
    inp$sample_missingness, 
    col_types = cols(
      .default = col_integer(), 
      INDV = col_character(), 
      F_MISS = col_double()
    ), 
    progress = FALSE
  ) %>%
  select(INDV, F_MISS)
sample_info <- left_join(metadata, coi, by = "Sample") %>%
  left_join(prop_miss, by = c("Sample" = "INDV")) %>%
  rename(sample_id = Sample)

# Filter samples based on QC, polyclonality, and longitudinal studies --
sample_info_initialfilter <- sample_info %>%
  filter(`QC pass`) %>%
  # Curiously, there are 41 samples without Fws values after the ones 
  # that failed QC are removed - the documentation says Fws was 
  # calculated for all QC pass samples
  filter(! is.na(Fws)) %>%
  # Remove polyclonal samples
  filter(Fws >= 0.95) %>%
  # Remove longitudinal samples
  filter(! str_detect(`All samples same individual`, ",")) %>%
  select(-`QC pass`, -`All samples same individual`, -`Exclusion reason`)

# Plot missingness -----------------------------------------------------
sample_info_initialfilter %>%
  ggplot(mapping = aes(x = F_MISS)) +
  geom_histogram(bins = 20) +
  labs(x = "Proportion of Missing Data", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/sample_properties_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This histogram shows the proportions of missing variants in each sample,
once samples that fail QC, polyclonal samples, and samples from
longitudinal studies are removed. The proportions of missing data in
these remaining samples are quite low, so no filter for missing data is
necessary.

# Sample Distribution by Year

``` r
# Plot sample size by year ---------------------------------------------
sample_info_initialfilter %>%
  ggplot(mapping = aes(x = Year)) +
  geom_bar() +
  labs(y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/sample_properties_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

This bar plot shows the sample count by year. Some samples in the
dataset are quite old. This is impressive, but we wish to focus on
recent samples in this study.

``` r
# Zoom in on recent data -----------------------------------------------
sample_info_initialfilter %>%
  filter(Year >= 2010) %>%
  ggplot(mapping = aes(x = Year)) +
  geom_bar() +
  labs(y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/sample_properties_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This plot shows the same thing, with samples from before 2010 removed.
We want to avoid analyzing a large span of years, to focus on genetic
relatedness between geographic locations, not through time. From this
plot, we can see that the years 2013 through 2016 have the highest
concentration of data. However, we judge four years to still be too
great a range of time, so we select the last two years from this period,
2015 and 2016, for further analysis.

``` r
sample_info_filterbyyear <- sample_info_initialfilter %>%
  filter(Year > 2014, Year < 2017)
```

# Sample Sizes

In the final analysis set, this leaves us with a sample size of 187. The
sample sizes for each site are given in the following table:

``` r
# Compute sample sizes by country and site -----------------------------
sample_info_filterbyyear %>%
  group_by(Country, Site) %>%
  summarize(n_samp = n(), .groups = "drop")
```

    ## # A tibble: 20 × 3
    ##    Country     Site                                  n_samp
    ##    <chr>       <chr>                                  <int>
    ##  1 Afghanistan Jalalabad                                 12
    ##  2 Afghanistan Laghman                                    4
    ##  3 Brazil      Manaus                                     3
    ##  4 Cambodia    Oddar Meanchey                            19
    ##  5 Colombia    Choco                                      3
    ##  6 Colombia    Santa Cecilia                              3
    ##  7 Colombia    Tierralta                                 17
    ##  8 Ethiopia    Amhara                                    15
    ##  9 Ethiopia    Gondar                                     4
    ## 10 Ethiopia    Jimma                                     17
    ## 11 Indonesia   Papua Indonesia                           18
    ## 12 Indonesia   Papua Indonesia (returning traveller)      2
    ## 13 Peru        Iquitos                                    6
    ## 14 Philippines Rio Tuba                                   2
    ## 15 Thailand    Mae Sot                                    4
    ## 16 Thailand    Umphang                                    5
    ## 17 Vietnam     Binh Phuoc                                 8
    ## 18 Vietnam     Dak O                                      9
    ## 19 Vietnam     Ho Chi Min                                18
    ## 20 Vietnam     Krong Pa                                  18

For the relatedness analysis, we will be merging some geographies and
removing others to create populations that each have an adequate sample
size.

``` r
# Write sample information to disk -------------------------------------
sample_info_filterbyyear %>%
  select(sample_id) %>%
  write_csv(out$filtered_samples, col_names = FALSE)
write_csv(sample_info_filterbyyear, out$filtered_metadata)
```
