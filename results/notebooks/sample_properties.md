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
  # Remove returning travelers
  filter(! `Is returning traveller`) %>%
  select(
    -`QC pass`, 
    -`All samples same individual`, 
    -`Exclusion reason`, 
    -`Is returning traveller`
  )

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

This bar plot shows the sample count by year.

``` r
# Zoom in on recent data -----------------------------------------------
sample_info_initialfilter %>%
  filter(Year >= 2010) %>%
  ggplot(mapping = aes(x = Year)) +
  geom_bar() +
  facet_wrap(vars(Population)) +
  labs(y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/sample_properties_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This plot shows the same thing, with samples from before 2010 removed.
In addition, the plots have been facetted by region.

From this figure, it seems that Africa, Southeast Asia, and Latin
America may all have suitably dense sample sets. These are visualized in
more detail below.

``` r
# Zoom in on recent data -----------------------------------------------
sample_info_initialfilter %>%
  filter(Year >= 2010, Population %in% c("AF", "WSEA", "ESEA", "LAM")) %>%
  ggplot(mapping = aes(x = Year, fill = Country)) +
  geom_bar(position = "stack") +
  facet_wrap(vars(Population)) +
  labs(y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/sample_properties_files/figure-gfm/unnamed-chunk-4-1.png" width="100%" />

As above, except limited to Africa, Southeast Asia, and Latin America.
Bars are now stacked, with colors scaled by country.

From this, three clusters in space and time present themselves as likely
populations for paneljudge analysis: Vietnam and Cambodia in 2015 and
2016; Ethiopia in 2013; and Brazil, Colombia, and Peru in 2013 and 2014.

``` r
sample_info_filterbyyear <- sample_info_initialfilter %>%
  filter(
    (Population == "ESEA" & Year %in% c(2015, 2016)) |
    (Population == "LAM" & Year %in% c(2013, 2014)) |
    (Country == "Ethiopia" & Year == 2013)
  )
```

# Sample Sizes

In the final analysis set, this leaves us with a sample size of 164. The
sample sizes for each site are given in the following table:

``` r
# Compute sample sizes by country and site -----------------------------
sample_info_filterbyyear %>%
  group_by(Country, Site) %>%
  summarize(n_samp = n(), .groups = "drop")
```

    ## # A tibble: 16 × 3
    ##    Country  Site                                            n_samp
    ##    <chr>    <chr>                                            <int>
    ##  1 Brazil   Manaus                                               8
    ##  2 Cambodia Oddar Meanchey                                      19
    ##  3 Colombia Antioquia                                            1
    ##  4 Colombia Buenaventura                                         3
    ##  5 Colombia Choco                                                7
    ##  6 Colombia Colombia                                             2
    ##  7 Colombia Cordoba                                              1
    ##  8 Colombia Tierralta                                           16
    ##  9 Colombia Tumaco                                               2
    ## 10 Ethiopia Oromia                                              26
    ## 11 Ethiopia South Nations Nationalities and Peoples' Region     23
    ## 12 Peru     Delta 1                                              3
    ## 13 Vietnam  Binh Phuoc                                           8
    ## 14 Vietnam  Dak O                                                9
    ## 15 Vietnam  Ho Chi Min                                          18
    ## 16 Vietnam  Krong Pa                                            18

# Study Information

The following table shows which studies generated the data in the final
analysis set.

``` r
# Print studies present in filtered data -------------------------------
sample_info_filterbyyear %>%
  group_by(Study, Site, Country, Year) %>%
  summarize(n_samp = n(), .groups = "drop")
```

    ## # A tibble: 20 × 5
    ##    Study                     Site                           Country  Year n_samp
    ##    <chr>                     <chr>                          <chr>   <int>  <int>
    ##  1 1098-PF-ET-GOLASSA        Oromia                         Ethiop…  2013     26
    ##  2 1128-PV-MULTI-GSK         Ho Chi Min                     Vietnam  2015     10
    ##  3 1128-PV-MULTI-GSK         Ho Chi Min                     Vietnam  2016      8
    ##  4 1128-PV-MULTI-GSK         Manaus                         Brazil   2014      8
    ##  5 1128-PV-MULTI-GSK         Oddar Meanchey                 Cambod…  2015     19
    ##  6 1157-PV-MULTI-PRICE       Antioquia                      Colomb…  2014      1
    ##  7 1157-PV-MULTI-PRICE       Binh Phuoc                     Vietnam  2015      8
    ##  8 1157-PV-MULTI-PRICE       Choco                          Colomb…  2014      5
    ##  9 1157-PV-MULTI-PRICE       Colombia                       Colomb…  2014      2
    ## 10 1157-PV-MULTI-PRICE       Cordoba                        Colomb…  2014      1
    ## 11 1157-PV-MULTI-PRICE       Dak O                          Vietnam  2015      4
    ## 12 1157-PV-MULTI-PRICE       Dak O                          Vietnam  2016      5
    ## 13 1157-PV-MULTI-PRICE       Krong Pa                       Vietnam  2015     12
    ## 14 1157-PV-MULTI-PRICE       Krong Pa                       Vietnam  2016      6
    ## 15 1157-PV-MULTI-PRICE       South Nations Nationalities a… Ethiop…  2013     23
    ## 16 X0001-PV-MULTI-HUPALO2016 Buenaventura                   Colomb…  2013      3
    ## 17 X0001-PV-MULTI-HUPALO2016 Choco                          Colomb…  2013      2
    ## 18 X0001-PV-MULTI-HUPALO2016 Delta 1                        Peru     2013      3
    ## 19 X0001-PV-MULTI-HUPALO2016 Tierralta                      Colomb…  2013     16
    ## 20 X0001-PV-MULTI-HUPALO2016 Tumaco                         Colomb…  2013      2

``` r
# Change spelling of Ho Chi Minh ---------------------------------------
sample_info_filterbyyear <- sample_info_filterbyyear %>%
  mutate(Site = str_replace(Site, "Ho Chi Min", "Ho Chi Minh"))

# Write sample information to disk -------------------------------------
sample_info_filterbyyear %>%
  select(sample_id) %>%
  write_csv(out$filtered_samples, col_names = FALSE)
write_csv(sample_info_filterbyyear, out$filtered_metadata)
```
