Filter Microhaplotype Data Prior to LD Analysis
================
Alfred Hubbard

The method that we are using to estimate linkage disequilibrium cannot
handle polyclonal loci. This notebook assesses the impact of various
methods of filtering out polyclonal loci, and saves filtered data from
the selected method.

# Impact of Discarding Polyclonal Loci

The simplest approach is simply to remove all polyclonal loci.

## Reduction in Data Volume

``` r
# Read and prepare microhaplotype data ---------------------------------
mh_data <- read_csv(
    inp$microhap_tidy,  
    col_types = cols(
      .default = col_character(), 
      n_read = col_integer(), 
      n_samp_wdata = col_integer(), 
      n_trg_wdata = col_integer(), 
      ParasiteDensity = col_double(), 
      CT_Pv = col_double(), 
      CT_Pf18s = col_double(), 
      CT_PfvarATS = col_double(), 
      Age = col_integer(), 
      n_hap = col_integer(), 
      moi = col_integer()
    ), 
    progress = FALSE
  ) %>%
  select(-n_samp_wdata)
```

    ## Warning: The following named parsers don't match the column names: Age

``` r
# Plot amount of data per sample prior to filtering --------------------
mh_data %>%
  select(sample_id, n_trg_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_trg_wdata)) +
  geom_histogram(bins = 20) +
  labs(x = "No. of Targets With Data", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/filter4ld_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This histogram gives the data volume prior to filtering, as the
distribution of the number of targets that have data for each sample.

``` r
# Filter out polyclonal loci and recompute n_trg_wdata -----------------
mh_data_filtered <- mh_data %>%
  filter(n_hap == 1) %>%
  select(-n_trg_wdata, -n_hap) %>%
  group_by(sample_id) %>%
  mutate(n_trg_wdata = n_distinct(target)) %>%
  ungroup()

# Plot amount of data per sample prior to filtering --------------------
mh_data_filtered %>%
  select(sample_id, n_trg_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_trg_wdata)) +
  geom_histogram(bins = 20) +
  labs(x = "No. of Targets With Data", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/filter4ld_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

The same histogram, after the filter has been applied. The reduction in
the overall data volume is not particularly substantial.

# Amount of Data by Marker

``` r
# Compute number of samples with data for each loci --------------------
mh_data_filtered <- mh_data_filtered %>%
  group_by(target) %>%
  mutate(n_samp_wdata = n_distinct(sample_id)) %>%
  ungroup()

# Plot data quantity by marker -----------------------------------------
mh_data_filtered %>%
  select(target, n_samp_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_samp_wdata)) +
  geom_histogram(bins = 25) +
  labs(x = "No. of Samples With Data", y = "No. of Markers")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/filter4ld_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This histogram gives the distribution of the number of samples that have
data for each marker. This shows that a few targets would have very
little data left if this method were applied (fewer than 10 samples with
data). It would be preferable to avoid this.

# Impact of Discarding Minor Alleles

An alternative is simply to take the major allele, as determined by read
count, for each sample/target combination.

``` r
# Filter out minor alleles ---------------------------------------------
mh_data_filtered <- mh_data %>%
  group_by(sample_id, target) %>%
  filter(n_read == max(n_read)) %>%
  ungroup()
```

A surprisingly large number of alleles are tied for maximum read count
in the dataset. Further investigation of this phenomenon could be
fruitful.  
Information on these alleles is printed below.

``` r
# Print instances where the read count is tied -------------------------
mh_data_filtered %>%
  group_by(sample_id, target) %>%
  filter(n() > 1) %>%
  select(target, sample_id, n_read, ASV)
```

    ## # A tibble: 354 × 4
    ## # Groups:   sample_id, target [177]
    ##    target                      sample_id     n_read ASV   
    ##    <chr>                       <chr>          <int> <chr> 
    ##  1 PvP01_05_v1_1352001_1352200 PV_SWGA_OLD1      89 ASV7  
    ##  2 PvP01_05_v1_1352001_1352200 PV_SWGA_OLD1      89 ASV177
    ##  3 PvP01_11_v1_1483801_1484000 PV_SWGA_OLD1       7 ASV83 
    ##  4 PvP01_11_v1_1483801_1484000 PV_SWGA_OLD1       7 ASV364
    ##  5 PvP01_05_v1_1352001_1352200 PV_SWGA_OLD10     79 ASV7  
    ##  6 PvP01_05_v1_1352001_1352200 PV_SWGA_OLD10     79 ASV177
    ##  7 PvP01_05_v1_1352001_1352200 PV_SWGA_OLD11    103 ASV7  
    ##  8 PvP01_05_v1_1352001_1352200 PV_SWGA_OLD11    103 ASV177
    ##  9 PvP01_11_v1_1483801_1484000 PV_SWGA_OLD11      8 ASV83 
    ## 10 PvP01_11_v1_1483801_1484000 PV_SWGA_OLD11      8 ASV364
    ## # ℹ 344 more rows

For now, these alleles are simply filtered out for the LD calculation.

``` r
# Filter out alleles where the read count is tied ----------------------
mh_data_filtered <- mh_data_filtered %>%
  group_by(sample_id, target) %>%
  filter(n() == 1) %>%
  ungroup()

# Compute number of samples with data for each loci --------------------
mh_data_filtered <- mh_data_filtered %>%
  group_by(target) %>%
  mutate(n_samp_wdata = n_distinct(sample_id)) %>%
  ungroup()

# Plot data quantity by marker -----------------------------------------
mh_data_filtered %>%
  select(target, n_samp_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_samp_wdata)) +
  geom_histogram(bins = 25) +
  labs(x = "No. of Samples With Data", y = "No. of Markers")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/filter4ld_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

This histogram gives the distribution of the number of samples with data
for each marker, as above, after all minor alleles have been discarded.
Unfortunately, removing the alleles with tied read counts has had a
substantial impact on the data available for some markers. As with the
other method of filtering, some markers have fewer than 10 samples with
data.

Additional information on the affected targets is printed below:

``` r
# Print targets with very few samples after filtering ------------------
mh_data_filtered %>%
  filter(n_samp_wdata < 15) %>%
  select(sample_id, target, n_read, ASV)
```

    ## # A tibble: 10 × 4
    ##    sample_id     target                      n_read ASV   
    ##    <chr>         <chr>                        <int> <chr> 
    ##  1 PV_SWGA_OLD12 PvP01_05_v1_1352001_1352200    800 ASV156
    ##  2 PV_SWGA_OLD12 PvP01_12_v1_1861301_1861500    268 ASV157
    ##  3 PV_SWGA_OLD31 PvP01_12_v1_1861301_1861500    232 ASV157
    ##  4 PV_SWGA_OLD33 PvP01_12_v1_1861301_1861500    248 ASV157
    ##  5 PV_SWGA_OLD44 PvP01_05_v1_1352001_1352200   3208 ASV84 
    ##  6 PV_SWGA_OLD47 PvP01_12_v1_1861301_1861500    651 ASV169
    ##  7 PV_SWGA_OLD50 PvP01_05_v1_1352001_1352200    415 ASV154
    ##  8 PV_SWGA_OLD6  PvP01_12_v1_1861301_1861500     51 ASV157
    ##  9 PV_SWGA_OLD60 PvP01_12_v1_1861301_1861500   2544 ASV103
    ## 10 PV_SWGA_OLD60 PvP01_05_v1_1352001_1352200    109 ASV256

These targets are removed, to avoid basing LD calculations on very
little data.

In the end, we decided to use this method of filtering (discarding minor
alleles). While in both cases some targets are reduced to very little
data, comparing the upper ends of both histograms shows that this method
preserves more data overall.

``` r
# Write filtered data to disk ------------------------------------------
mh_data_filtered %>%
  filter(n_samp_wdata >=15) %>%
  write_csv(out$microhap_filtered)
```
