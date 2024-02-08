Tidying of Microhaplotype Results
================
Alfred Hubbard

# Missing Data by Marker

``` r
# Read in final, filtered data -----------------------------------------
seqtab <- read_tsv(
    inp$seqtab, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  pivot_longer(-sample, names_to = "target_cigar", values_to = "n_read") %>%
  separate_wider_delim(target_cigar, ",", names = c("target", "cigar")) %>%
  rename(sample_id = sample)

# Read and join ASV IDs ------------------------------------------------
asv2cigar <- read_tsv(
  inp$asv2cigar, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)
mh_data <- seqtab %>%
  left_join(
    asv2cigar, 
    by = c("target" = "Amplicon", "cigar" = "CIGAR"), 
    relationship = "many-to-many"
  )

# Compute proportion of missing samples for each loci ------------------
mh_data <- mh_data %>%
  complete(target, sample_id, fill = list(n_read = 0))
n_samp_wdata <- mh_data %>%
  group_by(target, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(target) %>%
  summarize(n_samp_wdata = sum(n_read > 0), .groups = "drop")
mh_data <- mh_data %>%
  left_join(n_samp_wdata, by = "target")

# Plot data quantity by marker -----------------------------------------
n_samp_wdata_thres <- 20
mh_data %>%
  select(target, n_samp_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_samp_wdata)) +
  geom_histogram(bins = 25) +
  geom_vline(xintercept = n_samp_wdata_thres, color = "blue") +
  labs(x = "No. of Samples With Data", y = "No. of Markers")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This histogram shows the distribution of markers according to how many
samples have data for each marker. Given that the total number of
samples in the dataset is 68, these results are pretty good. Many
markers have data for all of the samples, and about half of the rest
have data for more than half of the samples.

A few markers, though, did not perform well. The vertical blue line
marks 20 samples with data. This threshold is used to filter markers for
subsequent analyses, removing those that fall below.

``` r
# Filter out loci with lots of missing data ----------------------------
mh_data_filtered <- mh_data %>%
  filter(n_samp_wdata >= n_samp_wdata_thres)
```

# Missing Data by Sample

``` r
# Compute proportion of missing loci for each sample -------------------
n_trg_wdata <- mh_data_filtered %>%
  group_by(target, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(sample_id) %>%
  summarize(n_trg_wdata = sum(n_read > 0), .groups = "drop")
mh_data_filtered <- mh_data_filtered %>%
  left_join(n_trg_wdata, by = "sample_id")

# Plot distribution of missing data proportion by marker ---------------
# This threshold guarantees ten markers worth of overlap between any 
# two samples
n_trg_wdata_thres <- round((n_distinct(mh_data_filtered$target)/2) + 10)
mh_data_filtered %>%
  select(sample_id, n_trg_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_trg_wdata)) +
  geom_histogram(bins = 20) +
  geom_vline(xintercept = n_trg_wdata_thres, color = "blue") +
  labs(x = "No. of Targets With Data", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This histogram shows the distribution of samples according to how many
markers have data for each sample. Note the marker filter of 20 samples
with data has already been applied, and this histogram is based on the
remaining data. There are 65 markers left in the dataset after the
filtering above. This histogram shows that almost all of the samples
have data for most of these markers.

The vertical blue line shows the chosen threshold of targets with data
that will be used for filtering the samples. This value, 42, was
selected with the following formula, rounding up:

$$
Thres. = \frac{No. Targets}{2} + 10
$$

This threshold ensures that every pair of samples has a minimum of 10
targets worth of overlap, hopefully ensuring sufficient data for
relatedness analysis. The blue line indicates only two samples will be
removed by applying this threshold.

``` r
# Filter out loci with lots of missing data ----------------------------
mh_data_filtered <- mh_data_filtered %>%
  filter(n_trg_wdata >= n_trg_wdata_thres) %>%
  # Now that the missing data filters have been applied, remove rows 
  # of missing data added by complete()
  filter(n_read > 0)
```

# Stats

Number of markers in the filtered dataset:

``` r
n_distinct(mh_data_filtered$target)
```

    ## [1] 65

Number of samples in the filtered dataset:

``` r
n_distinct(mh_data_filtered$sample_id)
```

    ## [1] 66

``` r
# Save filtered data ---------------------------------------------------
mh_data_filtered %>%
  # Calculate MOI. It would be more logical to do this in moi.Rmd, but 
  # that would complicate the pipeline.
  group_by(sample_id, target) %>%
  mutate(n_hap = n()) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(moi = max(n_hap)) %>%
  ungroup() %>%
  write_csv(out$microhap_tidy)
```
