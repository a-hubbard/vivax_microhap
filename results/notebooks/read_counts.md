Visualization of Read Counts from AmpSeQC
================
Alfred Hubbard

# Read Counts

``` r
# Read in read counts --------------------------------------------------
read_counts <- read_tsv(
    inp$read_counts, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  rename(sample_id = sample)
read_counts_failed <- read_counts %>%
  select(
    sample_id, 
    `__no_feature`, 
    `__ambiguous`, 
    `__too_low_aQual`, 
    `__not_aligned`, 
    `__alignment_not_unique`
  ) %>%
  pivot_longer(-sample_id, names_to = "reason", values_to = "n_read")
read_counts_bytrg <- read_counts %>%
  select(
    -`__no_feature`, 
    -`__ambiguous`, 
    -`__too_low_aQual`, 
    -`__not_aligned`, 
    -`__alignment_not_unique`
  ) %>%
  pivot_longer(-sample_id, names_to = "target", values_to = "n_read")

# Plot distribution of failed reads by reason --------------------------
read_counts_failed %>%
  ggplot(mapping = aes(x = n_read)) +
  geom_histogram(bins = 30) +
  facet_wrap(vars(reason)) +
  labs(x = "No. of Reads", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/read_counts_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

These histograms show the number of samples with a certain number of
reads that failed for each reason identified by AmpSeQC. For all of the
failure conditions, the read counts are at or near 0, indicating no
cause for concern.

The total number of reads that passed AmpSeQC’s filters is:

``` r
read_count_heatmap <- function(reads, 
                               fill_scale = scale_fill_fermenter(
                                 palette = "YlGnBu", 
                                 name = "Reads", 
                                 breaks = c(10, 100, 1000, 2000)
                               ), 
                               theme_tweaks = NULL) {
  ggplot(data = reads, mapping = aes(x = target, y = sample_id, fill = n_read)) +
    geom_tile() +
    fill_scale +
    labs(x = "Marker", y = "Sample") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    ) +
    theme_tweaks
}

# Print total read count -----------------------------------------------
sum(read_counts_bytrg$n_read)
```

    ## [1] 3063481

``` r
# Plot read counts by marker and replicate -----------------------------
read_counts_bytrg %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/read_counts_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

This heatmap shows the number of reads identified for each marker/sample
combination. The results vary considerably by marker. There are many
where almost every sample has more than 10 reads, but there are also
quite a few where most samples have *fewer* than 10 reads. The results
vary less by sample, but there are still a handful that seem to have
yielded overall better results than the others.

# Comparison with SeekDeep

``` r
# Read in data ---------------------------------------------------------
sd_read_counts <- read_tsv(
    inp$sd_read_counts, 
    col_names = c("target", "sample_id", "n_read"), 
    col_types = cols(.default = col_character(), n_read = col_integer()), 
    progress = FALSE
  ) %>%
  complete(target, sample_id) %>%
  replace_na(list(n_read = 0))

# Plot read counts by marker and replicate -----------------------------
sd_read_counts %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/read_counts_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This heatmap is the same as the one in the above section, except A) it
is using the results from SeekDeep, which was run in a previous version
of the project, and B) it shows both NEW and OLD samples. These
designations refer to two separate sequencing runs, in which essentially
the same set of samples was used, except that in the NEW set a group of
SNP markers was added to the panel.

This figure shows the following:

1.  The results from SeekDeep are qualitatively similar to those from
    AmpSeQC, in that the yield varies a lot by marker and to a lesser
    extent by sample. Visual comparison of the two figures suggests the
    specific markers and samples that did and did not amplify well are
    more or less the same between the two methods.
2.  Across the board, the SNP markers did not amplify well in the NEW
    samples. As stated above, they were not included in the OLD samples.
3.  While the overall level of amplification does not look terribly
    different between NEW and OLD samples, there does appear to be some
    variability in terms of which markers and samples performed well
    between NEW and OLD. This is examined more in the next section.

## OLD Versus NEW

``` r
# Compute difference between NEW and OLD -------------------------------
new_old_diff <- sd_read_counts %>%
  mutate(sample_id = str_remove(sample_id, "_S[0-9]+_L001")) %>%
  mutate(sample_num = str_extract(sample_id, "[0-9]+")) %>%
  mutate(new_or_old = str_extract(sample_id, "NEW|OLD")) %>%
  select(-sample_id) %>%
  pivot_wider(names_from = new_or_old, values_from = n_read) %>%
  # SNPs were not sequenced for the OLD samples
  filter(! str_detect(target, "SNP")) %>%
  filter(! is.na(NEW) & ! is.na(OLD)) %>%
  mutate(new_old_diff = NEW - OLD)

# Visualize difference -------------------------------------------------
new_old_diff %>%
  rename(sample_id = sample_num, n_read = new_old_diff) %>%
  read_count_heatmap(
    fill_scale = scale_fill_fermenter(
      palette = "BrBG", 
      breaks = c(-25000, -10000, -1000, -100, 100, 1000, 10000, 25000),
      name = "NEW - OLD"
    ), 
    theme_tweaks = theme(legend.text = element_text(angle = 90, hjust = 1))
  )
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/read_counts_files/figure-gfm/unnamed-chunk-4-1.png" width="100%" />

This heatmap shows the difference ($NEW - OLD$) in read counts for each
sample/marker combination between new and old samples. This figure shows
that, for most of the combinations, the difference is not too
substantial. In the rest of the cases, there is not a clear consensus as
to which run produced more reads. It varies a good bit by sample (and to
a lesser extent by marker) with perhaps slightly better overall results
obtained with NEW.

However, this is the number of samples included in the NEW set:

``` r
sd_read_counts %>%
  select(sample_id) %>%
  filter(str_detect(sample_id, "NEW")) %>%
  n_distinct()
```

    ## [1] 67

And this is the number in the OLD set:

``` r
sd_read_counts %>%
  select(sample_id) %>%
  filter(str_detect(sample_id, "OLD")) %>%
  n_distinct()
```

    ## [1] 68

We decided that the fact that OLD contains one more sample than NEW
leans the tradeoffs in favor of OLD. We decided to proceed with this
sample set.

Note that these analyses were done with the SeekDeep output. This means
that the $NEW - OLD$ results might be slightly different if computed
with AmpSeQC output. While we expect this would not change the overall
conclusion, we still plan to change the analysis to use AmpSeQC results
in the future, if time permits.
