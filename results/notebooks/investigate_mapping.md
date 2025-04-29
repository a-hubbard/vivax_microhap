Investigation of Primer Mapping for Issues
================
Alfred Hubbard

# Primer Mapping Success

The set of primers that was output by `refmt_primers.R`, and therefore
intended for subsequent use in the pipeline, but not included in the
mapping results, is below:

``` r
# Read in primers and mapping results ----------------------------------
primers_all <- read_tsv(
  inp$primers, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)
mapping_results <- read_tsv(
    inp$target_coords, 
    col_names = c("chrom", "start_pos", "end_pos", "target", "score", "strand"), 
    col_types = cols(
      .default = col_character(), 
      start_pos = col_integer(), 
      end_pos = col_integer()
    ), 
    progress = FALSE
  ) %>%
  mutate(target_length = end_pos - start_pos) %>%
  select(target, target_length)

# Identify primers that did not map successfully -----------------------
setdiff(primers_all$target, mapping_results$target)
```

    ## character(0)

This set is empty, which is what we want to see. This indicates mapping
was successful for each target.

The opposite, the difference between the mapping results and the input
primer set, is:

``` r
setdiff(mapping_results$target, primers_all$target)
```

    ## character(0)

This is empty, which is also expected. A non-empty set would indicate a
bug.

# Captured Sequence Length

``` r
# Plot lengths of captured sequences -----------------------------------
mapping_results %>%
  ggplot(mapping = aes(x = target_length)) +
  geom_histogram(bins = 30) +
  labs(x = "Length of Captured Sequence", y = "No. of Targets")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/investigate_mapping_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

All of the captured sequence lengths should be substantially less than
500 bp (because we are using 250 bp reads), and this is the case.
