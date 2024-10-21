Test for Selective Pressure on each Locus
================
Alfred Hubbard

The Tajima test was used to look for an impact of natural selection on
the loci selected for this panel. Tajima exploited the fact that two
metrics of genetic diversity, the number of segregating sites and the
average number of pairwise nucleotide differences, differ in their
response to selection to create a test for selection based on the
difference between these two metrics.

Note that the malaria system violates the assumptions of the Tajima
test, because mating is non-random and recombination is present. Tajima
(1989) briefly discusses the impact of recombination on the test
statistic, stating that the variance is reduced and the test therefore
becomes more conservative. Our goal in testing for selection is to
screen out loci that may not be neutral prior to performing downstream
population genetic analyses. A too-conservative test for selection would
potentially lead to data being removed unnecessarily. We are comfortable
with this possibility, as we would rather filter out too much data than
too little.

Given that the assumptions of the test are not met, we do not treat the
results as concrete evidence that a given loci IS under selection -
rather we consider it an indication that the loci MAY BE under
selection. We interpret the results accordingly.

The test was performed with the pegas R package. It was performed
separately for the data corresponding to each population, as defined by
MalariaGEN.

Some of the tests generated warnings - the set of warning messages
obtained is printed below:

``` r
# Read in data ---------------------------------------------------------
mh_popbin <- read_rds(inp$mh_popbin)

# Perform Tajima test --------------------------------------------------
quiet_tajima <- quietly(pegas::tajima.test)
selec_res <- mh_popbin %>%
  mutate(selec_out = map(seqs_dnabin, quiet_tajima)) %>%
  mutate(selec_res = map(selec_out, pluck, "result")) %>%
  mutate(selec_warn = map(selec_out, pluck, "warnings")) %>%
  select(-seqs_dnabin)

# Analyze warnings -----------------------------------------------------
selec_res %>%
  select(locus, selec_warn) %>%
  unnest(selec_warn) %$%
  unique(selec_warn)
```

    ## [1] "no segregating sites"

This warning is expected for some of the locus/population combinations,
and is not cause for concern.

The distribution of *p*-values obtained for these tests is given below.

``` r
# Extract test statistic and p-value -----------------------------------
selec_res <- selec_res %>%
  mutate(taj_d = map_dbl(selec_res, pluck, "D")) %>%
  # According to Tajima (1989), the beta distribution typically 
  # provides a better approximation for the D statistic distribution 
  # than the normal distribution does
  mutate(pval = map_dbl(selec_res, pluck, "Pval.beta")) %>%
  select(-selec_out, -selec_warn, -selec_res) %>%
  # Remove tests with NaN for results, which means there were no 
  # segregating sites. These are excluded from the Bonferroni threshold 
  # calculation, with the logic that if data conditions prevent the 
  # test statistic from being calculated, it doesn't really "count" as 
  # performing a hypothesis test.
  filter(! is.nan(taj_d))

# Visualize distribution of p-values -----------------------------------
selec_res %>%
  ggplot(mapping = aes(x = pval)) +
  geom_histogram(bins = 20) +
  labs(x = "p-value", y = "No. of Tests")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/selection_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

If the null hypothesis is true, a uniform distribution of *p*-values is
expected. This distribution is more normal than uniform, with a slight
right skew. This suggests that some of the tests may be significant, but
likely not many.

To assess significance of these results, the Bonferroni correction for
multiple testing is applied. This approach is sometimes considered too
conservative, so a generous *p*-value of 0.05 was used to help
compensate for that possibility. The Bonferroni-corrected significance
threshold for this case is:

``` r
# Apply Bonferroni correction ------------------------------------------
bf_thres <- 0.05 / nrow(selec_res)
bf_thres
```

    ## [1] 0.000154321

``` r
selec_res_signif <- selec_res %>%
  filter(pval < bf_thres)
selec_res_signif %>%
  select(-chrom)
```

    ## # A tibble: 3 × 4
    ##   locus                       pop_id taj_d        pval
    ##   <chr>                       <chr>  <dbl>       <dbl>
    ## 1 pvcrt_o.10k.indel           WAS    -2.90 0          
    ## 2 pvcrt_o.10k.indel           ESEA   -3.14 0          
    ## 3 PvP01_10_v1_1072001_1072200 WAS    -2.56 0.000000115

The table above shows the locus/population combinations that yielded a
significant result at this corrected threshold. These results suggest
that pvcrt_o.10k.indel may be under selective pressure in the WAS and
ESEA populations, and PvP01_10_v1_1072001_1072200 may be under selective
pressure in the WAS population. These loci will be removed for
downstream population genetic analyses that assume neutral markers.

``` r
# Save list of loci to filter ---------------------------------------
selec_res_signif %>%
  select(locus) %>%
  distinct() %>%
  write_csv(out$loci2filter)
```
