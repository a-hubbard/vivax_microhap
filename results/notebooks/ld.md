Linkage Disequilibrium Test for Microhaplotype Data
================
Alfred Hubbard

To evaluate linkage disequilibrium (LD) between pairs of loci, the
$\bar{r}_{d}$ measure was used in favor of the more common index of
association. The index of association is a ratio of the variance
observed between individuals at each locus and the variance expected
under conditions of linkage equilibrium. This measure increases with the
number of loci, so $\bar{r}_{d}$ was developed as an alternative that
approximates the index of association while avoiding this limitation.

LD was estimated using the populations defined by the MalariaGEN
authors. A correction for clonal samples was performed prior to making
the LD calculation, because having multiple copies of the same
multi-locus genotype can bias the calculation. Significance was
estimated using permutations.

The distribution of *p*-values obtained between all sample pairs is
shown below.

``` r
# Read and tidy LD results ---------------------------------------------
pairwise_ld_res <- read_rds(inp$ld_res) %>%
  as_tibble(rownames = "locus_pair") %>%
  mutate(across(-locus_pair, as.double)) %>%
  separate(locus_pair, c("locus_a", "locus_b"), sep = ":")

# Visualize distribution of p-values -----------------------------------
pairwise_ld_res %>%
  ggplot(mapping = aes(x = p.rD)) +
  geom_histogram(bins = 20) +
  labs(x = "p-value", y = "No. of Locus Pairs")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ld_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

A uniform distribution is expected if the null hypothesis of no linkage
is true. This histogram shows a distinct left skew, suggesting there may
be linkage in our dataset.

However, many hypothesis tests were performed to obtain these results,
so before concluding that any pairs are in significant LD, a correction
for multiple testing should be performed. This was done with the
Bonferroni method. This approach is sometimes considered too
conservative, so a generous *p*-value of 0.05 was used to help
compensate for that possibility. The Bonferroni-corrected significance
threshold for this case is:

``` r
# Compute and apply Bonferroni correction for multiple testing ---------
bf_thres <- 0.05 / nrow(pairwise_ld_res)
bf_thres
```

    ## [1] 1.434309e-05

``` r
pairwise_ld_res_signif <- pairwise_ld_res %>%
  filter(p.rD < bf_thres)
```

There are 0 loci pairs that have significant LD at this threshold.
Therefore, no filtering is necessary to address LD.
