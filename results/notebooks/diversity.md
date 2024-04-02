Marker and Sample Diversity
================
Alfred Hubbard

# Marker Diversity

Nucleotide diversity was computed for each marker with
`pegas::nuc.div()`. Since this was done with pegas `haplotype` objects,
the frequencies of each haplotype were used to perform this calculation.

Note that when identifying haplotypes, pegas treats gaps as separate
characters (meaning a Hamming distance of 1 from any other character)
but ambiguities are treated as having a Hamming distance of 0 from
nucleotides included in that ambiguity. Even so, sometimes it is not
possible to assign a sequence with an ambiguity at a given position to
other haplotypes that lack ambiguities at that position. In this
scenario, pegas assigns this sequence to a new haplotype. That means
that there tend to be a few low frequency haplotypes for each marker
that were assigned this way because of ambiguities, and would be merged
with one of the other haplotypes if the ambiguity could be resolved.

All of that is to say the diversity numbers computed by pegas may be
slightly inflated. Relative comparisons should be valid, but absolute
interpretation, or richness statistics, should be avoided.

Note also that only gaps in the middle of the sequence are treated as
separate characters. Gaps at the end of the sequences are treated as Ns,
and so will generally be assigned to other haplotypes.

``` r
# Compute and visualize nucleotide diversity ---------------------------
nuc_div_res <- read_rds(inp$mh_hap) %>%
  # Compute nucleotide diversity, telling pegas to remove sites with 
  # missing data in a pairwise manner to preserve the maximum amount of 
  # data in each individual distance calculation
  mutate(nuc_div = map_dbl(hap, pegas::nuc.div, pairwise.deletion = TRUE)) %>%
  select(-hap)
nuc_div_res %>%
  ggplot(mapping = aes(x = target, y = nuc_div)) +
  geom_col() +
  labs(x = "Target", y = "Nucleotide Diversity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/diversity_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This bar plot shows the nucleotide diversity for each marker. There is
considerable variability in the diversity of these markers, and a number
don’t have any segregating sites, leading to a diversity of 0.

# Sample Diversity

To complement marker diversity, sample-level diversity was also
estimated with Nei’s expected heterozygosity. For this analysis, samples
were grouped by geography. Samples from sites with little data were
either removed or merged with other, nearby sites.

``` r
# Read in microhaplotype genind object ---------------------------------
mh_data_gd <- read_rds(inp$mh_gd)

# Compute summary statistics, filtering to those of interest -----------
poppr::poppr(mh_data_gd) %>%
  select(Pop, N, MLG, Hexp)
```

    ##            Pop   N MLG  Hexp
    ## 1  Afghanistan  16  14 0.371
    ## 2       Amhara  19  16 0.357
    ## 3   Binh Phuoc  17  17 0.384
    ## 4     Cambodia  19  19 0.364
    ## 5     Colombia  23  13 0.290
    ## 6   Ho Chi Min  18  17 0.335
    ## 7    Indonesia  18  18 0.308
    ## 8        Jimma  17  17 0.360
    ## 9     Krong Pa  18  15 0.351
    ## 10       Total 165 146 0.523

This table gives the number of individuals, the number of multi-locus
genotypes (MLGs) observed in those individuals, and Nei’s expected
heterozygosity for each population.

This suggests that the overall genetic diversity is fairly low, likely
driven by low evenness for many of the markers (i.e., most markers are
dominated by one genotype). The geographies are fairly comparable, with
Colombia and Indonesia standing out as slightly less diverse than the
other locations.
