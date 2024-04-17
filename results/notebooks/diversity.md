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

# Genotype Accumulation Curve

A genotype accumulation curve was calculated to assess whether there are
sufficient markers in the dataset to distinguish individual samples.

First, loci with only one allele were filtered out of the dataset:

``` r
mh_data_gd <- read_rds(inp$mh_gd) %>%
  poppr::informloci()
```

    ## cutoff value: 1.08108108108108 % ( 2 samples ).

    ## MAF         : 0.01

    ## 
    ##  Found 15 uninformative loci 
    ##  ============================ 
    ##  15 loci found with a cutoff of 2 samples :
    ##  pvk12_124MI_151QK, pvk12_596KR, PvP01_01_v1_654301_654500,
    ## PvP01_02_v1_80401_80600, PvP01_03_v1_99001_99200,
    ## PvP01_04_v1_702801_703000, PvP01_04_v1_723401_723600,
    ## PvP01_04_v1_726601_726800, PvP01_04_v1_742901_743100,
    ## PvP01_04_v1_747401_747600, PvP01_06_v1_26301_26500,
    ## PvP01_07_v1_1467301_1467500, PvP01_07_v1_37801_38000,
    ## PvP01_10_v1_136201_136400, PvP01_14_v1_3091401_3091600 
    ##  15 loci found with MAF < 0.01 :
    ##  pvk12_124MI_151QK, pvk12_596KR, PvP01_01_v1_654301_654500,
    ## PvP01_02_v1_80401_80600, PvP01_03_v1_99001_99200,
    ## PvP01_04_v1_702801_703000, PvP01_04_v1_723401_723600,
    ## PvP01_04_v1_726601_726800, PvP01_04_v1_742901_743100,
    ## PvP01_04_v1_747401_747600, PvP01_06_v1_26301_26500,
    ## PvP01_07_v1_1467301_1467500, PvP01_07_v1_37801_38000,
    ## PvP01_10_v1_136201_136400, PvP01_14_v1_3091401_3091600

Then the curve was plotted:

``` r
poppr::genotype_curve(mh_data_gd, sample = 1000, quiet = TRUE)
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/diversity_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

From this figure, it is apparent that there are plenty of markers in the
dataset to distinguish different samples, even once loci with only one
allele were removed. In fact, it appears that only about 15 markers
would be sufficient for this purpose.

# Sample Diversity

To complement marker diversity, sample-level diversity was also
estimated with Nei’s expected heterozygosity and evenness was calculated
with the $E_5$ statistic. Two population strata were used for the entire
dataset: populations as defined by the MalariaGEN project (“Population”)
and countries. In addition, diversity was calculated at the site level
for Cambodia and Vietnam, to facilitate comparison with some of the
results in relatedness.Rmd.

## By Population

``` r
# Compute summary statistics, filtering to those of interest -----------
adegenet::setPop(mh_data_gd) <- ~Population
poppr::poppr(mh_data_gd) %>%
  select(Pop, N, MLG, Hexp, E.5)
```

    ##     Pop   N MLG  Hexp   E.5
    ## 1    AF  36  33 0.443 0.930
    ## 2   LAM  32  22 0.436 0.729
    ## 3  MSEA   2   2 0.478 1.000
    ## 4   WAS  16  14 0.451 0.947
    ## 5  ESEA  72  68 0.444 0.955
    ## 6  WSEA   9   9 0.477 1.000
    ## 7   OCE  18  18 0.375 1.000
    ## 8 Total 185 166 0.645 0.871

This table gives the number of individuals, the number of multi-locus
genotypes (MLGs) observed in those individuals, Nei’s expected
heterozygosity, and evenness for each population.

Diversity is generally fairly low and evenness is high, although Latin
America (LAM) has a lower evenness than other populations and Oceania
(OCE) has a particularly low diversity. The high evenness is not a
surprise given that the number of MLGs is high relative to the number of
individuals - most individuals constitute a unique MLG, meaning the
abundances of the MLGs tend to be low, and therefore even.

## By Country

``` r
# Compute summary statistics, filtering to those of interest -----------
adegenet::setPop(mh_data_gd) <- ~Country
poppr::poppr(mh_data_gd) %>%
  select(Pop, N, MLG, Hexp, E.5)
```

    ##            Pop   N MLG  Hexp   E.5
    ## 1     Ethiopia  36  33 0.443 0.930
    ## 2         Peru   6   6 0.495 1.000
    ## 3  Philippines   2   2 0.478 1.000
    ## 4  Afghanistan  16  14 0.451 0.947
    ## 5     Colombia  23  13 0.353 0.738
    ## 6      Vietnam  53  49 0.441 0.942
    ## 7     Thailand   9   9 0.477 1.000
    ## 8     Cambodia  19  19 0.441 1.000
    ## 9       Brazil   3   3 0.411 1.000
    ## 10   Indonesia  18  18 0.375 1.000
    ## 11       Total 185 166 0.645 0.871

Indonesia and Colombia stand out as having lower diversity than the
other countries, and Colombia also stands out as having lower evenness.
This makes sense in light of the population-level results.

## By Site

``` r
# Read metadata and make site/country key ------------------------------
site_country_key <- read_csv(
    inp$mh_csv, col_types = cols(
      .default = col_character(), 
      Lat = col_double(), 
      Long = col_double(), 
      Year = col_integer(), 
      `% callable` = col_double(), 
      Fws = col_double(), 
      F_MISS = col_double()
    ), 
    progress = FALSE
  ) %>%
  select(Site, Country) %>%
  distinct()

# Compute summary statistics, filtering to those of interest -----------
adegenet::setPop(mh_data_gd) <- ~Site
poppr::poppr(mh_data_gd) %>%
  left_join(site_country_key, by = c("Pop" = "Site")) %>%
  relocate(Country, .after = Pop) %>%
  filter(Country %in% c("Cambodia", "Vietnam")) %>%
  select(Pop, Country, N, MLG, Hexp, E.5)
```

    ##              Pop  Country  N MLG  Hexp  E.5
    ## 1       Krong Pa  Vietnam 18  15 0.428 0.89
    ## 2          Dak O  Vietnam  9   9 0.435 1.00
    ## 3     Binh Phuoc  Vietnam  8   8 0.496 1.00
    ## 4     Ho Chi Min  Vietnam 18  17 0.408 0.97
    ## 5 Oddar Meanchey Cambodia 19  19 0.441 1.00

As at the coarser levels of aggregation, we see that the sites in
Cambodia and Vietnam tend to have low diversity and high evenness.
