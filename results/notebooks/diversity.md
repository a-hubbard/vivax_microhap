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

    ## cutoff value: 1.06951871657754 % ( 2 samples ).

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
with the $E_5$ statistic. Two population strata were used: populations
as defined by the MalariaGEN project (“Population”) and “geographies”
defined in this project by either removing sites with little data or
merging them with other, nearby sites (“Geography”).

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
    ## 6   OCE  20  20 0.376 1.000
    ## 7  WSEA   9   9 0.477 1.000
    ## 8 Total 187 168 0.645 0.872

This table gives the number of individuals, the number of multi-locus
genotypes (MLGs) observed in those individuals, Nei’s expected
heterozygosity, and evenness for each population.

Diversity is generally fairly low and evenness is high, although Latin
America (LAM) has a lower evenness than other populations and Oceania
(OCE) has a particularly low diversity.

## By Geography

``` r
# Compute summary statistics, filtering to those of interest -----------
adegenet::setPop(mh_data_gd) <- ~Geography
# This throws a warning that it is filtering out individuals with 
# missing population information. This is the desired behavior, hence 
# why warnings are suppressed for this chunk.
poppr::poppr(mh_data_gd) %>%
  select(Pop, N, MLG, Hexp, E.5)
```

    ##            Pop   N MLG  Hexp   E.5
    ## 1        Jimma  17  17 0.438 1.000
    ## 2       Amhara  19  16 0.434 0.893
    ## 3  Afghanistan  16  14 0.451 0.947
    ## 4     Colombia  23  13 0.353 0.738
    ## 5     Krong Pa  18  15 0.428 0.890
    ## 6   Binh Phuoc  17  17 0.467 1.000
    ## 7   Ho Chi Min  18  17 0.408 0.970
    ## 8     Cambodia  19  19 0.441 1.000
    ## 9    Indonesia  18  18 0.375 1.000
    ## 10       Total 165 146 0.645 0.862

Again, overall genetic diversity is fairly low, likely driven by high
evenness for many of the markers (i.e., most markers are dominated by
one genotype). The geographies are fairly comparable, with Colombia and
Indonesia standing out as slightly less diverse than the other
locations. Evenness is generally high, with the exception of Colombia,
which is somewhat lower.
