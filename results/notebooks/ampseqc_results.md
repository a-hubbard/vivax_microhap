Visualization of Read Counts from AmpSeQC
================
Alfred Hubbard

# Read Counts

## May 2022

``` r
plot_failed_reads <- function(read_counts) {
  read_counts %>%
    select(
      sample_id, 
      `__no_feature`, 
      `__ambiguous`, 
      `__too_low_aQual`, 
      `__not_aligned`, 
      `__alignment_not_unique`
    ) %>%
    pivot_longer(-sample_id, names_to = "reason", values_to = "n_read") %>%
    ggplot(mapping = aes(x = n_read)) +
    geom_histogram(bins = 30) +
    facet_wrap(vars(reason)) +
    labs(x = "No. of Reads", y = "No. of Samples")
}

# Read in read counts --------------------------------------------------
may2022_read_counts <- read_tsv(
    inp$may2022_read_counts, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  rename(sample_id = sample)
may2022_read_counts_byloc <- may2022_read_counts %>%
  select(
    -`__no_feature`, 
    -`__ambiguous`, 
    -`__too_low_aQual`, 
    -`__not_aligned`, 
    -`__alignment_not_unique`
  ) %>%
  pivot_longer(-sample_id, names_to = "locus", values_to = "n_read")

# Plot distribution of failed reads by reason --------------------------
plot_failed_reads(may2022_read_counts)
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ampseqc_results_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

These histograms show the number of samples with a certain number of
reads that failed for each reason identified by AmpSeQC. There are quite
a few reads that failed because of the “no_feature” condition, which,
according to the `htseq-count` documentation, means that they aligned to
the reference but did not correspond to any of the provided features.

The total number of reads that passed AmpSeQC’s filters is:

``` r
read_count_heatmap <- function(reads, 
                               fill_scale = scale_fill_fermenter(
                                 palette = "YlGnBu", 
                                 name = "Reads", 
                                 breaks = c(10, 100, 1000, 2000)
                               ), 
                               theme_tweaks = NULL) {
  ggplot(
      data = reads, 
      mapping = aes(x = locus, y = sample_id, fill = n_read)
    ) +
    geom_tile() +
    fill_scale +
    labs(x = "Locus", y = "Replicate") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    ) +
    theme_tweaks
}

# Print total read count -----------------------------------------------
sum(may2022_read_counts_byloc$n_read)
```

    ## [1] 3289700

``` r
# Plot read counts by locus and replicate ------------------------------
may2022_read_counts_byloc %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ampseqc_results_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

This heatmap shows the number of reads identified for each locus/sample
combination. It is expected that some loci will have no reads, as the
primer set used for this run was reduced relative to the original panel
design. As for the remaining loci, read counts vary to a substantial
extent by replicate, but this is also expected to some extent because
several different library prep protocols are represented in this figure.
These are compared elsewhere.

## UCI 12/23

``` r
# Read in read counts --------------------------------------------------
uci1223_read_counts <- read_tsv(
    inp$uci1223_read_counts, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  rename(sample_id = sample)
uci1223_read_counts_byloc <- uci1223_read_counts %>%
  select(
    -`__no_feature`, 
    -`__ambiguous`, 
    -`__too_low_aQual`, 
    -`__not_aligned`, 
    -`__alignment_not_unique`
  ) %>%
  pivot_longer(-sample_id, names_to = "locus", values_to = "n_read")

# Plot distribution of failed reads by reason --------------------------
plot_failed_reads(uci1223_read_counts)
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ampseqc_results_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

These histograms show the number of samples with a certain number of
reads that failed for each reason identified by AmpSeQC. For all of the
failure conditions, the read counts are at or near 0, indicating no
cause for concern. Because the “no_feature” category is essentially
empty in this case, the reads falling into that category for the May
2022 data were not investigated further.

The total number of reads that passed AmpSeQC’s filters is:

``` r
# Print total read count -----------------------------------------------
sum(uci1223_read_counts_byloc$n_read)
```

    ## [1] 6707396

``` r
# Plot read counts by locus and replicate ------------------------------
uci1223_read_counts_byloc %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ampseqc_results_files/figure-gfm/unnamed-chunk-4-1.png" width="100%" />

This heatmap shows the number of reads identified for each locus/sample
combination. As above, some loci are expected to lack data because
primers were not included for those loci. Of the remaining loci, we
obtained adequate read counts for many replicates, but for many other
replicates we obtained almost nothing.

# Selected Loci

To select loci to include in the final panel, loci with good
amplification in both tests with field samples are identified. First,
only loci with data in both runs are considered.

``` r
# Identify loci present in both datasets -------------------------------
may2022_loci_wdata <- may2022_read_counts_byloc %>%
  filter(n_read > 0) %$%
  locus
uci1223_loci_wdata <- uci1223_read_counts_byloc %>%
  filter(n_read > 0) %$%
  locus
loci_good_amp <- union(may2022_loci_wdata, uci1223_loci_wdata)
```

Then, the read count heatmaps are recreated showing only these loci.

``` r
may2022_read_counts_byloc %>%
  filter(locus %in% loci_good_amp) %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ampseqc_results_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

``` r
uci1223_read_counts_byloc %>%
  filter(locus %in% loci_good_amp) %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ampseqc_results_files/figure-gfm/unnamed-chunk-7-1.png" width="100%" />

There are one or two loci in each run that consistently did not amplify
well, but they are not the same loci between runs. Therefore, we decided
not to remove any more loci. It has been demonstrated that all loci in
this panel have amplified well in at least one sequencing run, and we
can hope that under favorable conditions the set represented in these
figures will all yield good results.

There are 88 loci in this final set, and they are:

``` r
loci_good_amp
```

    ##  [1] "PvP01_01_v1_131101_131300"   "PvP01_01_v1_388001_388200"  
    ##  [3] "PvP01_01_v1_395301_395500"   "PvP01_01_v1_612101_612300"  
    ##  [5] "PvP01_01_v1_654301_654500"   "PvP01_01_v1_746501_746700"  
    ##  [7] "PvP01_02_v1_358801_359000"   "PvP01_02_v1_432501_432700"  
    ##  [9] "PvP01_02_v1_518101_518300"   "PvP01_02_v1_660501_660700"  
    ## [11] "PvP01_02_v1_660901_661100"   "PvP01_02_v1_80401_80600"    
    ## [13] "PvP01_03_v1_178801_179000"   "PvP01_03_v1_196901_197100"  
    ## [15] "PvP01_03_v1_334901_335100"   "PvP01_03_v1_717501_717700"  
    ## [17] "PvP01_03_v1_99001_99200"     "PvP01_04_v1_401001_401200"  
    ## [19] "PvP01_04_v1_401301_401500"   "PvP01_04_v1_702801_703000"  
    ## [21] "PvP01_04_v1_723401_723600"   "PvP01_04_v1_726601_726800"  
    ## [23] "PvP01_04_v1_742901_743100"   "PvP01_04_v1_747401_747600"  
    ## [25] "PvP01_05_v1_1352001_1352200" "PvP01_05_v1_601701_601900"  
    ## [27] "PvP01_06_v1_161001_161200"   "PvP01_06_v1_197901_198100"  
    ## [29] "PvP01_06_v1_220501_220700"   "PvP01_06_v1_26301_26500"    
    ## [31] "PvP01_06_v1_272601_272800"   "PvP01_06_v1_56701_56900"    
    ## [33] "PvP01_06_v1_596001_596200"   "PvP01_07_v1_1010001_1010200"
    ## [35] "PvP01_07_v1_1218801_1219000" "PvP01_07_v1_1417201_1417400"
    ## [37] "PvP01_07_v1_1467301_1467500" "PvP01_07_v1_37801_38000"    
    ## [39] "PvP01_07_v1_500001_500200"   "PvP01_07_v1_697401_697600"  
    ## [41] "PvP01_07_v1_784501_784700"   "PvP01_08_v1_1016101_1016300"
    ## [43] "PvP01_08_v1_1191801_1192000" "PvP01_08_v1_1279701_1279900"
    ## [45] "PvP01_08_v1_1473101_1473300" "PvP01_09_v1_1062901_1063100"
    ## [47] "PvP01_09_v1_1379201_1379400" "PvP01_09_v1_1564401_1564600"
    ## [49] "PvP01_09_v1_1814501_1814700" "PvP01_09_v1_282801_283000"  
    ## [51] "PvP01_09_v1_769601_769800"   "PvP01_10_v1_1072001_1072200"
    ## [53] "PvP01_10_v1_136201_136400"   "PvP01_10_v1_152201_152400"  
    ## [55] "PvP01_10_v1_874201_874400"   "PvP01_10_v1_878201_878400"  
    ## [57] "PvP01_10_v1_934401_934600"   "PvP01_11_v1_1068401_1068600"
    ## [59] "PvP01_11_v1_1483801_1484000" "PvP01_11_v1_2008301_2008500"
    ## [61] "PvP01_12_v1_1436301_1436500" "PvP01_12_v1_1861301_1861500"
    ## [63] "PvP01_12_v1_2032801_2033000" "PvP01_12_v1_3000601_3000800"
    ## [65] "PvP01_12_v1_327601_327800"   "PvP01_12_v1_461201_461400"  
    ## [67] "PvP01_12_v1_467301_467500"   "PvP01_12_v1_769701_769900"  
    ## [69] "PvP01_13_v1_1116901_1117100" "PvP01_13_v1_192401_192600"  
    ## [71] "PvP01_13_v1_258501_258700"   "PvP01_13_v1_396201_396400"  
    ## [73] "PvP01_14_v1_125901_126100"   "PvP01_14_v1_1887501_1887700"
    ## [75] "PvP01_14_v1_2261601_2261800" "PvP01_14_v1_2314101_2314300"
    ## [77] "PvP01_14_v1_2759601_2759800" "PvP01_14_v1_2861901_2862100"
    ## [79] "PvP01_14_v1_2957301_2957500" "PvP01_14_v1_3091401_3091600"
    ## [81] "pvcrt_o.10k.indel"           "pvdbp.503"                  
    ## [83] "pvdhfr.57FLI.58SR.61TM"      "pvdhps.382AG.383AG"         
    ## [85] "pvdhps.553AG"                "pvk12.124MI.151QK"          
    ## [87] "pvk12.596KR"                 "pvmdr1.1076FL"

``` r
tibble(locus = loci_good_amp) %>%
  write_csv(out$loci_good_amp)
```
