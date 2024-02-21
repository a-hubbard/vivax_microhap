Linkage Disequilibrium Test for Microhaplotype Data
================
Alfred Hubbard

``` r
# Read in data and covert to the loci class ----------------------------
mh_data_gd <- read_rds(inp$microhap_gd)
mh_data_l <- pegas::genind2loci(mh_data_gd)

# Create tibble of locus name/number pairs -----------------------------
locus_nums <- attr(mh_data_l, "locicol")
locus_names <- names(mh_data_l)[locus_nums]
locus_num_name_key <- tibble(locus_name = locus_names, locus_num = locus_nums)
ld_res_pegas <- combn(locus_nums, 2) %>%
  t() %>%
  as_tibble() %>%
  rename(locus_a_num = V1, locus_b_num = V2) %>%
  left_join(locus_num_name_key, by = c("locus_a_num" = "locus_num")) %>%
  rename(locus_a_name = locus_name) %>%
  left_join(locus_num_name_key, by = c("locus_b_num" = "locus_num")) %>%
  rename(locus_b_name = locus_name)

# Read in microhaplotype data table with metadata ----------------------
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
)

# Join and filter by chromosome information ----------------------------
chrom_info <- mh_data %>%
  select(target, chrom) %>%
  distinct()
ld_res_pegas <- ld_res_pegas %>%
  left_join(chrom_info, by = c("locus_a_name" = "target")) %>%
  rename(chrom_a = chrom) %>%
  left_join(chrom_info, by = c("locus_b_name" = "target")) %>%
  rename(chrom_b = chrom) %>%
  # Regardless of what the test shows, we know that targets on 
  # different chromosomes cannot be physically linked
  filter(chrom_a == chrom_b) %>%
  select(-chrom_a, -chrom_b)

# Estimate LD ----------------------------------------------------------
ld_res_pegas <- ld_res_pegas %>%
  # Note that pegas prints a warning message about unphased genotypes 
  # when an individual is missing data at a given locus, regardless of 
  # whether the genotypes actually are unphased (or, indeed, regardless 
  # of whether the ploidy is even greater than 1)
  mutate(
    ld_res = map2(
      locus_a_num, 
      locus_b_num, 
      # Note that while in the locicol attribute the pop column is 
      # included in locus numbering (meaning the first locus column is 
      # two), this function wants the locus numbering to start at 1...
      ~pegas::LD(mh_data_l, locus = c(.x - 1, .y - 1))
    )
  ) %>%
  mutate(T2_pval = map_dbl(ld_res, pluck, "T2", "P-val")) %>%
  # Remove NA p-values, which indicate at least one of the loci only 
  # has one allele
  filter(! is.na(T2_pval))

# Visualize distribution of p-values -----------------------------------
ld_res_pegas %>%
  ggplot(mapping = aes(x = T2_pval)) +
  geom_histogram(bins = 30) +
  labs(x = "T_2 Test p-value", y = "No. of Locus Pairs")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ld_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

The distribution of *p*-values is given above. It is clearly not
uniform, suggesting some linkage is present in the dataset.

``` r
bf_thres <- 0.05 / nrow(ld_res_pegas)
bf_thres
```

    ## [1] 0.0005813953

``` r
ld_res_pegas_signif <- ld_res_pegas %>%
  filter(T2_pval < bf_thres)
ld_res_pegas_signif %>%
  select(-locus_a_num, -locus_b_num, -ld_res)
```

    ## # A tibble: 6 × 3
    ##   locus_a_name                locus_b_name                 T2_pval
    ##   <chr>                       <chr>                          <dbl>
    ## 1 PvP01_04_v1_401001_401200   PvP01_04_v1_401301_401500   3.58e-16
    ## 2 PvP01_04_v1_401301_401500   PvP01_04_v1_747401_747600   4.44e- 8
    ## 3 PvP01_09_v1_1564401_1564600 PvP01_09_v1_282801_283000   3.18e- 4
    ## 4 PvP01_09_v1_1814501_1814700 PvP01_09_v1_769601_769800   1.21e- 4
    ## 5 PvP01_12_v1_1436301_1436500 PvP01_12_v1_2032801_2033000 1.46e- 5
    ## 6 PvP01_12_v1_2032801_2033000 PvP01_12_v1_3000601_3000800 1.21e- 6

``` r
for (i in seq(1:nrow(ld_res_pegas_signif))) {
  print(
    str_c(
      "Comparison between ", 
      ld_res_pegas_signif$locus_a_name[[i]], 
      " and ", 
      ld_res_pegas_signif$locus_b_name[[i]]
    )
  )
  print(ld_res_pegas_signif$ld_res[[i]])
}
```

    ## [1] "Comparison between PvP01_04_v1_401001_401200 and PvP01_04_v1_401301_401500"
    ## $`Observed frequencies`
    ##       ASV20 ASV54 ASV138 ASV89
    ## ASV14    52     1      4     0
    ## ASV59     0     3      0     1
    ## ASV72     0     2      0     0
    ## 
    ## $`Expected frequencies`
    ##           ASV20     ASV54    ASV138      ASV89
    ## ASV14 47.047619 5.4285714 3.6190476 0.90476190
    ## ASV59  3.301587 0.3809524 0.2539683 0.06349206
    ## ASV72  1.650794 0.1904762 0.1269841 0.03174603
    ## 
    ## $`Correlations among alleles`
    ##            ASV20      ASV54      ASV138      ASV89
    ## ASV14  0.7054131 -0.8157895  0.08447772 -0.3914407
    ## ASV59 -0.5661211  0.5807843 -0.06779661  0.4877532
    ## ASV72 -0.3936909  0.5580998 -0.04714700 -0.0229961
    ## 
    ## $`LRT (G-squared)`
    ## [1] NaN
    ## 
    ## $`Pearson's test (chi-squared)`
    ## [1] 118.9079
    ## 
    ## $T2
    ##           T2           df        P-val 
    ## 8.482966e+01 6.000000e+00 3.580548e-16 
    ## 
    ## [1] "Comparison between PvP01_04_v1_401301_401500 and PvP01_04_v1_747401_747600"
    ## $`Observed frequencies`
    ##        ASV25 ASV99 ASV153
    ## ASV20     50     2      0
    ## ASV54      6     0      0
    ## ASV138     4     0      0
    ## ASV89      0     0      1
    ## 
    ## $`Expected frequencies`
    ##            ASV25      ASV99     ASV153
    ## ASV20  49.523810 1.65079365 0.82539683
    ## ASV54   5.714286 0.19047619 0.09523810
    ## ASV138  3.809524 0.12698413 0.06349206
    ## ASV89   0.952381 0.03174603 0.01587302
    ## 
    ## $`Correlations among alleles`
    ##              ASV25       ASV99      ASV153
    ## ASV20   0.09349470  0.08328077 -0.27612739
    ## ASV54   0.07254763 -0.05874735 -0.04120428
    ## ASV138  0.05822225 -0.04714700 -0.03306802
    ## ASV89  -0.56796183 -0.02299610  1.00000000
    ## 
    ## $`LRT (G-squared)`
    ## [1] NaN
    ## 
    ## $`Pearson's test (chi-squared)`
    ## [1] 126.8077
    ## 
    ## $T2
    ##           T2           df        P-val 
    ## 4.511276e+01 6.000000e+00 4.444856e-08 
    ## 
    ## [1] "Comparison between PvP01_09_v1_1564401_1564600 and PvP01_09_v1_282801_283000"
    ## $`Observed frequencies`
    ##        ASV57 ASV197 ASV158
    ## ASV61     48      0      0
    ## ASV162     5      1      1
    ## ASV121     6      0      0
    ## 
    ## $`Expected frequencies`
    ##            ASV57     ASV197     ASV158
    ## ASV61  46.426230 0.78688525 0.78688525
    ## ASV162  6.770492 0.11475410 0.11475410
    ## ASV121  5.803279 0.09836066 0.09836066
    ## 
    ## $`Correlations among alleles`
    ##              ASV57      ASV197      ASV158
    ## ASV61   0.35378379 -0.24806947 -0.24806947
    ## ASV162 -0.51137189  0.35856858  0.35856858
    ## ASV121  0.06081116 -0.04264014 -0.04264014
    ## 
    ## $`LRT (G-squared)`
    ## [1] NaN
    ## 
    ## $`Pearson's test (chi-squared)`
    ## [1] 31.90315
    ## 
    ## $T2
    ##           T2           df        P-val 
    ## 2.098992e+01 4.000000e+00 3.181299e-04 
    ## 
    ## [1] "Comparison between PvP01_09_v1_1814501_1814700 and PvP01_09_v1_769601_769800"
    ## $`Observed frequencies`
    ##       ASV67 ASV126
    ## ASV48    46      6
    ## ASV98     1      4
    ## 
    ## $`Expected frequencies`
    ##           ASV67   ASV126
    ## ASV48 42.877193 9.122807
    ## ASV98  4.122807 0.877193
    ## 
    ## $`Correlations among alleles`
    ##            ASV67     ASV126
    ## ASV48  0.5091953 -0.5091953
    ## ASV98 -0.5091953  0.5091953
    ## 
    ## $`LRT (G-squared)`
    ## [1] 10.74502
    ## 
    ## $`Pearson's test (chi-squared)`
    ## [1] 29.55791
    ## 
    ## $T2
    ##           T2           df        P-val 
    ## 1.477895e+01 1.000000e+00 1.208772e-04 
    ## 
    ## [1] "Comparison between PvP01_12_v1_1436301_1436500 and PvP01_12_v1_2032801_2033000"
    ## $`Observed frequencies`
    ##        ASV13 ASV16 ASV55 ASV78 ASV106
    ## ASV19     45     2     1     2      0
    ## ASV113     1     1     1     0      1
    ## ASV68      4     0     1     2      1
    ## 
    ## $`Expected frequencies`
    ##            ASV13     ASV16     ASV55     ASV78    ASV106
    ## ASV19  40.322581 2.4193548 2.4193548 3.2258065 1.6129032
    ## ASV113  3.225806 0.1935484 0.1935484 0.2580645 0.1290323
    ## ASV68   6.451613 0.3870968 0.3870968 0.5161290 0.2580645
    ## 
    ## $`Correlations among alleles`
    ##             ASV13       ASV16      ASV55       ASV78     ASV106
    ## ASV19   0.4833333 -0.07978313 -0.2700352 -0.20370138 -0.3726780
    ## ASV113 -0.3698788  0.24673990  0.2467399 -0.06896552  0.3236377
    ## ASV68  -0.2985562 -0.08679261  0.1374216  0.29060425  0.2020344
    ## 
    ## $`LRT (G-squared)`
    ## [1] NaN
    ## 
    ## $`Pearson's test (chi-squared)`
    ## [1] 53.21667
    ## 
    ## $T2
    ##           T2           df        P-val 
    ## 3.644509e+01 8.000000e+00 1.455340e-05 
    ## 
    ## [1] "Comparison between PvP01_12_v1_2032801_2033000 and PvP01_12_v1_3000601_3000800"
    ## $`Observed frequencies`
    ##        ASV30 ASV50 ASV164
    ## ASV13     42     7      0
    ## ASV16      0     3      0
    ## ASV55      1     2      0
    ## ASV78      1     2      1
    ## ASV106     2     0      0
    ## 
    ## $`Expected frequencies`
    ##            ASV30      ASV50     ASV164
    ## ASV13  36.950820 11.2459016 0.80327869
    ## ASV16   2.262295  0.6885246 0.04918033
    ## ASV55   2.262295  0.6885246 0.04918033
    ## ASV78   3.016393  0.9180328 0.06557377
    ## ASV106  1.508197  0.4590164 0.03278689
    ## 
    ## $`Correlations among alleles`
    ##             ASV30      ASV50      ASV164
    ## ASV13   0.4835457 -0.4163880 -0.26087460
    ## ASV16  -0.3982721  0.4167077 -0.02936101
    ## ASV55  -0.2222243  0.2364299 -0.02936101
    ## ASV78  -0.3101081  0.1703976  0.48733972
    ## ASV106  0.1051370 -0.1004857 -0.02376913
    ## 
    ## $`LRT (G-squared)`
    ## [1] NaN
    ## 
    ## $`Pearson's test (chi-squared)`
    ## [1] 66.02019
    ## 
    ## $T2
    ##           T2           df        P-val 
    ## 4.225892e+01 8.000000e+00 1.210834e-06

# Overall

Linkage disequilibrium was tested using the $\bar{r}_{d}$ measure. This
was used in favor of the more common index of association because the
index of association is a ratio of the variance observed between
individuals at each locus and the variance expected under conditions of
linkage equilibrium. This measure increases with the number of loci, so
$\bar{r}_{d}$ was developed as an alternative that approximates the
index of association while avoiding this limitation.

Having multiple copies of a given multi-locus genotype can bias the
estimation of LD. To avoid this bias, these duplicate genotypes were
removed using `poppr::clonecorrect`. LD was then estimated using
`poppr::ia`.

``` r
ld_res_poppr <- read_rds(inp$ld_res_poppr)
ld_res_poppr
```

    ##        Ia      p.Ia     rbarD      p.rD 
    ## 5.3697483 0.0010000 0.1091568 0.0010000

Higher values of $\bar{r}_{d}$ indicate higher levels of linkage
disequilibrium.

# Pairwise

To understand which pairs of loci are responsible for this LD,
$\bar{r}_{d}$ was calculated and plotted separately for each pair of
loci.

``` r
pairwise_ld_res <- poppr::clonecorrect(mh_data_gd, strata = NA) %>%
  poppr::pair.ia()
```

    ## Warning in poppr::clonecorrect(mh_data_gd, strata = NA): Strata is not set for
    ## mh_data_gd. Clone correct will be performed without population information.

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ld_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

Targets with gray swatches only have one allele in the dataset, and so
$\bar{r}_{d}$ cannot be calculated.

``` r
trgs2filter <- read_csv(
    inp$trgs2filter, 
    col_types = cols(target = col_character()), 
    progress = FALSE
  ) %$%
  target
trgs2filter
```

    ## [1] "PvP01_04_v1_401001_401200"   "PvP01_12_v1_1436301_1436500"

``` r
# Prepare microhaplotype data for conversion to genind -----------------
mh_data <- mh_data %>%
  select(target, sample_id, ASV, pop) %>%
  filter(! is.na(pop)) %>%
  # Required by adegenet
  mutate(target = str_replace_all(target, "\\.", "_")) %>%
  rename(sample_ID = sample_id, locus = target, allele = ASV)

# Filter data and estimate updated LD results --------------------------
mh_data %>%
  filter(! locus %in% trgs2filter) %>%
  hubpopgen::tib2genind() %>%
  poppr::clonecorrect(strata = NA) %>%
  poppr::ia(sample = 999, plot = FALSE)
```

    ## Warning: Ploidy assumed to be 1.

    ## Warning in poppr::clonecorrect(., strata = NA): Strata is not set for .. Clone
    ## correct will be performed without population information.

    ##        Ia      p.Ia     rbarD      p.rD 
    ## 4.9385086 0.0010000 0.1045474 0.0010000

``` r
pairwise_ld_res_filtered <- mh_data %>%
  filter(! locus %in% trgs2filter) %>%
  hubpopgen::tib2genind() %>%
  poppr::clonecorrect(strata = NA) %>%
  poppr::pair.ia()
```

    ## Warning: Ploidy assumed to be 1.

    ## Warning: Strata is not set for .. Clone correct will be performed without
    ## population information.

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/ld_files/figure-gfm/unnamed-chunk-8-1.png" width="100%" />

Targets with gray swatches only have one allele in the dataset, and so
$\bar{r}_{d}$ cannot be calculated.
