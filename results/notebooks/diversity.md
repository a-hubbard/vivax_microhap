Marker Diversity, Sample Diversity, and MOI
================
Alfred Hubbard

# Marker Diversity

``` r
# Read in microhaplotype data ------------------------------------------
mh_data <- read_csv(
    inp$mh_csv, 
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
  ) %>%
  select(target, sample_id, ASV, pop, moi)

# Read and join ASV sequences ------------------------------------------
asv_seqs_seqfadna <- seqinr::read.fasta(inp$asv_seqs, as.string = TRUE)
asv_seqs <- tibble(
  asv_id = names(asv_seqs_seqfadna), 
  asv_seq = map_chr(asv_seqs_seqfadna, pluck, 1)
)
mh_data <- mh_data %>%
  left_join(asv_seqs, by = c("ASV" = "asv_id"))

# filter(mh_data, target == "PvP01_14_v1_2861901_2862100") %$%
#   as.matrix(asv_seq) %>%
#   pegas::haplotype() %>%
#   pegas::nuc.div()
```

# Sample Diversity

``` r
# Read in microhaplotype genind object ---------------------------------
mh_data_gd <- read_rds(inp$mh_gd)

# Compute summary statistics, filtering to those of interest -----------
poppr::poppr(mh_data_gd) %>%
  select(Pop, N, MLG, Hexp)
```

    ##     Pop  N MLG  Hexp
    ## 1  ABOB  5   5 0.480
    ## 2  ALPC 14  14 0.306
    ## 3  UKUN  4   4 0.281
    ## 4  VI89 40  40 0.280
    ## 5 Total 63  63 0.318

This table gives the number of individuals, the number of multi-locus
genotypes (MLGs) observed in those individuals, and Nei’s expected
heterozygosity for each population.

This suggests that the overall genetic diversity is fairly low, likely
driven by low evenness for most of the markers (i.e., most markers are
dominated by one genotype).

# MOI

## AmpSeq

``` r
# Plot distribution of MOI ---------------------------------------------
mh_data %>%
  select(sample_id, moi) %>%
  distinct() %>%
  ggplot(mapping = aes(x = moi)) +
  geom_bar() +
  labs(x = "MOI", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/diversity_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This bar plot shows the MOI values, calculated for each sample as the
maximum number of unique haplotypes across all loci, obtained from the
Neafsey lab malaria amplicon pipeline. Two and three are most prevalent,
and seemingly no infections are monoclonal.

## SeekDeep

``` r
# Read in SeekDeep output ----------------------------------------------
sd_moi <- read_tsv(
    inp$sd_moi, 
    col_names = c("filename", "n_hap"), 
    col_types = cols(.default = col_character(), n_hap = col_integer()), 
    progress = FALSE
  ) %>%
  separate(filename, c("target", "final", "file_basename"), sep = "/") %>%
  select(-final) %>%
  filter(str_detect(file_basename, "OLD")) %>%
  mutate(sample_id = str_extract(file_basename, "PV_SWGA_OLD[0-9]+")) %>%
  select(-file_basename) %>%
  group_by(sample_id) %>%
  summarize(moi = max(n_hap), .groups = "drop")

# Plot distribution of MOI ---------------------------------------------
sd_moi %>%
  ggplot(mapping = aes(x = moi)) +
  geom_bar() +
  labs(x = "MOI", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/diversity_files/figure-gfm/unnamed-chunk-4-1.png" width="100%" />

As above, but for the SeekDeep results. MOI values tend to be higher in
this case, with three and four being most prevalent. Again, no
monoclonal infections are present.

Without replicates, SeekDeep has been shown to identify false positives
(Early et al., 2019). This means it is expected that SeekDeep will
predict more unique haplotypes than DADA2. Therefore, A) the results of
this comparison are expected and B) we are inclined to interpret the
AmpSeq results as more representative of reality.
