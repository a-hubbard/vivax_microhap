Assessment of SNV and Indel Distances Found by DADA2
================
Alfred Hubbard

# SNV Distance

``` r
# Read in target lengths and ASV table ---------------------------------
target_lengths <- read_tsv(
    inp$target_coords, 
    col_names = c("chrom", "start_pos", "end_pos", "target"), 
    col_types = cols(
      .default = col_character(), 
      start_pos = col_integer(), 
      end_pos = col_integer()
    ), 
    progress = FALSE
  ) %>%
  mutate(target_length = end_pos - start_pos) %>%
  select(target, target_length)
asvtab <- read_tsv(
    inp$asvtab, 
    col_types = cols(
      .default = col_integer(), 
      hapid = col_character(), 
      strain = col_character(), 
      refid_PvP01 = col_character(), 
      indel_filter = col_character(), 
      bimera = col_logical()
    ), 
    progress = FALSE
  ) %>%
  # These will be removed by ASV_to_CIGAR.py, and so should not be 
  # analyzed here
  filter(indel_filter == "PASS", bimera == FALSE) %>%
  select(-indel_filter, -bimera) %>%
  left_join(target_lengths, by = c("refid_PvP01" = "target"))

# Plot distribution of SNV distances -----------------------------------
asvtab %>%
  ggplot(mapping = aes(x = snv_dist_from_PvP01)) +
  geom_histogram(bins = 30) +
  labs(x = "SNV Distance from PvP01", y = "No. of Haplotypes")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/check_snvs_indels_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This histogram shows the distribution of microhaplotypes according to
their SNV distance from the PvP01 reference.

``` r
# Plot distribution of SNV distances as a fraction of target length ----
asvtab %>%
  mutate(fractional_snv_dist = snv_dist_from_PvP01 / target_length) %>%
  ggplot(mapping = aes(x = fractional_snv_dist)) +
  geom_histogram(bins = 30) +
  labs(
    x = "SNV Distance from PvP01 as Fraction of Target Length", 
    y = "No. of Haplotypes"
  )
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/check_snvs_indels_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

This histogram shows the same thing as above, except expressed as a
fraction of the length of the target.

Considered together, these histograms show that there are two distinct
groups of amplicons: a larger group with a short distance to the
reference, and a smaller group with a much larger distance to the
reference. It is reasonable to assume that the second group consists of
chimeras and other erroneous genotypes. An SNV distance threshold of 50
will be used to remove these in `ASV_to_CIGAR.py`.

# Indel Distance

``` r
# Plot distribution of indel distances ---------------------------------
asvtab %>%
  ggplot(mapping = aes(x = indel_dist_from_PvP01)) +
  geom_histogram(bins = 30) +
  labs(x = "Indel Distance from PvP01", y = "No. of Haplotypes")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/check_snvs_indels_files/figure-gfm/unnamed-chunk-3-1.png" width="100%" />

This histogram shows the distribution of microhaplotypes according to
their indel distance from the PvP01 reference. An indel distance filter,
expressed as a fraction of target length, has already been used to label
the sequences for indel filtering by `postProc_dada2.R`. Amplicons with
this label were removed prior to the generation of this figure.

The remaining microhaplotypes are mostly quite close to PvP01. However,
there are a couple with somewhat large indel distances (19 and 24). For
now, an additional threshold will be applied in `ASV_to_CIGAR.py` to
remove these sequences, but if time permits it would be best to
investigate in more detail.
