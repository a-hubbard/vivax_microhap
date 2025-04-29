Quantification of Microhaplotypes Identified by AmpSeq
================
Alfred Hubbard

# Final Data Quantity

## Reads

### By Replicate

``` r
read_count_heatmap_byrep <- function(reads, loi) {
  lib_reads <- reads %>%
    filter(lib_name == loi)
  samples_reps <- lib_reads %>%
    distinct(sample_id, rep_id) %>%
    arrange(sample_id)
  lib_reads %>%
    arrange(sample_id) %>%
    # Make rep_id a factor with ordering corresponding to samples_reps
    mutate(rep_id = factor(rep_id, samples_reps$rep_id, ordered = TRUE)) %>%
    ggplot(mapping = aes(x = locus, y = rep_id, fill = n_read)) +
    geom_tile() +
    # Use sample IDs from samples_reps as y-axis labels to show which 
    # replicates came from the same sample
    scale_y_discrete(label = samples_reps$sample_id) +
    scale_fill_fermenter(
      palette = "YlGnBu", 
      name = "Reads", 
      breaks = c(10, 100, 1000, 2000)
    ) +
    labs(x = "Locus", y = "Replicate", title = loi) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    )
}

# Read in final, filtered data -----------------------------------------
variant_read_counts <- read_tsv(
    inp$cigar_table, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  pivot_longer(-sample, names_to = "locus_cigar", values_to = "n_read") %>%
  separate_wider_delim(locus_cigar, ",", names = c("locus_pos", "cigar")) %>%
  # Remove strand information
  mutate(locus_pos = str_remove(locus_pos, "\\(.\\)")) %>%
  rename(rep_id = sample)

# Check for implicit missing data --------------------------------------
n_loci_by_rep <- variant_read_counts %>%
  group_by(rep_id) %>%
  summarize(n_loci = n_distinct(locus_pos), .groups = "drop")
reps_wo_all_loci <- n_loci_by_rep %>%
  filter(n_loci < n_distinct(variant_read_counts$locus_pos))
if (nrow(reps_wo_all_loci) > 0) {
  stop(
    "Not all replicates have entries for all loci in AmpSeq output", 
    call. = FALSE
  )
}

# Read and join locus names --------------------------------------------
locus_metadata <- read_tsv(
    inp$locus_metadata, 
    col_types = cols(
      .default = col_character(), 
      start_pos = col_integer(), 
      end_pos = col_integer()
    ), 
    progress = FALSE
  ) %>%
  select(locus, chrom, start_pos, end_pos) %>%
  unite(start_end, start_pos, end_pos, sep = "-") %>%
  unite(locus_pos, chrom, start_end, sep = ":")
variant_read_counts <- variant_read_counts %>%
  left_join(locus_metadata, by = "locus_pos") %>%
  select(-locus_pos)

# Read replicate metadata ----------------------------------------------
rep_metadata <- read_xlsx(
    inp$rep_metadata, 
    sheet = "pv242"
  ) %>%
  select(
    SampleID, 
    `library name`, 
    MH_DNA_Well, 
    DARC_Phenotype, 
    Site, 
    Method, 
    `Ct...21`, 
    `Ct...23`
  ) %>%
  rename(
    rep_id = SampleID, 
    lib_name = `library name`, 
    sample_id = MH_DNA_Well, 
    ct1 = `Ct...21`, 
    ct2 = `Ct...23`
  ) %>%
  pivot_longer(c(ct1, ct2), names_to = "ct_rep", values_to = "ct") %>%
  mutate(ct = na_if(ct, "No Ct")) %>%
  mutate(ct = as.numeric(ct)) %>%
  group_by(rep_id, lib_name, sample_id, DARC_Phenotype, Site, Method) %>%
  summarize(ct = mean(ct), .groups = "drop")
# Make sure each sample ID has one and only one Ct value
sample_ct_values <- rep_metadata %>%
  distinct(sample_id, ct) %>%
  group_by(sample_id) %>%
  mutate(n_ct = n()) %>%
  ungroup()
if (max(sample_ct_values$n_ct) > 1) {
  stop("Some sample IDs have more than one Ct value", call. = FALSE)
}

# Summarize read counts and join metadata ------------------------------
rl_read_counts <- variant_read_counts %>%
  # Calculate total reads for each replicate/locus
  group_by(rep_id, locus) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  left_join(rep_metadata, by = "rep_id")

# Plot reads by replicate and locus ------------------------------------
for (l in unique(rl_read_counts$lib_name)) {
  rl_read_counts %>%
    read_count_heatmap_byrep(l) %>%
    print()
}
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" /><img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-1-2.png" width="100%" /><img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-1-3.png" width="100%" />

These heatmaps show the number of reads in the final dataset for each
locus-replicate combination. The plots have been grouped according to
library. LB2 contained samples from Duffy positives, sequenced in
duplicate, while LB5 and LB7 contained samples from Duffy negatives that
were sequenced 5x.

Sequencing yield was reasonably consistent across replicates, though
there were a few exceptions to this. Sequencing yield was a lot better
with samples from Duffy positive individuals, which is presumably
related to parasitemia.

### Sample Means

``` r
read_count_heatmap_bysample <- function(reads, loi) {
  reads %>%
    filter(lib_name == loi) %>%
    ggplot(mapping = aes(x = locus, y = sample_id, fill = n_read)) +
    geom_tile() +
    scale_fill_fermenter(
      palette = "YlGnBu", 
      name = "Reads", 
      breaks = c(10, 100, 1000, 2000)
    ) +
    labs(x = "Locus", y = "Sample", title = loi) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    )
}

# Compute mean read counts for each sample/locus -----------------------
sl_mean_read_counts <- rl_read_counts %>%
  group_by(sample_id, locus, lib_name) %>%
  summarize(n_read = mean(n_read), .groups = "drop")

# Plot reads by sample and locus ---------------------------------------
for (l in unique(sl_mean_read_counts$lib_name)) {
  sl_mean_read_counts %>%
    read_count_heatmap_bysample(l) %>%
    print()
}
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" /><img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-2-2.png" width="100%" /><img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-2-3.png" width="100%" />

These heatmaps show the mean number of reads for each sample-locus
combination. Once again, a separate plot was made for each library.

# Filtering at Each Step

``` r
# Read in read filtering information -----------------------------------
reads_summary <- read_table(
    inp$reads_summary, 
    col_names = c(
      "sample", 
      "input", 
      "filtered", 
      "denoisedF", 
      "denoisedR", 
      "merged", 
      "nonchim"
    ), 
    col_types = cols(.default = "i", sample = "c"), 
    skip = 1, 
    progress = FALSE
  ) %>%
  rename(rep_id = sample)

# Join final read counts -----------------------------------------------
reads_by_step <- rl_read_counts %>%
  group_by(rep_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  right_join(reads_summary, by = "rep_id") %>%
  rename(final = n_read) %>%
  relocate(final, .after = last_col()) %>%
  pivot_longer(-rep_id, names_to = "workflow_step", values_to = "n_read") %>%
  group_by(workflow_step) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  mutate(
    workflow_step = factor(
      workflow_step, 
      c(
        "input", 
        "filtered", 
        "denoisedF", 
        "denoisedR", 
        "merged", 
        "nonchim", 
        "final"
      )
    )
  ) %>%
  arrange(workflow_step)
reads_by_step
```

    ## # A tibble: 7 × 2
    ##   workflow_step  n_read
    ##   <fct>           <int>
    ## 1 input         7118376
    ## 2 filtered      6441139
    ## 3 denoisedF     6417124
    ## 4 denoisedR     6407467
    ## 5 merged        6340656
    ## 6 nonchim       6329246
    ## 7 final         6220699

This table shows the total number of reads after each step in the DADA2
part of the pipeline plus the final filtering done by `ASV_to_CIGAR.py`.
At no point is there an order(s) of magnitude reduction in the data
volume, which would indicate cause for concern.

# Relationship with Parasitemia

This table shows mean total sample read counts and parasitemia. The two
Ct replicates have been averaged.

``` r
# Compute and print sample read count and parasitemia ------------------
sample_mean_total_read_counts <- rl_read_counts %>%
  group_by(rep_id, sample_id, lib_name, ct, DARC_Phenotype) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(sample_id, lib_name, ct, DARC_Phenotype) %>%
  summarize(n_read = mean(n_read), .groups = "drop") %>%
  filter(! is.na(ct)) %>%
  # Two samples have a suspiciously low Ct
  filter(ct > 15) %>%
  mutate(parasitemia = (10^((ct - 41.663)/-3.289)))
sample_mean_total_read_counts %>%
  arrange(parasitemia)
```

    ## # A tibble: 54 × 6
    ##    sample_id lib_name    ct DARC_Phenotype n_read parasitemia
    ##    <chr>     <chr>    <dbl> <chr>           <dbl>       <dbl>
    ##  1 MH2_C01   LB2       38.4 positive        1188.        9.97
    ##  2 MH1_C02   LB2       36.0 heterozygous    1027        51.0 
    ##  3 MH4_A11   LB2       36.0 positive        2046.       51.2 
    ##  4 MH2_D04   LB2       34.7 positive        6337       131.  
    ##  5 MH1_A04   LB2       34.7 heterozygous    1014.      135.  
    ##  6 MH2_D11   LB2       34.2 positive       68071       187.  
    ##  7 MH4_B04   LB2       33.6 positive        3241       286.  
    ##  8 MH1_B10   LB2       33.5 heterozygous    3251       303.  
    ##  9 MH4_D12   LB2       33.1 positive       16656       411.  
    ## 10 MH2_B08   LB2       32.8 positive       43558       499.  
    ## # ℹ 44 more rows

Even at parasitemias of 10-50 parasites/$\mu$L, more than 1000 reads
were obtained.

``` r
# Plot parasitemia versus sample read count ----------------------------
sample_mean_total_read_counts %>%
  ggplot(mapping = aes(x = parasitemia, y = n_read)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = expression("Parasite Density (parasites" ~ "/" ~ mu ~ "L)"), 
    y = "Mean Reads per Replicate"
  )
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

This scatterplot shows the relationship between parasitemia and mean
total read count obtained for each sample. Both axes are $log_{10}$
scaled.

The fitted linear relationship suggests that at a parasitemia of 100
parasites/$\mu$L one can expect to obtain about 10,000 reads. This is
encouraging, as this would yield more than 100 reads per locus if they
are evenly distributed across the panel.

``` r
# Save data for publication figures ------------------------------------
sl_mean_read_counts %>%
  write_csv(out$sl_read_counts)
sample_mean_total_read_counts %>%
  write_csv(out$s_read_counts_parasitemia)
```
