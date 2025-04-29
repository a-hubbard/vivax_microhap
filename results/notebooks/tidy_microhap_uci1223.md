Tidying of Microhaplotype Data from UCI 12/23 Run
================
Alfred Hubbard

# Missing Data By Sample

For each sample, the replicate with the most loci that have 100 reads or
more is kept. At the same time, any samples where no replicates have 100
reads or more are removed. Duffy negative samples are also removed.

``` r
# Read in final, filtered data -----------------------------------------
allele_table <- read_tsv(
    inp$cigar_table, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  pivot_longer(-sample, names_to = "locus_cigar", values_to = "n_read") %>%
  separate_wider_delim(locus_cigar, ",", names = c("locus_pos", "cigar")) %>%
  rename(rep_id = sample) %>%
  # Remove explicit missing data
  filter(n_read > 0)

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
allele_table <- allele_table %>%
  left_join(locus_metadata, by = "locus_pos") %>%
  select(-locus_pos)

# Read and join replicate metadata -------------------------------------
rep_metadata <- read_xlsx(
    inp$rep_metadata, 
    sheet = "pv242"
  ) %>%
  select(SampleID, MH_DNA_Well, DARC_Phenotype) %>%
  rename(rep_id = SampleID, sample_id = MH_DNA_Well)
```

    ## New names:
    ## • `Well` -> `Well...20`
    ## • `Ct` -> `Ct...21`
    ## • `Well` -> `Well...22`
    ## • `Ct` -> `Ct...23`
    ## • `` -> `...25`
    ## • `` -> `...27`

``` r
allele_table <- allele_table %>%
  left_join(rep_metadata, by = "rep_id")

# Filter out non-diversity loci ----------------------------------------
allele_table <- allele_table %>%
  filter(str_detect(locus, "PvP01"))

# Filter out Duffy negatives -------------------------------------------
allele_table <- allele_table %>%
  filter(DARC_Phenotype != "negative") %>%
  select(-DARC_Phenotype)

# Get number of samples before additional filtering --------------------
n_Duffy_pos <- allele_table %$%
  n_distinct(sample_id)

# For each sample, keep replicate with most loci with >= 100 reads -----
reps_w_highreads <- allele_table %>%
  group_by(rep_id, locus, sample_id) %>%
  summarize(rl_total_reads = sum(n_read), .groups = "drop") %>%
  filter(rl_total_reads >= 100) %>%
  group_by(rep_id, sample_id) %>%
  summarize(n_loc_highreads = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  filter(n_loc_highreads == max(n_loc_highreads)) %>%
  ungroup()
```

    ## Warning: There was 1 warning in `filter()`.
    ## ℹ In argument: `n_loc_highreads == max(n_loc_highreads)`.
    ## Caused by warning in `max()`:
    ## ! no non-missing arguments to max; returning -Inf

``` r
allele_table <- allele_table %>%
  filter(rep_id %in% reps_w_highreads$rep_id) %>%
  select(-rep_id)

# Compute number of loci with data for each sample ---------------------
n_loc_wdata <- allele_table %>%
  group_by(locus, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(sample_id) %>%
  summarize(n_loc_wdata = sum(n_read > 0), .groups = "drop")
allele_table <- allele_table %>%
  left_join(n_loc_wdata, by = "sample_id")

# Plot distribution of missing data proportion by locus ----------------
n_loc_wdata_thres <- 60
allele_table %>%
  select(sample_id, n_loc_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_loc_wdata)) +
  geom_histogram(bins = 20) +
  geom_vline(xintercept = n_loc_wdata_thres, color = "blue") +
  labs(x = "No. of Loci With Data", y = "No. of Samples")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_uci1223_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

This histogram shows the distribution of samples according to how many
loci have data for each sample. The vertical blue line shows the
threshold of 60 loci with data that will be used for filtering the
samples.

``` r
# Filter out samples with lots of missing data -------------------------
allele_table <- allele_table %>%
  filter(n_loc_wdata >= n_loc_wdata_thres) %>%
  select(-n_loc_wdata)
```

# Missing Data by Locus

``` r
# Compute number of samples with data for each loci --------------------
n_samp_wdata <- allele_table %>%
  group_by(locus, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(locus) %>%
  summarize(n_samp_wdata = sum(n_read > 0), .groups = "drop")
allele_table <- allele_table %>%
  left_join(n_samp_wdata, by = "locus")

# Plot data quantity by locus ------------------------------------------
n_samp_wdata_thres <- 20
allele_table %>%
  select(locus, n_samp_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_samp_wdata)) +
  geom_histogram(bins = 25) +
  geom_vline(xintercept = n_samp_wdata_thres, color = "blue") +
  labs(x = "No. of Samples With Data", y = "No. of Loci")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_uci1223_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

This histogram shows the distribution of loci according to how many
samples have data for each locus. Note the sample filter of 60 loci with
data has already been applied, and this histogram is based on the
remaining data.

The vertical blue line marks 20 samples with data. This threshold is
used to filter loci for subsequent analyses, removing those that fall
below.

``` r
# Filter out loci with lots of missing data ----------------------------
allele_table <- allele_table %>%
  filter(n_samp_wdata >= n_samp_wdata_thres) %>%
  select(-n_samp_wdata)
```

# Stats

Number of population genetics loci in the filtered dataset:

``` r
allele_table %$%
  n_distinct(locus)
```

    ## [1] 0

Number of samples in the filtered dataset:

``` r
n_distinct(allele_table$sample_id)
```

    ## [1] 0

…out of a total of 0 Duffy positive samples.

``` r
# Save filtered data
allele_table %>%
  # Rename to column names expected by PGEcore tools
  rename(
    specimen_id = sample_id, 
    target_id = locus, 
    # This is not a sequence, obviously, but that won't matter for the 
    # tools we are using
    seq = cigar, 
    read_count = n_read
  ) %>%
  write_tsv(out$allele_table)
```
