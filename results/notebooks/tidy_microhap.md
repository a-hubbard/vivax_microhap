Tidying of Microhaplotype Results
================
Alfred Hubbard

# Final Data Quantity

## Unfiltered

``` r
read_count_heatmap <- function(reads) {
  ggplot(
      data = reads, 
      mapping = aes(x = target, y = sample_id, fill = n_read)
    ) +
    geom_tile() +
    scale_fill_fermenter(
      palette = "YlGnBu", 
      name = "Reads", 
      breaks = c(10, 100, 1000, 2000)
    ) +
    labs(x = "Marker", y = "Sample") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1), 
      legend.position = "bottom", 
      legend.key.width = unit(0.35, "in")
    )
}

# Read in final, filtered data -----------------------------------------
seqtab <- read_tsv(
    inp$seqtab, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  pivot_longer(-sample, names_to = "target_cigar", values_to = "n_read") %>%
  separate_wider_delim(target_cigar, ",", names = c("target", "cigar")) %>%
  rename(sample_id = sample)

# Read and join ASV IDs ------------------------------------------------
asv2cigar <- read_tsv(
  inp$asv2cigar, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)
mh_data <- seqtab %>%
  left_join(
    asv2cigar, 
    by = c("target" = "Amplicon", "cigar" = "CIGAR"), 
    relationship = "many-to-many"
  )

# Read in sample metadata ----------------------------------------------
vera_lane_summary <- read_csv(
    inp$vera_lane_summary, 
    col_types = cols(.default = col_character()), 
    progress = FALSE
  ) %>%
  filter(str_detect(Sample, "OLD")) %>%
  separate(`Barcode sequence`, c("i7_seq", "i5_seq"), sep = "\\+")
vera_library <- read_csv(
    inp$vera_library, 
    col_names = c(
      "tube", 
      "lib_name", 
      "well", 
      "i7_name", 
      "i7_seq", 
      "i5_name", 
      "i5_seq"
    ), 
    col_types = cols(.default = col_character()), 
    skip = 1, 
    progress = FALSE
  ) %>%
  filter(str_detect(lib_name, "old")) %>%
  select(well, i7_seq, i5_seq)
pv_samples_idaho <- read_csv(
    inp$pv_samples_idaho, 
    col_types = cols(
      .default = col_double(), 
      SeqPlate_Well = col_character(), 
      `Sample Name` = col_character(), 
      Dilution = col_character()
    ), 
    progress = FALSE
  ) %>%
  mutate(SeqPlate_Well = str_sub(SeqPlate_Well, start = 4L)) %>%
  rename(UCI_SID = `Sample Name`)
gbpcd_meta <- read_csv(
  inp$gbpcd_meta, 
  col_types = cols(.default = col_character(), Age = col_integer()), 
  progress = FALSE
)

# Join sample metadata into one tibble, then to microhaplotypes --------
sample_metadata <- left_join(
    vera_lane_summary, 
    vera_library, 
    by = c("i7_seq", "i5_seq")
  ) %>%
  select(-i7_seq, -i5_seq) %>%
  left_join(pv_samples_idaho, by = c("well" = "SeqPlate_Well")) %>%
  left_join(gbpcd_meta, by = "UCI_SID")
mh_data <- mh_data %>%
  left_join(sample_metadata, by = c("sample_id" = "Sample")) %>%
  mutate(pop = str_split_i(DBS_ID, "-", 1))

# Plot reads by sample and target --------------------------------------
undiluted_reads <- mh_data %>%
  filter(Dilution == "Neat") %>%
  group_by(sample_id, target) %>%
  summarize(n_read = sum(n_read), .groups = "drop")
undiluted_reads %>%
  # Filter out controls
  filter(
    ! sample_id %in% c("PV_SWGA_OLD42", "PV_SWGA_OLD48", "PV_SWGA_OLD66")
  ) %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This heatmap shows the number of reads in the final dataset for each
target-sample combination. Controls and dilutions have been removed. The
results are quite similar to those found by AmpSeQC (analyzed in
`read_counts.Rmd`), which is expected (note that some markers present in
the AmpSeQC output have been filtered out here).

## Negative Controls

``` r
# Plot reads for negative controls -------------------------------------
mh_data %>%
  filter(sample_id %in% c("PV_SWGA_OLD42", "PV_SWGA_OLD48")) %>%
  group_by(sample_id, target) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-2-1.png" width="100%" />

This heatmap shows the number of reads in the final dataset for the two
negative controls. Clearly, there has been quite a bit of contamination,
albeit at fairly low read counts.

This table shows the genotypes that were found in the negative controls
at a relatively high read count:

``` r
# Print instances of substantial contamination -------------------------
neg_controls <- mh_data %>%
  filter(sample_id %in% c("PV_SWGA_OLD42", "PV_SWGA_OLD48")) %>%
  select(sample_id, target, ASV, n_read, cigar) %>%
  relocate(ASV, .before = cigar) %>%
  relocate(n_read, .before = cigar) %>%
  filter(n_read > 0)
neg_controls %>%
  filter(n_read > 100)
```

    ## # A tibble: 14 × 5
    ##    sample_id     target                      ASV    n_read cigar    
    ##    <chr>         <chr>                       <chr>   <int> <chr>    
    ##  1 PV_SWGA_OLD42 PvP01_01_v1_654301_654500   ASV1      238 63A      
    ##  2 PV_SWGA_OLD42 pvcrt_o.10k.indel           ASV2      170 85C      
    ##  3 PV_SWGA_OLD42 PvP01_13_v1_1116901_1117100 ASV3      165 20C      
    ##  4 PV_SWGA_OLD42 PvP01_07_v1_784501_784700   ASV4      142 .        
    ##  5 PV_SWGA_OLD42 pvk12.596KR                 ASV6      115 .        
    ##  6 PV_SWGA_OLD48 PvP01_01_v1_654301_654500   ASV1      263 63A      
    ##  7 PV_SWGA_OLD48 pvcrt_o.10k.indel           ASV2      204 85C      
    ##  8 PV_SWGA_OLD48 PvP01_13_v1_1116901_1117100 ASV3      161 20C      
    ##  9 PV_SWGA_OLD48 PvP01_07_v1_784501_784700   ASV4      144 .        
    ## 10 PV_SWGA_OLD48 PvP01_12_v1_769701_769900   ASV5      121 .        
    ## 11 PV_SWGA_OLD48 pvk12.596KR                 ASV6      108 .        
    ## 12 PV_SWGA_OLD48 PvP01_05_v1_1352001_1352200 ASV7      112 30C69C99G
    ## 13 PV_SWGA_OLD48 PvP01_05_v1_1352001_1352200 ASV177    112 30C69C99G
    ## 14 PV_SWGA_OLD48 PvP01_12_v1_2032801_2033000 ASV13     111 6G56G113G

There are not too many of these, and they do not exceed 250 reads. This
is not a trivial amount of sequence, but it also isn’t enormous.

The following two tables show 1) instances where a genotype was found in
one negative control and not the other and 2) instances where a genotype
was found in both negative controls but there is a substantial
difference in the read count.

``` r
# Assess similarity of two negative controls ---------------------------
neg_control_diff <- neg_controls %>%
  pivot_wider(names_from = sample_id, values_from = n_read) %>%
  mutate(control_diff = PV_SWGA_OLD48 - PV_SWGA_OLD42) %>%
  relocate(cigar, .after = last_col())
neg_control_diff %>%
  filter(is.na(control_diff)) %>%
  select(-control_diff)
```

    ## # A tibble: 16 × 5
    ##    target                      ASV    PV_SWGA_OLD42 PV_SWGA_OLD48 cigar         
    ##    <chr>                       <chr>          <int>         <int> <chr>         
    ##  1 PvP01_06_v1_197901_198100   ASV82              9            NA .             
    ##  2 PvP01_14_v1_1887501_1887700 ASV88              6            NA 126G132A      
    ##  3 PvP01_11_v1_2008301_2008500 ASV93              4            NA .             
    ##  4 PvP01_07_v1_500001_500200   ASV94              4            NA .             
    ##  5 PvP01_07_v1_37801_38000     ASV107             5            NA 76A86C135A    
    ##  6 PvP01_01_v1_746501_746700   ASV127             3            NA 160A          
    ##  7 PvP01_12_v1_1861301_1861500 ASV35             NA            33 .             
    ##  8 PvP01_12_v1_1861301_1861500 ASV37             NA            33 .             
    ##  9 PvP01_11_v1_1068401_1068600 ASV51             NA            25 58A60G174T177…
    ## 10 PvP01_11_v1_1068401_1068600 ASV160            NA            25 58A60G174T177…
    ## 11 PvP01_04_v1_702801_703000   ASV69             NA            10 7T18A46T67G74…
    ## 12 PvP01_07_v1_1417201_1417400 ASV108            NA             3 28G           
    ## 13 PvP01_01_v1_395301_395500   ASV114            NA             5 142C157C      
    ## 14 PvP01_11_v1_2008301_2008500 ASV140            NA             5 89G           
    ## 15 PvP01_07_v1_1218801_1219000 ASV144            NA             7 55A138G179A   
    ## 16 PvP01_10_v1_136201_136400   ASV202            NA             4 83T109A134T17…

``` r
neg_control_diff %>%
  filter(abs(control_diff) > 10)
```

    ## # A tibble: 12 × 6
    ##    target                   ASV   PV_SWGA_OLD42 PV_SWGA_OLD48 control_diff cigar
    ##    <chr>                    <chr>         <int>         <int>        <int> <chr>
    ##  1 PvP01_01_v1_654301_6545… ASV1            238           263           25 63A  
    ##  2 pvcrt_o.10k.indel        ASV2            170           204           34 85C  
    ##  3 PvP01_12_v1_769701_7699… ASV5             97           121           24 .    
    ##  4 PvP01_05_v1_1352001_135… ASV7             68           112           44 30C6…
    ##  5 PvP01_05_v1_1352001_135… ASV1…            68           112           44 30C6…
    ##  6 PvP01_06_v1_596001_5962… ASV8             75            87           12 .    
    ##  7 PvP01_03_v1_178801_1790… ASV9             78            63          -15 23C6…
    ##  8 PvP01_08_v1_1191801_119… ASV10            50            64           14 37C  
    ##  9 PvP01_14_v1_2861901_286… ASV12            70            86           16 32G1…
    ## 10 PvP01_12_v1_2032801_203… ASV13            78           111           33 6G56…
    ## 11 PvP01_08_v1_1473101_147… ASV21            41            53           12 .    
    ## 12 PvP01_04_v1_723401_7236… ASV41             8            20           12 23T4…

Taken together, these tables suggest concurrence between the negative
controls, with the possible exception of the four haplotypes (from two
loci) that show up in PV_SWGA_OLD48 at a read count greater than 10 but
do not appear in PV_SWGA_OLD42. This may warrant further investigation.

## Filtered for Contamination

To address this contamination, genotypes that match genotypes found in
the negative controls are removed from the dataset unless the read count
exceeds that found in the negative controls by at least a factor of two.
Note that no attempt to adjust the read counts of genotypes that match
the contaminating genotypes but exceed this threshold was made - it is
not clear how this could be done, mathematically speaking.

``` r
# Filter genotype read counts according to contamination in controls ---
neg_control_filterthres <- neg_controls %>%
  # Take maximum number of reads between two negative controls 
  group_by(target, cigar, ASV) %>%
  summarize(n_read_negcontrol = max(n_read), .groups = "drop_last")
mh_data_filteredbycontrol <- mh_data %>%
  # Filter out negative controls
  filter(! sample_id %in% c("PV_SWGA_OLD42", "PV_SWGA_OLD48")) %>%
  left_join(neg_control_filterthres, by = c("target", "cigar", "ASV")) %>%
  replace_na(list(n_read_negcontrol = 0)) %>%
  # Set read count to 0 for a genotype if it does not exceed the read 
  # count for any contamination for that genotype by a factor of 2
  mutate(n_read = if_else(n_read > n_read_negcontrol * 2, n_read, 0)) %>%
  select(-n_read_negcontrol)

# Plot filtered data ---------------------------------------------------
mh_data_filteredbycontrol %>%
  filter(Dilution == "Neat") %>%
  group_by(target, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-5-1.png" width="100%" />

This heatmap shows the read count for each marker/sample combination
after this filter was applied. By comparing with the read count heatmap
before this filter was applied, it can be seen that there were low
levels of data for some markers that seemed to be present in almost all
samples prior to the filter, and that this apparent data has been
removed.

This pattern, along with inspection of the contaminating genotypes in
the tables above, is consistent with the idea that one sample
contaminated all of the wells at a relatively low concentration. The
filtered heatmap suggests that about 20 samples still have a substantial
amount of data after the effects of this contamination were filtered
out.

## Positive Control

``` r
# Plot reads for positive control --------------------------------------
pos_control <- mh_data_filteredbycontrol %>%
  filter(sample_id == "PV_SWGA_OLD66") %>%
  select(sample_id, target, ASV, n_read, cigar)
pos_control %>%
  filter(n_read > 0)
```

    ## # A tibble: 4 × 5
    ##   sample_id     target                      ASV    n_read cigar   
    ##   <chr>         <chr>                       <chr>   <dbl> <chr>   
    ## 1 PV_SWGA_OLD66 PvP01_02_v1_80401_80600     ASV131      3 2G21T   
    ## 2 PV_SWGA_OLD66 PvP01_14_v1_2759601_2759800 ASV132      4 55A157G 
    ## 3 PV_SWGA_OLD66 PvP01_10_v1_136201_136400   ASV147      2 109A172C
    ## 4 PV_SWGA_OLD66 PvP01_06_v1_197901_198100   ASV149     10 190C

``` r
# pos_control %>%
#   group_by(target, sample_id) %>%
#   summarize(n_read = sum(n_read), .groups = "drop") %>%
#   read_count_heatmap()
```

After filtering out contaminated reads, there is almost nothing left in
the positive control.

``` r
# Check to make sure it is monoclonal ----------------------------------
# pos_control %>%
#   group_by(sample_id, target) %>%
#   filter(n() > 2)

# Filter out positive control ------------------------------------------
mh_data_filteredbycontrol <- mh_data_filteredbycontrol %>%
  filter(sample_id != "PV_SWGA_OLD66")

# Separate and save dilutions ------------------------------------------
diluted_samples <- mh_data_filteredbycontrol %>%
  filter(Dilution != "Neat") %$%
  unique(UCI_SID)
mh_data_filteredbycontrol %>%
  filter(UCI_SID %in% diluted_samples) %>%
  # Remove explicit missing data
  filter(n_read > 0) %>%
  select(-Age, -Sex, -pop) %>%
  write_csv(out$microhap_dilutions)

# Filter out dilutions for remainder of analysis -----------------------
mh_data_undiluted <- mh_data_filteredbycontrol %>%
  filter(Dilution == "Neat") %>%
  select(-Dilution)
```

# Missing Data by Sample

``` r
# Compute number of loci with data for each sample ---------------------
n_trg_wdata <- mh_data_undiluted %>%
  group_by(target, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(sample_id) %>%
  summarize(n_trg_wdata = sum(n_read > 0), .groups = "drop")
mh_data_undiluted <- mh_data_undiluted %>%
  left_join(n_trg_wdata, by = "sample_id")

# Plot distribution of missing data proportion by marker ---------------
n_trgs_filtered <- mh_data_undiluted %>%
  filter(n_read > 0) %$%
  n_distinct(target)
n_trg_wdata_thres <- 20
mh_data_undiluted %>%
  select(sample_id, n_trg_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_trg_wdata)) +
  geom_histogram(bins = 20) +
  geom_vline(xintercept = n_trg_wdata_thres, color = "blue") +
  labs(x = "No. of Targets With Data", y = "No. of Samples")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-8-1.png" width="100%" />

This histogram shows the distribution of samples according to how many
markers have data for each sample. There are 78 markers in the dataset.

The vertical blue line shows the chosen threshold of targets with data
that will be used for filtering the samples.

``` r
# Filter out samples with lots of missing data -------------------------
mh_data_filteredbysamp <- mh_data_undiluted %>%
  filter(n_trg_wdata >= n_trg_wdata_thres) %>%
  select(-n_trg_wdata)
```

# Missing Data by Marker

``` r
# Compute number of samples with data for each loci --------------------
n_samp_wdata <- mh_data_filteredbysamp %>%
  group_by(target, sample_id) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  group_by(target) %>%
  summarize(n_samp_wdata = sum(n_read > 0), .groups = "drop")
mh_data_filteredbysamp <- mh_data_filteredbysamp %>%
  left_join(n_samp_wdata, by = "target")

# Plot data quantity by marker -----------------------------------------
n_samp_wdata_thres <- 10
mh_data_filteredbysamp %>%
  select(target, n_samp_wdata) %>%
  distinct() %>%
  ggplot(mapping = aes(x = n_samp_wdata)) +
  geom_histogram(bins = 20) +
  geom_vline(xintercept = n_samp_wdata_thres, color = "blue") +
  labs(x = "No. of Samples With Data", y = "No. of Markers")
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/tidy_microhap_files/figure-gfm/unnamed-chunk-10-1.png" width="100%" />

This histogram shows the distribution of markers according to how many
samples have data for each marker. Note that the threshold for number of
targets with data chosen above has already been applied. After this
filter, 20 samples remain in the dataset.

The vertical blue line marks 10 samples with data. This threshold is
used to filter markers for subsequent analyses, removing those that fall
below.

``` r
# Filter out loci with lots of missing data ----------------------------
mh_data_filteredbyloci <- mh_data_filteredbysamp %>%
  filter(n_samp_wdata >= n_samp_wdata_thres) %>%
  select(-n_samp_wdata) %>%
  # Remove explicit missing data
  filter(n_read > 0)
```

# Stats

Number of markers in the filtered dataset:

``` r
n_distinct(mh_data_filteredbyloci$target)
```

    ## [1] 63

Number of samples in the filtered dataset:

``` r
n_distinct(mh_data_filteredbyloci$sample_id)
```

    ## [1] 20

``` r
# Read in and join chromosome information ------------------------------
chrom_info <- read_tsv(
    inp$trg_coords, 
    col_names = c("chrom", "start_pos", "end_pos", "target"), 
    col_types = cols(
      .default = col_character(), 
      start_pos = col_integer(), 
      end_pos = col_integer()
    ), 
    progress = FALSE
  ) %>%
  select(chrom, target)
mh_data_filteredbyloci <- mh_data_filteredbyloci %>%
  left_join(chrom_info, by = "target")

# Save filtered data ---------------------------------------------------
mh_data_filteredbyloci %>%
  # Calculate MOI. It would be more logical to do this in 
  # diversity.Rmd, but that would complicate the pipeline.
  group_by(sample_id, target) %>%
  mutate(n_hap = n()) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(moi = max(n_hap)) %>%
  ungroup() %>%
  write_csv(out$microhap_undiluted)
```
