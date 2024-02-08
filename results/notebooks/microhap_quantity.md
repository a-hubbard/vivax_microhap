Quantification of Microhaplotypes Identified by AmpSeq
================
Alfred Hubbard

# Final Data Quantity

``` r
read_count_heatmap <- function(reads) {
  ggplot(data = reads, mapping = aes(x = target, y = sample, fill = n_read)) +
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
final_data <- read_tsv(
    inp$cigar_table, 
    col_types = cols(.default = col_integer(), sample = col_character()), 
    progress = FALSE
  ) %>%
  pivot_longer(-sample, names_to = "target_cigar", values_to = "n_read") %>%
  separate_wider_delim(target_cigar, ",", names = c("target", "cigar"))

# Plot reads by sample and target --------------------------------------
final_reads <- final_data %>%
  group_by(sample, target) %>%
  summarize(n_read = sum(n_read), .groups = "drop")
final_reads %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This heatmap shows the number of reads in the final dataset for each
target-sample combination. The results are quite similar to those found
by AmpSeQC (analyzed in `read_counts.Rmd`), which is expected (note that
some markers present in the AmpSeQC output have been filtered out here.

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
)

# Join final read counts -----------------------------------------------
reads_by_step <- final_reads %>%
  group_by(sample) %>%
  summarize(n_read = sum(n_read), .groups = "drop") %>%
  right_join(reads_summary, by = "sample") %>%
  rename(final = n_read) %>%
  relocate(final, .after = last_col()) %>%
  pivot_longer(-sample, names_to = "workflow_step", values_to = "n_read") %>%
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
    ## 1 input         3067656
    ## 2 filtered      3046526
    ## 3 denoisedF     3022110
    ## 4 denoisedR     3027157
    ## 5 merged        3000987
    ## 6 nonchim       2994741
    ## 7 final         2763299

This table shows the total number of reads after each step in the DADA2
part of the pipeline plus the final filtering done by `ASV_to_CIGAR.py`.
