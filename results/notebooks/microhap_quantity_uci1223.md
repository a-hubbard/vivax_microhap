Quantification of Microhaplotypes Identified by AmpSeq
================
Alfred Hubbard

# Final Data Quantity

## Reads

``` r
read_count_heatmap <- function(reads) {
  ggplot(data = reads, mapping = aes(x = locus, y = sample, fill = n_read)) +
    geom_tile() +
    scale_fill_fermenter(
      palette = "YlGnBu", 
      name = "Reads", 
      breaks = c(10, 100, 1000, 2000)
    ) +
    labs(x = "Locus", y = "Sample") +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
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
  pivot_longer(-sample, names_to = "locus_cigar", values_to = "n_read") %>%
  separate_wider_delim(locus_cigar, ",", names = c("locus", "cigar"))

# Plot reads by sample and locus ---------------------------------------
final_reads <- final_data %>%
  group_by(sample, locus) %>%
  summarize(n_read = sum(n_read), .groups = "drop")
# Save for publication figures
final_reads %>%
  write_csv(out$read_counts)
final_reads %>%
  # filter(! str_detect(sample, "Control")) %>%
  read_count_heatmap()
```

<img src="/users/ahubba16/projects/vivax_microhap/results/notebooks/microhap_quantity_uci1223_files/figure-gfm/unnamed-chunk-1-1.png" width="100%" />

This heatmap shows the number of reads in the final dataset for each
locus-sample combination. Replicates and controls have been removed.
Some loci have a decent amount of data, but most lack data for most
samples. Also, many samples have essentially no data.

This broadly mirrors the patterns seen in `read_counts.Rmd`, which shows
the same information but with an alignment based method and less
filtering. `read_counts.Rmd` shows somewhat more data in many cases,
which makes sense because, again, there was overall less filtering
applied in that analysis.

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
    ## 1 input         7118376
    ## 2 filtered      6441139
    ## 3 denoisedF     6417124
    ## 4 denoisedR     6407467
    ## 5 merged        6340656
    ## 6 nonchim       6329246
    ## 7 final         5438391

This table shows the total number of reads after each step in the DADA2
part of the pipeline plus the final filtering done by `ASV_to_CIGAR.py`.
