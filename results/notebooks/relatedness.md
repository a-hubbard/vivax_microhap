Dcifer Relatedness for Microhaplotype Data
================
Alfred Hubbard

# Relatedness Histograms

``` r
# Convert the values in an array to a tidy tibble with columns for the 
# row names, column names, and each layer in the third dimension of the 
# input array
array2tib <- function(arr, rowname_lbl, colname_lbl) {
  # Convert each layer of the 3D array to a tidy tibble
  map(1:dim(arr)[3], matrix2tib, arr, rowname_lbl, colname_lbl) %>%
    # Join all tidy tibbles
    reduce(left_join, by = c(rowname_lbl, colname_lbl))
}

# Convert a 2D matrix to a tidy tibble, with columns for the rownames, 
# column names, and values
matrix2tib <- function(stack_num, arr, rowname_lbl, colname_lbl) {
  arr[,,stack_num] %>%
    as_tibble(rownames = rowname_lbl) %>%
    pivot_longer(
      -all_of(rowname_lbl), 
      names_to = colname_lbl, 
      values_to = dimnames(arr)[[3]][stack_num]
    )
}

# Read data ------------------------------------------------------------
rel_res_arr <- read_rds(inp$rel_mat)
rel_res <- rel_res_arr %>%
  array2tib("sample_a", "sample_b") %>%
  filter(! is.na(estimate))

# Save tidy version of sample-level relatedness ------------------------
write_csv(rel_res, out$sample_rel)
```

``` r
# Plot histogram of sample-level relatedness ---------------------------
ggplot(data = rel_res, mapping = aes(x = estimate)) +
  geom_histogram(bins = 50) +
  labs(x = "Relatedness", y = "No. of Sample Pairs")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

This histogram shows the individual-level relatedness values present in
the dataset. As expected, there are a great many sample pairs with 0
relatedness. These are removed to better visualize the remaining pairs.

``` r
# Plot histogram of sample-level relatedness ---------------------------
filter(rel_res, estimate != 0) %>%
  ggplot(mapping = aes(x = estimate)) +
  geom_histogram(bins = 50) +
  labs(x = "Relatedness", y = "No. of Sample Pairs")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

This histogram shows the same thing, except pairs with a relatedness of
0 have been removed to allow closer inspection of the remaining values.
This figure shows that most of the non-zero relatedness pairs have
values between 0 and 0.25, after which point the frequency declines
fairly continuously up to 0.5. A few pairs have very high relatedness,
with values near 1.

Overall, these results make sense. In a geographically-disparate dataset
such as this one, most sample pairs should have a fairly low
relatedness. Also, while we expect most pairs are not clonally related,
this does occur in low transmission settings, so the fact that a few
pairs have relatedness at or near 1 is not surprising.

# Relatedness Heatmaps

In this section, sample-level relatedness is aggregated to the group
level, starting with MalariaGEN populations.

## Population-level

``` r
# Order sample pairs according to population pairs in a tibble
order_pairs_by_pop <- function(pair_data, label_a, label_b) {
  pair_data %>%
    # This temporary column is necessary because group_modify doesn't 
    # want the grouping variables to be passed back from the internal 
    # function
    unite(pops, all_of(label_a), all_of(label_b), remove = FALSE) %>%
    group_by(pops) %>%
    group_modify(~ order_pairs(.x, label_a, label_b)) %>%
    ungroup() %>%
    select(-pops)
}

# Order sample pairs for one population pair in a tibble
order_pairs <- function(pair_data, label_a, label_b) {
  if (pair_data[[1,label_a]] > pair_data[[1,label_b]]) {
    old_pop_a <- pair_data[[label_a]]
    old_pop_b <- pair_data[[label_b]]
    pair_data[[label_a]] <- old_pop_b
    pair_data[[label_b]] <- old_pop_a
    old_sample_a <- pair_data[["sample_a"]]
    old_sample_b <- pair_data[["sample_b"]]
    pair_data[["sample_a"]] <- old_sample_b
    pair_data[["sample_b"]] <- old_sample_a
  } 
  pair_data
}

relatedness_heatmap <- function(rel_tib, label_a, label_b, metric) {
  scale_label_key <- list(
    "mean_r" = "Mean r", 
    "frac_signif" = "Frac. Signif.", 
    "mean_signif_r" = "Mean r if Signif."
  )
  ggplot(
      data = rel_tib, 
      mapping = aes(
        x = .data[[label_a]], 
        y = .data[[label_b]], 
        fill = .data[[metric]]
      )
    ) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_distiller(palette = "YlGnBu") +
    labs(x = "Pop. A", y = "Pop. B", fill = scale_label_key[[metric]]) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

# Read and join metadata -----------------------------------------------
mh_meta <- read_csv(
    inp$mh_data, 
    col_types = cols(
      .default = col_character(), 
      Lat = col_double(), 
      Long = col_double(), 
      Year = col_integer(), 
      `% callable` = col_double(), 
      Fws = col_double(), 
      F_MISS = col_double()
    ), 
    progress = FALSE
  ) %>%
  select(sample_id, Population, geog) %>%
  distinct()
rel_res <- rel_res %>%
  left_join(mh_meta, by = c("sample_a" = "sample_id")) %>%
  rename(pop_a = Population, geog_a = geog) %>%
  left_join(mh_meta, by = c("sample_b" = "sample_id")) %>%
  rename(pop_b = Population, geog_b = geog)

# Compute population-level relatedness metrics -------------------------
low_n_pops <- c("MSEA", "WSEA")
pop_rel <- rel_res %>%
  # Remove populations with few samples
  filter((! pop_a %in% low_n_pops) & (! pop_b %in% low_n_pops)) %>%
  # Reorder columns to ensure pop_a comes before pop_b
  order_pairs_by_pop("pop_a", "pop_b") %>%
  group_by(pop_a, pop_b) %>%
  summarize(mean_r = mean(estimate), .groups = "drop")
# As a sanity check, plot heatmap before reverse combinations are 
# added. There should only be swatches in the lower triangle.
relatedness_heatmap(pop_rel, "pop_a", "pop_b", "mean_r")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

This heatmap should only have color swatches in the lower triangle. This
is the case.

``` r
# Add reverse combinations and duplicate values ------------------------
rev_combs <- filter(pop_rel, pop_a != pop_b) %>%
  mutate(
    pop_a_temp = pop_b, 
    pop_b = pop_a, 
    pop_a = pop_a_temp
  ) %>%
  select(-pop_a_temp)
pop_rel <- bind_rows(pop_rel, rev_combs)

# Convert population fields to factors to control sorting --------------
pop_levels <- c("LAM", "AF", "WAS", "WSEA", "ESEA", "MSEA", "OCE")
pop_rel <- pop_rel %>%
  mutate(pop_a = factor(pop_a, pop_levels)) %>%
  mutate(pop_b = factor(pop_b, pop_levels))

# Plot mean relatedness by population ----------------------------------
# GOAL: If more involved heatmaps are needed, rework this code to use 
# base R matrices and plotting, which has better support for this sort 
# of thing than the tidyverse. If I really need something advanced, 
# the ComplexHeatmap package might be worth a look.
relatedness_heatmap(pop_rel, "pop_a", "pop_b", "mean_r")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

This heatmap shows the mean relatedness within and among all sample
pairs corresponding to pairs of MalariaGEN populations. The key for
these populations is given below:

``` r
# Read and print key for MalariaGEN population abbreviations -----------
pop_key <- read_csv(
  inp$mg_pop_key, 
  col_types = cols(.default = col_character()), 
  progress = FALSE
)
pop_key %>%
  filter(! pop_abbrev %in% low_n_pops)
```

    ## # A tibble: 5 × 2
    ##   pop_abbrev pop_name            
    ##   <chr>      <chr>               
    ## 1 LAM        Latin America       
    ## 2 AF         Africa              
    ## 3 WAS        West Asia           
    ## 4 ESEA       East south-east Asia
    ## 5 OCE        Oceania

This analysis shows, unsurprisingly, that within-population relatedness
is much higher than between-population relatedness. It is particularly
high in Latin America and Oceania. This makes particular sense in the
case of Latin America, which is considered to be isolated compared to
other *P. vivax* populations.

The between-population relatedness patterns also conform to
expectations: relatedness is higher between Africa and West Asia and
between East south-east Asia and Oceania than between other pairs. These
pairs are also in the closest geographic proximity to each other
relative to the others.

## Geography-Level

In this section, relatedness is analyzed at the level of more localized
geographies, defined manually by merging some sample sites and removing
others to create a set of distinct geographies that each has an
acceptable sample size.

``` r
# Compute population-level relatedness metrics -------------------------
geog_rel <- rel_res %>%
  # Remove sample pairs without geography defined, which are intended 
  # to be omitted from this level of analysis
  filter((! is.na(geog_a)) & (! is.na(geog_b))) %>%
  # Reorder columns to ensure geog_a comes before geog_b
  order_pairs_by_pop("geog_a", "geog_b") %>%
  mutate(is_signif = p_value < 0.05) %>%
  group_by(geog_a, geog_b) %>%
  summarize(
    mean_r = mean(estimate), 
    frac_signif = sum(is_signif) / n(), 
    .groups = "drop"
  )

# Add reverse combinations and duplicate values ------------------------
rev_combs <- filter(geog_rel, geog_a != geog_b) %>%
  mutate(
    geog_a_temp = geog_b, 
    geog_b = geog_a, 
    geog_a = geog_a_temp
  ) %>%
  select(-geog_a_temp)
geog_rel <- bind_rows(geog_rel, rev_combs)

# Convert population fields to factors to control sorting --------------
geog_levels <- c(
  "Colombia", 
  "Jimma", 
  "Amhara", 
  "Afghanistan", 
  "Cambodia", 
  "Ho Chi Min", 
  "Binh Phuoc", 
  "Krong Pa", 
  "Indonesia"
)
geog_rel <- geog_rel %>%
  mutate(geog_a = factor(geog_a, geog_levels)) %>%
  mutate(geog_b = factor(geog_b, geog_levels))

# Plot mean relatedness by population ----------------------------------
# GOAL: If more involved heatmaps are needed, rework this code to use 
# base R matrices and plotting, which has better support for this sort 
# of thing than the tidyverse. If I really need something advanced, 
# the ComplexHeatmap package might be worth a look.
relatedness_heatmap(geog_rel, "geog_a", "geog_b", "mean_r")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

This heatmap shows the relatedness within and between geographies,
sorted on a rough east-west axis. Exactly as one would expect,
relatedness is higher within versus between populations. Indeed,
relatedness between geographies from the same populations is comparable
to relatedness within populations, which indicates a high degree of
transmission and gene flow within these populations (i.e., these results
support the population definitions created by MalariaGEN).

One can also observe higher relatedness between Afghanistan and the
locations in Ethiopia, and between Indonesia and the Southeast Asia
locations, than between other geography pairs. This makes sense, as
these locations are all from the same regions, and is also consistent
with the population-level results, above.

### Southeast Asia

One of the desired capabilities of this panel is to be able to identify
within-region population structure, so in this section the geographies
from Southeast Asia are examined in more detail.

``` r
# Make heatmap just for Southeast Asia ---------------------------------
se_asia_geogs <- c("Binh Phuoc", "Cambodia", "Ho Chi Min", "Krong Pa")
se_asia_geog_rel <- geog_rel %>%
  filter(geog_a %in% se_asia_geogs & geog_b %in% se_asia_geogs)
se_asia_geog_rel %>%
  relatedness_heatmap("geog_a", "geog_b", "mean_r")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

This heatmap shows the mean relatedness within and between the
geographies from Southeast Asia. The spread of mean relatedness is
fairly low within this region, providing limited ability to identify
patterns. For this reason, another metric, the fraction of
significantly-related sample pairs, was examined as a complement.

``` r
se_asia_geog_rel %>%
  relatedness_heatmap("geog_a", "geog_b", "frac_signif")
```

![](/users/ahubba16/projects/vivax_microhap/results/notebooks/relatedness_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

This heatmap shows the fraction of sample pairs that are significantly
related within and between the geographies from Southeast Asia. Ranging
from roughly 0.5 to 0.7, this metric provides a somewhat improved
ability to discriminate between geography pairs than mean relatedness.

The patterns evident in this figure do not conform to a simple pattern
of isolation-by-distance. Cambodia and Ho Chi Min are as closely-related
as Ho Chi Min and Krong Pa, which are much closer. In general, Ho Chi
Min is closely-related to the other locations in the region, in some
cases with higher relatedness than was found within that location. This
might be explained by the fact that Ho Chi Min is a large city with lots
of travel in and out.

However, with that hypothesis, one would expect the within-population
relatedness to be lower in Ho Chi Min than in other geographies, because
greater human movement in and out of the city should increase parasite
diversity. The opposite pattern was observed, however - Ho Chi Min has
the highest within-population relatedness of any of these locations.
