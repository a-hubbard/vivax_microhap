# Filter and reformat gene annotations to create locus metadata file

# Load required libraries ----------------------------------------------
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--annot", help = "TSV file containing annotations"), 
  make_option("--out", help = "TSV to contain output")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$annot <- "../../results/panel_annot.tsv"
}

extract_locus_prop <- function(locus_matches) {
  # Separate out coding sequences and genes that overlap with this 
  # locus
  cds_matches <- locus_matches %>%
    filter(feature_type == "CDS")
  gene_matches <- locus_matches %>%
    filter(feature_type %in% c("gene", "pseudogene"))
  # Check for multiple gene or coding sequence matches, which would 
  # violate assumptions of logic below
  if (nrow(gene_matches) > 1) {
    warning("Multiple genes detected for one locus:", call. = FALSE)
    print(gene_matches)
  }
  if (nrow(cds_matches) > 1) {
    warning("Multiple CDSes overlap one locus:", call. = FALSE)
    print(cds_matches)
  }
  # Use CDS and gene matches, plus positional information, to classify 
  # loci as exons, introns, a mix of the two, or noncoding
  if (identical(nrow(gene_matches), 0L)) {
    seq_type <- "noncoding"
    gene_type <- NA
  } else {
    gene_type <- gene_matches[[1,"feature_type"]]
    if (nrow(cds_matches) > 0) {
      if (
        (cds_matches[[1,"locus_start_pos"]] >= cds_matches[[1,"ref_start_pos"]]) 
        && (cds_matches[[1,"locus_end_pos"]] <= cds_matches[[1,"ref_end_pos"]])
      ) {
        seq_type <- "exon"
      } else {
        seq_type <- "mix"
      }
    } else {
      seq_type <- "intron"
    }
  }
  # Use positional information to define type of overlap with the gene
  if (identical(seq_type, "noncoding")) {
    gene_overlap <- NA
  } else {
    if (
      (gene_matches[[1,"locus_start_pos"]] >= gene_matches[[1,"ref_start_pos"]]) 
      && (gene_matches[[1,"locus_end_pos"]] <= gene_matches[[1,"ref_end_pos"]])
    ) {
      gene_overlap <- "within"
    } else if (
      (gene_matches[[1,"locus_start_pos"]] < gene_matches[[1,"ref_start_pos"]]) 
      && (gene_matches[[1,"locus_end_pos"]] < gene_matches[[1,"ref_end_pos"]])
    ) {
      gene_overlap <- "left_hang"
    } else if (
      (gene_matches[[1,"locus_start_pos"]] > gene_matches[[1,"ref_start_pos"]]) 
      && (gene_matches[[1,"locus_end_pos"]] > gene_matches[[1,"ref_end_pos"]])
    ) {
      gene_overlap <- "right_hang"
    } else if (
      (gene_matches[[1,"locus_start_pos"]] <= gene_matches[[1,"ref_start_pos"]]) 
      && (gene_matches[[1,"locus_end_pos"]] >= gene_matches[[1,"ref_end_pos"]])
    ) {
      gene_overlap <- "contains"
    } else {
      stop("Gene overlap logic did not find a matching condition for :")
      print(locus_matches)
    }
  }
  # Extract gene symbol, if applicable
  if (identical(seq_type, "noncoding")) {
    gene_sym <- NA
  } else {
    gene_sym <- gene_matches[[1,"attributes"]] %>%
      str_extract("Name=PVP01_[0-9]+") %>%
      str_remove("Name=")
  }
  # Extract protein product, if applicable
  if (identical(seq_type, "noncoding")) {
    product <- NA
  } else {
    if (identical(nrow(cds_matches), 0L)) {
      warning(
        "No CDS matches for gene ", 
        gene_sym, 
        ". Cannot identify protein product", 
        call. = FALSE
      )
    }
    split_attr <- cds_matches[[1,"attributes"]] %>%
      str_split_1(";")
    product <- split_attr %>%
      str_subset("product=") %>%
      str_remove("product=")
    # For pseudogenes, need to use Note attribute
    if (identical(length(product), 0L)) {
      product <- split_attr %>%
        str_subset("Note=") %>%
        str_remove("Note=")
    }
    product <- product %>%
      str_replace("%2C", ",")
  }
  # Assemble output tibble
  locus_matches %>%
    select(-ref_start_pos, -ref_end_pos, -feature_type, -attributes) %>%
    rename(
      start_pos = locus_start_pos, 
      end_pos = locus_end_pos, 
      strand = locus_strand
    ) %>%
    distinct() %>%
    mutate(
      gene_type = gene_type, 
      seq_type = seq_type, 
      gene_overlap = gene_overlap, 
      gene_sym = gene_sym, 
      product = product
    )
}

# Read in annotations --------------------------------------------------
annot <- read_tsv(
    arg$annot, 
    col_names = c(
      "locus_chrom", 
      "locus_start_pos", 
      "locus_end_pos", 
      "locus", 
      "locus_score", 
      "locus_strand", 
      "ref_chrom", 
      "feature_source", 
      "feature_type", 
      "ref_start_pos", 
      "ref_end_pos", 
      "ref_score", 
      "ref_strand", 
      "phase", 
      "attributes"
    ), 
    col_types = cols(
      .default = col_character(), 
      locus_start_pos = col_integer(), 
      locus_end_pos = col_integer(), 
      ref_start_pos = col_integer(), 
      ref_end_pos = col_integer()
    ), 
    progress = FALSE
  ) %>%
  select(
    locus_chrom, 
    locus_start_pos, 
    locus_end_pos, 
    locus_strand, 
    locus, 
    feature_type, 
    ref_start_pos, 
    ref_end_pos, 
    attributes
  ) %>%
  rename(chrom = locus_chrom)

# Extract locus properties and write to output -------------------------
annot %>%
  nest(locus_matches = ! c(chrom, locus)) %>%
  mutate(locus_matches = map(locus_matches, extract_locus_prop)) %>%
  unnest(locus_matches) %>%
  relocate(locus, .before = 1) %>%
  write_tsv(arg$out)
