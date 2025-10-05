# Tidy MalariaGEN microhaplotypes and join metadata

# Load required libraries
library(Biostrings)
library(msa)
library(seqinr)
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(magrittr)
library(optparse)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--mh_fasta_dir", 
    help = "Directory containing FASTAs with microhaplotyes"
  ), 
  make_option(
    "--sample_metadata", 
    help = "CSV file containing sample metadata"
  ), 
  make_option("--trg_coords", help = "BED file containing target coordinates"), 
  make_option(
    "--aligned_haps", 
    help = "RDS file to contain tibble with aligned haplotypes"
  ), 
  make_option("--out_csv", help = "CSV file to contain tidied data"), 
  make_option(
    "--alignments_fasta_dir", 
    help = "Directory to contain FASTAs with alignments for each locus"
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   mh_fasta_dir = "../../results/microhap/MalariaGEN/fasta", 
#   trg_coords = 
#     "../../results/primer_mapping/captured_seq_coords_good_amp.bed", 
#   sample_metadata = 
#     "../../results/MalariaGEN/preprocessed/sample_metadata.csv", 
#   aligned_haps = "../../results/microhap/MalariaGEN/mg_microhap_aligned.rds", 
#   out_csv = "../../results/microhap/MalariaGEN/mg_microhap.csv", 
#   alignments_fasta_dir = "../../results/microhap/MalariaGEN/aligned_fastas"
# )

seqlist2tib <- function(seqlist) {
  tibble(seqloc = names(seqlist), hapseq = map_chr(seqlist, pluck, 1))
}

align_seqtib <- function(seqtib) {
  seq_vec <- seqtib$hapseq
  names(seq_vec) <- seqtib$sample_id
  seq_vec %>%
    Biostrings::DNAStringSet() %>%
    # ClustalW does not produce a good alignment for pvcrt_o.10k.indel
    msa::msa(method = "Muscle")
}

alignment2tib <- function(alignment) {
  alignment_stringset <- Biostrings::unmasked(alignment)
  tibble(
    sample_id = names(alignment_stringset), 
    hapseq = as.character(alignment_stringset)
  )
}

alignedseq2fasta <- function(locus, seqs) {
  fasta <- file.path(arg$alignments_fasta_dir, str_c(locus, ".fa"))
  seqinr::write.fasta(as.list(seqs$hapseq), seqs$sample_id, fasta)
}

# Read haplotype sequences into a tibble -------------------------------
mh_data <- tibble(
    filename = list.files(arg$mh_fasta_dir, full.names = TRUE)
  ) %>%
  mutate(sample_id = str_replace(basename(filename), "_hapseq.fa", "")) %>%
  mutate(
    seqlist = map(
      filename, 
      seqinr::read.fasta, 
      seqtype = "DNA", 
      as.string = TRUE
    )
  ) %>%
  mutate(seqtib = map(seqlist, seqlist2tib)) %>%
  unnest(seqtib) %>%
  select(-filename, -seqlist)

# Read in target coordinates and use to join in locus names ------------
trg_loc_name_key <- read_tsv(
    arg$trg_coords, 
    col_types = cols(.default = col_character()), 
    col_names = c("chrom", "start_pos", "end_pos", "locus", "score", "strand"), 
    progress = FALSE
  ) %>%
  select(-score, -strand) %>%
  unite(coords, start_pos, end_pos, sep = "-") %>%
  unite(trgloc, chrom, coords, sep = ":")
mh_data <- mh_data %>%
  left_join(trg_loc_name_key, by = c("seqloc" = "trgloc")) %>%
  separate(seqloc, "chrom", sep = ":", extra = "drop")

# Align sequences for each locus ---------------------------------------
# Used for development
# mh_data <- filter(mh_data, locus == "pvcrt_o.10k.indel")
mh_data <- mh_data %>%
  nest(locus_seqs = c(sample_id, hapseq)) %>%
  mutate(locus_seqs_aligned = map(locus_seqs, align_seqtib)) %>%
  select(-locus_seqs)
# Save aligned sequences for diversity and selection analyses
write_rds(mh_data, arg$aligned_haps)

# Convert aligned sequences back to character vectors ------------------
mh_data <- mh_data %>%
  mutate(seqtib = map(locus_seqs_aligned, alignment2tib)) %>%
  select(-locus_seqs_aligned) %>%
  unnest(seqtib)
# Print whether any loci have gaps
loci_w_gaps <- mh_data %>%
  filter(str_detect(hapseq, "-")) %$%
  unique(locus)
if (length(loci_w_gaps > 0)) {
  print("The following loci have gaps in their alignments:")
  print(loci_w_gaps)
}

# Read and join sample metadata ----------------------------------------
sample_metadata <- read_csv(
    arg$sample_metadata, 
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
  )
mh_data <- mh_data %>%
  left_join(sample_metadata, by = "sample_id")

# Write to disk --------------------------------------------------------
# CSV containing tidy data
write_csv(mh_data, arg$out_csv)
# FASTA files with alignments, in case visualization is desired
mh_data %>%
  select(locus, sample_id, hapseq) %>%
  nest(locus_seqs = c(sample_id, hapseq)) %$%
  walk2(locus, locus_seqs, alignedseq2fasta)
