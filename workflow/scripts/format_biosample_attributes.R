# Reformat primer list and filter out targets with known issues

# Load required libraries
# These libraries are referenced without the :: construction, and so 
# are loaded second to avoid masking
library(dplyr)
library(magrittr)
library(optparse)
library(readr)
library(readxl)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--uci1223_rep_metadata", 
    help = "Excel file containing replicate metadata"
  ), 
  make_option(
    "--biosample_attributes", 
    help = "Path of TSV file to contain SRA biosample attributes"
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg <- list(
    uci1223_rep_metadata = 
      "../../resources/MH_template__LB_EHS_dz07242024_to_liz.xlsx", 
    biosample_attributes = "../../results/uci1223_biosample_attributes.tsv"
  )
}

# Read in replicate metadata -------------------------------------------
rep_metadata <- read_xlsx(
    arg$uci1223_rep_metadata, 
    sheet = "pv242"
  ) %>%
  select(
    SampleID, 
    `library name`, 
    MH_DNA_Well, 
    Site, 
    Method
  )

# Condense to sample level and reformat to SRA standards ---------------
rep_metadata %>%
  filter(`library name` == "LB2") %>%
  select(-`library name`) %>%
  distinct(MH_DNA_Well, Site, Method) %>%
  rename(sample_name = MH_DNA_Well) %>%
  rename(geo_loc_name = Site) %>%
  mutate(
    geo_loc_name = if_else(
      geo_loc_name == "ET_AD", 
      "Ethiopia:Arjo-Didessa sugarcane plantation", 
      "Ethiopia:Gambella rice development areas"
    )
  ) %>%
  rename(sampling_method = Method) %>%
  mutate(
    sampling_method = if_else(
      sampling_method == "PCD", 
      "Passive case detection", 
      "Mass blood survey"
    )
  ) %>%
  mutate(organism = "Plasmodium vivax") %>%
  mutate(isolate = sample_name) %>%
  mutate(collected_by = "Jimma University, Ethiopia") %>%
  mutate(collection_date = "2018-10") %>%
  mutate(
    lat_lon = if_else(
      geo_loc_name == "Ethiopia:Arjo-Didessa sugarcane plantation", 
      "8.69498 N 36.43214 E", 
      "7.87972 N 34.55132 E"
    )
  ) %>%
  mutate(host = "Homo sapiens") %>%
  mutate(sample_type = "Dried blood spot") %>%
  write_tsv(arg$biosample_attributes)
