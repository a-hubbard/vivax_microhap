# Box plots showing mean marker heterozygosity across populations

# Load required libraries ----------------------------------------------
library(adegenet)
library(poppr)
library(rstatix)
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(patchwork)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--mh_data_gd", 
    help = "RDS file containing genind with microhaplotype data"
  ), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   mh_data_gd = 
#     "../../results/microhap/MalariaGEN/mg_microhap_filtered4popgen_gd.rds"
# )

mean_het_boxplot <- function(gendata, pop_level) {

  # Instantiate tibble with populations of interest --------------------
  if (pop_level == "Country") {
    adegenet::setPop(gendata) <- ~Country
    mean_het_by_pop <- poppr::poppr(gendata, quiet = TRUE) %>%
      as_tibble() %>%
      filter(Pop != "Total") %>%
      filter(N >= 15) %>%
      select(Pop)
  } else if (pop_level == "Site") {
    sites_of_interest <- c(
      "Oddar Meanchey", 
      "Ho Chi Min", 
      "Binh Phuoc", 
      "Dak O", 
      "Krong Pa"
    )
    adegenet::setPop(gendata) <- ~Site
    mean_het_by_pop <- poppr::poppr(gendata, quiet = TRUE) %>%
      as_tibble() %>%
      filter(Pop %in% sites_of_interest) %>%
      select(Pop)
  }

  # Compute mean marker heterozygosity ---------------------------------
  mean_het_by_pop <- mean_het_by_pop %>%
    mutate(
      locus_stats = map(
        Pop, ~poppr::locus_table(gendata, pop = .x, information = FALSE)
      )
    ) %>%
    mutate(locus_stats = map(locus_stats, as_tibble)) %>%
    unnest(locus_stats) %>%
    select(-allele, -`1-D`, -Evenness) %>%
    mutate(Hexp = as.double(Hexp))

  # Perform t-tests between population pairs ---------------------------
  pw_ttests <- rstatix::pairwise_t_test(
      mean_het_by_pop, 
      Hexp ~ Pop, 
      p.adjust.method = "bonferroni"
    ) %>%
    filter(p.adj.signif != "ns")
  if (nrow(pw_ttests > 0)) {
    print("Significant differences between populations were found!")
  }

  # Make box plot ------------------------------------------------------
  if (pop_level == "Country") {
    mean_het_by_pop <- mean_het_by_pop %>%
      mutate(
        Pop = factor(
          Pop, 
          c(
            "Colombia", 
            "Ethiopia", 
            "Afghanistan", 
            "Cambodia", 
            "Vietnam", 
            "Indonesia"
          )
        )
      )
    plot_labs <- labs(x = "Country", y = "Mean Expected Heterozygosity")
  } else if (pop_level == "Site") {
    plot_labs <- labs(x = "Site", y = "Mean Expected Heterozygosity")
    mean_het_by_pop <- mean_het_by_pop %>%
      mutate(Pop = factor(Pop, sites_of_interest))
  }
  mean_het_by_pop %>%
    ggplot(mapping = aes(x = Pop, y = Hexp)) +
    geom_boxplot() +
    plot_labs +
    ylim(0, 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

}

# Read in data ---------------------------------------------------------
mh_data_gd <- read_rds(arg$mh_data_gd) %>%
  # Remove loci with only one allele
  poppr::informloci()

# Make figure and save -------------------------------------------------
country_plot <- mean_het_boxplot(mh_data_gd, "Country")
site_plot <- mean_het_boxplot(mh_data_gd, "Site")
fig <- (country_plot | site_plot) +
  plot_annotation(tag_levels = "A")
w <- 8
h <- 5
ggsave(
  str_c(arg$out_base, ".pdf"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".png"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
