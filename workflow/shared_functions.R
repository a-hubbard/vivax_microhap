# Convert allele frequency from Dcifer format to wide format
af_dcifer2wide <- function(af_long) {
  af_wide <- af_long %>%
    group_by(target_id) %>%
    mutate(allele_id = str_c("Allele_", row_number())) %>%
    ungroup() %>%
    select(-seq) %>%
    pivot_wider(names_from = "allele_id", values_from = "freq") %>%
    mutate(across(-target_id, ~replace_na(.x, 0)))

  af_wide
}