# scripts/pdl_data.R

# libraries
library(here)
library(dplyr)
library(stringr)
library(lubridate)
library(janitor)
library(vroom)

# Load Texas PDL data

tx_formulary_drug <- vroom(
  here("raw_data/PDL/TX/formulary_drug.txt"),
  delim = "|",
  trim_ws = TRUE,
  col_types = cols(.default = "c")
) |>
  clean_names() |>
  select(
    drug_ndc,
    drug_med_effdate,
    drug_med_end_date,
    drug_med_code,
    drug_pdl_pa_required,
    drug_pdl_eff_date,
    drug_clinical_pa_required,
    drug_tnp
  ) |>
  rename(
    ndc = drug_ndc,
    medicaid_coverage_start = drug_med_effdate,
    medicaid_coverage_end = drug_med_end_date,
    medicaid_covered = drug_med_code
  ) |>
  mutate(
    medicaid_coverage_start = mdy(medicaid_coverage_start),
    medicaid_coverage_end = mdy(medicaid_coverage_end))

# Develop indicators for passing PDL criterion

tx_formulary_drug <- tx_formulary_drug |>
  mutate(
    price_stage_pass = ifelse(drug_pdl_pa_required == "Yes", 0, 1),
    pt_stage_pass = ifelse(drug_clinical_pa_required == "Yes", 0, 1),
    remain_in_list = ifelse(
      !is.na(medicaid_coverage_start) & is.na(medicaid_coverage_end) |
      year(medicaid_coverage_end) >= 2025, 
      1, 0)
  )

# Develop preferred indicator

tx_formulary_drug <- tx_formulary_drug |>
  mutate(
    preferred = ifelse(
      medicaid_covered == "Yes" & price_stage_pass == 1 & pt_stage_pass == 1,
      1, 0),
    covered_not_preferred = ifelse(
      medicaid_covered == "Yes" & (price_stage_pass == 0 | pt_stage_pass == 0),
      1, 0),
    not_covered = ifelse(medicaid_covered == "No", 1, 0)
  )

# Obtain NDC codes

tx_formulary_drug <- tx_formulary_drug |>
  mutate(
    labeler_ndc = substr(ndc, 1, 5),
    product_ndc = substr(ndc, 6, 9),
    package_ndc = substr(ndc, 10, 11)
  )

# Select relevant variables

tx_formulary_drug <- tx_formulary_drug |>
  select(
    ndc,
    labeler_ndc,
    product_ndc,
    package_ndc,
    medicaid_coverage_start,
    medicaid_coverage_end,
    medicaid_covered,
    price_stage_pass,
    pt_stage_pass,
    remain_in_list
  )

# Save cleaned data

write.csv(
  tx_formulary_drug,
  "processed_data/PDL/texas_PDL_data.csv",
  row.names = FALSE
)

# Remove ndc to remove duplication of variable in merge

tx_formulary_drug <- tx_formulary_drug |>
  select(-ndc)

# Import SDUD data cleaned

SDUD_firm_mapped <- readRDS(
  file = "processed_data/state_drug_utilisation_data/SDUD_firm_mapped.rds"
)

# Merge PDL data with SDUD data

SDUD_PDL <- SDUD_firm_mapped |>
  mutate(package_code = substr(ndc, 10, 11)) |>
  left_join(
    tx_formulary_drug,
    by = c("labeler_code" = "labeler_ndc",
           "product_code" = "product_ndc",
           "package_code" = "package_ndc")
  )

# Save merged data

saveRDS(
  SDUD_PDL,
  file = "processed_data/state_drug_utilisation_data/SDUD_firm_PDL.rds"
)