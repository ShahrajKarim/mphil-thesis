# exploration/pdl_overlap_analysis.r

# libraries
library(here)
library(dplyr)
library(stringr)
library(lubridate)
library(janitor)
library(vroom)
library(readr)
library(readxl)
library(knitr)
library(tibble)

# Load Texas PDL data

texas_pdl <- read_csv("processed_data/PDL/texas_PDL_data.csv") |>
    filter(price_stage_pass == 1 & pt_stage_pass == 1) |>
    select(labeler_ndc, product_ndc) |>
    mutate(texas_row = 1) |>
    distinct(labeler_ndc, product_ndc, .keep_all = TRUE)

# Load New York PDL data

ny_pdl <- read_csv("raw_data/PDL/NY/MedReimbDrugsFormulary.csv") |>
  clean_names() |>
  mutate(
    labeler_ndc = substr(ndc, 1, 5),
    product_ndc = substr(ndc, 7, 10),
    package_ndc = substr(ndc, 10, 11),
    ny_row = 1
  )

# Load California PDL data

ca_pdl <- read_xlsx("raw_data/PDL/CA/Medi-Cal_Rx_Approved_NDC_List.xlsx",
                    sheet = "NDC",
                    skip = 2) |>
  clean_names() |>
  filter(prior_authorization == "Yes") |>
  select(product_id) |>
  rename(ndc = product_id) |>
  mutate(
    labeler_ndc = substr(ndc, 1, 5),
    product_ndc = substr(ndc, 6, 9),
  ) |>
  distinct(labeler_ndc, product_ndc, .keep_all = TRUE) |>
  mutate(ca_row = 1)

# Load Washington PDL data

wa_pdl <- read_xlsx("raw_data/PDL/WA/apple-health-preferred-drug-list.xlsx",
                    sheet = "PDL",
                    skip = 3) |>
  clean_names() |>
  filter(identifier_type == "NDC") |>
  filter(preferred_status == "P") |>
  select(product_identifier) |>
  rename(ndc = product_identifier) |>
  mutate(
    labeler_ndc = substr(ndc, 1, 5),
    product_ndc = substr(ndc, 6, 9)
  ) |>
  distinct(labeler_ndc, product_ndc, .keep_all = TRUE) |>
  mutate(wa_row = 1)


# Note:
# The NY PDL data only contains the latest cross section of preferred drugs
# The Texas PDL data contains the historical list
# Representativeness test:
# Observe how many of the drugs in the Texas PDL can be found in the NY PDL
# Left join with NY PDL as base and TX PDL as the one being joined to

ny_texas_pdl_matches <- ny_pdl |>
  left_join(texas_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(texas_row = ifelse(is.na(texas_row), 0, 1)) |>
  summarise(mean_match = mean(texas_row)) |>
  pull()

texas_ny_pdl_matches <- texas_pdl |>
  left_join(ny_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(ny_row = ifelse(is.na(ny_row), 0, 1)) |>
  summarise(mean_match = mean(ny_row)) |>
  pull()

ca_texas_pdl_matches <- ca_pdl |>
  left_join(texas_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(texas_row = ifelse(is.na(texas_row), 0, 1)) |>
  summarise(mean_match = mean(texas_row)) |>
  pull()

texas_ca_pdl_matches <- texas_pdl |>
  left_join(ca_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(ca_row = ifelse(is.na(ca_row), 0, 1)) |>
  summarise(mean_match = mean(ca_row)) |>
  pull()

wa_texas_pdl_matches <- wa_pdl |>
  left_join(texas_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(texas_row = ifelse(is.na(texas_row), 0, 1)) |>
  summarise(mean_match = mean(texas_row)) |>
  pull()

texas_wa_pdl_matches <- texas_pdl |>
  left_join(wa_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(wa_row = ifelse(is.na(wa_row), 0, 1)) |>
  summarise(mean_match = mean(wa_row)) |>
  pull()

ca_ny_pdl_matches <- ca_pdl |>
  left_join(ny_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(ny_row = ifelse(is.na(ny_row), 0, 1)) |>
  summarise(mean_match = mean(ny_row)) |>
  pull()

ny_ca_pdl_matches <- ny_pdl |>
  left_join(ca_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(ca_row = ifelse(is.na(ca_row), 0, 1)) |>
  summarise(mean_match = mean(ca_row)) |>
  pull()

ca_wa_pdl_matches <- ca_pdl |>
  left_join(wa_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(wa_row = ifelse(is.na(wa_row), 0, 1)) |>
  summarise(mean_match = mean(wa_row)) |>
  pull()

wa_ca_pdl_matches <- wa_pdl |>
  left_join(ca_pdl, by = c("labeler_ndc", "product_ndc")) |>
  mutate(ca_row = ifelse(is.na(ca_row), 0, 1)) |>
  summarise(mean_match = mean(ca_row)) |>
  pull()

overlap_table <- tibble(
  state_1 = c("NY", "CA", "WA", "CA", "CA"),
  state_2 = c("TX", "TX", "TX", "NY", "WA"),
  overlap_1 = c(
    ny_texas_pdl_matches,
    ca_texas_pdl_matches,
    wa_texas_pdl_matches,
    ca_ny_pdl_matches,
    ca_wa_pdl_matches
  ),
  overlap_2 = c(
    texas_ny_pdl_matches,
    texas_ca_pdl_matches,
    texas_wa_pdl_matches,
    ny_ca_pdl_matches,
    wa_ca_pdl_matches
  )
) |>
  mutate(max_overlap = pmax(overlap_1, overlap_2))

kable(
  overlap_table,
  digits = 3,
  col.names = c("State 1", "State 2", "Overlap 1", "Overlap 2", "Max overlap"),
  caption = "Pairwise overlap of Medicaid Preferred Drug Lists (PDLs). Overlap 1 denotes the share of drugs in State 1's PDL that also appear in State 2's PDL. Overlap 2 denotes the share of drugs in State 2's PDL that also appear in State 1's PDL."
)