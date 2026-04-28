# regressions/RQ3.r
#
# Final RQ3 regressions: triple-difference allowing class-specific slopes
# on MedicaidShare x Post. Reference (omitted) class is V (Various) to
# match the theoretical model's choice of baseline class.
#
# Outputs (regressions/output/RQ3/):
#   - fda_models.tex,    patent_models.tex (modelsummary tables of beta_c)

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(fixest)
library(knitr)
library(modelsummary)
library(kableExtra)
library(here)

# --- Load data ------------------------------------------------------------

fda_approvals <- read_csv(
  here("processed_data/fda_approvals/fda_sdud_pdl_merged.csv")
) |>
  rename(approvals = n_approvals, firm = sponsor_name)

patents <- read_csv(
  here("processed_data/patents/patents_sdud_pdl_merged.csv")
) |>
  rename(atc_level1 = ATC1, patents = num_patents)

# --- Build ex-ante demand --- #

fda_ex_ante_stats <- fda_approvals |>
  filter(approval_year <= 2013 & approval_year >= 2010) |>
  distinct(atc_level1, approval_year, approval_quarter, prescriptions)

fda_ex_ante_total_prescriptions <- fda_ex_ante_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

fda_ex_ante_demand <- fda_ex_ante_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_ante = sum(prescriptions, na.rm = TRUE) /
                        fda_ex_ante_total_prescriptions,
            .groups = "drop")

patents_ex_ante_stats <- patents |>
  filter(patent_year <= 2013 & patent_year >= 2010) |>
  distinct(atc_level1, patent_year, patent_quarter, prescriptions)

patents_ex_ante_total_prescriptions <- patents_ex_ante_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

patents_ex_ante_demand <- patents_ex_ante_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_ante = sum(prescriptions, na.rm = TRUE) /
                        patents_ex_ante_total_prescriptions,
            .groups = "drop")

# --- Build ex-post demand --- #

fda_ex_post_stats <- fda_approvals |>
  filter(approval_year %in% 2014:2015) |>
  distinct(atc_level1, approval_year, approval_quarter, prescriptions)

fda_ex_post_total_prescriptions <- fda_ex_post_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

fda_ex_post_demand <- fda_ex_post_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_post = sum(prescriptions, na.rm = TRUE) /
                        fda_ex_post_total_prescriptions,
            .groups = "drop")

patents_ex_post_stats <- patents |>
  filter(patent_year %in% 2014:2015) |>
  distinct(atc_level1, patent_year, patent_quarter, prescriptions)

patents_ex_post_total_prescriptions <- patents_ex_post_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

patents_ex_post_demand <- patents_ex_post_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_post = sum(prescriptions, na.rm = TRUE) /
                        patents_ex_post_total_prescriptions,
            .groups = "drop")

# --- Merge demand measures + post indicator --- #

fda_approvals <- fda_approvals |>
  left_join(fda_ex_ante_demand, by = "atc_level1") |>
  left_join(fda_ex_post_demand, by = "atc_level1") |>
  filter(approval_year >= 2010) |>
  mutate(post = ifelse(approval_year >= 2014, 1, 0))

patents <- patents |>
  left_join(patents_ex_ante_demand, by = "atc_level1") |>
  left_join(patents_ex_post_demand, by = "atc_level1") |>
  filter(patent_year >= 2010) |>
  mutate(post = ifelse(patent_year >= 2014, 1, 0))

# --- FDA regressions ------------------------------------------------------
# Baseline (omitted) class is V.

fda_ex_ante_data <- fda_approvals |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

fda_ex_post_data <- fda_approvals |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)

fda_rq3_feols_exante <- feols(
  approvals ~ medicaid_market_share : post : i(atc_level1, ref = "V") |
              firm + atc_level1 + approval_year,
  data = fda_ex_ante_data, cluster = ~ firm
)

fda_rq3_feols_expost <- feols(
  approvals ~ medicaid_market_share : post : i(atc_level1, ref = "V") |
              firm + atc_level1 + approval_year,
  data = fda_ex_post_data, cluster = ~ firm
)

fda_rq3_model_list <- list(
  "FEOLS: Ex Ante" = fda_rq3_feols_exante,
  "FEOLS: Ex Post" = fda_rq3_feols_expost
)

modelsummary(
  fda_rq3_model_list,
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  output = here("output", "RQ3", "fda_models.tex")
)

# --- Patent regressions ---------------------------------------------------

patents_ex_ante_data <- patents |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

patents_ex_post_data <- patents |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)

patent_rq3_feols_exante <- feols(
  patents ~ medicaid_market_share : post : i(atc_level1, ref = "V") |
            firm + atc_level1 + patent_year,
  data = patents_ex_ante_data, cluster = ~ firm
)

patent_rq3_feols_expost <- feols(
  patents ~ medicaid_market_share : post : i(atc_level1, ref = "V") |
            firm + atc_level1 + patent_year,
  data = patents_ex_post_data, cluster = ~ firm
)

patent_rq3_model_list <- list(
  "FEOLS: Ex Ante" = patent_rq3_feols_exante,
  "FEOLS: Ex Post" = patent_rq3_feols_expost
)

modelsummary(
  patent_rq3_model_list,
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  output = here("output", "RQ3", "patent_models.tex")
)

cat("\nRQ3 outputs written to: output/RQ3/\n")
