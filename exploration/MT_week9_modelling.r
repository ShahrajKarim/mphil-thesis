#exploration/MT_week9_modelling.r


# Load necessary libraries
library(tidyverse)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(fixest)
library(knitr)
library(modelsummary)
library(kableExtra)
library(fixest)
library(TwoWayFEWeights)
library(DIDmultiplegtDYN)
library(did)
library(here)

# Load in the relevant datasets

fda_approvals <- read_csv(
    here(
  "processed_data/fda_approvals/fda_sdud_merged.csv"
    )
)

patents <- read_csv(
    here(
  "processed_data/patents/patents_sdud_merged.csv"
    )
)

expansion_prop_quarterly <- read_csv(
    here(
  "aux_data/medicaid_expansion_proportion_by_quarter.csv"
    )
)

expansion_prop_annual <- read_csv(
    here(
  "aux_data/medicaid_expansion_proportion_by_year.csv"
    )
)


# Merge rollout_data with fda_approvals and patents datasets

last_prop <- expansion_prop_quarterly |>
  arrange(year, quarter) |>
  slice_tail(n = 1) |>
  pull(prop_expanded)

fda_approvals <- fda_approvals |>
  left_join(expansion_prop_quarterly, by = c("approval_year" = "year", "approval_quarter" = "quarter")) |>
  mutate(
    prop_expanded = ifelse(approval_year < 2014, 0, prop_expanded),
    prop_expanded = ifelse(approval_year > 2023, last_prop, prop_expanded)
    )

patents <- patents |>  
  left_join(expansion_prop_quarterly, by = c("patent_year" = "year", "patent_quarter" = "quarter") ) |>
  mutate(
    prop_expanded = ifelse(patent_year < 2014, 0, prop_expanded),
    prop_expanded = ifelse(patent_year > 2023, last_prop, prop_expanded)
    )


# Create national-level summaries 

fda_national <- fda_approvals |>
  group_by(
    approval_year,
    approval_quarter,
    sponsor_name,
    labeler_code,
    atc_level1,
    atc_level2
  ) |>
  mutate(
    prescriptions_total = sum(prescriptions, na.rm = TRUE),
    medicaid_reimbursed_total = sum(medicaid_reimbursed, na.rm = TRUE),
    non_medicaid_reimbursed_total = sum(non_medicaid_reimbursed, na.rm = TRUE)
  ) |>
  ungroup() |>
  distinct(
    approval_year,
    approval_quarter,
    sponsor_name,
    labeler_code,
    atc_level1,
    atc_level2,
    prescriptions_total,
    medicaid_reimbursed_total,
    non_medicaid_reimbursed_total,
    .keep_all = TRUE
  ) |>
  select(-state)
  

patents_national <- patents |>
  group_by(
    patent_year,
    patent_quarter,
    ndc_firm_id,
    firm,
    ATC1
  ) |>
  mutate(
    prescriptions_total = sum(prescriptions, na.rm = TRUE),
    medicaid_reimbursed_total = sum(medicaid_reimbursed, na.rm = TRUE),
    non_medicaid_reimbursed_total = sum(non_medicaid_reimbursed, na.rm = TRUE),

  ) |>
  ungroup() |>
  select(-state)

# Research Question 1: Effect of Medicaid Expansion on Innovation i.e. market size

fda_atc_panel <- fda_national |>
  group_by(atc_level1, approval_year, approval_quarter) |>
  summarise(
    approvals = sum(n_approvals, na.rm = TRUE),
    .groups = "drop"
  )

patents_atc_panel <- patents_national |>
  group_by(ATC1, patent_year, patent_quarter) |>
  summarise(
    patents = sum(num_patents, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(atc_level1 = ATC1)

# Build ex-ante demand

pre_period <- fda_national |>
  filter(approval_year <= 2013)

ex_ante_window <- pre_period |>
  group_by(atc_level1) |>
  summarise(
    pre_mcaid_prescriptions = sum(ATC1_prescriptions, na.rm = TRUE),
    .groups = "drop"
  )

pre_total_prescriptions <- sum(ex_ante_window$pre_mcaid_prescriptions, na.rm = TRUE)

demand_ex_ante <- ex_ante_window |>
  mutate(
    d_ex_ante = pre_mcaid_prescriptions / pre_total_prescriptions
  ) |>
  select(atc_level1, d_ex_ante)

# Build ex-post demand

ex_post_window <- fda_national |>
  filter(approval_year %in% 2014:2015) |>
  group_by(atc_level1) |>
  summarise(
    post_mcaid_prescriptions = sum(ATC1_prescriptions, na.rm = TRUE),
    .groups = "drop"
  )

post_total_prescriptions <- sum(ex_post_window$post_mcaid_prescriptions, na.rm = TRUE)

demand_ex_post <- ex_post_window |>
  mutate(
    d_ex_post = post_mcaid_prescriptions / post_total_prescriptions
  ) |>
  select(atc_level1, d_ex_post)

fda_atc_panel <- fda_atc_panel |>
  left_join(demand_ex_ante, by = "atc_level1") |>
  left_join(demand_ex_post, by = "atc_level1")

patents_atc_panel <- patents_atc_panel |>
  left_join(demand_ex_ante, by = "atc_level1") |>
  left_join(demand_ex_post, by = "atc_level1")

# Run regression for FDA data

fda_rq1_feols_exante <- feols(
  approvals ~ i(approval_year, d_ex_ante, ref = 2013) | 
              atc_level1 + approval_year,
  data = fda_atc_panel
)

fda_rq1_feols_expost <- feols(
  approvals ~ i(approval_year, d_ex_post, ref = 2013) | 
              atc_level1 + approval_year,
  data = fda_atc_panel
)

fda_rq1_nb_exante <- fenegbin(
  approvals ~ i(approval_year, d_ex_ante, ref = 2013) | 
              atc_level1 + approval_year,
  data = fda_atc_panel
)

fda_rq1_nb_expost <- fenegbin(
  approvals ~ i(approval_year, d_ex_post, ref = 2013) | 
              atc_level1 + approval_year,
  data = fda_atc_panel
)

fda_rq1_model_list <- list(
  "FEOLS: Ex Ante"  = fda_rq1_feols_exante,
  "FEOLS: Ex Post"  = fda_rq1_feols_expost,
  "NB: Ex Ante"     = fda_rq1_nb_exante,
  "NB: Ex Post"     = fda_rq1_nb_expost
)

fda_rq1_models <- modelsummary(
             fda_rq1_model_list,
             stars = c(`*` = 0.05, `**` = 0.01, `***` = 0.001)
             )


# Run regression for patent data

patent_rq1_feols_exante <- feols(
  patents ~ i(patent_year, d_ex_ante, ref = 2013) | 
              atc_level1 + patent_year,
  data = patents_atc_panel
)

patent_rq1_feols_expost <- feols(
  patents ~ i(patent_year, d_ex_post, ref = 2013) | 
              atc_level1 + patent_year,
  data = patents_atc_panel
)

patent_rq1_nb_exante <- fenegbin(
  patents ~ i(patent_year, d_ex_ante, ref = 2013) | 
              atc_level1 + patent_year,
  data = patents_atc_panel
)

patent_rq1_nb_expost <- fenegbin(
  patents ~ i(patent_year, d_ex_post, ref = 2013) | 
              atc_level1 + patent_year,
  data = patents_atc_panel
)

patent_rq1_model_list <- list(
  "FEOLS: Ex Ante"  = patent_rq1_feols_exante,
  "FEOLS: Ex Post"  = patent_rq1_feols_expost,
  "NB: Ex Ante"     = patent_rq1_nb_exante,
  "NB: Ex Post"     = patent_rq1_nb_expost
)

patent_rq1_models <- modelsummary(
             patent_rq1_model_list,
             stars = c(`*` = 0.05, `**` = 0.01, `***` = 0.001)
             )


# --- Research Question 2: Reallocation across therapeutic areas --- #

# FDA side of RQ2

fda_national <- fda_national |>
  mutate(post = as.integer(approval_year >= 2014))

fda_atc_agg <- fda_national |>
  group_by(
    labeler_code,
    sponsor_name,
    atc_level1,
    post
  ) |>
  summarise(
    approvals_atc_post = sum(n_approvals),
    .groups = "drop"
  )

fda_df_total_firm_post <- fda_atc_agg |>
  group_by(
    labeler_code,
    sponsor_name,
    post
  ) |>
  summarise(
    total_approvals_firm_post = sum(approvals_atc_post),
    .groups = "drop"
  )


fda_rq2_df <- fda_atc_agg |>
  inner_join(fda_df_total_firm_post, by = c("labeler_code", "sponsor_name", "post")) |>
  mutate(
    share_approvals = ifelse(
      total_approvals_firm_post > 0,
      approvals_atc_post / total_approvals_firm_post,
      0
    )
  )

fda_names <- fda_national |>
  select(sponsor_name, prescriptions) |>
  group_by(sponsor_name) |>
  summarise(
    medicaid_present = as.integer(any(!is.na(prescriptions))),
    .groups = "drop"
  ) |>
  ungroup()

fda_rq2_df <- fda_rq2_df |>
    left_join(fda_names, by = "sponsor_name")

fda_rq2_reg <- feols(
  share_approvals ~ medicaid_present * post,
  cluster = ~ sponsor_name,
  data = fda_rq2_df
)

# Patent side of RQ2

patents_national <- patents_national |>
  mutate(post = as.integer(patent_year >= 2014))

patent_atc_agg <- patents_national |>
  group_by(
    ndc_firm_id,
    firm,
    ATC1,
    post
  ) |>
  summarise(
    patents_atc_post = sum(num_patents),
    .groups = "drop"
  )

patent_df_total_firm_post <- patent_atc_agg |>
  group_by(
    ndc_firm_id,
    firm,
    post
  ) |>
  summarise(
    total_patents_firm_post = sum(patents_atc_post),
    .groups = "drop"
  )


patent_rq2_df <- patent_atc_agg |>
  inner_join(patent_df_total_firm_post, by = c("ndc_firm_id", "firm", "post")) |>
  mutate(
    share_patents = ifelse(
      total_patents_firm_post > 0,
      patents_atc_post / total_patents_firm_post,
      0
    )
  )

patent_names <- patents_national |>
  select(firm, prescriptions) |>
  group_by(firm) |>
  summarise(
    medicaid_present = as.integer(any(!is.na(prescriptions))),
    .groups = "drop"
  ) |>
  ungroup()

patent_rq2_df <- patent_rq2_df |>
    left_join(patent_names, by = "firm")

patent_rq2_reg <- feols(
  share_patents ~ medicaid_present * post,
  cluster = ~ firm,
  data = patent_rq2_df
)

rq2_model_list <- list(
  "Patents"  = patent_rq2_reg,
  "FDA Approvals"  = fda_rq2_reg
)

rq2_models <- modelsummary(
             rq2_model_list,
             stars = c(`*` = 0.05, `**` = 0.01, `***` = 0.001)
             )

# --- Research Question 3: Separating Market Size and Pricing Pressure --- #

fda_rq3_df <- fda_atc_panel |>
    mutate(post = ifelse(approval_year >= 2014, 1, 0)) |>
  left_join(expansion_prop_quarterly, by = c("approval_year" = "year", "approval_quarter" = "quarter")) |>
  mutate(
    prop_expanded = ifelse(approval_year < 2014, 0, prop_expanded),
    prop_expanded = ifelse(approval_year > 2023, last_prop, prop_expanded)
    )

fda_rq3_ex_ante <- feols(
    approvals ~ prop_expanded * post + d_ex_ante * post,
    data = fda_rq3_df
)