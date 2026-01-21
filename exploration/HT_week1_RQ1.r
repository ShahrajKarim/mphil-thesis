#exploration/HT_week1_RQ1.r


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
  mutate(firm = ifelse(is.na(firm), sponsor_name, firm))
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
  rename(atc_level1 = ATC1) |>
  mutate(
    prescriptions_total = sum(prescriptions, na.rm = TRUE),
    medicaid_reimbursed_total = sum(medicaid_reimbursed, na.rm = TRUE),
    non_medicaid_reimbursed_total = sum(non_medicaid_reimbursed, na.rm = TRUE),
  ) |>
  ungroup() |>
  select(-state)

# Research Question 1: Effect of Medicaid Expansion on Innovation i.e. market size

fda_atc_panel <- fda_national |>
  group_by(atc_level1, firm, approval_year, approval_quarter) |>
  distinct(atc_level1, firm, approval_year, approval_quarter, n_approvals) |>
  summarise(
    approvals = sum(n_approvals, na.rm = TRUE),
    .groups = "drop"
  )

patents_atc_panel <- patents_national |>
  group_by(atc_level1, firm, patent_year, patent_quarter) |>
  distinct(atc_level1, firm, patent_year, patent_quarter, num_patents) |>
  summarise(
    patents = sum(num_patents, na.rm = TRUE),
    .groups = "drop"
  )

# Build ex-ante demand

pre_period_fda <- fda_national |>
  filter(approval_year <= 2013 & approval_year >= 2010)

ex_ante_window_fda <- pre_period_fda |>
  group_by(atc_level1) |>
  summarise(
    pre_mcaid_prescriptions = sum(ATC1_prescriptions, na.rm = TRUE),
    .groups = "drop"
  )

pre_total_prescriptions_fda <- sum(ex_ante_window_fda$pre_mcaid_prescriptions, na.rm = TRUE)

demand_ex_ante_fda <- ex_ante_window_fda |>
  mutate(
    d_ex_ante = pre_mcaid_prescriptions / pre_total_prescriptions_fda
  ) |>
  select(atc_level1, d_ex_ante)

pre_period_patents <- patents_national |>
  filter(patent_year <= 2013 & patent_year >= 2010)

ex_ante_window_patents <- pre_period_patents |>
  group_by(atc_level1) |>
  summarise(
    pre_mcaid_prescriptions = sum(ATC1_prescriptions, na.rm = TRUE),
    .groups = "drop"
  )

pre_total_prescriptions_patents <- sum(ex_ante_window_patents$pre_mcaid_prescriptions)

demand_ex_ante_patents <- ex_ante_window_patents |>
  mutate(
    d_ex_ante = pre_mcaid_prescriptions / pre_total_prescriptions_patents
  ) |>
  select(atc_level1, d_ex_ante)



# Build ex-post demand

ex_post_window_fda <- fda_national |>
  filter(approval_year %in% 2014:2015) |>
  group_by(atc_level1) |>
  summarise(
    post_mcaid_prescriptions = sum(ATC1_prescriptions, na.rm = TRUE),
    .groups = "drop"
  )

post_total_prescriptions_fda <- sum(ex_post_window_fda$post_mcaid_prescriptions, na.rm = TRUE)

demand_ex_post_fda <- ex_post_window_fda |>
  mutate(
    d_ex_post = post_mcaid_prescriptions / post_total_prescriptions_fda
  ) |>
  select(atc_level1, d_ex_post)

ex_post_window_patents <- patents_national |>
  filter(patent_year %in% 2014:2015) |>
  group_by(atc_level1) |>
  summarise(
    post_mcaid_prescriptions = sum(ATC1_prescriptions, na.rm = TRUE),
    .groups = "drop"
  )

post_total_prescriptions_patents <- sum(ex_post_window_patents$post_mcaid_prescriptions, na.rm = TRUE)

demand_ex_post_patents <- ex_post_window_patents |>
  mutate(
    d_ex_post = post_mcaid_prescriptions / post_total_prescriptions_patents
  ) |>
  select(atc_level1, d_ex_post)

fda_atc_panel <- fda_atc_panel |>
  left_join(demand_ex_ante_fda, by = "atc_level1") |>
  left_join(demand_ex_post_fda, by = "atc_level1") |>
  filter(approval_year >= 2010)

patents_atc_panel <- patents_atc_panel |>
  left_join(demand_ex_ante_patents, by = "atc_level1") |>
  left_join(demand_ex_post_patents, by = "atc_level1") |>
  filter(patent_year >= 2010)

# Run regression for FDA data

fda_ex_ante_data <- fda_atc_panel |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

fda_ex_post_data <- fda_atc_panel |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)


fda_rq1_feols_exante <- feols(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + approval_year,
              data = fda_ex_ante_data,
              cluster = ~ firm
)

fda_rq1_feols_expost <- feols(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + approval_year,
              data = fda_ex_post_data,
              cluster = ~ firm
)

fda_rq1_nb_exante <- fenegbin(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + approval_year,
              data = fda_ex_ante_data,
              cluster = ~ firm
)

fda_rq1_nb_expost <- fenegbin(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + approval_year,
              data = fda_ex_post_data,
              cluster = ~ firm
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

patents_ex_ante_data <- patents_atc_panel |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

patents_ex_post_data <- patents_atc_panel |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)


patent_rq1_feols_exante <- feols(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + patent_year,
              data = patents_ex_ante_data,
              cluster = ~ firm
)

patent_rq1_feols_expost <- feols(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + patent_year,
              data = patents_ex_post_data,
              cluster = ~ firm
)

patent_rq1_nb_exante <- fenegbin(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + patent_year,
              data = patents_ex_ante_data,
              cluster = ~ firm
)

patent_rq1_nb_expost <- fenegbin(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) | 
              firm^atc_level1 + patent_year,
              data = patents_ex_post_data,
              cluster = ~ firm
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

# Plotting functions for parallel trends test:

iplot(
  patent_rq1_nb_exante,
  ci.level = 0.95,
  ref.line = 0,
  xlab = "Year",
  ylab = "Coefficient (relative to 2013)",
  main = "Event study: patents × ex-ante exposure"
)