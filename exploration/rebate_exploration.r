# exploration/rebate_analysis.r

# This script tests whether a correlation exists between market size and rebates
# A positive correlation suggests increasing monopsony power with state size
# Results can be used to justify assumptions made within theoretical model.


# libraries
library(here)
library(dplyr)
library(readr)
library(broom)
library(modelsummary)

# Generate lists for Multi-state PDL alliances

ssdc <- c("DE", "IA", "KY", "ME", "MS", "ND", "OH", "OK",
          "OR", "PA", "SD", "UT", "VT", "WV", "WY")

nmpi <- c("AK", "AR", "DC", "MI", "MN", "MT", "NV",
          "NH", "NY", "NC", "RI", "SC", "VA")

top_dollar <- c("CT", "GA", "KS", "LA", "MD", "NE", "WI")

# Load SDUD data

SDUD <- readRDS(here("processed_data/state_drug_utilisation_data/SDUD_full.rds")) |>
  mutate(
    state = ifelse(state %in% ssdc, "SSDC", state),
    state = ifelse(state %in% nmpi, "NMPI", state),
    state = ifelse(state %in% top_dollar, "Top$", state)
  ) |>
  group_by(state, ndc, year) |>
  summarise(
    units = sum(number_of_prescriptions, na.rm = TRUE),
    reimbursement_amount = sum(total_amount_reimbursed, na.rm = TRUE),
    rebate = ifelse(units > 0, reimbursement_amount / units, 0),
    .groups = "drop"
  )

market_sizes <- SDUD |>
  group_by(state, year) |>
  summarise(market_size = sum(units, na.rm = TRUE), .groups = "drop")


# Need to include multi-state alliances
# Market size would be combination of states

SDUD <- SDUD |>
  left_join(market_sizes, by = c("state", "year")) |>
  filter(between(
    rebate,
    quantile(rebate, 0.025, na.rm = TRUE),
    quantile(rebate, 0.975, na.rm = TRUE)
  ))

# Map NDCs to ATC mapping file

ndc_atc <- read_csv(here("aux_data/ndc_atc_mapping.csv")) |>
  select(ndc, atc = ATC_Level1)

ndc_firm <- read_csv(here("aux_data/ndc_fda_firm_mapping.csv"))

SDUD <- SDUD |>
  mutate(labeler_code = substr(ndc, 1, 5)) |>
  left_join(ndc_atc, by = "ndc") |>
  left_join(ndc_firm, by = "labeler_code")

# Regression for rebate against market size

m1 <- lm(rebate ~ log(market_size) + factor(atc) + factor(year), data = SDUD)

tidy(m1)

# Coefficient on log(market_size): 5.51
# Standard error: 0.0348
# t-statistic: 158
# p-value: < 0.001
# Interpretation:
# There is a statistically significant positive relation between market size and rebates
# Larger markets may have greater bargaining power to negotiate higher rebates.