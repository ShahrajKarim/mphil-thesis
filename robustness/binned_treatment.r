# robustness/binned_treatment.r
#
# Robustness: binned (tercile) treatment event study. Partitions
# d_ex_ante into terciles at the ATC class level and replaces the
# continuous interaction with two dummy interactions (Medium × year,
# High × year), using the Low tercile as the omitted group.
#
# This tests the linearity-in-dose assumption imposed by the main
# FEOLS specification: if effects are proportional to Medicaid exposure,
# the High-bin coefficients should be roughly twice the Medium-bin
# coefficients. Non-monotone or non-proportional patterns indicate
# nonlinearity in the dose-response.
#
# Note on scale: d_ex_ante ranges 0-0.36 since class shares sum to 1
# across all ATC Level 1 categories. Tercile cutpoints therefore sit
# at approximately 0.12 and 0.24, capturing near-zero, moderate, and
# high-exposure classes respectively.
#
# Outputs (output/robustness/):
#   - binned_patent_medium_es.png  (event study: patents, medium bin)
#   - binned_patent_high_es.png    (event study: patents, high bin)
#   - binned_fda_medium_es.png     (event study: FDA approvals, medium bin)
#   - binned_fda_high_es.png       (event study: FDA approvals, high bin)
#   - binned_patent_models.tex    (coefficient table: patents)
#   - binned_fda_models.tex       (coefficient table: FDA approvals)

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readr)
library(ggplot2)
library(fixest)
library(modelsummary)
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

# --- Build ex-ante demand (identical to RQ1.r) ----------------------------

fda_ex_ante_stats <- fda_approvals |>
  filter(approval_year <= 2013 & approval_year >= 2010) |>
  distinct(atc_level1, approval_year, approval_quarter, prescriptions)

fda_ex_ante_total <- fda_ex_ante_stats |>
  summarise(total = sum(prescriptions, na.rm = TRUE)) |>
  pull()

fda_ex_ante_demand <- fda_ex_ante_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_ante = sum(prescriptions, na.rm = TRUE) / fda_ex_ante_total,
            .groups = "drop")

patents_ex_ante_stats <- patents |>
  filter(patent_year <= 2013 & patent_year >= 2010) |>
  distinct(atc_level1, patent_year, patent_quarter, prescriptions)

patents_ex_ante_total <- patents_ex_ante_stats |>
  summarise(total = sum(prescriptions, na.rm = TRUE)) |>
  pull()

patents_ex_ante_demand <- patents_ex_ante_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_ante = sum(prescriptions, na.rm = TRUE) / patents_ex_ante_total,
            .groups = "drop")

fda_approvals <- fda_approvals |>
  left_join(fda_ex_ante_demand, by = "atc_level1") |>
  filter(approval_year >= 2010)

patents <- patents |>
  left_join(patents_ex_ante_demand, by = "atc_level1") |>
  filter(patent_year >= 2010)

# --- Collapse patents to yearly panel -------------------------------------
# The two i() interaction terms double the model matrix relative to a single
# interaction spec and exhaust vector memory at the quarterly grain.
# patent_year is the only time index in the model so aggregation is lossless.

patents_year <- patents |>
  group_by(firm, atc_level1, patent_year) |>
  summarise(
    patents   = sum(patents,   na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    .groups   = "drop"
  ) |>
  mutate(firm = as.integer(as.factor(firm)))

# --- Assign tercile bins --------------------------------------------------
# Bins are computed at the ATC class level (13 classes) so that
# assignment is not contaminated by unequal firm counts across classes.

assign_bins <- function(df) {
  class_terciles <- df |>
    distinct(atc_level1, d_ex_ante) |>
    mutate(
      dose_bin   = ntile(d_ex_ante, 3),  # 1 = Low, 2 = Medium, 3 = High
      bin_medium = as.integer(dose_bin == 2),
      bin_high   = as.integer(dose_bin == 3)
      # Low (dose_bin == 1) is the omitted reference group
    ) |>
    select(atc_level1, dose_bin, bin_medium, bin_high)

  df |> left_join(class_terciles, by = "atc_level1")
}

patents_year  <- assign_bins(patents_year)
fda_approvals <- assign_bins(fda_approvals)

# --- Event-study specifications -------------------------------------------
# Y_{fct} = alpha_{fc} + lambda_t
#   + sum_tau beta_tau^M (1{t=tau} x bin_medium_c)
#   + sum_tau beta_tau^H (1{t=tau} x bin_high_c)
#   + eps_{fct}
#
# Low bin is omitted; medium and high coefficients are relative to Low.

# Patents (yearly panel)
patent_binned <- feols(
  patents ~ i(patent_year, bin_medium, ref = 2013) +
            i(patent_year, bin_high,   ref = 2013) |
            firm^atc_level1 + patent_year,
  data    = patents_year,
  cluster = ~ firm
)

# FDA approvals (quarterly — smaller so no memory issue)
fda_binned <- feols(
  approvals ~ i(approval_year, bin_medium, ref = 2013) +
              i(approval_year, bin_high,   ref = 2013) |
              firm^atc_level1 + approval_year,
  data    = fda_approvals,
  cluster = ~ firm
)

# --- Plots ----------------------------------------------------------------

save_iplot <- function(model, file, main, i.select = 1) {
  png(here("output", "robustness", "binned_treatment", file), width = 1600, height = 1000, res = 200)
  iplot(
    model,
    i.select = i.select,
    ci.level = 0.95,
    ref.line = 0,
    xlab     = "Year",
    ylab     = "Coefficient (relative to Low bin, 2013)",
    main     = main
  )
  dev.off()
}

# i.select = 1: bin_medium; i.select = 2: bin_high
save_iplot(patent_binned, "binned_patent_medium_es.png", "Binned Event Study: Patents — Medium Bin", i.select = 1)
save_iplot(patent_binned, "binned_patent_high_es.png",   "Binned Event Study: Patents — High Bin",   i.select = 2)
save_iplot(fda_binned,    "binned_fda_medium_es.png",    "Binned Event Study: FDA Approvals — Medium Bin", i.select = 1)
save_iplot(fda_binned,    "binned_fda_high_es.png",      "Binned Event Study: FDA Approvals — High Bin",   i.select = 2)

# --- Coefficient tables ---------------------------------------------------

modelsummary(
  list("Patents: Binned" = patent_binned),
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_treatment", "binned_patent_models.tex")
)

modelsummary(
  list("FDA Approvals: Binned" = fda_binned),
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_treatment", "binned_fda_models.tex")
)

# --- End of file ----------------------------------------------------------
