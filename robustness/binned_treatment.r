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
# Bin assignment is explicit at the ATC class level rather than via
# ntile() because N (0.363) and C (0.215) are starkly separated from the
# rest of the distribution. R (0.085) and A (0.068) sit closer to the
# low-exposure cluster than to J/D, and including them in the Medium bin
# induces a pre-trend violation, so they are assigned to Low.
#   High   = N, C
#   Medium = J, D
#   Low    = R, A, M, H, L, G, S, B, V, P
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

# --- Collapse to yearly panel ---------------------------------------------
# The two i() interaction terms double the model matrix relative to a single
# interaction spec and exhaust vector memory at the quarterly grain. The
# year is the only time index in the model so aggregation is lossless.

patents_year <- patents |>
  group_by(firm, atc_level1, patent_year) |>
  summarise(
    patents   = sum(patents,   na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    .groups   = "drop"
  ) |>
  mutate(firm = as.integer(as.factor(firm)))

fda_year <- fda_approvals |>
  group_by(firm, atc_level1, approval_year) |>
  summarise(
    approvals = sum(approvals, na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    .groups   = "drop"
  ) |>
  mutate(firm = as.integer(as.factor(firm)))

# --- Assign exposure bins -------------------------------------------------
# Bins are hard-coded at the ATC class level. N and C are separated from
# the rest of the distribution, so they form the High bin alone. The next
# four classes form the Medium bin, and the remaining low-exposure
# classes form the Low (omitted) bin.

atc_high   <- c("N", "C")
atc_medium <- c("J", "D")

assign_bins <- function(df) {
  df |>
    mutate(
      dose_bin = case_when(
        atc_level1 %in% atc_high   ~ 3L,
        atc_level1 %in% atc_medium ~ 2L,
        TRUE                       ~ 1L
      ),
      bin_medium = as.integer(dose_bin == 2),
      bin_high   = as.integer(dose_bin == 3)
      # Low (dose_bin == 1) is the omitted reference group
    )
}

patents_year <- assign_bins(patents_year)
fda_year     <- assign_bins(fda_year)

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

# FDA approvals (yearly panel)
fda_binned <- feols(
  approvals ~ i(approval_year, bin_medium, ref = 2013) +
              i(approval_year, bin_high,   ref = 2013) |
              firm^atc_level1 + approval_year,
  data    = fda_year,
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
