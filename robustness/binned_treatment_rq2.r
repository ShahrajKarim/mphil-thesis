# robustness/binned_treatment_rq2.r
#
# Robustness: binned (tercile) treatment event study for RQ2. Extends
# binned_treatment.r by simultaneously binning both d_ex_ante (class-
# level Medicaid exposure) and in_pdl (firm-class PDL exposure) into
# terciles, then replacing both continuous interactions from the main
# RQ2 FEOLS with dummy × year interactions.
#
# Medicaid bins are computed at the ATC class level (13 classes).
# PDL bins are computed at the firm-class level using the pre-2014
# time-average of in_pdl so that the assignment is predetermined.
# Low bin is the omitted reference group for both sets of interactions.
#
# Two separate event-study plots are produced: one showing how the
# response varies across Medicaid exposure bins (holding PDL fixed),
# and one showing how it varies across PDL bins (holding Medicaid fixed).
#
# Note: contdid and didimputation are not extended to RQ2 due to the
# computational cost of handling two continuous treatments simultaneously.
#
# Outputs (output/robustness/):
#   - binned_rq2_patent_med_medicaid_es.png
#   - binned_rq2_patent_high_medicaid_es.png
#   - binned_rq2_patent_med_pdl_es.png
#   - binned_rq2_patent_high_pdl_es.png
#   - binned_rq2_fda_med_medicaid_es.png
#   - binned_rq2_fda_high_medicaid_es.png
#   - binned_rq2_fda_med_pdl_es.png
#   - binned_rq2_fda_high_pdl_es.png
#   - binned_rq2_patent_models.tex
#   - binned_rq2_fda_models.tex

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

# --- Build ex-ante demand (identical to RQ1.r / RQ2.r) -------------------

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

# --- Collapse patents to yearly panel (identical to RQ2.r) ----------------

patents_year <- patents |>
  group_by(firm, atc_level1, patent_year) |>
  summarise(
    patents   = sum(patents,   na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    in_pdl    = max(in_pdl,   na.rm = TRUE),
    .groups   = "drop"
  ) |>
  mutate(firm = as.integer(as.factor(firm)))

fda_year <- fda_approvals |>
  group_by(firm, atc_level1, approval_year) |>
  summarise(
    approvals = sum(approvals, na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    in_pdl    = max(in_pdl,   na.rm = TRUE),
    .groups   = "drop"
  )

# --- Assign Medicaid tercile bins -----------------------------------------
# Computed at the ATC class level so bin assignment is the same for all
# firms within a class, mirroring the class-level variation in d_ex_ante.

assign_medicaid_bins <- function(df) {
  class_bins <- df |>
    distinct(atc_level1, d_ex_ante) |>
    mutate(
      medicaid_bin      = ntile(d_ex_ante, 3),
      bin_med_medicaid  = as.integer(medicaid_bin == 2),
      bin_high_medicaid = as.integer(medicaid_bin == 3)
      # Low (medicaid_bin == 1) is the omitted reference group
    ) |>
    select(atc_level1, medicaid_bin, bin_med_medicaid, bin_high_medicaid)

  df |> left_join(class_bins, by = "atc_level1")
}

# --- Assign PDL tercile bins ----------------------------------------------
# in_pdl varies at the firm-class level. Bins are computed using the
# pre-2014 average of in_pdl for each firm-class pair so that assignment
# is predetermined with respect to the policy shock.

assign_pdl_bins <- function(df, tvar) {
  pre_pdl <- df |>
    filter(.data[[tvar]] < 2014) |>
    group_by(firm, atc_level1) |>
    summarise(in_pdl_pre = mean(in_pdl, na.rm = TRUE), .groups = "drop")

  pdl_bins <- pre_pdl |>
    mutate(
      pdl_bin      = ntile(in_pdl_pre, 3),
      bin_med_pdl  = as.integer(pdl_bin == 2),
      bin_high_pdl = as.integer(pdl_bin == 3)
      # Low (pdl_bin == 1) is the omitted reference group
    ) |>
    select(firm, atc_level1, pdl_bin, bin_med_pdl, bin_high_pdl)

  df |> left_join(pdl_bins, by = c("firm", "atc_level1"))
}

patents_year <- patents_year |>
  assign_medicaid_bins() |>
  assign_pdl_bins("patent_year") |>
  select(firm, atc_level1, patent_year, patents,
         bin_med_medicaid, bin_high_medicaid, bin_med_pdl, bin_high_pdl)

fda_year <- fda_year |>
  assign_medicaid_bins() |>
  assign_pdl_bins("approval_year") |>
  select(firm, atc_level1, approval_year, approvals,
         bin_med_medicaid, bin_high_medicaid, bin_med_pdl, bin_high_pdl)

# --- Event-study specifications -------------------------------------------
# Single model with all four dummy interactions. Binary dummies keep Low as
# the natural omitted group. Formula order:
#   i=1 bin_med_medicaid, i=2 bin_high_medicaid,
#   i=3 bin_med_pdl,      i=4 bin_high_pdl

patent_binned_rq2 <- feols(
  patents ~ i(patent_year, bin_med_medicaid,  ref = 2013) +
            i(patent_year, bin_high_medicaid, ref = 2013) +
            i(patent_year, bin_med_pdl,       ref = 2013) +
            i(patent_year, bin_high_pdl,      ref = 2013) |
            firm^atc_level1 + patent_year,
  data    = patents_year,
  cluster = ~ firm
)

fda_binned_rq2 <- feols(
  approvals ~ i(approval_year, bin_med_medicaid,  ref = 2013) +
              i(approval_year, bin_high_medicaid, ref = 2013) +
              i(approval_year, bin_med_pdl,       ref = 2013) +
              i(approval_year, bin_high_pdl,      ref = 2013) |
              firm^atc_level1 + approval_year,
  data    = fda_year,
  cluster = ~ firm
)

# --- Plot helper ----------------------------------------------------------
# i.select = c(1,2) shows Medicaid Med + High on one plot;
# i.select = c(3,4) shows PDL Med + High on one plot.

save_iplot_rq2 <- function(model, file, main, i.select) {
  png(here("output", "robustness", "binned_treatment_rq2", file), width = 1600, height = 1000, res = 200)
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

# i=1 bin_med_medicaid, i=2 bin_high_medicaid,
# i=3 bin_med_pdl,      i=4 bin_high_pdl

# Patents
save_iplot_rq2(patent_binned_rq2, "binned_rq2_patent_med_medicaid_es.png",
               "Binned Event Study: Patents — Medicaid Medium", i.select = 1)
save_iplot_rq2(patent_binned_rq2, "binned_rq2_patent_high_medicaid_es.png",
               "Binned Event Study: Patents — Medicaid High",   i.select = 2)
save_iplot_rq2(patent_binned_rq2, "binned_rq2_patent_med_pdl_es.png",
               "Binned Event Study: Patents — PDL Medium",      i.select = 3)
save_iplot_rq2(patent_binned_rq2, "binned_rq2_patent_high_pdl_es.png",
               "Binned Event Study: Patents — PDL High",        i.select = 4)

# FDA approvals
save_iplot_rq2(fda_binned_rq2, "binned_rq2_fda_med_medicaid_es.png",
               "Binned Event Study: FDA Approvals — Medicaid Medium", i.select = 1)
save_iplot_rq2(fda_binned_rq2, "binned_rq2_fda_high_medicaid_es.png",
               "Binned Event Study: FDA Approvals — Medicaid High",   i.select = 2)
save_iplot_rq2(fda_binned_rq2, "binned_rq2_fda_med_pdl_es.png",
               "Binned Event Study: FDA Approvals — PDL Medium",      i.select = 3)
save_iplot_rq2(fda_binned_rq2, "binned_rq2_fda_high_pdl_es.png",
               "Binned Event Study: FDA Approvals — PDL High",        i.select = 4)

# --- Coefficient tables ---------------------------------------------------

modelsummary(
  list("Patents: Binned RQ2" = patent_binned_rq2),
  stars   = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw,        ~clean,      ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_treatment_rq2", "binned_rq2_patent_models.tex")
)

modelsummary(
  list("FDA Approvals: Binned RQ2" = fda_binned_rq2),
  stars   = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw,        ~clean,      ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_treatment_rq2", "binned_rq2_fda_models.tex")
)

# --- End of file ----------------------------------------------------------
