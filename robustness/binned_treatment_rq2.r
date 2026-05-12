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
#   - binned_rq2_patent_medicaid_es.png
#   - binned_rq2_patent_pdl_es.png
#   - binned_rq2_fda_medicaid_es.png
#   - binned_rq2_fda_pdl_es.png
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
  )

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
      medicaid_bin    = ntile(d_ex_ante, 3),
      bin_med_medicaid  = as.integer(medicaid_bin == 2),
      bin_high_medicaid = as.integer(medicaid_bin == 3)
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
      pdl_bin        = ntile(in_pdl_pre, 3),
      bin_med_pdl    = as.integer(pdl_bin == 2),
      bin_high_pdl   = as.integer(pdl_bin == 3)
    ) |>
    select(firm, atc_level1, pdl_bin, bin_med_pdl, bin_high_pdl)

  df |> left_join(pdl_bins, by = c("firm", "atc_level1"))
}

patents_year  <- patents_year  |> assign_medicaid_bins() |> assign_pdl_bins("patent_year")
fda_year      <- fda_year      |> assign_medicaid_bins() |> assign_pdl_bins("approval_year")

# --- Event-study specifications -------------------------------------------
# Y_{fct} = alpha_{fc} + lambda_t
#   + sum_tau beta_tau^{M,med}  (1{t=tau} x bin_med_medicaid_c)
#   + sum_tau beta_tau^{M,high} (1{t=tau} x bin_high_medicaid_c)
#   + sum_tau beta_tau^{P,med}  (1{t=tau} x bin_med_pdl_{fc})
#   + sum_tau beta_tau^{P,high} (1{t=tau} x bin_high_pdl_{fc})
#   + eps_{fct}
#
# Low bin is omitted for both sets of interactions.

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

save_binned_rq2_plot <- function(model, pattern_med, pattern_high,
                                  file, title_str) {
  ct        <- as.data.frame(model$coeftable)
  ct$term   <- rownames(ct)

  extract_bin <- function(pattern, label) {
    ct |>
      filter(str_detect(term, pattern)) |>
      mutate(
        year  = as.integer(str_extract(term, "\\d{4}")),
        bin   = label,
        ci_lo = Estimate - 1.96 * `Std. Error`,
        ci_hi = Estimate + 1.96 * `Std. Error`
      ) |>
      select(year, bin, est = Estimate, ci_lo, ci_hi)
  }

  df <- bind_rows(
    extract_bin(pattern_med,  "Medium"),
    extract_bin(pattern_high, "High")
  ) |>
    bind_rows(
      tibble(year = 2013L, bin = c("Medium", "High"),
             est = 0, ci_lo = 0, ci_hi = 0)
    ) |>
    arrange(bin, year)

  ggplot(df, aes(x = year, y = est,
                 colour = bin, fill = bin,
                 linetype = bin, shape = bin)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 2013.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.10,
                colour = NA) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_colour_manual(values = c("Medium" = "#6e6e6e", "High" = "#1a1a1a")) +
    scale_fill_manual(  values = c("Medium" = "#6e6e6e", "High" = "#1a1a1a")) +
    scale_linetype_manual(values = c("Medium" = "dashed", "High" = "solid")) +
    scale_shape_manual(   values = c("Medium" = 17,       "High" = 16)) +
    labs(
      x        = "Year",
      y        = "Coefficient (relative to Low bin, 2013)",
      title    = title_str,
      colour   = "Dose bin", fill     = "Dose bin",
      linetype = "Dose bin", shape    = "Dose bin"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(size = 11),
      legend.position = "bottom",
      axis.line       = element_line(colour = "grey40"),
      axis.ticks      = element_line(colour = "grey40"),
      plot.margin     = margin(8, 12, 8, 8)
    ) |>
    (\(p) ggsave(here("output", "robustness", file),
                 plot = p, width = 8, height = 4, dpi = 200))()
}

# Patents: Medicaid bins
save_binned_rq2_plot(
  patent_binned_rq2,
  pattern_med  = "bin_med_medicaid",
  pattern_high = "bin_high_medicaid",
  file         = "binned_rq2_patent_medicaid_es.png",
  title_str    = "Binned RQ2 event study: patents — Medicaid exposure bins (PDL controlled)"
)

# Patents: PDL bins
save_binned_rq2_plot(
  patent_binned_rq2,
  pattern_med  = "bin_med_pdl",
  pattern_high = "bin_high_pdl",
  file         = "binned_rq2_patent_pdl_es.png",
  title_str    = "Binned RQ2 event study: patents — PDL exposure bins (Medicaid controlled)"
)

# FDA: Medicaid bins
save_binned_rq2_plot(
  fda_binned_rq2,
  pattern_med  = "bin_med_medicaid",
  pattern_high = "bin_high_medicaid",
  file         = "binned_rq2_fda_medicaid_es.png",
  title_str    = "Binned RQ2 event study: FDA approvals — Medicaid exposure bins (PDL controlled)"
)

# FDA: PDL bins
save_binned_rq2_plot(
  fda_binned_rq2,
  pattern_med  = "bin_med_pdl",
  pattern_high = "bin_high_pdl",
  file         = "binned_rq2_fda_pdl_es.png",
  title_str    = "Binned RQ2 event study: FDA approvals — PDL exposure bins (Medicaid controlled)"
)

# --- Coefficient tables ---------------------------------------------------

modelsummary(
  list("Patents: Binned RQ2" = patent_binned_rq2),
  stars   = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_rq2_patent_models.tex")
)

modelsummary(
  list("FDA Approvals: Binned RQ2" = fda_binned_rq2),
  stars   = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_rq2_fda_models.tex")
)

# --- End of file ----------------------------------------------------------
