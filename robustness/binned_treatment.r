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
#   - binned_patent_es.png        (event study: patents, ex-ante)
#   - binned_fda_es.png           (event study: FDA approvals, ex-ante)
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

# --- Assign tercile bins --------------------------------------------------
# Bins are computed at the ATC class level (13 classes) so that
# assignment is not contaminated by unequal firm counts across classes.

assign_bins <- function(df) {
  class_terciles <- df |>
    distinct(atc_level1, d_ex_ante) |>
    mutate(
      dose_bin = ntile(d_ex_ante, 3),  # 1 = Low, 2 = Medium, 3 = High
      bin_medium = as.integer(dose_bin == 2),
      bin_high   = as.integer(dose_bin == 3)
      # Low (dose_bin == 1) is the omitted reference group
    ) |>
    select(atc_level1, dose_bin, bin_medium, bin_high)

  df |> left_join(class_terciles, by = "atc_level1")
}

patents      <- assign_bins(patents)
fda_approvals <- assign_bins(fda_approvals)

# Report tercile cutpoints for reference
cat("Patent panel tercile cutpoints for d_ex_ante:\n")
patents |>
  distinct(atc_level1, d_ex_ante, dose_bin) |>
  arrange(dose_bin, d_ex_ante) |>
  print()

# --- Event-study specifications -------------------------------------------
# Y_{fct} = alpha_{fc} + lambda_t
#   + sum_tau beta_tau^M (1{t=tau} x bin_medium_c)
#   + sum_tau beta_tau^H (1{t=tau} x bin_high_c)
#   + eps_{fct}
#
# Low bin is omitted; medium and high coefficients are relative to Low.

# Patents
patent_binned <- feols(
  patents ~ i(patent_year, bin_medium, ref = 2013) +
            i(patent_year, bin_high,   ref = 2013) |
            firm^atc_level1 + patent_year,
  data    = patents,
  cluster = ~ firm
)

# FDA approvals
fda_binned <- feols(
  approvals ~ i(approval_year, bin_medium, ref = 2013) +
              i(approval_year, bin_high,   ref = 2013) |
              firm^atc_level1 + approval_year,
  data    = fda_approvals,
  cluster = ~ firm
)

# --- Plot helper ----------------------------------------------------------
# Extracts both sets of event-study coefficients and plots on one panel.

save_binned_plot <- function(model, tvar, file, title_str) {

  ct <- as.data.frame(model$coeftable)
  ct$term <- rownames(ct)

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
    extract_bin("bin_medium", "Medium"),
    extract_bin("bin_high",   "High")
  )

  # Add 2013 reference rows at zero
  ref_rows <- tibble(
    year  = 2013L,
    bin   = c("Medium", "High"),
    est   = 0, ci_lo = 0, ci_hi = 0
  )
  df <- bind_rows(df, ref_rows) |> arrange(bin, year)

  p <- ggplot(df, aes(x = year, y = est,
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
    scale_shape_manual(   values = c("Medium" = 17, "High" = 16)) +
    labs(
      x        = "Year",
      y        = "Coefficient (relative to Low bin, 2013)",
      title    = title_str,
      colour   = "Dose bin",
      fill     = "Dose bin",
      linetype = "Dose bin",
      shape    = "Dose bin"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(size = 11),
      legend.position = "bottom",
      axis.line       = element_line(colour = "grey40"),
      axis.ticks      = element_line(colour = "grey40"),
      plot.margin     = margin(8, 12, 8, 8)
    )

  ggsave(here("output", "robustness", file),
         plot = p, width = 8, height = 4, dpi = 200)
}

save_binned_plot(patent_binned, "patent_year",
                 "binned_patent_es.png",
                 "Binned event study: patents (ex-ante Medicaid exposure, Low bin omitted)")

save_binned_plot(fda_binned, "approval_year",
                 "binned_fda_es.png",
                 "Binned event study: FDA approvals (ex-ante Medicaid exposure, Low bin omitted)")

# --- Coefficient tables ---------------------------------------------------

modelsummary(
  list("Patents: Binned" = patent_binned),
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_patent_models.tex")
)

modelsummary(
  list("FDA Approvals: Binned" = fda_binned),
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs",      "Num. obs.", 0,
    "r.squared", "R2",        3
  ),
  output = here("output", "robustness", "binned_fda_models.tex")
)

# --- End of file ----------------------------------------------------------
