# robustness/didimputation.r
#
# Robustness: Borusyak, Jaravel & Spiess (2024, ReStud) imputation
# estimator. Estimates untreated potential outcomes using only
# pre-2014 data and low-exposure classes, then constructs post-2014
# counterfactuals from those fitted values. Avoids negative weighting
# by never using post-2014 control-group variation to identify
# the baseline trend.
#
# Treatment is binarised: units with d_ex_ante above the sample median
# are treated in 2014 (G = 2014); units at or below the median serve
# as the never-treated comparison group (G = Inf). This removes the
# linearity-in-dose assumption from the main FEOLS specification.
#
# Install (if needed):
#   install.packages("didimputation")
#   # or: devtools::install_github("kylebutts/didimputation")
#
# Outputs (output/robustness/):
#   - didimputation_patent_es.png   (event study: patents)
#   - didimputation_fda_es.png      (event study: FDA approvals)

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readr)
library(ggplot2)
library(didimputation)
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

patents_year <- patents |>
  group_by(firm, atc_level1, patent_year) |>
  summarise(
    patents   = sum(patents, na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    .groups   = "drop"
  )

fda_year <- fda_approvals |>
  group_by(firm, atc_level1, approval_year) |>
  summarise(
    approvals = sum(approvals, na.rm = TRUE),
    d_ex_ante = first(d_ex_ante),
    .groups   = "drop"
  )

# --- Binarise treatment ---------------------------------------------------
# G = 2014 for units with d_ex_ante above the class-level median;
# G = Inf  for units at or below the median (never-treated comparison).
# The median is computed at the ATC class level to avoid artefacts from
# unbalanced firm counts across classes.

compute_binary_G <- function(df) {
  class_doses <- df |>
    distinct(atc_level1, d_ex_ante)
  med <- median(class_doses$d_ex_ante, na.rm = TRUE)
  df |>
    mutate(
      unit_id = as.integer(as.factor(paste(firm, atc_level1, sep = "_"))),
      G       = if_else(d_ex_ante > med, 2014, NA_real_)
    )
}

patents_di <- compute_binary_G(patents_year) |>
  complete(unit_id, patent_year, fill = list(patents = 0)) |>
  group_by(unit_id) |>
  fill(G, d_ex_ante, firm, atc_level1, .direction = "downup") |>
  ungroup()

fda_di <- compute_binary_G(fda_year) |>
  complete(unit_id, approval_year, fill = list(approvals = 0)) |>
  group_by(unit_id) |>
  fill(G, d_ex_ante, firm, atc_level1, .direction = "downup") |>
  ungroup()

# --- didimputation estimation ---------------------------------------------
# did_imputation() fits the untreated model on pre-treatment observations
# and never-treated units, imputes counterfactuals, then regresses
# the gap on event-time indicators.

patent_dimp <- did_imputation(
  data      = patents_di,
  yname     = "patents",
  gname     = "G",
  tname     = "patent_year",
  idname    = "unit_id",
  horizon   = 0:11,
  pretrends = -3:-1,
  cluster_var = "unit_id"
)

fda_dimp <- did_imputation(
  data      = fda_di,
  yname     = "approvals",
  gname     = "G",
  tname     = "approval_year",
  idname    = "unit_id",
  horizon   = 0:11,
  pretrends = -3:-1,
  cluster_var = "unit_id"
)

# --- Plot helper ----------------------------------------------------------

save_dimp_plot <- function(res, file, title_str) {
  df <- res |>
    mutate(
      ci_lo = estimate - 1.96 * std.error,
      ci_hi = estimate + 1.96 * std.error,
      horizon = as.integer(term)
    )

  p <- ggplot(df, aes(x = horizon, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "#404040") +
    geom_line(colour = "#404040") +
    geom_point(colour = "#404040", size = 1.5) +
    labs(
      x     = "Years relative to 2014",
      y     = "Estimate (relative to 2013)",
      title = title_str
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title  = element_text(size = 11),
      axis.line   = element_line(colour = "grey40"),
      axis.ticks  = element_line(colour = "grey40"),
      plot.margin = margin(8, 12, 8, 8)
    )

  ggsave(here("output", "robustness", file),
         plot = p, width = 8, height = 4, dpi = 200)
}

# --- Save plots -----------------------------------------------------------

save_dimp_plot(patent_dimp,
               "didimputation_patent_es.png",
               "BJS imputation event study: patents (above-median Medicaid exposure)")

save_dimp_plot(fda_dimp,
               "didimputation_fda_es.png",
               "BJS imputation event study: FDA approvals (above-median Medicaid exposure)")

# --- End of file ----------------------------------------------------------
