# robustness/contdid.r
#
# Robustness: Callaway, Goodman-Bacon & Sant'Anna (2024) continuous
# treatment DiD estimator. Relaxes the linearity-in-dose assumption
# imposed by FEOLS by estimating ATT(d, t) nonparametrically across
# the distribution of pre-expansion Medicaid exposure.
#
# The dose variable is d_ex_ante (class-level ex-ante Medicaid share,
# range 0-0.36). Units in the lowest-dose tercile serve as the
# implicit comparison group as there are no true never-treated units.
#
# Install (if needed):
#   devtools::install_github("bcallaway11/contdid")
#
# Outputs (output/robustness/):
#   - contdid_patent_es.png   (event study: patents)
#   - contdid_fda_es.png      (event study: FDA approvals)

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readr)
library(ggplot2)
library(contdid)
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
# contdid requires a balanced panel; quarterly aggregation is lossless
# for this specification since patent_year is the only time index.

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

# --- Prepare contdid panel ------------------------------------------------
# contdid requires:
#   idname  : numeric unit identifier
#   tname   : time variable
#   yname   : outcome
#   dname   : dose variable — 0 in pre-period, d_ex_ante post-2014
#   gname   : first period of treatment (2014 for all units)
#
# Since no class has d_ex_ante = 0, the lowest-dose tercile is used
# as the implicit comparison group by contdid internally.

prep_contdid <- function(df, tvar) {
  df |>
    mutate(
      unit_id = as.integer(as.factor(paste(firm, atc_level1, sep = "_"))),
      # Dose turns on in 2014: 0 before, d_ex_ante after
      dose    = if_else(.data[[tvar]] >= 2014, d_ex_ante, 0),
      # All units first treated in 2014
      G       = 2014L
    )
}

patents_cd <- prep_contdid(patents_year, "patent_year")
fda_cd     <- prep_contdid(fda_year,    "approval_year")

# --- contdid estimation ---------------------------------------------------

# Patents
patent_contdid <- att_d(
  yname            = "patents",
  tname            = "patent_year",
  idname           = "unit_id",
  gname            = "G",
  dname            = "dose",
  data             = patents_cd,
  target_parameter = "level"   # ATT level effect; use "slope" for ACRT
)

patent_contdid_es <- aggte(patent_contdid, type = "dynamic")

# FDA approvals
fda_contdid <- att_d(
  yname            = "approvals",
  tname            = "approval_year",
  idname           = "unit_id",
  gname            = "G",
  dname            = "dose",
  data             = fda_cd,
  target_parameter = "level"
)

fda_contdid_es <- aggte(fda_contdid, type = "dynamic")

# --- Plot helper ----------------------------------------------------------

save_contdid_plot <- function(es_obj, file, title_str) {
  df <- data.frame(
    t        = es_obj$egt,
    att      = es_obj$att.egt,
    se       = es_obj$se.egt,
    crit_val = es_obj$crit.val.egt
  ) |>
    mutate(
      ci_lo = att - crit_val * se,
      ci_hi = att + crit_val * se
    )

  p <- ggplot(df, aes(x = t, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, fill = "#404040") +
    geom_line(colour = "#404040") +
    geom_point(colour = "#404040", size = 1.5) +
    labs(
      x     = "Years relative to 2014",
      y     = "ATT (relative to 2013)",
      title = title_str
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(size = 11),
      axis.line    = element_line(colour = "grey40"),
      axis.ticks   = element_line(colour = "grey40"),
      plot.margin  = margin(8, 12, 8, 8)
    )

  ggsave(here("output", "robustness", file),
         plot = p, width = 8, height = 4, dpi = 200)
}

# --- Save plots -----------------------------------------------------------

save_contdid_plot(patent_contdid_es,
                  "contdid_patent_es.png",
                  "contdid event study: patents (ex-ante Medicaid exposure)")

save_contdid_plot(fda_contdid_es,
                  "contdid_fda_es.png",
                  "contdid event study: FDA approvals (ex-ante Medicaid exposure)")

# --- End of file ----------------------------------------------------------
