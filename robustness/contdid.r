# robustness/contdid.r
#
# Robustness: Callaway, Goodman-Bacon & Sant'Anna (2024) continuous
# treatment DiD estimator via cont_did(). Relaxes the linearity-in-dose
# assumption imposed by FEOLS by estimating ATT(d, t) using a
# nonparametric B-spline model in the dose dimension.
#
# Three key departures from the old att_d() script:
#   1. Function is cont_did(), not att_d() — the API changed in the
#      alpha release.
#   2. The dose variable must be TIME-INVARIANT across all periods.
#      The old code set dose = 0 pre-2014, which the README explicitly
#      says not to do.
#   3. A comparison group is required. Since there are no true
#      never-treated units, the lowest-dose tercile is designated as
#      the comparison group (G = 0). Middle and upper terciles are
#      treated in 2014 (G = 2014). This mirrors the approach in
#      binned_treatment.r.
#
# Install (if needed):
#   devtools::install_github("bcallaway11/contdid")
#
# Note: biters = 199 bootstrap iterations. Expect runtime of several
# minutes per outcome on a large panel.
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

# --- Collapse to yearly panel ---------------------------------------------

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

# --- Prepare cont_did panel -----------------------------------------------
# cont_did() requires:
#   - dname  : time-invariant dose; same value in ALL periods (pre and post)
#   - gname  : 0 for never-treated comparison units, 2014 for treated units
#   - Balanced panel: every unit must appear in every time period
#
# Tercile assignment mirrors binned_treatment.r (ATC class level).
# Lowest tercile = G = 0 (comparison group).
# Middle + upper terciles = G = 2014 (treated in 2014).

prep_contdid <- function(df, tvar, yvar) {
  # Terciles computed at ATC class level (13 classes)
  class_terciles <- df |>
    distinct(atc_level1, d_ex_ante) |>
    mutate(tercile = ntile(d_ex_ante, 3))

  df |>
    left_join(class_terciles |> select(atc_level1, tercile),
              by = "atc_level1") |>
    mutate(
      unit_id = as.integer(as.factor(paste(firm, atc_level1, sep = "_"))),
      dose    = d_ex_ante,                           # time-invariant dose
      G       = if_else(tercile == 1, 0L, 2014L)    # 0 = comparison group
    ) |>
    # Balance the panel: insert zero-outcome rows for missing unit-year combos
    complete(unit_id, !!sym(tvar),
             fill = setNames(list(0), yvar)) |>
    group_by(unit_id) |>
    fill(G, dose, firm, atc_level1, tercile, .direction = "downup") |>
    ungroup() |>
    as.data.frame()
}

patents_cd <- prep_contdid(patents_year, "patent_year",   "patents")
fda_cd     <- prep_contdid(fda_year,    "approval_year",  "approvals")

# --- cont_did estimation --------------------------------------------------

patent_contdid <- cont_did(
  yname            = "patents",
  tname            = "patent_year",
  idname           = "unit_id",
  dname            = "dose",
  data             = patents_cd,
  gname            = "G",
  target_parameter = "level",       # ATT-type; use "slope" for ACRT
  aggregation      = "eventstudy",
  treatment_type   = "continuous",
  control_group    = "nevertreated",
  biters           = 199,
  cband            = TRUE
)

fda_contdid <- cont_did(
  yname            = "approvals",
  tname            = "approval_year",
  idname           = "unit_id",
  dname            = "dose",
  data             = fda_cd,
  gname            = "G",
  target_parameter = "level",
  aggregation      = "eventstudy",
  treatment_type   = "continuous",
  control_group    = "nevertreated",
  biters           = 199,
  cband            = TRUE
)

# --- Extract att_gt results -----------------------------------------------
# ggcont_did() and the default aggregation rely on $egt fields which are NaN
# when there is only one cohort (G = 2014). The group-time ATTs in $att_gt
# are fine and are exactly the event study for a single-cohort design
# (event_time = t - 2014). We extract them directly.

extract_att_gt <- function(cd_res) {
  ag <- cd_res$att_gt
  data.frame(
    group = ag$group,
    time  = ag$t,
    att   = ag$att,
    se    = ag$se,
    ci_lo = ag$att - ag$crit_val * ag$se,
    ci_hi = ag$att + ag$crit_val * ag$se
  ) |>
    mutate(event_time = time - group)
}

patent_att_gt <- extract_att_gt(patent_contdid)
fda_att_gt    <- extract_att_gt(fda_contdid)

# --- Export LaTeX table ---------------------------------------------------
# Combined talltblr table matching the modelsummary house style.
# Patents and FDA side by side; SE in parentheses on second row.
# hline separates pre-treatment from post-treatment periods.

save_contdid_tex <- function(pat_df, fda_df, file) {
  fmt_num <- function(x) sprintf("\\num{%.4f}", x)

  fmt_col <- function(df) {
    df |>
      arrange(event_time) |>
      mutate(
        stars   = if_else(ci_lo > 0 | ci_hi < 0, "*", ""),
        att_fmt = paste0(fmt_num(att), stars),
        se_fmt  = paste0("(", fmt_num(se), ")")
      )
  }

  pat <- fmt_col(pat_df)
  fda <- fmt_col(fda_df)

  combined <- full_join(
    pat |> select(event_time, pat_att = att_fmt, pat_se = se_fmt),
    fda |> select(event_time, fda_att = att_fmt, fda_se = se_fmt),
    by = "event_time"
  ) |> arrange(event_time)

  # hline{n} sits above row n. Header = row 1; each event time = 2 rows.
  # Place separator between last pre-treatment SE row and first post row.
  n_pre      <- sum(combined$event_time < 0)
  hline_pos  <- 1 + n_pre * 2 + 1

  lines <- c(
    "\\begin{table}",
    "\\centering",
    "\\begin{talltblr}[         %% tabularray outer open",
    "entry=none,label=none,",
    "note{}={* 95\\% uniform confidence band does not cover 0},",
    "]                     %% tabularray outer close",
    "{                     %% tabularray inner open",
    "colspec={Q[]Q[]Q[]},",
    "column{2,3}={}{halign=c,},",
    "column{1}={}{halign=l,},",
    sprintf("hline{%d}={1,2,3}{solid, black, 0.05em},", hline_pos),
    "}                     %% tabularray inner close",
    "\\toprule",
    "Event Time & Patents & FDA Approvals \\\\ \\midrule %% TinyTableHeader"
  )

  for (i in seq_len(nrow(combined))) {
    lines <- c(lines,
      sprintf("%d & %s & %s \\\\",
              combined$event_time[i],
              combined$pat_att[i], combined$fda_att[i]),
      sprintf("   & %s & %s \\\\",
              combined$pat_se[i],  combined$fda_se[i])
    )
  }

  lines <- c(lines,
    "\\bottomrule",
    "\\end{talltblr}",
    "\\end{table}"
  )

  writeLines(lines, here("output", "robustness", "contdid", file))
}

save_contdid_tex(patent_att_gt, fda_att_gt, "contdid_es.tex")

# --- Plot helper ----------------------------------------------------------

save_contdid_plot <- function(df, file, title_str) {
  p <- ggplot(df, aes(x = event_time, y = att)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", colour = "grey50") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.2, colour = "#404040") +
    geom_point(colour = "#404040", size = 1.5) +
    labs(
      x     = "Years relative to 2014",
      y     = "ATT (relative to 2013)",
      title = title_str
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.margin  = margin(8, 12, 8, 8)
    )

  ggsave(here("output", "robustness", "contdid", file),
         plot = p, width = 8, height = 4, dpi = 200)
}

# --- Save plots -----------------------------------------------------------

save_contdid_plot(patent_att_gt,
                  "contdid_patent_es.png",
                  "Continuous Treatment DiD: Patents")

save_contdid_plot(fda_att_gt,
                  "contdid_fda_es.png",
                  "Continuous Treatment DiD: FDA Approvals")

# --- End of file ----------------------------------------------------------
