# robustness/patent_matching.r
#
# Robustness: assess whether the low patent-assignee to NDC firm match
# rate (27.1%) threatens the analysis. The match rate is measured at the
# assignee level, but the analysis is patent-weighted, so the question is
# whether the matched assignees account for a disproportionately large
# share of patents and whether the unmatched group differs systematically
# on assignee type or therapeutic footprint.
#
# Inputs (saved by scripts/patent_data.r — script not modified):
#   - processed_data/patents/pharma_patents_clean.csv
#   - aux_data/patent_firm_to_ndc_firm_mapping.csv
#
# Outputs (output/robustness/patent_matching/):
#   - patent_matching_summary.tex   (headline comparison table)
#   - patent_matching_atc.tex       (ATC1 share comparison)

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(here)

# --- Load saved artefacts -------------------------------------------------

patents <- read_csv(
  here("processed_data/patents/pharma_patents_clean.csv"),
  col_types = cols(patent_id = col_character())
)

firm_map <- read_csv(
  here("aux_data/patent_firm_to_ndc_firm_mapping.csv")
)

# --- Flag matched vs unmatched assignees ----------------------------------
# pharma_patents_clean.csv is exploded across cpc_group and ATC1 mappings,
# so patent_id is repeated. n_distinct(patent_id) is the correct count.

patents <- patents |>
  mutate(matched = firm %in% firm_map$firm)

# --- Pre/post-2014 patent share -------------------------------------------

period_stats <- patents |>
  distinct(patent_id, patent_year, matched) |>
  filter(patent_year >= 2005, patent_year <= 2024) |>
  mutate(period = if_else(patent_year < 2014, "Pre-2014", "Post-2014")) |>
  count(matched, period, name = "n_patents") |>
  group_by(matched) |>
  mutate(share = n_patents / sum(n_patents)) |>
  ungroup()

# --- ATC1 distribution ----------------------------------------------------

atc_stats <- patents |>
  distinct(patent_id, ATC1, matched) |>
  filter(!is.na(ATC1)) |>
  count(matched, ATC1, name = "n_patents") |>
  group_by(matched) |>
  mutate(share = n_patents / sum(n_patents)) |>
  ungroup()

# --- Build headline comparison table --------------------------------------
# Single booktabs tabular matching the style of class_substitutability.tex.

pct  <- function(x) sprintf("%.1f\\%%", 100 * x)

period_share <- function(group, m) {
  v <- period_stats$share[period_stats$matched == m & period_stats$period == group]
  if (length(v) == 0) 0 else v
}

rows <- tribble(
  ~stat, ~matched, ~unmatched,
  "Pre-2014 patent share",
    pct(period_share("Pre-2014", TRUE)),
    pct(period_share("Pre-2014", FALSE)),
  "Post-2014 patent share",
    pct(period_share("Post-2014", TRUE)),
    pct(period_share("Post-2014", FALSE))
)

lines <- c(
  "\\begin{table}[!h]",
  "\\centering",
  "\\caption{\\label{tab:patent_matching_summary}Pre-2014 patent share: matched vs unmatched assignees}",
  "\\centering",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rr}",
  "\\toprule",
  " & Matched & Unmatched \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(rows))) {
  lines <- c(lines,
    sprintf("%s & %s & %s \\\\", rows$stat[i], rows$matched[i], rows$unmatched[i])
  )
}
lines <- c(lines, "\\bottomrule", "\\end{tabular*}", "\\end{table}", "")

writeLines(
  lines,
  here("output", "robustness", "patent_matching", "patent_matching_summary.tex")
)

# --- Build ATC distribution table -----------------------------------------

atc_wide <- atc_stats |>
  select(ATC1, matched, share) |>
  pivot_wider(names_from = matched, values_from = share, names_prefix = "m_") |>
  mutate(across(starts_with("m_"), \(x) replace_na(x, 0))) |>
  arrange(ATC1)

atc_lines <- c(
  "\\begin{table}[!h]",
  "\\centering",
  "\\caption{\\label{tab:patent_matching_atc}ATC1 share of A61P patents: matched vs unmatched}",
  "\\centering",
  "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}rr}",
  "\\toprule",
  "ATC1 & Matched & Unmatched \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(atc_wide))) {
  atc_lines <- c(atc_lines,
    sprintf("%s & %s & %s \\\\",
            atc_wide$ATC1[i],
            pct(atc_wide$m_TRUE[i]),
            pct(atc_wide$m_FALSE[i]))
  )
}
atc_lines <- c(atc_lines, "\\bottomrule", "\\end{tabular*}", "\\end{table}", "")

writeLines(
  atc_lines,
  here("output", "robustness", "patent_matching", "patent_matching_atc.tex")
)

# --- End of file ----------------------------------------------------------
