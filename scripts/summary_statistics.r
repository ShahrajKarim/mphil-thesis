# scripts/summary_statistics.r

# Builds a source-level descriptive statistics table followed by
# summary statistics for each balanced panel.
# Outputs output/summary_statistics/summary_statistics.tex

# libraries

library(here)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(knitr)
library(kableExtra)
library(tibble)

# --- SDUD --- #

sdud <- readRDS(here("processed_data/state_drug_utilisation_data/SDUD_full.rds"))

sdud_counts <- tibble(
  Variable = c("Unique labeler codes", "Unique NDCs", "States"),
  N = c(n_distinct(sdud$labeler_code), n_distinct(sdud$ndc), n_distinct(sdud$state)),
  Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_
)

sdud_dist <- sdud |>
  mutate(
    medicaid_amount_reimbursed = medicaid_amount_reimbursed / 1000,
    non_medicaid_amount_reimbursed = non_medicaid_amount_reimbursed / 1000
  ) |>
  summarise(across(
    c(number_of_prescriptions, medicaid_amount_reimbursed, non_medicaid_amount_reimbursed),
    list(N = ~sum(!is.na(.)), Mean = ~mean(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE), Min = ~min(., na.rm = TRUE), Max = ~max(., na.rm = TRUE))
  )) |>
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(N|Mean|SD|Min|Max)"
  ) |>
  mutate(Variable = recode(Variable,
    "number_of_prescriptions" = "Prescriptions",
    "medicaid_amount_reimbursed" = "Medicaid reimbursement (\\$000s)",
    "non_medicaid_amount_reimbursed" = "Non-Medicaid reimbursement (\\$000s)"
  ))

sdud_stats <- bind_rows(sdud_counts, sdud_dist)

# --- PATENTS --- #

patents <- read_csv(
  here("processed_data/patents/pharma_patents_clean.csv"),
  show_col_types = FALSE
)

patent_counts <- tibble(
  Variable = c("Unique assignees", "Unique ATC classes"),
  N = c(n_distinct(patents$firm), n_distinct(patents$ATC1)),
  Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_
)

patent_dist <- patents |>
  filter(!is.na(patent_year), patent_year >= 2005) |>
  group_by(firm, patent_year, ATC1) |>
  summarise(num_patents = n(), .groups = "drop") |>
  summarise(
    Variable = "Patents per firm-year-ATC",
    N = n(),
    Mean = mean(num_patents, na.rm = TRUE),
    SD = sd(num_patents, na.rm = TRUE),
    Min = min(num_patents, na.rm = TRUE),
    Max = max(num_patents, na.rm = TRUE)
  )

patent_stats <- bind_rows(patent_counts, patent_dist)

# --- FDA APPROVALS --- #

fda_mapped <- read_csv(
  here("processed_data/fda_approvals/fda_approvals_atc_firm_mapped.csv"),
  show_col_types = FALSE,
  col_types = cols(labeler_code = col_character())
)

fda_counts <- tibble(
  Variable = c("Unique sponsors", "Unique ATC classes"),
  N = c(n_distinct(fda_mapped$sponsor_name), n_distinct(fda_mapped$atc_level1)),
  Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_
)

fda_dist <- fda_mapped |>
  filter(!is.na(approval_year)) |>
  group_by(sponsor_name, approval_year, atc_level1) |>
  summarise(n_approvals = n(), .groups = "drop") |>
  summarise(
    Variable = "Approvals per sponsor-year-ATC",
    N = n(),
    Mean = mean(n_approvals, na.rm = TRUE),
    SD = sd(n_approvals, na.rm = TRUE),
    Min = min(n_approvals, na.rm = TRUE),
    Max = max(n_approvals, na.rm = TRUE)
  )

fda_stats <- bind_rows(fda_counts, fda_dist)

# --- TEXAS PDL --- #

pdl <- read_csv(
  here("processed_data/PDL/texas_PDL_data.csv"),
  show_col_types = FALSE
) |>
  mutate(
    preferred = as.integer(medicaid_covered == "Yes" & price_stage_pass == 1 & pt_stage_pass == 1),
    covered_not_pref = as.integer(medicaid_covered == "Yes" & (price_stage_pass == 0 | pt_stage_pass == 0)),
    not_covered = as.integer(medicaid_covered == "No")
  )

pdl_counts <- tibble(
  Variable = "Unique NDCs",
  N = n_distinct(pdl$ndc),
  Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_
)

pdl_dist <- pdl |>
  summarise(across(
    c(preferred, pt_stage_pass, price_stage_pass),
    list(N = ~sum(!is.na(.)), Mean = ~mean(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE), Min = ~min(., na.rm = TRUE), Max = ~max(., na.rm = TRUE))
  )) |>
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(N|Mean|SD|Min|Max)"
  ) |>
  mutate(Variable = recode(Variable,
    "preferred" = "Preferred (PDL)",
    "pt_stage_pass" = "PT stage pass",
    "price_stage_pass" = "Price stage pass"
  ))

pdl_stats <- bind_rows(pdl_counts, pdl_dist)

# --- PATENTS PANEL --- #

patents_panel <- read_csv(
  here("processed_data/patents/patents_sdud_pdl_merged.csv"),
  show_col_types = FALSE
) |>
  mutate(medicaid_share = ifelse(
    ATC1_medicaid_reimbursed + ATC1_non_medicaid_reimbursed > 0,
    ATC1_medicaid_reimbursed / (ATC1_medicaid_reimbursed + ATC1_non_medicaid_reimbursed),
    NA_real_
  ))

patents_panel_counts <- tibble(
  Variable = c("Unique firms", "Unique ATC classes", "Years covered"),
  N = c(n_distinct(patents_panel$firm), n_distinct(patents_panel$ATC1), n_distinct(patents_panel$patent_year)),
  Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_
)

patents_panel_dist <- patents_panel |>
  summarise(across(
    c(num_patents, in_pdl, medicaid_share),
    list(N = ~sum(!is.na(.)), Mean = ~mean(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE), Min = ~min(., na.rm = TRUE), Max = ~max(., na.rm = TRUE))
  )) |>
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(N|Mean|SD|Min|Max)"
  ) |>
  mutate(Variable = recode(Variable,
    "num_patents" = "Patents",
    "in_pdl" = "PDL exposure",
    "medicaid_share" = "Medicaid share (class-level)"
  ))

patents_panel_stats <- bind_rows(patents_panel_counts, patents_panel_dist)

# --- FDA PANEL --- #

fda_panel <- read_csv(
  here("processed_data/fda_approvals/fda_sdud_pdl_merged.csv"),
  show_col_types = FALSE
) |>
  mutate(medicaid_share = ifelse(
    ATC1_medicaid_reimbursed + ATC1_non_medicaid_reimbursed > 0,
    ATC1_medicaid_reimbursed / (ATC1_medicaid_reimbursed + ATC1_non_medicaid_reimbursed),
    NA_real_
  ))

fda_panel_counts <- tibble(
  Variable = c("Unique firms", "Unique ATC classes", "Years covered"),
  N = c(n_distinct(fda_panel$sponsor_name), n_distinct(fda_panel$atc_level1), n_distinct(fda_panel$approval_year)),
  Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_
)

fda_panel_dist <- fda_panel |>
  summarise(across(
    c(n_approvals, in_pdl, medicaid_share),
    list(N = ~sum(!is.na(.)), Mean = ~mean(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE), Min = ~min(., na.rm = TRUE), Max = ~max(., na.rm = TRUE))
  )) |>
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(N|Mean|SD|Min|Max)"
  ) |>
  mutate(Variable = recode(Variable,
    "n_approvals" = "FDA approvals",
    "in_pdl" = "PDL exposure",
    "medicaid_share" = "Medicaid share (class-level)"
  ))

fda_panel_stats <- bind_rows(fda_panel_counts, fda_panel_dist)

# --- Build source-level table --- #

table_format <- function(df) {
  df |> mutate(
    N = format(as.integer(N), big.mark = ","),
    Mean = ifelse(is.na(Mean), "-", format(round(Mean, 2), nsmall = 2, big.mark = ",")),
    SD = ifelse(is.na(SD), "-", format(round(SD, 2), nsmall = 2, big.mark = ",")),
    Min = ifelse(is.na(Min), "-", format(round(Min, 2), big.mark = ",")),
    Max = ifelse(is.na(Max), "-", format(round(Max, 2), big.mark = ","))
  )
}

source_combined <- bind_rows(
  sdud_stats, patent_stats, fda_stats, pdl_stats
  ) |> table_format()

s1 <- 1; e1 <- nrow(sdud_stats)
s2 <- e1 + 1; e2 <- e1 + nrow(patent_stats)
s3 <- e2 + 1; e3 <- e2 + nrow(fda_stats)
s4 <- e3 + 1; e4 <- e3 + nrow(pdl_stats)

tbl_source <- source_combined |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    col.names = c("Variable", "N", "Mean", "SD", "Min", "Max"),
    align = c("l", "r", "r", "r", "r", "r"),
    escape = FALSE,
    caption = "Source-level descriptive statistics",
    label = "source_descriptives"
  ) |>
  kable_styling(latex_options = c("hold_position", "scale_down")) |>
  pack_rows("State Drug Utilisation Data", s1, e1, bold = FALSE, italic = TRUE, latex_gap_space = "0.5em") |>
  pack_rows("USPTO PatentsView", s2, e2, bold = FALSE, italic = TRUE, latex_gap_space = "0.5em") |>
  pack_rows("Drugs@FDA", s3, e3, bold = FALSE, italic = TRUE, latex_gap_space = "0.5em") |>
  pack_rows("Texas Vendor Drug Program", s4, e4, bold = FALSE, italic = TRUE, latex_gap_space = "0.5em")

writeLines(as.character(tbl_source), here("output/summary_statistics/summary_statistics.tex"))

# --- Build panel table --- #

panel_combined <- bind_rows(
  patents_panel_stats, fda_panel_stats
  ) |> table_format()

p1 <- 1; q1 <- nrow(patents_panel_stats)
p2 <- q1 + 1; q2 <- q1 + nrow(fda_panel_stats)

tbl_panel <- panel_combined |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    col.names = c("Variable", "N", "Mean", "SD", "Min", "Max"),
    align = c("l", "r", "r", "r", "r", "r"),
    escape = FALSE,
    caption = "Balanced panel summary statistics",
    label = "panel_descriptives"
  ) |>
  kable_styling(latex_options = c("hold_position", "scale_down")) |>
  pack_rows("SDUD Patents panel", p1, q1, bold = FALSE, italic = TRUE, latex_gap_space = "0.5em", escape = FALSE) |>
  pack_rows("SDUD FDA approvals panel", p2, q2, bold = FALSE, italic = TRUE, latex_gap_space = "0.5em", escape = FALSE)

writeLines(as.character(tbl_panel), here("output/summary_statistics/panel_statistics.tex"))