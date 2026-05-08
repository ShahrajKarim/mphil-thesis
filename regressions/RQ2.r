# regressions/RQ2.r
#
# Final RQ2 regressions: PDL counterfactual analysis. Adds a year x in_pdl
# event-study term to the RQ1 specification to isolate the PDL-winner
# advantage. Firm x class and year fixed effects, clustered at the firm
# level.
#
# Outputs (regressions/output/RQ2/):
#   - fda_models.tex,    patent_models.tex          (modelsummary tables)
#   - fda_feols_{exante,expost}_{medicaid,pdl}.png  (event-study iplots,
#     one per i() term)
#   - patent_feols_{exante,expost}_{medicaid,pdl}.png

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(fixest)
library(knitr)
library(modelsummary)
library(kableExtra)
library(here)

save_iplot <- function(model, file, main, i.select = 1) {
  png(here("output", "RQ2", file), width = 1600, height = 1000, res = 200)
  iplot(
    model,
    i.select = i.select,
    ci.level = 0.95,
    ref.line = 0,
    xlab = "Year",
    ylab = "Coefficient (relative to 2013)",
    main = main
  )
  dev.off()
}

# --- Load data ------------------------------------------------------------

fda_approvals <- read_csv(
  here("processed_data/fda_approvals/fda_sdud_pdl_merged.csv")
) |>
  rename(approvals = n_approvals, firm = sponsor_name)

patents <- read_csv(
  here("processed_data/patents/patents_sdud_pdl_merged.csv")
) |>
  rename(atc_level1 = ATC1, patents = num_patents)

# --- Build ex-ante demand --- #

fda_ex_ante_stats <- fda_approvals |>
  filter(approval_year <= 2013 & approval_year >= 2010) |>
  distinct(atc_level1, approval_year, approval_quarter, prescriptions)

fda_ex_ante_total_prescriptions <- fda_ex_ante_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

fda_ex_ante_demand <- fda_ex_ante_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_ante = sum(prescriptions, na.rm = TRUE) /
                        fda_ex_ante_total_prescriptions,
            .groups = "drop")

patents_ex_ante_stats <- patents |>
  filter(patent_year <= 2013 & patent_year >= 2010) |>
  distinct(atc_level1, patent_year, patent_quarter, prescriptions)

patents_ex_ante_total_prescriptions <- patents_ex_ante_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

patents_ex_ante_demand <- patents_ex_ante_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_ante = sum(prescriptions, na.rm = TRUE) /
                        patents_ex_ante_total_prescriptions,
            .groups = "drop")

# --- Build ex-post demand --- #

fda_ex_post_stats <- fda_approvals |>
  filter(approval_year %in% 2014:2015) |>
  distinct(atc_level1, approval_year, approval_quarter, prescriptions)

fda_ex_post_total_prescriptions <- fda_ex_post_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

fda_ex_post_demand <- fda_ex_post_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_post = sum(prescriptions, na.rm = TRUE) /
                        fda_ex_post_total_prescriptions,
            .groups = "drop")

patents_ex_post_stats <- patents |>
  filter(patent_year %in% 2014:2015) |>
  distinct(atc_level1, patent_year, patent_quarter, prescriptions)

patents_ex_post_total_prescriptions <- patents_ex_post_stats |>
  summarise(total_prescriptions = sum(prescriptions, na.rm = TRUE),
            .groups = "drop") |>
  pull()

patents_ex_post_demand <- patents_ex_post_stats |>
  group_by(atc_level1) |>
  summarise(d_ex_post = sum(prescriptions, na.rm = TRUE) /
                        patents_ex_post_total_prescriptions,
            .groups = "drop")

# --- Merge demand measures --- #

fda_approvals <- fda_approvals |>
  left_join(fda_ex_ante_demand, by = "atc_level1") |>
  left_join(fda_ex_post_demand, by = "atc_level1") |>
  filter(approval_year >= 2010)

patents <- patents |>
  left_join(patents_ex_ante_demand, by = "atc_level1") |>
  left_join(patents_ex_post_demand, by = "atc_level1") |>
  filter(patent_year >= 2010)

# --- FDA regressions ------------------------------------------------------

fda_ex_ante_data <- fda_approvals |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

fda_ex_post_data <- fda_approvals |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)

fda_rq2_feols_exante <- feols(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013)
            + i(approval_year, in_pdl, ref = 2013) |
              firm^atc_level1 + approval_year,
  data = fda_ex_ante_data, cluster = ~ firm
)

fda_rq2_feols_expost <- feols(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013)
            + i(approval_year, in_pdl, ref = 2013) |
              firm^atc_level1 + approval_year,
  data = fda_ex_post_data, cluster = ~ firm
)

fda_rq2_model_list <- list(
  "FEOLS: Ex Ante" = fda_rq2_feols_exante,
  "FEOLS: Ex Post" = fda_rq2_feols_expost
)

modelsummary(
  fda_rq2_model_list,
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  output = here("output", "RQ2", "fda_models.tex")
)

save_iplot(fda_rq2_feols_exante, "fda_feols_exante_medicaid.png",
           "Event study: Medicaid market size (ex ante)", i.select = 1)
save_iplot(fda_rq2_feols_exante, "fda_feols_exante_pdl.png",
           "Event study: PDL exposure (ex ante)", i.select = 2)
save_iplot(fda_rq2_feols_expost, "fda_feols_expost_medicaid.png",
           "Event study: Medicaid market size (ex post)", i.select = 1)
save_iplot(fda_rq2_feols_expost, "fda_feols_expost_pdl.png",
           "Event study: PDL exposure (ex post)", i.select = 2)

# --- Patent regressions ---------------------------------------------------

patents_ex_ante_data <- patents |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

patents_ex_post_data <- patents |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)

# Collapse quarterly patents to yearly before estimating: the second i() term
# in RQ2 doubles the model matrix relative to RQ1 and exhausts vector memory
# at the quarterly grain. patent_year is the only time index in the model, so
# aggregation is lossless for this specification.
patents_ex_ante_year <- patents_ex_ante_data |>
  group_by(firm, atc_level1, patent_year) |>
  summarise(
    patents = sum(patents, na.rm = TRUE),
    medicaid_market_share = first(medicaid_market_share),
    in_pdl = max(in_pdl, na.rm = TRUE),
    .groups = "drop"
  )

patents_ex_post_year <- patents_ex_post_data |>
  group_by(firm, atc_level1, patent_year) |>
  summarise(
    patents = sum(patents, na.rm = TRUE),
    medicaid_market_share = first(medicaid_market_share),
    in_pdl = max(in_pdl, na.rm = TRUE),
    .groups = "drop"
  )

patent_rq2_feols_exante <- feols(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013)
          + i(patent_year, in_pdl, ref = 2013) |
            firm^atc_level1 + patent_year,
  data = patents_ex_ante_year, cluster = ~ firm
)

patent_rq2_feols_expost <- feols(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013)
          + i(patent_year, in_pdl, ref = 2013) |
            firm^atc_level1 + patent_year,
  data = patents_ex_post_year, cluster = ~ firm
)

patent_rq2_model_list <- list(
  "FEOLS: Ex Ante" = patent_rq2_feols_exante,
  "FEOLS: Ex Post" = patent_rq2_feols_expost
)

modelsummary(
  patent_rq2_model_list,
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map = tribble(
    ~raw, ~clean, ~fmt,
    "nobs", "Num. obs.", 0,
    "r.squared", "R2", 3,
    "FE", "FE", 0
  ),
  output = here("output", "RQ2", "patent_models.tex")
)

save_iplot(patent_rq2_feols_exante, "patent_feols_exante_medicaid.png",
           "Event study: Medicaid market size (ex ante)", i.select = 1)
save_iplot(patent_rq2_feols_exante, "patent_feols_exante_pdl.png",
           "Event study: PDL exposure (ex ante)", i.select = 2)
save_iplot(patent_rq2_feols_expost, "patent_feols_expost_medicaid.png",
           "Event study: Medicaid market size (ex post)", i.select = 1)
save_iplot(patent_rq2_feols_expost, "patent_feols_expost_pdl.png",
           "Event study: PDL exposure (ex post)", i.select = 2)

# --- End of file --- #