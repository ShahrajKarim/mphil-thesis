# regressions/RQ1.r
#
# Final RQ1 regressions: aggregate effect of MedicaidShare x Year on
# innovation outcomes (FDA approvals, patents). Event-study specification
# with firm x class and year fixed effects, clustered at the firm level.
#
# Outputs (regressions/output/RQ1/):
#   - fda_models.tex,    patent_models.tex          (modelsummary tables)
#   - fda_feols_exante.png, fda_feols_expost.png     (event-study iplots)
#   - fda_poi_exante.png,   fda_poi_expost.png
#   - patent_feols_exante.png, patent_feols_expost.png
#   - patent_poi_exante.png,   patent_poi_expost.png

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

save_iplot <- function(model, file, main) {
  png(here("output", "RQ1", file), width = 1600, height = 1000, res = 200)
  iplot(
    model,
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

fda_rq1_feols_exante <- feols(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) |
              firm^atc_level1 + approval_year,
  data = fda_ex_ante_data, cluster = ~ firm
)

fda_rq1_feols_expost <- feols(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) |
              firm^atc_level1 + approval_year,
  data = fda_ex_post_data, cluster = ~ firm
)

fda_rq1_poi_exante <- fepois(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) |
              firm^atc_level1 + approval_year,
  data = fda_ex_ante_data, cluster = ~ firm
)

fda_rq1_poi_expost <- fepois(
  approvals ~ i(approval_year, medicaid_market_share, ref = 2013) |
              firm^atc_level1 + approval_year,
  data = fda_ex_post_data, cluster = ~ firm
)

fda_rq1_model_list <- list(
  "FEOLS: Ex Ante"  = fda_rq1_feols_exante,
  "FEOLS: Ex Post"  = fda_rq1_feols_expost,
  "FEPOIS: Ex Ante" = fda_rq1_poi_exante,
  "FEPOIS: Ex Post" = fda_rq1_poi_expost
)

modelsummary(
  fda_rq1_model_list,
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  output = here("output", "RQ1", "fda_models.tex")
)

save_iplot(fda_rq1_feols_exante, "fda_feols_exante.png",
           "FEOLS: FDA approvals (ex-ante exposure)")
save_iplot(fda_rq1_feols_expost, "fda_feols_expost.png",
           "FEOLS: FDA approvals (ex-post exposure)")
save_iplot(fda_rq1_poi_exante,   "fda_poi_exante.png",
           "FE Poisson: FDA approvals (ex-ante exposure)")
save_iplot(fda_rq1_poi_expost,   "fda_poi_expost.png",
           "FE Poisson: FDA approvals (ex-post exposure)")

# --- Patent regressions ---------------------------------------------------

patents_ex_ante_data <- patents |>
  select(-d_ex_post) |>
  rename(medicaid_market_share = d_ex_ante)

patents_ex_post_data <- patents |>
  select(-d_ex_ante) |>
  rename(medicaid_market_share = d_ex_post)

patent_rq1_feols_exante <- feols(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) |
            firm^atc_level1 + patent_year,
  data = patents_ex_ante_data, cluster = ~ firm
)

patent_rq1_feols_expost <- feols(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) |
            firm^atc_level1 + patent_year,
  data = patents_ex_post_data, cluster = ~ firm
)

patent_rq1_poi_exante <- fepois(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) |
            firm^atc_level1 + patent_year,
  data = patents_ex_ante_data, cluster = ~ firm
)

patent_rq1_poi_expost <- fepois(
  patents ~ i(patent_year, medicaid_market_share, ref = 2013) |
            firm^atc_level1 + patent_year,
  data = patents_ex_post_data, cluster = ~ firm
)

patent_rq1_model_list <- list(
  "FEOLS: Ex Ante"  = patent_rq1_feols_exante,
  "FEOLS: Ex Post"  = patent_rq1_feols_expost,
  "FEPOIS: Ex Ante" = patent_rq1_poi_exante,
  "FEPOIS: Ex Post" = patent_rq1_poi_expost
)

modelsummary(
  patent_rq1_model_list,
  stars  = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  output = here("output", "RQ1", "patent_models.tex")
)

save_iplot(patent_rq1_feols_exante, "patent_feols_exante.png",
           "FEOLS: Patents (ex-ante exposure)")
save_iplot(patent_rq1_feols_expost, "patent_feols_expost.png",
           "FEOLS: Patents (ex-post exposure)")
save_iplot(patent_rq1_poi_exante,   "patent_poi_exante.png",
           "FE Poisson: Patents (ex-ante exposure)")
save_iplot(patent_rq1_poi_expost,   "patent_poi_expost.png",
           "FE Poisson: Patents (ex-post exposure)")

# --- End of file --- #