# regressions/pdl_counterfactual_analysis.r
#
# PDL counterfactual analysis. Decomposes the post-2014 innovation response
# by where a firm-class pair sits in the PDL selection process:
#
#   beta_1 (pt_stage_pass x Post)    — passed P&T review only (clinical quality
#                                      without price competitiveness)
#   beta_2 (price_stage_pass x Post) — passed price review only (cheap but not
#                                      clinically preferred)
#   beta_3 (in_pdl x Post)           — passed both stages (on the PDL)
#
# Since in_pdl ≈ pt_stage_pass × price_stage_pass, the specification is
# effectively a triple-difference that isolates the innovation response to each
# counterfactual PDL outcome, relative to the baseline of failing both stages.
#
# The regression is estimated separately for patents and FDA approvals;
# results are collected into a single combined table for comparison.
#
# Outputs (output/PDL_counterfactual/):
#   - pdl_counterfactual_models.tex   (combined patents + FDA table)

# --- Libraries ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readr)
library(fixest)
library(modelsummary)
library(here)

# --- Load data ------------------------------------------------------------

fda_approvals <- read_csv(
  here("processed_data/fda_approvals/fda_sdud_pdl_merged.csv")
) |>
  rename(approvals = n_approvals, firm = sponsor_name) |>
  filter(approval_year >= 2010)

patents <- read_csv(
  here("processed_data/patents/patents_sdud_pdl_merged.csv")
) |>
  rename(atc_level1 = ATC1, patents = num_patents) |>
  filter(patent_year >= 2010)

# --- Construct post indicator ---------------------------------------------

patents <- patents |>
  mutate(post = as.integer(patent_year >= 2014))

fda_approvals <- fda_approvals |>
  mutate(post = as.integer(approval_year >= 2014))

# --- Collapse patents to yearly panel -------------------------------------
# Avoids memory exhaustion when interacting post with three continuous
# variables; patent_year is the only time index in the model so
# aggregation is lossless for this specification.

patents_year <- patents |>
  group_by(firm, atc_level1, patent_year, post) |>
  summarise(
    patents        = sum(patents,        na.rm = TRUE),
    pt_stage_pass  = max(pt_stage_pass,  na.rm = TRUE),
    price_stage_pass = max(price_stage_pass, na.rm = TRUE),
    in_pdl         = max(in_pdl,         na.rm = TRUE),
    .groups        = "drop"
  )

# FDA approvals are already at the firm-class-year level in the source;
# collapse to be safe and to keep structure consistent.

fda_year <- fda_approvals |>
  group_by(firm, atc_level1, approval_year, post) |>
  summarise(
    approvals      = sum(approvals,      na.rm = TRUE),
    pt_stage_pass  = max(pt_stage_pass,  na.rm = TRUE),
    price_stage_pass = max(price_stage_pass, na.rm = TRUE),
    in_pdl         = max(in_pdl,         na.rm = TRUE),
    .groups        = "drop"
  )

# --- Regressions ----------------------------------------------------------
#
# Specification (same for both outcomes):
#
#   Y_{f,c,t} = alpha_{f,c} + delta_t
#             + beta_1 (pt_stage_pass_{f,c} x post_t)
#             + beta_2 (price_stage_pass_{f,c} x post_t)
#             + beta_3 (in_pdl_{f,c} x post_t)
#             + eps_{f,c,t}
#
# Fixed effects: firm x class (alpha_{f,c}), year (delta_t).
# Standard errors clustered at the firm level.

patent_pdl_cf <- feols(
  patents ~ post:pt_stage_pass + post:price_stage_pass + post:in_pdl |
            firm^atc_level1 + patent_year,
  data    = patents_year,
  cluster = ~ firm
)

fda_pdl_cf <- feols(
  approvals ~ post:pt_stage_pass + post:price_stage_pass + post:in_pdl |
              firm^atc_level1 + approval_year,
  data    = fda_year,
  cluster = ~ firm
)

# --- Combined coefficient table -------------------------------------------
#
# Coefficient labels map to the three counterfactual firm types:
#   pass P&T only | pass price only | pass both (on PDL)

coef_map <- c(
  "post:pt_stage_pass"    = "Post $\\times$ P\\&T Pass (Clinical Only)",
  "post:price_stage_pass" = "Post $\\times$ Price Pass (Price Only)",
  "post:in_pdl"           = "Post $\\times$ In PDL (Both Stages)"
)

modelsummary(
  list(
    "Patents"       = patent_pdl_cf,
    "FDA Approvals" = fda_pdl_cf
  ),
  coef_map = coef_map,
  escape   = FALSE,
  stars    = c(`*` = 0.05, `**` = 0.01, `***` = 0.001),
  gof_map  = tribble(
    ~raw,        ~clean,       ~fmt,
    "nobs",      "Num. obs.",  0,
    "r.squared", "R2",         3
  ),
  output = here("output", "PDL_counterfactual", "pdl_counterfactual_models.tex")
)

# --- End of file ----------------------------------------------------------
