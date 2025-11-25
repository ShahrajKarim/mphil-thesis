# scripts/patent_data.r

# ================================================================
# patent_data.r
# Full pipeline to build:
#   - Patent × CPC × firm × disease mapping
#   - Innovation panel: year × disease_area × firm
# ================================================================

library(vroom)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(lubridate)

# ----------------------------
# 1. File paths
# ----------------------------
dir <- "raw_data/patents/"

path_patent <- paste0(dir, "g_patent.tsv")
path_cpc <- paste0(dir, "g_cpc_current.tsv")
path_assignee <- paste0(dir, "g_assignee_disambiguated.tsv")
path_title <- paste0(dir, "g_cpc_title.tsv")
path_tech_field <- paste0(dir, "g_wipo_technology.tsv")
# path_abstract <- paste0(dir, "g_patent_abstract.tsv")


# ----------------------------
# 2. Load base patent data
# ----------------------------
patent <- vroom(
  path_patent,
  col_select = c(
    patent_id,
    patent_date,
    patent_title
  )
)

# ----------------------------
# 3. Load CPC classifications
# ----------------------------
cpc <- vroom(
  path_cpc,
  col_select = c(
    patent_id,
    cpc_subclass,
    cpc_group
  )
) |>
  mutate(patent_id = as.character(patent_id))

# ----------------------------
# 4. Load assignee data
# ----------------------------
assignee <- vroom(
  path_assignee,
  col_select = c(
    patent_id,
    assignee_id,
    disambig_assignee_organization,
    assignee_type
  )
)

# ----------------------------
# 5. Load cpc_class data
# ----------------------------
titles <- vroom(
  path_title,
  col_select = c(
    cpc_class,
    cpc_subclass,
    cpc_group,
    cpc_group_title,
  )
)

# ----------------------------
# 6. Load tech field data
# ----------------------------

pharma_fields <- c("Pharmaceuticals")

tech_field <- vroom(
  path_tech_field,
  col_select = c(
    patent_id,
    wipo_field_title,
  )
) |>
  filter(wipo_field_title %in% pharma_fields) |>
  mutate(patent_id = as.character(patent_id))

# ================================================================
# 7. Make a clean dataset based on these pharma fields
# ================================================================

pharma_patents_full <- tech_field |>
  left_join(cpc, by = "patent_id", relationship = "many-to-many") |>
  left_join(titles, by = c("cpc_subclass", "cpc_group")) |>
  left_join(assignee, by = "patent_id", relationship = "many-to-many") |>
  left_join(patent, by = "patent_id") |>
  mutate(
    patent_year = year(patent_date),
    firm = case_when(
      assignee_type == 2 ~ disambig_assignee_organization,
      TRUE ~ NA_character_
    )
  )

# Filter for only A61P CPC subclass (therapeutic activity of drugs)
# https://www.uspto.gov/web/patents/classification/cpc/html/cpc-A61P.html

pharma_patents_full <- pharma_patents_full %>%
  filter(cpc_subclass == "A61P")

# ================================================================
# 8. Find a method to clean CPC group titles! - TBC 
# ================================================================


