# scripts/patent_data.r

# libraries
library(vroom)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(lubridate)

# File paths for patent data

dir <- "raw_data/patents/"

path_patent <- paste0(dir, "g_patent.tsv")
path_cpc <- paste0(dir, "g_cpc_current.tsv")
path_assignee <- paste0(dir, "g_assignee_disambiguated.tsv")
path_title <- paste0(dir, "g_cpc_title.tsv")
path_tech_field <- paste0(dir, "g_wipo_technology.tsv")

# path_abstract <- paste0(dir, "g_patent_abstract.tsv") # (use abstracts using NLP if required)


# Load base patent data

patent <- vroom(
  path_patent,
  col_select = c(
    patent_id,
    patent_date,
    patent_title
  )
)

# Load CPC classifications

cpc <- vroom(
  path_cpc,
  col_select = c(
    patent_id,
    cpc_subclass,
    cpc_group
  )
) |>
  mutate(patent_id = as.character(patent_id))

# Load assignee data

assignee <- vroom(
  path_assignee,
  col_select = c(
    patent_id,
    assignee_id,
    disambig_assignee_organization,
    assignee_type
  )
)

# Load cpc_class data

titles <- vroom(
  path_title,
  col_select = c(
    cpc_class,
    cpc_subclass,
    cpc_group,
    cpc_group_title,
  )
)

# Load tech field data

tech_field <- vroom(
  path_tech_field,
  col_select = c(
    patent_id,
    wipo_field_title,
  )
) |>
  filter(wipo_field_title == "Pharmaceuticals") |>
  mutate(patent_id = as.character(patent_id))

# Make a clean dataset based on these pharma fields

pharma_patents_full <- tech_field |>
  left_join(cpc, by = "patent_id", relationship = "many-to-many") |>
  left_join(titles, by = c("cpc_subclass", "cpc_group")) |>
  left_join(assignee, by = "patent_id", relationship = "many-to-many") |>
  left_join(patent, by = "patent_id") |>
  mutate(patent_year = year(patent_date)) |>
  rename(firm = disambig_assignee_organization)

# Filter for only A61P CPC subclass (therapeutic activity of drugs)
# https://www.uspto.gov/web/patents/classification/cpc/html/cpc-A61P.html

pharma_patents_full <- pharma_patents_full %>%
  filter(cpc_subclass == "A61P")

# Remove variables which are same across all rows
pharma_patents_full <- pharma_patents_full |>
select(
  -cpc_subclass,
  -cpc_class,
  -wipo_field_title  
)

# save the full pharma patents data
write_csv(
  pharma_patents_full,
  file = paste0("processed_data/patents/", "pharma_patents_full.csv")
)

