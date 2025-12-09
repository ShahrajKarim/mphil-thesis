# scripts/clinical_trials_data.r

# libraries
library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(purrr)
library(jsonlite)
library(httr)
library(qs)

# Read Clinical Trials data
ctg_raw <- read_csv("raw_data/clinical_trials/ctg-studies.csv",
                    show_col_types = FALSE) |>
  clean_names()

# --- Clean Clinical Trials Data --- #

# Keep interventional studies
ctg_cleaned <- ctg_raw |>
  filter(!is.na(interventions)) |>
  filter(study_type == "INTERVENTIONAL")

# Extract drug/biological names from the 'interventions' text field
ctg_cleaned <- ctg_cleaned |>
  mutate(
    extracted_names = str_extract_all(
      interventions,
      pattern = "(DRUG|BIOLOGICAL):\\s*[^\\|]+"
    )
  )

# Begin rowwise operations and remove DRUG:/BIOLOGICAL: prefixes
ctg_cleaned <- ctg_cleaned |>
  rowwise() |>
  mutate(
    clean_list = list(
      str_remove(extracted_names, pattern = "^(DRUG|BIOLOGICAL):\\s*")
    )
  )

# Filter out placebo, comparator, saline, procedure, and non-drug entries
ctg_cleaned <- ctg_cleaned |>
  mutate(
    clean_list = list(
      str_subset(
        clean_list,
        pattern = "(?i)(placebo|comparator|unspecified|saline|vehicle|procedure|radiation)",
        negate = TRUE
      )
    )
  )

# Split combination therapies into individual components (+, /, and, parentheses)
ctg_cleaned <- ctg_cleaned |>
  mutate(
    split_ingredients = list(
      clean_list |>
        map(~ unlist(str_split(.x, pattern = "\\(|\\+|\\/| and "))) |>
        unlist() |>
        str_squish()
    )
  )


ctg_cleaned <- ctg_cleaned |>
  mutate(
    final_ingredients = list(
      split_ingredients |>
        str_subset(pattern = "^$|^\\s*$", negate = TRUE) |>
        str_subset(pattern = "(?i)(placebo|comparator|unspecified|saline|vehicle)", negate = TRUE) |>
        unique()
    ),

    clean_drug_names = ifelse(
      length(final_ingredients) == 0,
      NA_character_,
      paste(final_ingredients, collapse = ", ")
    )
  )

# Ungroup and drop intermediate helper columns
ctg_cleaned <- ctg_cleaned |>
  ungroup() |>
  select(
    -extracted_names,
    -clean_list,
    -split_ingredients,
    -final_ingredients
  ) |>
  # Keep only trials with at least one cleaned drug name extracted
  filter(!is.na(clean_drug_names))


# Check the result with a few examples
print(head(select(ctg_cleaned, interventions, clean_drug_names), 5))

# save cleaned data
write.csv(ctg_cleaned, "processed_data/clinical_trials/clinical_trials_drug_cleaned.csv", row.names = FALSE)

ctg_cleaned <- read.csv(
  "processed_data/clinical_trials/clinical_trials_drug_cleaned.csv"
)

