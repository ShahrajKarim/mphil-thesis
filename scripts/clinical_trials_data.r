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

# --- RxClass API Mapping (Draft) -- #

unique_ings <- ctg_cleaned |>
  select(clean_drug_names) |>
  separate_rows(clean_drug_names, sep = ",\\s*") |>
  pull(clean_drug_names) |>
  unique() |>
  str_squish() |> # Final cleanup
  str_subset(pattern = "^$|^\\s*$", negate = TRUE)


# RxClass lookup function

rxclass_lookup <- function(drug_name) {
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxclass/class/byDrugName.json?drugName=",
    URLencode(drug_name)
  )

  res <- tryCatch(httr::GET(url), error = function(e) NULL)
  if (is.null(res) || res$status_code != 200) return(NULL)

  json <- tryCatch(
    jsonlite::fromJSON(content(res, "text", encoding = "UTF-8"), flatten = TRUE),
    error = function(e) NULL
  )
  if (is.null(json)) return(NULL)

  info <- json$rxclassDrugInfoList$rxclassDrugInfo
  if (is.null(info)) return(NULL)

  # Extract rxclass block
  rx_info <- info$rxclass
  if (is.null(rx_info)) return(NULL)

  # turn rx_info into a proper tibble safely
  rx_info <- tryCatch(as_tibble(rx_info), error = function(e) NULL)
  if (is.null(rx_info)) return(NULL)

  # If "source" column is missing or nested, flatten again
  if (!"source" %in% names(rx_info)) {
    rx_info <- tryCatch(
      jsonlite::flatten(rx_info),
      error = function(e) rx_info
    )
  }

  # Ensure `source` column is atomic
  if (is.list(rx_info$source)) {
    rx_info$source <- sapply(rx_info$source, function(x) if (length(x) == 0) NA else as.character(x))
  }

  # Filter ATC
  atc <- rx_info |>
    filter(!is.na(source), source == "ATC") |>
    select(classId, className) # nolint: object_usage_linter.

  if (nrow(atc) == 0) return(NULL)
  return(atc)
}

cache_file <- "rxclass_cache.qs"
if (file.exists(cache_file)) {
  rxclass_cache <- qs::qread(cache_file)
} else {
  rxclass_cache <- list()
}

rxclass_results <- list()

# Process unique ingredients using the cache

for (ing in unique_ings) {
  
  # Use cache if exists
  if (!is.null(rxclass_cache[[ing]])) {
    rxclass_results[[ing]] <- rxclass_cache[[ing]]
    next
  }

  # Call the robust lookup function (Handles API logic and error checking)
  df_result <- rxclass_lookup(ing)
  
  # NOTE: Using a delay (Sys.sleep) is crucial here to avoid API rate-limiting!
  Sys.sleep(0.1) 

  # Save result to the final list and the cache
  rxclass_results[[ing]] <- df_result
  rxclass_cache[[ing]] <- df_result
  
  # Print progress
  cat("Processed:", ing, " | Classes found:", ifelse(is.null(df_result), 0, nrow(df_result)), "\n")
}

# Save the updated cache
qs::qsave(rxclass_cache, cache_file)

# Combine all list elements into a single tibble
atc_crosswalk_raw <- bind_rows(rxclass_results, .id = "clean_drug_names")

# Filter for only ATC codes and clean up the code structure
atc_crosswalk <- atc_crosswalk_raw |>
  # Ensure the 'source' column exists and is not nested before filtering
  mutate(source = unlist(source)) |>
  filter(!is.na(source), source == "ATC") |>
  
  # Select the drug name, the full ATC code (classId), and the 4th level
  select(
    clean_drug_names,
    atc_code_full = classId,
    atc_name = className
  ) |>
  mutate(
    # ATC-4 is the therapeutic subgroup—perfect for mapping to the disease level
    atc_4 = str_sub(atc_code_full, 1, 4)
  ) |>
  # Ensure the drug name is lower-case for merging reliability
  mutate(clean_drug_names = str_to_lower(clean_drug_names)) |>
  # Keep only unique ATC codes per ingredient (as one ingredient may have multiple rows)
  distinct(clean_drug_names, atc_4, .keep_all = TRUE)


# Merge ATC Codes back into Clinical Trials Data 

# Split the original data into long format for joining
ctg_long <- ctg_cleaned |>
  select(nct_number, clean_drug_names, sponsor, phases, primary_completion_date) |>
  separate_rows(clean_drug_names, sep = ",\\s*") |>
  mutate(clean_drug_names = str_to_lower(clean_drug_names)) |>
  filter(!is.na(clean_drug_names))

# Join with the ATC crosswalk
ctg_atc_mapped <- ctg_long |>
  left_join(
    atc_crosswalk,
    by = "clean_drug_names",
    relationship = "many-to-many" # One drug name can map to multiple ATC codes
  )

# Aggregate back to the trial level (NCT Number)
ctg_final_data <- ctg_atc_mapped |>
  group_by(nct_number) |>
  summarise(
    # Create the key variable for merging with Medicaid data
    atc_4_list = paste(unique(atc_4[!is.na(atc_4)]), collapse = ", "),
    
    # Keep the core study variables
    across(
      c(sponsor, phases, primary_completion_date),
      first
    ),
    .groups = "drop"
  ) |>
  # Filter out trials that could not be mapped to any ATC code
  filter(atc_4_list != "")