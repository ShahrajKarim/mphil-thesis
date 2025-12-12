# scripts/clinical_trials_data.r

# libraries
library(tidyverse)
library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(purrr)
library(DBI)
library(RSQLite)

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


# --- Clean Drug-NCT & Condition-NCT dataframes --- #

ctg_cleaned <- read.csv(
  "processed_data/clinical_trials/clinical_trials_drug_cleaned.csv"
)

final_drug_agg <- ctg_cleaned |>
  select(nct_number, interventions) |>

  separate_rows(interventions, sep = "\\|") |>
  mutate(interventions = str_trim(interventions)) |>
  filter(!is.na(interventions), interventions != "") |>
  filter(str_detect(interventions, "^(DRUG|BIOLOGICAL):")) |>
  mutate(
    clean_drug_name = interventions |>
      # Remove "DRUG: " or "BIOLOGICAL: " prefix
      str_replace_all("^(DRUG|BIOLOGICAL):\\s*", "") |> 
      # Remove dosage, parenthetical info, and trailing non-name characters
      str_replace_all("\\s*\\(.*|\\s*\\;.*|\\s*\\[.*|\\s*\\<.*|\\s*\\>.*", "") |>
      str_trim()
  ) |>
  filter(
    !str_detect(clean_drug_name, regex("Placebo|Saline|Vehicle|Control|Observation", ignore_case = TRUE)),
    nchar(clean_drug_name) > 0
  ) |>
  group_by(nct_number) |>
  summarise(
    clean_drug_names = paste(unique(clean_drug_name), collapse = "; "),
    .groups = 'drop'
  ) |>
  separate_rows(clean_drug_names, sep = ";") |>
  mutate(clean_drug_names = trimws(clean_drug_names))

conditions <- ctg_cleaned |>
  select(nct_number, conditions) |>
  separate_rows(conditions, sep = "\\|") |>
  mutate(conditions = str_squish(conditions)) |>
  filter(!is.na(conditions), conditions != "")

unique_conditions <- conditions |>
  select(conditions) |>
  distinct()

# --- RSQL UMLS Mapping --- #

DB_FILE <- "aux_data/umls.db"

classify_conditions_from_db <- function(conditions_df, db_file) {
    
    con <- dbConnect(RSQLite::SQLite(), db_file)
    
    dbWriteTable(con, "RAW_CONDITIONS", conditions_df, overwrite = TRUE)

    message("Executing complex SQL join against UMLS Metathesaurus for classification...")
    
    # --- UPDATED SQL QUERY ---
    sql_query <- "
    SELECT
      T1.nct_number,
      T1.conditions,
      T4.STY AS final_condition_ta,
      -- NEW: Pull TTY and SAB to aid in R-side granular mapping
      T2.TTY AS source_term_type,
      T2.SAB AS source_vocabulary
    FROM RAW_CONDITIONS AS T1
    
    -- Join raw condition to the UMLS concept table
    INNER JOIN MRCONSO AS T2
      ON T2.STR LIKE '%' || T1.conditions || '%' 

    -- Join concept to Semantic Type table
    INNER JOIN MRSTY AS T4
      ON T2.CUI = T4.CUI
      
    -- Filter for high-quality sources and term types
    WHERE
      T2.SAB IN ('MSH', 'SNOMEDCT_US', 'ICD10CM', 'ICD9CM') 
      AND T2.TTY IN ('PT', 'HT', 'SY')
    
    -- Group to select the best match
    GROUP BY
      T1.nct_number, T1.conditions;
    "
    
    classified_results <- dbGetQuery(con, sql_query)
    
    dbRemoveTable(con, "RAW_CONDITIONS")
    dbDisconnect(con)
    return(classified_results)
}

# --- Execute Classification ---
classified_data <- classify_conditions_from_db(conditions, DB_FILE)

# Aggregate the final TA back to the trial level (NCT)
final_condition_ta_agg <- classified_data |>
  group_by(nct_number) |>
  summarise(
    final_condition_ta = paste(unique(final_condition_ta), collapse = "; "),
    .groups = 'drop'
  )




















































# --- API draft - new method --- #

nct_ids <- ctg_cleaned$nct_number |> unique()

nct_ids <- nct_ids[1:100]

fetch_nih_therapeutic_areas <- function(nct_id) {
  
  url <- paste0(
    "https://clinicaltrials.gov/api/v2/studies/",
    nct_id
  )
  
  # ... (HTTP GET and status check logic remains here, skipping for brevity) ...
  res <- tryCatch(httr::GET(url), error = function(e) NULL) 
  if (is.null(res) || res$status_code != 200) return(NULL)
  
  raw_content <- content(res, "text", encoding = "UTF-8")
  json <- tryCatch(
    jsonlite::fromJSON(raw_content, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(json)) return(NULL)
  
  # --- CRITICAL FIX: Check the 'derivedSection' FIRST for MeSH data ---
  
  study_data <- json[["protocolSection"]]
  derived_data <- json[["derivedSection"]] # New key section to check
  
  # 1. ATTEMPT 1: Get MeSH terms from the DERIVED section (most reliable)
  if (!is.null(derived_data)) {
    path_to_branch <- derived_data[["conditionBrowse"]][["conditionBrowseBranch"]]
    
    if (!is.null(path_to_branch)) {
      df <- tryCatch(
        dplyr::bind_rows(path_to_branch) |> 
          dplyr::select(branchAbbrev, branchName) |> 
          dplyr::rename(code = branchAbbrev, name = branchName),
        error = function(e) NULL
      )
      if (!is.null(df) && nrow(df) > 0) {
        cat("\nSUCCESS: Data extracted from MeSH (Derived Section).")
        df$nct_id <- nct_id
        return(df)
      }
    }
  }
  
  # 2. ATTEMPT 2: Fallback to Primary Study Conditions (from protocolSection)
  
  # Conditions are in 'conditions' under 'conditionsModule'
  if (!is.null(study_data)) {
    path_to_conditions <- study_data[["conditionsModule"]][["conditions"]]
    
    if (!is.null(path_to_conditions)) {
      # Convert the list of conditions into a dataframe
      df_conditions <- data.frame(
        code = NA_character_, 
        name = unlist(path_to_conditions),
        nct_id = nct_id
      )
      cat("\nSUCCESS: Data extracted from Fallback Conditions (protocolSection).")
      return(df_conditions)
    }
  }
  
  cat("\nERROR: No MeSH or Conditions data found for this trial.")
  return(NULL)
}

cache_file_nih <- "ctgov_therapeutic_area_cache.qs"

if (file.exists(cache_file_nih)) {
  nih_cache <- qs::qread(cache_file_nih)
} else {
  nih_cache <- list()
}


nih_results <- list()

for (nct in nct_ids) {
  
  # If cached, reuse
  if (!is.null(nih_cache[[nct]])) {
    nih_results[[nct]] <- nih_cache[[nct]]
    next
  }
  
  df <- fetch_nih_therapeutic_areas(nct)
  
  Sys.sleep(0.15)  # avoid rate-limiting
  
  nih_results[[nct]] <- df
  nih_cache[[nct]] <- df
  
  cat("Processed:", nct, "\n")
}

qs::qsave(nih_cache, cache_file_nih)

nih_raw <- bind_rows(nih_results, .id = "nct_number")

# Only keep rows that have successfully extracted data
nih_therapeutic <- nih_raw %>%
  filter(!is.na(BranchAbbrev)) %>% 
  group_by(nct_number) %>%
  summarise(
    # Use the MeSH code (BranchAbbrev) as the main grouping variable
    nih_mesh_codes = paste(unique(BranchAbbrev), collapse = "; "),
    nih_mesh_names = paste(unique(BranchName), collapse = "; "),
    .groups = "drop"
  )

# Join the NIH/MeSH therapeutic areas back to your cleaned data
ctg_with_ta <- ctg_cleaned %>%
  left_join(nih_therapeutic, by = "nct_number")
































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