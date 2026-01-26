# scripts/fda_approval_data.r

# Pipeline overview:
# 1. FDA approval extraction:
#    - Identify first approval dates (NDA/ANDA/BLA) from Submissions.txt.
#    - Merge with Applications.txt for sponsor information.
#    - Add drug names and active ingredients from Products.txt.
#
# 2. Ingredient processing:
#    - Standardise and split multi-ingredient entries.
#    - Create a normalised ingredient identifier.

# 3. Therapeutic classification:
#    - Use RxClass API (cached) to map ingredients to ATC classes.
#    - Derive ATC Level 1–2 categories from ATC codes.
#
# 4. Final approvals dataset:
#    - Collapse ingredient-level ATC results to the application level.
#    - Merge firm, approval date, and therapeutic area information.
#    - Export clean dataset for downstream merging with SDUD.

# libraries
library(readr)
library(dplyr)
library(janitor)
library(stringr)
library(tidyr)
library(httr)
library(httr2)
library(jsonlite)
library(purrr)
library(qs)
library(lubridate)
library(rxnorm)

# --- Load FDA approval data --- #

fda_path <- "raw_data/Drugs@FDA"

applications <- read_delim(
  file.path(fda_path, "Applications.txt"),
  delim = "\t",
  col_types = cols(.default = "c")
) |> clean_names()

products <- read_delim(
  file.path(fda_path, "Products.txt"),
  delim = "\t",
  col_types = cols(.default = "c")
) |> clean_names()

submissions <- read_delim(
  file.path(fda_path, "Submissions.txt"),
  delim = "\t",
  col_types = cols(.default = "c")
) |> clean_names()

marketing_status <- read_delim(
  file.path(fda_path, "MarketingStatus.txt"),
  delim = "\t",
  col_types = cols(.default = "c")
) |> clean_names()

te_codes <- read_delim(
  file.path(fda_path, "TE.txt"),
  delim = "\t",
  col_types = cols(.default = "c")
) |> clean_names()

# Extract approval events

approvals <- submissions |>
  filter(submission_status == "AP") |>
  rename(submission_date = submission_status_date) |>
  mutate(
    submission_date = str_sub(submission_date, 1, 10),
    submission_date = as.Date(submission_date, format = "%Y-%m-%d")
  ) |>
  group_by(appl_no) |>
  summarise(
    first_approval_date = min(submission_date, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    first_approval_date = ifelse(
                                 is.infinite(first_approval_date),
                                 NA,
                                 first_approval_date)
  ) |>
  mutate(
    first_approval_date = as.Date(first_approval_date, format = "%Y-%m-%d")
  )

# Merge approvals with applications and products

fda_drugs <- applications |>
  left_join(approvals, by = "appl_no") |>
  left_join(products,  by = "appl_no") |>
  select(
         -c(appl_public_notes,
           form,
           strength,
           reference_drug,
           reference_standard
         ))


# Split active ingredients into separate strings

active_ingredients <- fda_drugs |>
  mutate(
    active_ing = gsub("\\|\\|", ";", active_ingredient),
    active_ing = gsub("\\(.*?\\)", "", active_ing),
    active_ing = trimws(active_ing),
    active_ing_split = strsplit(active_ing, ";")
  ) |>
  unnest(active_ing_split) |>
  mutate(active_ing_split = trimws(active_ing_split)) |>
  filter(active_ing_split != "")

# --- RXCLASS ATC MAPPING --- #

# Extract unique ingredients
unique_ings <- active_ingredients$active_ing_split |>
  toupper() |>
  unique()

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

# Caching setup
cache_file <- "rxclass_cache.qs"
if (file.exists(cache_file)) {
  rxclass_cache <- qs::qread(cache_file)
} else {
  rxclass_cache <- list()
}

# Batch API processing
# Rebuild rxclass_results (fixing ATC extraction bug) 

rxclass_results <- list()

for (ing in unique_ings) {
  
  # Use cache if exists
  if (!is.null(rxclass_cache[[ing]])) {
    rxclass_results[[ing]] <- rxclass_cache[[ing]]
    next
  }

  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxclass/class/byDrugName.json?drugName=",
    URLencode(ing)
  )

  res <- tryCatch(httr::GET(url), error = function(e) NULL)
  if (is.null(res) || res$status_code != 200) {
    rxclass_results[[ing]] <- NULL
    next
  }

  json <- tryCatch(
    jsonlite::fromJSON(content(res, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(json)) {
    rxclass_results[[ing]] <- NULL
    next
  }

  df <- tryCatch(
    as.data.frame(json$rxclassDrugInfoList$rxclassDrugInfo$rxclass),
    error = function(e) NULL
  )

  rxclass_results[[ing]] <- df
  rxclass_cache[[ing]] <- df

  cat("Processed:", ing, "\n")
}

qs::qsave(rxclass_cache, "rxclass_cache.qs")

# Build ATC mapping table
rxclass_atc_map <- bind_rows(
  lapply(names(rxclass_results), function(nm) {

    df <- rxclass_results[[nm]]

    if (is.null(df)) return(NULL)
    if (!is.list(df)) return(NULL)

    df <- tryCatch(
      as.data.frame(df, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(df)) return(NULL)

    # Required columns
    needed <- c("classId", "className", "classType")
    if (!all(needed %in% names(df))) return(NULL)

    # Use ATC1-4 instead of ATC
    df <- df[df$classType == "ATC1-4", ]
    if (nrow(df) == 0) return(NULL)

    df$ingredient_norm <- nm

    df[, c("ingredient_norm", "classId", "className")]
  })
)

names(rxclass_atc_map)[names(rxclass_atc_map) == "classId"] <- "ATC_CODE"
names(rxclass_atc_map)[names(rxclass_atc_map) == "className"] <- "ATC_NAME"

rxclass_atc_map <- distinct(rxclass_atc_map)

# Join ATC back to FDA data + derive ATC levels
active_ingredients_atc <- active_ingredients |>
  mutate(ingredient_norm = toupper(active_ing_split)) |>
  left_join(rxclass_atc_map, by = "ingredient_norm")

# Derive ATC hierarchy levels 
active_ingredients_atc <- active_ingredients_atc |>
  mutate(
    ATC1 = ifelse(!is.na(ATC_CODE), substr(ATC_CODE, 1, 1), NA),
    ATC2 = ifelse(!is.na(ATC_CODE), substr(ATC_CODE, 1, 3), NA),
    ATC3 = ifelse(!is.na(ATC_CODE), substr(ATC_CODE, 1, 4), NA)
  )

# Collapse to application-level mapping
fda_with_atc <- active_ingredients_atc |>
  group_by(appl_no) |>
  summarise(
    atc_codes = paste(sort(unique(na.omit(ATC_CODE))), collapse = "; "),
    atc_level1 = paste(sort(unique(na.omit(ATC1))), collapse = "; "),
    atc_level2 = paste(sort(unique(na.omit(ATC2))), collapse = "; "),
    atc_level3 = paste(sort(unique(na.omit(ATC3))), collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(fda_drugs, by = "appl_no") |>
  mutate(
    atc_matched = atc_codes != "" & !is.na(atc_codes)
  )

# Merge diagnostics
message("Total rows: ", nrow(fda_with_atc)) # 50,368
message("Matched rows: ", sum(fda_with_atc$atc_matched)) # 48,143
message("Unmatched rows: ", sum(!fda_with_atc$atc_matched)) # 2,225
message("Match rate: ", round(mean(fda_with_atc$atc_matched) * 100, 2), "%") # 95.58%


# Save cleaned file

write.csv(fda_with_atc, "processed_data/fda_approvals/fda_approvals_cleaned.csv", row.names = FALSE)

# Conduct ATC mapping

atc_map <- read.csv("aux_data/atc_level_mapping.csv")

atc_level_1 <- atc_map |>
  distinct(ATC1_CODE, ATC1_NAME)

atc_level_2 <- atc_map |>
  distinct(ATC2_CODE, ATC2_NAME)

fda_atc_expanded <- fda_with_atc |>
  separate_rows(atc_level1, sep = ";\\s*") |>
  separate_rows(atc_level2, sep = ";\\s*") |>
  distinct() |>
  filter(atc_matched) |>
  select(
         -c(
           atc_codes,
           atc_level3,
           appl_type,
           atc_matched,
           product_no
         ))

fda_atc_expanded <- fda_atc_expanded |>
  left_join(atc_level_1, by = c("atc_level1" = "ATC1_CODE")) |>
  left_join(atc_level_2, by = c("atc_level2" = "ATC2_CODE"))


# Save final FDA approvals with ATC mapping

write.csv(fda_atc_expanded, "processed_data/fda_approvals/fda_approvals_atc_mapped.csv", row.names = FALSE)

# --- Firm - NDC MAPPING --- #

fda_atc_expanded <- read.csv("processed_data/fda_approvals/fda_approvals_atc_mapped.csv")

# Extract firm names if date >= 2005-01-01
fda_firms <- fda_atc_expanded |>
  filter(as.Date(first_approval_date) >= as.Date("2005-01-01")) |>
  distinct(sponsor_name) |>
  rename(firm = sponsor_name)

# Import NDC firm mapping
ndc_firms <- read_csv("aux_data/ndc_firm_mapping.csv")

# Clean FDA firm names
fda_firms <- fda_firms |>
  mutate(
    firm_clean = firm |> 
      toupper() |>
      gsub("[^A-Z0-9 ]", " ", x = _) |>
      gsub("\\b(CORP|CORPORATION|INC|CO|COMPANY|LLC|LTD|LP|PLC|AG|SA|PHARMA|PHARMACEUTICALS)\\b", "", x = _) |>
      gsub("\\s+", " ", x = _) |>
      trimws()
  )

# Clean NDC firm names
ndc_firms <- ndc_firms |>
  mutate(
    firm_utf8 = iconv(firm, from = "", to = "UTF-8", sub = ""),   # drop invalid bytes
    firm_clean = firm_utf8 |>
      toupper() |>
      gsub("[^A-Z0-9 ]", " ", x = _) |>
      gsub("\\b(CORP|CORPORATION|INC|CO|COMPANY|LLC|LTD|LP|PLC|AG|SA|PHARMA|PHARMACEUTICALS)\\b", "", x = _) |>
      gsub("\\s+", " ", x = _) |>
      trimws()
  ) |>
  filter(!is.na(firm_clean), firm_clean != "") |>
  arrange(firm_clean) |>
  group_by(firm_clean) |>
  summarise(
    labeler_code = first(labeler_code_pad),
    .groups = "drop"
  )

# Exact cleaned match
stage1 <- fda_firms |>
  left_join(ndc_firms, by = "firm_clean", suffix = c("_fda", "_ndc"))

# Identify unmatched firms
unmatched <- stage1 |>
  filter(is.na(labeler_code)) |>
  mutate(
    root = str_match(firm_clean, "\\b([A-Z0-9]{4,})\\b")[, 2],
    root = ifelse(
      is.na(root),
      str_extract(firm_clean, "^[A-Z0-9]+"),
      root
    )
  ) |>
  select(-labeler_code)

# Build NDC firm root table - first word with > 3 chars

ndc_roots <- ndc_firms |>
  mutate(
    root = str_match(firm_clean, "\\b([A-Z0-9]{4,})\\b")[, 2],
    root = ifelse(is.na(root),
                  str_extract(firm_clean, "^[A-Z0-9]+"),
                  root)
  ) |>
  filter(!is.na(root)) |>

  # Pick the most frequent labeler per root
  # Some misidentification will occur but numbers are small
  group_by(root, labeler_code) |>
  summarise(n = n(), .groups = "drop") |>
  group_by(root) |>
  slice_max(n, with_ties = FALSE) |>  # ensures exactly ONE labeler per root
  ungroup() |>
  select(root, labeler_code)

# Merge unmatched with root-based mapping
fallback <- unmatched |>
  left_join(ndc_roots, by = "root")

# Combine exact matches with fallbacks and non-matches
firm_matches <- stage1 |>
  filter(!is.na(labeler_code)) |>
  bind_rows(
    fallback
  ) |>
  select(firm, labeler_code)

# Merge diagnostics
message("Total rows: ", nrow(firm_matches)) # 1,465
message("Matched rows: ", sum(!is.na(firm_matches$labeler_code))) # 1015
message("Unmatched rows: ", sum(is.na(firm_matches$labeler_code))) # 450
message("Match rate: ", round(mean(!is.na(firm_matches$labeler_code)) * 100, 2), "%") # 69.28%

# Save firm - NDC mapping

write.csv(
  firm_matches,
  file = "aux_data/ndc_fda_firm_mapping.csv",
  row.names = FALSE
)

# Merge firm mapping with FDA data

fda_atc_expanded <- fda_atc_expanded |>
  left_join(firm_matches, by = c("sponsor_name" = "firm")) |>
  filter(as.Date(first_approval_date) >= as.Date("2005-01-01"))

message("Match rate: ", round(mean(!is.na(fda_atc_expanded$labeler_code)) * 100, 2), "%") # 82.26%


# Generate Year & Quarter variables

fda_atc_expanded <- fda_atc_expanded |>
  mutate(
    approval_year    = year(first_approval_date),
    approval_quarter = quarter(first_approval_date)
  )

# Save final FDA approvals with firm mapping

write.csv(
  fda_atc_expanded,
  "processed_data/fda_approvals/fda_approvals_atc_firm_mapped.csv",
  row.names = FALSE
)

# --- SDUD MAPPING for FDA data --- #

SDUD_firm_mapped <- readRDS("processed_data/state_drug_utilisation_data/SDUD_firm_mapped.rds")

ndc_list <- SDUD_firm_mapped |> 
  distinct(ndc) |> 
  pull(ndc)

# Define function mapping NDC to RXCUI

rxnorm_ndc_to_rxcui <- function(ndc_vec) {

  ndc_vec <- ndc_vec |> unique()

  safely_get <- safely(function(ndc) {
    url <- paste0("https://rxnav.nlm.nih.gov/REST/ndcstatus.json?ndc=", ndc)
    resp <- request(url) |> req_perform()
    json <- resp_body_json(resp)

    rxcui <- json$ndcStatus$rxcui
    if (is.null(rxcui) || rxcui == "") return(NA)
    return(rxcui)
  })

  results <- map_chr(ndc_vec, ~ safely_get(.x)$result %||% NA)

  tibble(
    ndc = ndc_vec,
    rxcui = results
  )
}

# Define function mapping RXCUI to ATC

rxcui_to_atc <- function(rxcui_df) {

  valid_rxcui <- rxcui_df |> 
    filter(!is.na(rxcui)) |>
    distinct(rxcui) |>
    pull(rxcui)

  atc_mapping <- map_dfr(valid_rxcui, function(rxcui_id) {

    # Call the API for the single CUI
    atc_codes <- get_atc(rxcui_id)

    # If result is NA (no ATC found)
    if (all(is.na(atc_codes))) {
      return(tibble(rxcui = rxcui_id, atc = NA_character_))
    }

    # Convert the vector result into a two-column tibble: CUI and ATC
    return(tibble(
      rxcui = rxcui_id,
      atc = atc_codes
    ))
  }) |>
    group_by(rxcui) |>
    summarise(
      atc = paste(unique(atc[!is.na(atc)]), collapse = "; "),
      .groups = 'drop'
    ) |>
    mutate(atc = ifelse(atc == "", NA_character_, atc))

  return(atc_mapping)
}

# Define function to generate ATC levels 1 and 2

generate_atc_levels <- function(atc_mapped_df) {

  atc_level_df <- atc_mapped_df |>
    filter(!is.na(atc), atc != "Unclassified") |>
    separate_rows(atc, sep = ";\\s*") |>
    mutate(
      ATC_Level1 = case_when(
        nchar(atc) >= 1 ~ stringr::str_sub(atc, 1, 1),
        TRUE ~ NA_character_
      )
    ) |>
    mutate(
      ATC_Level2 = case_when(
        nchar(atc) >= 3 ~ stringr::str_sub(atc, 1, 3),
        TRUE ~ NA_character_
      )
    ) |>
    select(rxcui, atc, ATC_Level1, ATC_Level2) |>
    distinct()

  return(atc_level_df)

}

# Define function that combines the three steps together

map_ndc_to_atc_full <- function(ndc_vec) {

  step1_df <- rxnorm_ndc_to_rxcui(ndc_vec)

  step2_df <- rxcui_to_atc(step1_df)

  step3_df <- generate_atc_levels(step2_df)

  final_mapping <- step1_df |>
    left_join(step3_df, by = "rxcui", relationship = "many-to-many")

  return(final_mapping)
}

ndc_atc_mapping <- map_ndc_to_atc_full(ndc_list)

write.csv(
  ndc_atc_mapping,
  file = "aux_data/ndc_atc_mapping.csv",
  row.names = FALSE
)

# Now merge the mapping with the SDUD data

SDUD_firm_mapped <- SDUD_firm_mapped |>
  left_join(ndc_atc_mapping, by = "ndc", relationship = "many-to-many") |>
  select(- rxcui)

# Merge diagnostics

message("Total rows: ", nrow(SDUD_firm_mapped)) # 42373338
message("Matched rows: ", sum(!is.na(SDUD_firm_mapped$atc))) # 40722005
message("Unmatched rows: ", sum(is.na(SDUD_firm_mapped$atc))) # 1651333
message("Match rate: ", round(mean(!is.na(SDUD_firm_mapped$atc)) * 100, 2), "%") # 96.10%

SDUD_firm_mapped <- SDUD_firm_mapped |>
  filter(!is.na(atc))

# Load in FDA approvals data with firm mapping

fda_atc_expanded <- read_csv(
  "processed_data/fda_approvals/fda_approvals_atc_firm_mapped.csv",
  col_types = cols(
    labeler_code = col_character()
  )
)

# Final name cleaning for FDA firms

sponsor_final <- fda_atc_expanded |>
  group_by(labeler_code) |>
  summarise(
    sponsor_canonical = first(sponsor_name),
    .groups = "drop"
  ) |>
  select(labeler_code, sponsor_canonical) |>
  rename(sponsor_name = sponsor_canonical)

fda_atc_expanded <- fda_atc_expanded |>
  select(-sponsor_name) |>
  left_join(sponsor_final, by = "labeler_code")

# Collapse FDA data to firm-year-quarter-ATC level

fda_approvals_panel <- fda_atc_expanded |>
  group_by(
    sponsor_name,
    labeler_code,
    approval_year,
    approval_quarter,
    atc_level1,
    atc_level2
  ) |>
  summarise(
    n_approvals = n(),
    .groups = "drop"
  )

# Collapse SDUD data to firm-year-quarter-ATC level

SDUD_firm_ATC <- SDUD_firm_mapped |>
  group_by(
    state,
    firm,
    labeler_code,
    year,
    quarter,
    ATC_Level1,
    ATC_Level2
  ) |>
  summarise(
    prescriptions = sum(number_of_prescriptions, na.rm = TRUE),
    medicaid_reimbursed = sum(medicaid_amount_reimbursed, na.rm = TRUE),
    non_medicaid_reimbursed = sum(non_medicaid_amount_reimbursed, na.rm = TRUE),
    .groups = "drop"
  )

# SDUD ATC level 1 aggregate indicators

SDUD_ATC1 <- SDUD_firm_ATC |>
  group_by(year, quarter, ATC_Level1) |>
  summarise(
    ATC1_prescriptions = sum(prescriptions, na.rm = TRUE),
    ATC1_medicaid_reimbursed = sum(medicaid_reimbursed, na.rm = TRUE),
    ATC1_non_medicaid_reimbursed = sum(non_medicaid_reimbursed, na.rm = TRUE),
    .groups = "drop"
  )

# SDUD ATC level 2 aggregate indicators

SDUD_ATC2 <- SDUD_firm_ATC |>
  group_by(year, quarter, ATC_Level2) |>
  summarise(
    ATC2_prescriptions = sum(prescriptions, na.rm = TRUE),
    ATC2_medicaid_reimbursed = sum(medicaid_reimbursed, na.rm = TRUE),
    ATC2_non_medicaid_reimbursed = sum(non_medicaid_reimbursed, na.rm = TRUE),
    .groups = "drop"
  )

# Merge SDUD datasets to obtain resulting SDUD master file

SDUD_master <- SDUD_firm_ATC |>
  left_join(
    SDUD_ATC1,
    by = c("year", "quarter", "ATC_Level1")
  ) |>
  left_join(
    SDUD_ATC2,
    by = c("year", "quarter", "ATC_Level2")
  )


# Merge FDA approvals panel with SDUD data

fda_sdud_merged <- fda_approvals_panel |>
  left_join(
    SDUD_master,
    by = c(
      "labeler_code" = "labeler_code",
      "approval_year" = "year",
      "approval_quarter" = "quarter",
      "atc_level1" = "ATC_Level1",
      "atc_level2" = "ATC_Level2"
    )
  )

# Save final merged dataset

write.csv(
  fda_sdud_merged,
  "processed_data/fda_approvals/fda_sdud_merged.csv",
  row.names = FALSE
)

# Repeat merging with dataset for SDUD + PDL data

# Load in NDC ATC mapping

ndc_atc_mapping <- read.csv(
  "aux_data/ndc_atc_mapping.csv",
  colClasses = c(ndc = "character")
)

# Load in SDUD PDL data

SDUD_PDL <- readRDS(
  file = "processed_data/state_drug_utilisation_data/SDUD_firm_PDL.rds"
) |>
select(
  state, ndc, firm, labeler_code, product_code,
  year, quarter, number_of_prescriptions,
  medicaid_amount_reimbursed, non_medicaid_amount_reimbursed,
  medicaid_coverage_start, medicaid_coverage_end, medicaid_covered,
  price_stage_pass, pt_stage_pass, remain_in_list
) |>
  left_join(
    ndc_atc_mapping,
    by = "ndc",
    relationship = "many-to-many"
  ) |>
  filter(!is.na(atc)) |>
  select(- c(rxcui, atc, ATC_Level2)) |>
  mutate(
    in_pdl = ifelse(price_stage_pass == 1 & pt_stage_pass == 1, 1, 0)
  )

# Generate Medicaid exposure

PDL_exposure <- SDUD_PDL |>
  select(labeler_code, year, price_stage_pass, pt_stage_pass, in_pdl, ATC_Level1) |>
  filter(year >= 2010 & year <= 2014) |>
  mutate(
    price_stage_pass = ifelse(is.na(price_stage_pass), 0, price_stage_pass),
    pt_stage_pass = ifelse(is.na(pt_stage_pass), 0, pt_stage_pass),
    in_pdl = ifelse(is.na(in_pdl), 0, in_pdl)
  ) |>
  group_by(labeler_code, ATC_Level1) |>
  summarise(
    pt_stage_pass_prop = mean(pt_stage_pass),
    price_stage_pass_prop = mean(price_stage_pass),
    in_pdl_prop = mean(in_pdl),
    .groups = "drop"
  )

# Collapse SDUD-PDL data to firm-year-quarter-ATC level

SDUD_PDL <- SDUD_PDL |>
  left_join(
    PDL_exposure,
    by = c("labeler_code", "ATC_Level1")
  ) |>
  group_by(
    labeler_code,
    year,
    quarter,
    ATC_Level1,
  ) |>
  summarise(
    prescriptions = sum(number_of_prescriptions, na.rm = TRUE),
    medicaid_reimbursed = sum(medicaid_amount_reimbursed, na.rm = TRUE),
    non_medicaid_reimbursed = sum(non_medicaid_amount_reimbursed, na.rm = TRUE),
    price_stage_pass = first(price_stage_pass_prop),
    pt_stage_pass = first(pt_stage_pass_prop),
    in_pdl = first(in_pdl_prop),
    .groups = "drop"
  ) |>
  mutate(
    price_stage_pass = ifelse(is.na(price_stage_pass), 0, price_stage_pass),
    pt_stage_pass = ifelse(is.na(pt_stage_pass), 0, pt_stage_pass),
    in_pdl = ifelse(is.na(in_pdl), 0, in_pdl)
  )

SDUD_ATC1 <- SDUD_PDL |>
  group_by(year, quarter, ATC_Level1) |>
  summarise(
    ATC1_prescriptions = sum(prescriptions, na.rm = TRUE),
    ATC1_medicaid_reimbursed = sum(medicaid_reimbursed, na.rm = TRUE),
    ATC1_non_medicaid_reimbursed = sum(non_medicaid_reimbursed, na.rm = TRUE),
    .groups = "drop"
  )

SDUD_PDL <- SDUD_PDL |>
  left_join(
    SDUD_ATC1,
    by = c("year", "quarter", "ATC_Level1")
  )


# Merge FDA approvals panel with SDUD PDL data

fda_sdud_pdl_merged <- fda_approvals_panel |>
  left_join(
    SDUD_PDL,
    by = c(
      "labeler_code" = "labeler_code",
      "approval_year" = "year",
      "approval_quarter" = "quarter",
      "atc_level1" = "ATC_Level1"
    )
  ) |>
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# Save final merged dataset

write.csv(
  fda_sdud_pdl_merged,
  "processed_data/fda_approvals/fda_sdud_pdl_merged.csv",
  row.names = FALSE
)