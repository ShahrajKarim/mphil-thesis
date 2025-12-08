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
#    - Derive ATC Level 1–3 categories from ATC codes.
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
library(jsonlite)
library(purrr)
library(qs)
library(lubridate)

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


# -- Draft work for ATC mapping: RxClass API version (cached) -- #

#====================== RXCLASS ATC MAPPING (CACHED) ======================#

# 1. Extract unique ingredients
unique_ings <- active_ingredients$active_ing_split |>
  toupper() |>
  unique()

# 2. RxClass lookup function

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

  # ---- Key fix: turn rx_info into a proper tibble safely ----
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
    select(classId, className)

  if (nrow(atc) == 0) return(NULL)
  return(atc)
}

# 3. Caching setup
cache_file <- "rxclass_cache.qs"
if (file.exists(cache_file)) {
  rxclass_cache <- qs::qread(cache_file)
} else {
  rxclass_cache <- list()
}

# 4. Batch API processing
# --- Rebuild rxclass_results (fixing ATC extraction bug) ---

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

  # FIXED: correct path and no "source" filtering
  df <- tryCatch(
    as.data.frame(json$rxclassDrugInfoList$rxclassDrugInfo$rxclass),
    error = function(e) NULL
  )

  rxclass_results[[ing]] <- df
  rxclass_cache[[ing]] <- df

  cat("Processed:", ing, "\n")
}

qs::qsave(rxclass_cache, "rxclass_cache.qs")

# 5. Build ATC mapping table (robust version)
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

    # THE FIX: ATC1-4 instead of ATC
    df <- df[df$classType == "ATC1-4", ]
    if (nrow(df) == 0) return(NULL)

    df$ingredient_norm <- nm

    df[, c("ingredient_norm", "classId", "className")]
  })
)

names(rxclass_atc_map)[names(rxclass_atc_map) == "classId"] <- "ATC_CODE"
names(rxclass_atc_map)[names(rxclass_atc_map) == "className"] <- "ATC_NAME"

rxclass_atc_map <- distinct(rxclass_atc_map)

# 6. Join ATC back to FDA data + derive ATC levels
active_ingredients_atc <- active_ingredients |>
  mutate(ingredient_norm = toupper(active_ing_split)) |>
  left_join(rxclass_atc_map, by = "ingredient_norm")

# --- Derive ATC hierarchy levels ---
active_ingredients_atc <- active_ingredients_atc |>
  mutate(
    ATC1 = ifelse(!is.na(ATC_CODE), substr(ATC_CODE, 1, 1), NA),
    ATC2 = ifelse(!is.na(ATC_CODE), substr(ATC_CODE, 1, 3), NA),
    ATC3 = ifelse(!is.na(ATC_CODE), substr(ATC_CODE, 1, 4), NA)
  )

# --- Collapse to application-level mapping ---
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

# Diagnostics
message("Total rows: ", nrow(fda_with_atc)) # 50,368
message("Matched rows: ", sum(fda_with_atc$atc_matched)) # 48,143
message("Unmatched rows: ", sum(!fda_with_atc$atc_matched)) # 2,225
message("Match rate: ", round(mean(fda_with_atc$atc_matched) * 100, 2), "%") # 95.58%


# Save cleaned file

write.csv(fda_with_atc, "processed_data/fda_approvals/fda_approvals_cleaned.csv", row.names = FALSE)

# --- Conduct ATC mapping ---

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























# -- Draft work for ATC mapping -- #

# --- Load RxNorm data --- #

rxnorm_path <- "raw_data/RxNorm/rrf"


## ----------------------------- RXNCONSO --------------------------------- ##
# RXNCONSO: drug names + ATC vocabulary
conso_col_names <- c(
  "RXCUI","LAT","TS","LUI","STT","SUI","ISPREF","RXAUI","SAUI",
  "SCUI","SDUI","SAB","TTY","CODE","STR","SRL","SUPPRESS","CVF","EMPTY"
)

rxnconso <- read_delim(
  file.path(rxnorm_path, "RXNCONSO.RRF"),
  delim = "|",
  col_names = conso_col_names,
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = FALSE
) |>
  select(-EMPTY)

## ----------------------------- RXNREL (PIN → IN) ------------------------- ##
rel_col_names <- c(
  "RXCUI1","RXAUI1","STYPE1",
  "REL","RELA",
  "RXCUI2","RXAUI2","STYPE2",
  "RELB","RUI","SRUI","SRUI2",
  "SUPPRESS","CVF","EMPTY"
)

rxnrel <- read_delim(
  file.path(rxnorm_path, "RXNREL.RRF"),
  delim = "|",
  col_names = rel_col_names,
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = FALSE
) |>
  select(-EMPTY)

pin_to_in <- rxnrel |>
  filter(RELA == "ingredient_of") |>
  select(
    RXCUI_PIN = RXCUI1,
    RXCUI_IN  = RXCUI2
  ) |>
  distinct()

rxnconso_all <- rxnconso |>
  filter(SAB == "RXNORM") |>
  mutate(str_upper = toupper(STR)) |>
  select(RXCUI, STR, TTY, str_upper)

rxn_mapped <- rxnconso_all |>
  left_join(pin_to_in, by = c("RXCUI" = "RXCUI_PIN")) |>
  mutate(
    RXCUI_FINAL = ifelse(!is.na(RXCUI_IN), RXCUI_IN, RXCUI)
  )

rxnorm_ingredients_final <- rxn_mapped |>
  filter(TTY %in% c("IN", "PIN", "MIN", "PN")) |>
  mutate(ingredient_norm = str_upper) |>
  select(ingredient_norm, RXCUI_FINAL) |>
  distinct()


## ----------------------------- RXNSAT ------------------------------------ ##
# RXNSAT: ATC codes + ATC hierarchy levels
sat_col_names <- c(
  "RXCUI","LUI","SUI","RXAUI","STYPE","CODE","ATUI","SATUI",
  "SAB","ATN","ATV","SUPPRESS","CVF","EMPTY_COL"
)

rxnsat <- read_delim(
  file.path(rxnorm_path, "RXNSAT.RRF"),
  delim = "|",
  col_names = sat_col_names,
  col_types = cols(.default = "c"),
  escape_double = FALSE,
  trim_ws = FALSE
)

# ✔ Correct ATC extraction logic for *your* RxNorm version:
#   - SAB == "ATC_LEVEL"   → rows that include ATC codes
#   - CODE contains ATC_CODE (e.g. "R05CB05")
#   - ATV holds ATC_LEVEL (1–5) → not useful for mapping

ingredient_atc <- rxnsat |>
  filter(SAB == "ATC_LEVEL") |>
  filter(nchar(CODE) >= 5) |>            # ATC codes are 5–7 chars
  select(
    RXCUI,
    ATC_CODE = CODE
  ) |>
  distinct()

# ---------------------- Map FDA Ingredients → ATC ------------------------- #

active_ingredients_mapped <- active_ingredients |>
  mutate(ingredient_norm = toupper(active_ing_split)) |>
  left_join(
    rxnorm_ingredients_final,
    by = "ingredient_norm",
    relationship = "many-to-one"
  ) |>
  left_join(
    ingredient_atc,
    by = c("RXCUI_FINAL" = "RXCUI"),
    relationship = "many-to-many"
  )

fda_with_atc <- active_ingredients_mapped |>
  group_by(appl_no) |>
  summarise(
    atc_codes = paste(sort(unique(na.omit(ATC_CODE))), collapse = "; "),
    .groups = "drop"
  ) |>
  left_join(fda_drugs, by = "appl_no")


# Merge diagnostics
fda_with_atc <- fda_with_atc |>
  mutate(
    atc_matched = atc_codes != "" & !is.na(atc_codes)
  )

message("Total rows: ", nrow(fda_with_atc)) # 50368
message("Matched rows: ", sum(fda_with_atc$atc_matched)) # 25026
message("Unmatched rows: ", sum(!fda_with_atc$atc_matched)) # 25342
message("Match rate: ", round(mean(fda_with_atc$atc_matched) * 100, 2), "%") # 49.69%
