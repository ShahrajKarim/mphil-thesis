# scripts/fda_approval_data.r

# Pipeline overview:
# 1. Extract FDA approval data from Drugs@FDA:
#    - Identify first approval dates (NDA/ANDA/BLA) from Submissions.txt.
#    - Merge with Applications.txt for sponsor (firm) and drug information.
#    - Add active ingredients from Products.txt.
#
# 2. Clean and standardise firm identifiers:
#    - Normalise sponsor names (case, punctuation, corp. suffixes).
#    - Optionally apply fuzzy matching for consistent firm-level grouping.
#
# 3. Add therapeutic classification:
#    - Use RxNorm/RxClass API to map active ingredients to ATC classes.
#
# 4. Construct final approvals dataset

# libraries
library(readr)
library(dplyr)
library(janitor)
library(stringr)
library(tidyr)

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
message("Matched rows: ", sum(fda_with_atc$atc_matched)) # 23877
message("Unmatched rows: ", sum(!fda_with_atc$atc_matched)) # 26491
message("Match rate: ", round(mean(fda_with_atc$atc_matched) * 100, 2), "%") # 47.41%
