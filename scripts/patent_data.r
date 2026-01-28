# scripts/patent_data.r

# libraries
library(vroom)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(lubridate)
library(readxl)
library(stringdist)

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

pharma_patents_full <- pharma_patents_full |>
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

# Create a broad category mapping for CPC groups
# Source: https://www.uspto.gov/web/patents/classification/cpc/pdf/cpc-scheme-A61P.pdf

pharma_patents_full <- pharma_patents_full |>
  mutate(
    therapeutic_class = case_when(
      str_detect(cpc_group, "^A61P1")  ~ "Gastrointestinal",
      str_detect(cpc_group, "^A61P3")  ~ "Metabolic",
      str_detect(cpc_group, "^A61P5")  ~ "Endocrine",
      str_detect(cpc_group, "^A61P7")  ~ "Blood",
      str_detect(cpc_group, "^A61P9")  ~ "Cardiovascular",
      str_detect(cpc_group, "^A61P11") ~ "Respiratory",
      str_detect(cpc_group, "^A61P13") ~ "Urinary",
      str_detect(cpc_group, "^A61P15") ~ "Genital/Sexual",
      str_detect(cpc_group, "^A61P17") ~ "Dermatology",
      str_detect(cpc_group, "^A61P19") ~ "Skeletal",
      str_detect(cpc_group, "^A61P21") ~ "Muscular/Neuromuscular",
      str_detect(cpc_group, "^A61P23") ~ "Anaesthetics",
      str_detect(cpc_group, "^A61P25") ~ "CNS/Neurology",
      str_detect(cpc_group, "^A61P27") ~ "Ophthalmic/Otological",
      str_detect(cpc_group, "^A61P29") ~ "Anti-inflammatory/Analgesics",
      str_detect(cpc_group, "^A61P31") ~ "Anti-infectives",
      str_detect(cpc_group, "^A61P33") ~ "Anti-parasitics",
      str_detect(cpc_group, "^A61P35") ~ "Anti-neoplastics",
      str_detect(cpc_group, "^A61P37") ~ "Allergy/Immune disorders",
      str_detect(cpc_group, "^A61P39") ~ "Protective/Detoxification",
      str_detect(cpc_group, "^A61P41") ~ "Surgical",
      str_detect(cpc_group, "^A61P43") ~ "Other",
    )
  )

# Create a sub-level category mapping within therapeutic categories
# First load mapping file and then join together

cpc_mappings <- read.csv("aux_data/cpc_atc_mapped.csv") |>
  select(
    cpc_group,
    sub_class,
    ATC1
  ) |>
  mutate(sub_class = ifelse(sub_class == "-", NA, sub_class))

pharma_patents_full <- pharma_patents_full |>
  left_join(cpc_mappings, by = "cpc_group")

# Generate quarter variable

pharma_patents_full <- pharma_patents_full |>
  mutate(
    patent_quarter = quarter(patent_date)
  )

# save the cleaned pharmaceutical patents data

write_csv(
  pharma_patents_full,
  file = paste0("processed_data/patents/", "pharma_patents_clean.csv")
)

# Create firm name mapping based on NDC labeler names

ndc_firm_mapping <- read_csv(
  "aux_data/ndc_firm_mapping.csv",
  col_types = cols(
    labeler_code_pad = col_character()
  )
)

patent_firms <- pharma_patents_full |>
  select(firm) |>
  distinct()

# Define the set of suffixes to remove
SUFFIXES <- "\\b(CORP|CORPORATION|INC|CO|COMPANY|LLC|LTD|LP|PLC|AG|SA|PHARMA|PHARMACEUTICALS|LLP|S A|GMBH|S P A|S R L)\\b"

# Clean and distinct patent firm names, PRESERVING the original name
patent_firms_full_clean <- patent_firms |>
  rename(original_firm_name = firm) |>
  mutate(
    firm_clean = original_firm_name |>
      toupper() |>
      gsub("[^A-Z0-9 ]", " ", x = _) |>
      gsub(SUFFIXES, "", x = _) |>
      gsub("\\s+", " ", x = _) |>
      trimws()
  ) |>
  filter(!is.na(firm_clean), firm_clean != "") |>
  select(original_firm_name, firm_clean) |>
  distinct()

# Clean NDC firm names
ndc_firms_clean <- ndc_firm_mapping |>
  rename(ndc_firm_id = labeler_code_pad) |>
  mutate(
    firm_utf8 = iconv(firm, from = "", to = "UTF-8", sub = ""),
    firm_clean = firm_utf8 |>
      toupper() |>
      gsub("[^A-Z0-9 ]", " ", x = _) |>
      gsub(SUFFIXES, "", x = _) |>
      gsub("\\s+", " ", x = _) |>
      trimws()
  ) |>
  filter(!is.na(firm_clean), firm_clean != "") |>
  # Group and summarize to get one representative NDC firm ID per clean name
  group_by(firm_clean) |>
  summarise(
    ndc_firm_id = first(ndc_firm_id),
    .groups = "drop"
  ) |>
  select(firm_clean, ndc_firm_id)

# Initialize the final map list
final_map_list <- list()

# First use exact clean match

stage1_matches <- patent_firms_full_clean |>
  inner_join(ndc_firms_clean, by = "firm_clean") |>
  mutate(match_type = "Exact_Clean")

final_map_list[[1]] <- stage1_matches

# Firms remaining for Stage 2
unmatched_after_1 <- patent_firms_full_clean |>
  anti_join(stage1_matches, by = "firm_clean")


# Use root fallback matching if stage 1 fails

# Build NDC firm root table (needs to be calculated once)
ndc_roots <- ndc_firms_clean |>
  mutate(
    root = str_match(firm_clean, "\\b([A-Z0-9]{4,})\\b")[, 2],
    root = ifelse(is.na(root), str_extract(firm_clean, "^[A-Z0-9]+"), root)
  ) |>
  filter(!is.na(root)) |>
  group_by(root) |>
  summarise(
    ndc_firm_id_root = first(ndc_firm_id),
    .groups = "drop"
  ) |>
  select(root, ndc_firm_id = ndc_firm_id_root)

# Merge unmatched with root-based mapping
stage2_matches <- unmatched_after_1 |>
  mutate(
    root = str_match(firm_clean, "\\b([A-Z0-9]{4,})\\b")[, 2],
    root = ifelse(is.na(root), str_extract(firm_clean, "^[A-Z0-9]+"), root)
  ) |>
  inner_join(ndc_roots, by = "root") |>
  mutate(match_type = "Root_Fallback") |>
  select(original_firm_name, firm_clean, ndc_firm_id, match_type)

final_map_list[[2]] <- stage2_matches

# Firms remaining for Stage 3
unmatched_after_2 <- unmatched_after_1 |>
  anti_join(stage2_matches, by = "firm_clean")

# Conduct fuzzy matching on remaining firms

# Prepare data for cross-join
unmatched_names <- unmatched_after_2 |> select(firm_clean)
ndc_names <- ndc_firms_clean |> select(ndc_clean = firm_clean, ndc_firm_id)

# Create a data frame containing all possible match pairs
all_pairs <- unmatched_names |>
  mutate(temp_key = 1) |>
  left_join(
    ndc_names |> mutate(temp_key = 1),
    by = "temp_key",
    relationship = "many-to-many"
  ) |>
  select(-temp_key)

# Calculate Jaro-Winkler Similarity (>= 0.90 threshold)
stage3_matches <- all_pairs |>
  mutate(
    jw_score = stringsim(firm_clean, ndc_clean, method = "jw")
  ) |>
  filter(jw_score >= 0.90) |>
  group_by(firm_clean) |>
  slice_max(jw_score, with_ties = FALSE) |> # Select best score
  ungroup() |>
  mutate(match_type = paste0("Fuzzy_", round(jw_score, 2))) |>
  # Re-join original_firm_name from unmatched_after_2
  left_join(
    select(unmatched_after_2, original_firm_name, firm_clean) |> distinct(),
    by = "firm_clean"
  ) |>
  select(original_firm_name, firm_clean, ndc_firm_id, match_type)

final_map_list[[3]] <- stage3_matches


# Final combination and reporting

# Combine all stages into the final mapping table
patent_firm_ndc_firm_map <- bind_rows(final_map_list) |>
  rename(firm = original_firm_name) |>
  select(firm, ndc_firm_id)

# Final stats calculation
total_patents <- nrow(patent_firms)
matches <- nrow(patent_firm_ndc_firm_map)

# Calculate final merge rate
final_merge_rate_percent <- (num_total_unique_matches_final / num_patent_firms_unique) * 100

# Merge diagnostics
message("Total patent firms: ", total_patents) # 48,143
message("Matched rows: ", matches) # 7,556
message("Match rate: ", round(matches / total_patents * 100, 2), "%") # 27.1%

# Save final mapping file

write_csv(
  patent_firm_ndc_firm_map,
  "aux_data/patent_firm_to_ndc_firm_mapping.csv"
)

# Merge NDC-firm mapping with patent data

pharma_patents_full <- pharma_patents_full |>
  left_join(
    patent_firm_ndc_firm_map,
    by = "firm"
  )

# Collapse patents at firm-year-quarter-ATC level

pharma_patents_collapsed <- pharma_patents_full |>
  filter(patent_year >= 2005) |>
  group_by(
    ndc_firm_id,
    firm,
    patent_year,
    patent_quarter,
    ATC1
  ) |>
  summarise(
    num_patents = n(),
    .groups = "drop"
  )


# Now introduce SDUD data

SDUD_firm_mapped <- readRDS(
  file = "processed_data/state_drug_utilisation_data/SDUD_firm_mapped.rds"
)


# Merge ATC mapping into SDUD data
ndc_atc_mapping <- read_csv(
  "aux_data/ndc_atc_mapping.csv",
  col_types = cols(
    ndc = col_character(),
  )
)

SDUD_firm_mapped <- SDUD_firm_mapped |>
  left_join(ndc_atc_mapping, by = "ndc", relationship = "many-to-many") |>
  select(- rxcui, -ATC_Level2)

# Collapse SDUD data to firm-year-quarter-ATC level

SDUD_firm_ATC <- SDUD_firm_mapped |>
  group_by(
    state,
    firm,
    labeler_code,
    year,
    quarter,
    ATC_Level1
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

# Merge SDUD datasets to obtain resulting SDUD master file

SDUD_master <- SDUD_firm_ATC |>
  left_join(
    SDUD_ATC1,
    by = c("year", "quarter", "ATC_Level1")
  ) |>
  select(-firm)

# Merge with collapsed patents data

patent_sdud_merged <- pharma_patents_collapsed |>
  left_join(
    SDUD_master,
    by = c(
      "ndc_firm_id" = "labeler_code",
      "patent_year" = "year",
      "patent_quarter" = "quarter",
      "ATC1" = "ATC_Level1"
    ),
    relationship = "many-to-many"
  )

# Save final merged dataset

write.csv(
  patent_sdud_merged,
  "processed_data/patents/patents_sdud_merged.csv",
  row.names = FALSE
)


# Now merge patent data with SDUD - PDL data

# Read the pharma patents clean data and collapse

pharma_patents_full <- read_csv(
  file = "processed_data/patents/pharma_patents_clean.csv",
  col_types = cols(
    patent_id = col_character()
  )
)

patent_firm_ndc_firm_map <- read_csv(
  "aux_data/patent_firm_to_ndc_firm_mapping.csv"
)

pharma_patents_full <- pharma_patents_full |>
  left_join(
    patent_firm_ndc_firm_map,
    by = "firm"
  )

pharma_patents_collapsed <- pharma_patents_full |>
  filter(patent_year >= 2005) |>
  group_by(
    ndc_firm_id,
    firm,
    patent_year,
    patent_quarter,
    ATC1
  ) |>
  summarise(
    num_patents = n(),
    .groups = "drop"
  )


# Import ATC mapping for SDUD data

ndc_atc_mapping <- read_csv(
  "aux_data/ndc_atc_mapping.csv",
  col_types = cols(
    ndc = col_character(),
  )
)

# Import SDUD - PDL data

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
  left_join(ndc_atc_mapping, by = "ndc", relationship = "many-to-many") |>
  select(- rxcui, -ATC_Level2) |>
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

# Merge with collapsed patents data

patent_sdud_pdl_merged <- pharma_patents_collapsed |>
  left_join(
    SDUD_PDL,
    by = c(
      "ndc_firm_id" = "labeler_code",
      "patent_year" = "year",
      "patent_quarter" = "quarter",
      "ATC1" = "ATC_Level1"
    ),
    relationship = "many-to-many"
  ) |>
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# Save final merged dataset

write.csv(
  patent_sdud_pdl_merged,
  "processed_data/patents/patents_sdud_pdl_merged.csv",
  row.names = FALSE
)