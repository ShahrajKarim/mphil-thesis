# scripts/state_drug_utilisation_data.R

# libraries
library(httr2)
library(jsonlite)
library(purrr)
library(tibble)
library(dplyr)
library(readr)
library(janitor)
library(tidyr)

# Query Medicaid API for State Drug Utilization Data (SDUD) datasets

medicaid_url <- "https://data.medicaid.gov/api/1/search"

query_object <- request(medicaid_url) |>
  req_url_query(
    theme = "State Drug Utilization",
    "page-size" = 100
  ) |>
  req_perform()

medicaid_json <- resp_body_json(query_object, simplifyVector = FALSE)
medicaid_list <- medicaid_json$results

# Build a tibble with identifier, title, theme, download_url

datasets <- map_dfr(
  medicaid_list,
  function(x) {
    tibble(
      identifier = x$identifier,
      title = x$title,
      theme = paste(x$theme, collapse = "; "),
      download_url = x$distribution[[1]]$downloadURL
    )
  }
)

datasets <- datasets |>
  filter(grepl("State Drug Utilization Data", title)) |>
  arrange(title) |>
  mutate(
    year = as.integer(substr(title, nchar(title) - 3, nchar(title)))
    ) |>
    filter(year >= 2005) |>
    mutate(title = str_replace_all(title, "State Drug Utilization Data ", "SDUD_"))

# Download all CSVs

for (i in seq_len(nrow(datasets))) {

  url <- datasets$download_url[i]
  title <- datasets$title[i]

  message("Downloading: ", title, " from ", url)

  temp <- tempfile(fileext = ".csv")
  
  download.file(url, temp, mode = "wb", quiet = TRUE)

  df <- read_csv(temp, show_col_types = FALSE)


  write.csv(df, paste0("raw_data/state_drug_utilisation_data/", title, ".csv"), row.names = FALSE)
}

# Clean and merge the datasets

for (i in 2005:2024) {

  file_path <- paste0("raw_data/state_drug_utilisation_data/SDUD_", i, ".csv")
  df <- read_csv(file_path, show_col_types = FALSE)

  # Select relevant columns and rename
  df <- df |>
    clean_names() |>
    mutate(
      year = as.integer(i),
      ndc = as.character(ndc),
      package_size = as.numeric(package_size),
      units_reimbursed = as.numeric(units_reimbursed),
      number_of_prescriptions = as.numeric(number_of_prescriptions),
      total_amount_reimbursed = as.numeric(total_amount_reimbursed),
      medicaid_amount_reimbursed = as.numeric(medicaid_amount_reimbursed),
      non_medicaid_amount_reimbursed = as.numeric(non_medicaid_amount_reimbursed)
    ) |>
    filter(!suppression_used) |>
    select(
      -suppression_used,
    )

    if (i == 2005) {

      sdud_full <- df

    } else {
      
      sdud_full <- bind_rows(sdud_full, df)

    }

}

# Save the cleaned full dataset
saveRDS(
  sdud_full,
  file = "processed_data/state_drug_utilisation_data/SDUD_full.rds"
)

# Load NDC directory file

ndc_raw <- fromJSON("raw_data/FDA/drug-ndc-0001-of-0001.json", flatten = TRUE)

ndc <- ndc_raw$results |>
  as_tibble() |>
 select(
    -marketing_start_date,
    -marketing_end_date
  ) |>
  unnest(packaging) |>
  filter(nchar(package_ndc) >= 9) |>
  mutate(
    parts = strsplit(package_ndc, "-"),
    labeler  = str_pad(map_chr(parts, 1), width = 5, pad = "0"),
    product  = str_pad(map_chr(parts, 2), width = 4, pad = "0"),
    package  = str_pad(map_chr(parts, 3), width = 2, pad = "0"),
    ndc_clean = paste0(labeler, product, package)
  ) |>
  distinct(ndc_clean, .keep_all = TRUE) |>
  select(
    ndc_clean,
    generic_name,
    labeler_name,
    brand_name
  )

# Merge NDC directory with SDUD data

sdud_full <- readRDS("processed_data/state_drug_utilisation_data/SDUD_full.rds")

sdud_full <- sdud_full |>
  left_join(ndc, by = c("ndc" = "ndc_clean"))

# Merge diagnostics
sdud_full <- sdud_full |>
  mutate(ndc_matched = !is.na(generic_name))

message("Total rows: ", nrow(sdud_full))
message("Matched rows: ", sum(sdud_full$ndc_matched))
message("Unmatched rows: ", sum(!sdud_full$ndc_matched))
message("Match rate: ", round(mean(sdud_full$ndc_matched) * 100, 2), "%")

# Optional: save unmatched NDCs
unmatched_ndcs <- sdud_full |>
  filter(!ndc_matched) |>
  distinct(ndc)

write_csv(unmatched_ndcs, "processed_data/unmatched_ndcs.csv")