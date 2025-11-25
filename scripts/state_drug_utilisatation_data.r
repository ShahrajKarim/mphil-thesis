# scripts/state_drug_utilisation_data.R

# libraries
library(httr2)
library(jsonlite)
library(purrr)
library(tibble)
library(dplyr)
library(readr)

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
  results_list,
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
    filter(year >= 2005)

# Download all CSVs

for (i in seq_len(nrow(datasets))) {

  url <- datasets$download_url[i]
  title <- datasets$title[i]

  message("Downloading: ", title, " from ", url)

  df <- read_csv(url, show_col_types = FALSE)

  write.csv(df, paste0("raw_data/", title, ".csv"), row.names = FALSE)
}