# scripts/medicaid_expansion_data.r

# libraries
library(dplyr)
library(readr)
library(lubridate)
library(stringr)

# Tidy up rollout data

rollout_data <- read_csv(
  ("raw_data/medicaid_expansion_status.csv"),
  skip = 2,
  col_select = 1:5,
  show_col_types = FALSE
) |>
  mutate(
    `Expansion Implementation Date` = na_if(`Expansion Implementation Date`, ""),
    expansion_date = as.Date(`Expansion Implementation Date`, format = "%m/%d/%Y")
  )

rollout_data <- rollout_data |>
  mutate(
    `Expansion Implementation Date` = ifelse(
      `Expansion Implementation Date` %in% c("", "NA", "N/A", "—"), NA, `Expansion Implementation Date`),
    `Expansion Adopted Through Ballot Initiative` = ifelse(
      `Expansion Adopted Through Ballot Initiative` %in% c("", "NA", "N/A", "—"), NA, `Expansion Adopted Through Ballot Initiative`),
    expansion_date = coalesce(
      mdy(`Expansion Implementation Date`, quiet = TRUE),
      as.Date(suppressWarnings(as.numeric(`Expansion Implementation Date`)), origin = "1899-12-30")
    ),
    Location = gsub("\\s+", " ", trimws(Location))
  ) |>
  filter(Location %in% c(state.name, "District of Columbia")) |>
  mutate(
    state_abbrev = case_when(
      Location == "District of Columbia" ~ "DC",
      TRUE ~ state.abb[match(Location, state.name)]
    )
  ) |>
  mutate(
    expansion_year = year(expansion_date),
    expansion_quarter = quarter(expansion_date)
  ) |>
  select(state_abbrev, expansion_date, expansion_year, expansion_quarter)

write.csv(
  rollout_data,
  "aux_data/medicaid_expansion_dates_cleaned.csv",
  row.names = FALSE
)

# Build a dataset showing intensity of expansion over time

quarters <- seq(
  from = as.Date("2014-01-01"),
  to   = as.Date("2023-07-01"),
  by   = "quarter"
)

quarter_df <- tibble(
  quarter_date = quarters,
  year = year(quarter_date),
  quarter  = quarter(quarter_date)
)

panel <- expand_grid(
  state_abbrev = rollout_data$state_abbrev,
  quarter_df
)

panel <- panel |>
  left_join(
    rollout_data |> select(state_abbrev, expansion_date),
    by = "state_abbrev"
  )

panel <- panel |>
  mutate(
    expanded = ifelse(!is.na(expansion_date) & expansion_date <= quarter_date, 1, 0)
  )

prop_expanded_df <- panel |>
  group_by(year, quarter) |>
  summarise(
    prop_expanded = mean(expanded, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
          prop_expanded_df,
          "aux_data/medicaid_expansion_proportion_by_quarter.csv",
          row.names = FALSE)

prop_expanded_df <- panel |>
  group_by(year) |>
  summarise(
    prop_expanded = mean(expanded, na.rm = TRUE)
  )

write.csv(
          prop_expanded_df,
          "aux_data/medicaid_expansion_proportion_by_year.csv",
          row.names = FALSE)