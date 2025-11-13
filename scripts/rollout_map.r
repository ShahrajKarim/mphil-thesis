# scripts/rollout_map.R

# libraries
    library(dplyr)
    library(readr)
    library(sf)
    library(leaflet)
    library(tigris)
    library(lubridate)
    library(viridisLite)
    library(here)
    library(htmlwidgets)
    library(rmarkdown) 

# global options / params
options(tigris_use_cache = TRUE, tigris_class = "sf")

# Import rollout data

rollout_data <- read_csv(
  here("raw_data", "medicaid_expansion_status.csv"),
  skip = 2,
  col_select = 1:5,
  show_col_types = FALSE
) |>
  mutate(
    `Expansion Implementation Date` = na_if(`Expansion Implementation Date`, ""),
    expansion_date = as.Date(`Expansion Implementation Date`, format = "%m/%d/%Y")
  )

# Clean dates and drop non-state rows (keep 50 states + DC)
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
  filter(Location %in% c(state.name, "District of Columbia"))

  # Shapefiles for States
  states_sf <- states(cb = TRUE, year = 2019, progress_bar = FALSE) |>
    st_as_sf()

  states_sf <- st_transform(states_sf, 4326)

  # Restrict shapes to only those states present in rollout_data, then join attributes
  states_keep <- states_sf |>
    semi_join(rollout_data, by = c("NAME" = "Location"))

  mapdf <- states_keep |>
    left_join(rollout_data, by = c("NAME" = "Location"))

  pal <- colorNumeric(
  palette = cividis(256),
  domain = range(as.numeric(mapdf$expansion_date), na.rm = TRUE)
)

# Plot map
map <- leaflet(mapdf, options = leafletOptions(zoomControl = TRUE), width = "100%") |>
    addProviderTiles("CartoDB.Positron") |>
    addPolygons(
      fillColor   = ~ifelse(is.na(expansion_date), "#FF0000", pal(as.numeric(expansion_date))),
      color       = "#666",
      weight      = 1,
      opacity     = 1,
      fillOpacity = 0.3,
      label = ~paste0(
        NAME,
        ifelse(!is.na(expansion_date), paste0(": ", format(expansion_date, "%b %d, %Y")), ": No expansion")
      ),
      highlightOptions = highlightOptions(weight = 2, color = "#FFFF00", bringToFront = TRUE)
    ) |>
    addLegend(
      "bottomleft",
      pal = pal,
      values = ~as.numeric(expansion_date[!is.na(expansion_date)]),
      title = "Medicaid expansion date",
      labFormat = function(type, cuts, p) {
        format(as.Date(cuts, origin = "1970-01-01"), "%b %Y")
      }
    ) |>
    setView(lng = -98.0, lat = 38.5, zoom = 3)

# Save the widget as a self-contained HTML file under the project figures directory
out_dir <- here("figures")

saveWidget(
  map,
  file = file.path(out_dir, "rollout_map.html"),
  selfcontained = TRUE
)