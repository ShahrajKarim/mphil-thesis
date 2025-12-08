# scripts/fda_approval_data.r

# Pipeline overview:
# 1. Extract FDA approval data from Drugs@FDA:
#    - Pull NDA/ANDA/BLA approvals with sponsor (firm), active ingredients, and approval year.
#
# 2. Link approvals to Orange Book:
#    - Join on NDA/ANDA numbers to obtain product numbers and NDC labeler codes.
#
# 3. Map to Medicaid SDUD:
#    - Use NDC labeler codes; apply existing fuzzy-matching to harmonise firm names.
#
# 4. Add therapeutic classification:
#    - Map active ingredients to ATC classes via RxNorm (RxClass API).

# libraries
library(readr)
library(dplyr)
library(janitor)
library(stringr)

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


# --- Load Orange Book data --- #

# --------------------------------------------------------------------------- #
# NOTE TO SELF — Orange Book components and their purpose
#
# • ob_products (products.txt)
#     - Links NDA/ANDA (appl_no) → product_no → NDC labeler codes.
#     - This is the only OB dataset needed to merge FDA approvals into SDUD.
#
# • ob_patents (patent.txt)
#     - Optional.
#     - Lists Orange Book–reported patents for small-molecule drugs.
#     - Useful for analysing patent protection, innovation incentives, or
#       modelling expected exclusivity durations.
#     - NOT needed for SDUD linkage.
#
# • ob_exclusivity (exclusivity.txt)
#     - Optional.
#     - Contains FDA-granted exclusivity periods.
#     - Useful for research on market protection, competition timing, and
#       estimating effective post-approval monopoly periods.
#     - NOT needed for SDUD linkage.
#
# Summary:
#     For now: only load and use ob_products.
#     Later: load patents/exclusivity if investigating protection mechanisms.
# --------------------------------------------------------------------------- #


orange_book_path <- "raw_data/orange_book"

ob_products <- read_delim(
  file.path(orange_book_path, "products.txt"),
  delim = "~",
  col_types = cols(.default = "c")
) |> clean_names()

ob_exclusivity <- read_delim(
  file.path(orange_book_path, "exclusivity.txt"),
  delim = "~",
  col_types = cols(.default = "c")
) |> clean_names()

ob_patents <- read_delim(
  file.path(orange_book_path, "patent.txt"),
  delim = "~",
  col_types = cols(.default = "c")
) |> clean_names()

# Merge FDA approval data with Orange Book products

fda_ob <- fda_drugs |>
  left_join(ob_products, by = "appl_no", relationship = "many-to-many")
