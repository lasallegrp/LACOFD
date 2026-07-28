library(openxlsx)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
# ============================================================
# 1. File paths
# ============================================================

input_file <- paste0(
  "/quobyte/lasallegrp/projects/LACOFD/",
  "data/raw/wet_lab/dna_isolation/la_county/",
  "LACoFD_LaSalle_Sullivan_DNA_Isolation29June2026_GK.xlsx"
)

output_dir <- paste0(
  "/quobyte/lasallegrp/projects/LACOFD/",
  "data/processed/wet_lab/dna_isolation/la_county"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Read raw workbook
# ============================================================

la_raw <- read.xlsx(
  input_file,
  check.names = FALSE
)

cat("Raw dimensions:", nrow(la_raw), "rows x", ncol(la_raw), "columns\n")

head(la_raw, 10)
# ============================================================
# 3. Helper function
# ============================================================

first_non_missing <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA)
  }

  x[[1]]
}

first_positive <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x) & x > 0]

  if (length(x) == 0) {
    return(NA_real_)
  }

  x[[1]]
}

# ============================================================
# 4. Rename and standardise variables
# ============================================================

la_clean_long <- la_raw %>%
  rename(
    study_name = `SHL.Study.Name`,
    sample_id = `Globally.Unique.Sample.ID`,
    study_site_original = `SHL.Study.Site`,
    participant_id = `SHL.Participant.ID`,
    timepoint_original = `SHL.Time.Point`,
    aliquot_id = `Unique.Aliquot.ID`,
    aliquot_type = `SHL.Aliquot.Type`,
    current_amount_ml = `Current.Amount`,
    aliquot_status = `Aliquot.Status`,
    thaw_count = `Thaws`,
    plasma_volume_a_ul = `Plamsa.volume.A(even)(ul)`,
    plasma_volume_b_ul = `Plamsa.volume.B(odd)(ul)`,
    starting_plasma_ul = `Total.Plasma.volume(ul)`,
    pbs_added_ul = `1X.PBS.needeed`,
    extraction_date_excel = `Date.of.Isolation`,
    extraction_batch = `Batch.of.isolation`,
    dna_concentration_ng_ul = `DNA.conc(ng/ul)`,
    dna_elution_volume_ul = `Vol(ul)`,
    reported_total_dna_yield_ng = `Total.DNA`
  ) %>%
  mutate(
    site = "la_county",
    extraction_operator = "george_kuodza",

    sample_id = as.character(sample_id),
    participant_id = as.character(participant_id),
    aliquot_id = as.character(aliquot_id),

    timepoint_original = str_trim(timepoint_original),

    timepoint = case_when(
      str_to_lower(timepoint_original) == "baseline" ~ "pre",
      str_to_lower(timepoint_original) == "pre-exposure" ~ "pre",
      str_to_lower(timepoint_original) == "post-exposure" ~ "post",
      TRUE ~ NA_character_
    ),

    extraction_date = as.Date(
      suppressWarnings(as.numeric(extraction_date_excel)),
      origin = "1899-12-30"
    ),

    across(
      c(
        current_amount_ml,
        thaw_count,
        plasma_volume_a_ul,
        plasma_volume_b_ul,
        starting_plasma_ul,
        pbs_added_ul,
        extraction_batch,
        dna_concentration_ng_ul,
        dna_elution_volume_ul,
        reported_total_dna_yield_ng
      ),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

  head(la_clean_long)

# Save the cleaned aliquot-level file for provenance
write_tsv(
  la_clean_long,
  file.path(
    output_dir,
    "la_county_dna_isolation_aliquot_level_cleaned.tsv"
  ),
  na = ""
)

# ============================================================
# 6. Calculate final analysis variables
# ============================================================

la_extractions <- la_clean_long %>%
  mutate(
    starting_plasma_ml = starting_plasma_ul / 1000,

    final_processing_volume_ul =
      starting_plasma_ul + pbs_added_ul,

    calculated_total_dna_yield_ng =
      dna_concentration_ng_ul * dna_elution_volume_ul,

    total_dna_difference_ng =
      reported_total_dna_yield_ng -
      calculated_total_dna_yield_ng,

    dna_yield_ng_per_ml_plasma = if_else(
      !is.na(starting_plasma_ml) &
        starting_plasma_ml > 0 &
        !is.na(calculated_total_dna_yield_ng),
      calculated_total_dna_yield_ng / starting_plasma_ml,
      NA_real_
    )

  ) %>%
  select(
    sample_id,
    participant_id,
    site,
    timepoint,
    extraction_batch,
    extraction_date,
    extraction_operator,
    study_name,
    aliquot_type,
    aliquot_id,
    plasma_volume_a_ul,
    plasma_volume_b_ul,
    starting_plasma_ul,
    starting_plasma_ml,
    pbs_added_ul,
    final_processing_volume_ul,
    dna_concentration_ng_ul,
    dna_elution_volume_ul,
    reported_total_dna_yield_ng,
    calculated_total_dna_yield_ng,
    total_dna_difference_ng,
    dna_yield_ng_per_ml_plasma
  )

head(la_extractions)
dim(la_extractions)
str(la_extractions)


# ============================================================
# 7. Save analysis-ready LA County extraction dataset
# ============================================================

analysis_output_file <- file.path(
  output_dir,
  "la_county_dna_isolation_extraction_level.tsv"
)

write_tsv(
  la_extractions,
  analysis_output_file,
  na = ""
)

cat(
  "\nSaved analysis-ready dataset to:\n",
  analysis_output_file,
  "\n"
)

write_tsv(
  la_clean_long,
  file.path(
    output_dir,
    "la_county_dna_isolation_cleaned.tsv"
  ),
  na = ""
)

write_tsv(
  la_extractions,
  file.path(
    output_dir,
    "la_county_dna_isolation_analysis_ready.tsv"
  ),
  na = ""
)