#!/usr/bin/env Rscript

library(openxlsx)
library(dplyr)
library(stringr)
library(readr)

# ============================================================
# DNA isolation cleaning script
#
# Required arguments:
#   --input
#   --output-dir
#   --site
#   --batch
#   --operator
#
# Example:
#
# Rscript 01_clean_dna_isolation.R \
#   --input data/raw/wet_lab/dna_isolation/la_county/DNAEXT_LA_01/file.xlsx \
#   --output-dir data/processed/wet_lab/dna_isolation/la_county/DNAEXT_LA_01 \
#   --site la_county \
#   --batch DNAEXT_LA_01 \
#   --operator george_kuodza
# ============================================================


# ============================================================
# 1. Parse command-line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, required = TRUE, default = NULL) {

  index <- match(flag, args)

  if (is.na(index)) {
    if (required) {
      stop("Missing required argument: ", flag)
    } else {
      return(default)
    }
  }

  if (index == length(args)) {
    stop("No value supplied for argument: ", flag)
  }

  args[index + 1]
}


input_file <- get_arg("--input")
output_dir <- get_arg("--output-dir")
site <- get_arg("--site")
dna_extraction_batch <- get_arg("--batch")
extraction_operator <- get_arg("--operator")

# Optional output prefix
output_prefix <- get_arg(
  "--prefix",
  required = FALSE,
  default = dna_extraction_batch
)


# ============================================================
# 2. Validate arguments
# ============================================================

if (!file.exists(input_file)) {
  stop("Input file does not exist:\n", input_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n============================================\n")
cat("DNA isolation cleaning\n")
cat("============================================\n")
cat("Input file:       ", input_file, "\n")
cat("Output directory: ", output_dir, "\n")
cat("Site:             ", site, "\n")
cat("DNA extraction:   ", dna_extraction_batch, "\n")
cat("Operator:         ", extraction_operator, "\n")
cat("============================================\n\n")


# ============================================================
# 3. Read raw workbook
# ============================================================

raw_data <- read.xlsx(
  input_file,
  check.names = FALSE
)

cat(
  "Raw dimensions:",
  nrow(raw_data),
  "rows x",
  ncol(raw_data),
  "columns\n\n"
)


# ============================================================
# 4. Check required source columns
# ============================================================

required_columns <- c(
  "SHL.Study.Name",
  "Globally.Unique.Sample.ID",
  "SHL.Participant.ID",
  "SHL.Time.Point",
  "Unique.Aliquot.ID",
  "SHL.Aliquot.Type",
  "Current.Amount",
  "Thaws",
  "Plamsa.volume.A(even)(ul)",
  "Plamsa.volume.B(odd)(ul)",
  "Total.Plasma.volume(ul)",
  "1X.PBS.needeed",
  "Date.of.Isolation",
  "Batch.of.isolation",
  "DNA.conc(ng/ul)",
  "Vol(ul)",
  "Total.DNA"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (length(missing_columns) > 0) {

  stop(
    paste0(
      "The following required columns are missing:\n",
      paste(missing_columns, collapse = "\n")
    )
  )
}


# ============================================================
# 5. Rename and standardize variables
# ============================================================

clean_data <- raw_data %>%
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

    # Preserve the numeric/raw batch recorded in workbook
    source_extraction_batch = `Batch.of.isolation`,

    dna_concentration_ng_ul = `DNA.conc(ng/ul)`,
    dna_elution_volume_ul = `Vol(ul)`,
    reported_total_dna_yield_ng = `Total.DNA`
  ) %>%

  mutate(

    # Project-level metadata supplied as arguments
    site = site,

    dna_extraction_batch =
      dna_extraction_batch,

    extraction_operator =
      extraction_operator,

    # IDs
    sample_id =
      as.character(sample_id),

    participant_id =
      as.character(participant_id),

    aliquot_id =
      as.character(aliquot_id),

    # Standardize timepoint
    timepoint_original =
      str_trim(timepoint_original),

    timepoint = case_when(

      str_to_lower(timepoint_original) %in%
        c("baseline", "pre-exposure", "pre exposure") ~
        "pre",

      str_to_lower(timepoint_original) %in%
        c("post-exposure", "post exposure") ~
        "post",

      TRUE ~ NA_character_
    ),

    # Convert Excel date
    extraction_date = as.Date(
      suppressWarnings(
        as.numeric(extraction_date_excel)
      ),
      origin = "1899-12-30"
    ),

    # Numeric variables
    across(
      c(
        current_amount_ml,
        thaw_count,
        plasma_volume_a_ul,
        plasma_volume_b_ul,
        starting_plasma_ul,
        pbs_added_ul,
        source_extraction_batch,
        dna_concentration_ng_ul,
        dna_elution_volume_ul,
        reported_total_dna_yield_ng
      ),
      ~ suppressWarnings(as.numeric(.x))
    )
  )


# ============================================================
# 6. Calculate analysis variables
# ============================================================

analysis_data <- clean_data %>%

  mutate(

    starting_plasma_ml =
      starting_plasma_ul / 1000,

    final_processing_volume_ul =
      starting_plasma_ul +
      pbs_added_ul,

    calculated_total_dna_yield_ng =
      dna_concentration_ng_ul *
      dna_elution_volume_ul,

    total_dna_difference_ng =
      reported_total_dna_yield_ng -
      calculated_total_dna_yield_ng,

    dna_yield_ng_per_ml_plasma =
      if_else(
        !is.na(starting_plasma_ml) &
          starting_plasma_ml > 0 &
          !is.na(calculated_total_dna_yield_ng),

        calculated_total_dna_yield_ng /
          starting_plasma_ml,

        NA_real_
      )
  ) %>%

  select(
    sample_id,
    participant_id,

    site,
    dna_extraction_batch,
    source_extraction_batch,
    extraction_date,
    extraction_operator,

    timepoint,

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


# ============================================================
# 7. QC checks
# ============================================================

cat("\nAnalysis dimensions:\n")
cat(
  nrow(analysis_data),
  "rows x",
  ncol(analysis_data),
  "columns\n"
)

cat("\nTimepoint counts:\n")

print(
  analysis_data %>%
    count(timepoint)
)


# Unknown timepoints
unknown_timepoints <- clean_data %>%
  filter(is.na(timepoint)) %>%
  distinct(timepoint_original)

if (nrow(unknown_timepoints) > 0) {

  warning(
    "Unrecognized timepoint values were found."
  )

  print(unknown_timepoints)
}


# Check final processing volume
cat("\nFinal processing volumes:\n")

print(
  analysis_data %>%
    count(final_processing_volume_ul)
)


# Check reported vs calculated yield
yield_check <- analysis_data %>%
  filter(
    !is.na(total_dna_difference_ng),
    abs(total_dna_difference_ng) > 0.01
  )

cat(
  "\nSamples where reported and calculated ",
  "DNA yield differ by >0.01 ng: ",
  nrow(yield_check),
  "\n",
  sep = ""
)


# ============================================================
# 8. Save files
# ============================================================

clean_output <- file.path(
  output_dir,
  paste0(
    output_prefix,
    "_dna_isolation_cleaned.tsv"
  )
)

analysis_output <- file.path(
  output_dir,
  paste0(
    output_prefix,
    "_dna_isolation_analysis_ready.tsv"
  )
)


write_tsv(
  clean_data,
  clean_output,
  na = ""
)

write_tsv(
  analysis_data,
  analysis_output,
  na = ""
)


cat("\n============================================\n")
cat("Files successfully written\n")
cat("============================================\n")
cat("Cleaned data:\n", clean_output, "\n\n")
cat("Analysis-ready data:\n", analysis_output, "\n")
cat("============================================\n")