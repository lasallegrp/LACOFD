#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# DNA isolation paired pre/post plotting script
#
# Required:
#   --input
#   --output-dir
#   --label
#   --prefix
#
# Optional:
#   --concentration-threshold
#   --yield-threshold
#   --width
#   --height
# ============================================================


# ============================================================
# 1. Parse command-line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(
    flag,
    required = TRUE,
    default = NULL
) {

  index <- match(flag, args)

  if (is.na(index)) {

    if (required) {
      stop(
        "Missing required argument: ",
        flag
      )
    }

    return(default)
  }

  if (index == length(args)) {
    stop(
      "No value supplied for argument: ",
      flag
    )
  }

  args[index + 1]
}


input_file <- get_arg("--input")

output_dir <- get_arg(
  "--output-dir"
)

dataset_label <- get_arg(
  "--label"
)

output_prefix <- get_arg(
  "--prefix"
)


# Optional settings
concentration_threshold <- as.numeric(
  get_arg(
    "--concentration-threshold",
    required = FALSE,
    default = "0.01"
  )
)

yield_threshold <- as.numeric(
  get_arg(
    "--yield-threshold",
    required = FALSE,
    default = "0.1"
  )
)

plot_width <- as.numeric(
  get_arg(
    "--width",
    required = FALSE,
    default = "14"
  )
)

plot_height <- as.numeric(
  get_arg(
    "--height",
    required = FALSE,
    default = "7"
  )
)


# ============================================================
# 2. Validate files
# ============================================================

if (!file.exists(input_file)) {
  stop(
    "Input file does not exist:\n",
    input_file
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("\n============================================\n")
cat("DNA isolation plotting\n")
cat("============================================\n")
cat("Input:  ", input_file, "\n")
cat("Label:  ", dataset_label, "\n")
cat("Prefix: ", output_prefix, "\n")
cat("Output: ", output_dir, "\n")
cat("============================================\n\n")


# ============================================================
# 3. Read data
# ============================================================

extractions <- read_tsv(
  input_file,
  show_col_types = FALSE
)

cat(
  "Input dimensions:",
  nrow(extractions),
  "rows x",
  ncol(extractions),
  "columns\n"
)


# ============================================================
# 4. Prepare plotting data
# ============================================================

plot_data <- extractions %>%

  filter(
    !is.na(participant_id),
    !is.na(timepoint)
  ) %>%

  mutate(

    timepoint = case_when(
      timepoint == "pre" ~
        "Pre-exposure",

      timepoint == "post" ~
        "Post-exposure",

      TRUE ~
        as.character(timepoint)
    ),

    timepoint = factor(
      timepoint,
      levels = c(
        "Pre-exposure",
        "Post-exposure"
      )
    )
  )


# ============================================================
# 5. Check duplicate participant/timepoint records
# ============================================================

duplicate_check <- plot_data %>%
  count(
    participant_id,
    timepoint
  ) %>%
  filter(n > 1)


if (nrow(duplicate_check) > 0) {

  print(duplicate_check)

  stop(
    "Some participants have more than ",
    "one record per timepoint."
  )
}


# ============================================================
# 6. DNA concentration plot
# ============================================================

paired_concentration_ids <-
  plot_data %>%

  filter(
    !is.na(dna_concentration_ng_ul)
  ) %>%

  distinct(
    participant_id,
    timepoint
  ) %>%

  count(
    participant_id,
    name = "n_timepoints"
  ) %>%

  filter(
    n_timepoints == 2
  ) %>%

  pull(participant_id)


paired_concentration_data <-
  plot_data %>%

  filter(
    participant_id %in%
      paired_concentration_ids
  )


concentration_change <-
  paired_concentration_data %>%

  select(
    participant_id,
    timepoint,
    dna_concentration_ng_ul
  ) %>%

  pivot_wider(
    names_from = timepoint,
    values_from =
      dna_concentration_ng_ul
  ) %>%

  mutate(

    concentration_change_ng_ul =
      `Post-exposure` -
      `Pre-exposure`,

    change_direction =
      case_when(

        concentration_change_ng_ul >
          concentration_threshold ~
          "Increase",

        concentration_change_ng_ul <
          -concentration_threshold ~
          "Decrease",

        TRUE ~
          "No meaningful change"
      )
  )


paired_concentration_data <-
  paired_concentration_data %>%

  left_join(
    concentration_change %>%
      select(
        participant_id,
        concentration_change_ng_ul,
        change_direction
      ),
    by = "participant_id"
  )


concentration_summary <-
  concentration_change %>%

  count(change_direction) %>%

  mutate(
    percent =
      round(
        100 * n / sum(n),
        1
      )
  )


cat(
  "\nPaired concentration participants: ",
  length(paired_concentration_ids),
  "\n",
  sep = ""
)

print(concentration_summary)


p1 <- ggplot(
  paired_concentration_data,
  aes(
    x = timepoint,
    y = dna_concentration_ng_ul,
    group = participant_id
  )
) +

  geom_boxplot(
    aes(group = timepoint),
    width = 0.45,
    outlier.shape = NA,
    alpha = 0.15,
    colour = "black"
  ) +

  geom_line(
    aes(
      colour = change_direction
    ),
    alpha = 0.60,
    linewidth = 0.7
  ) +

  geom_point(
    aes(
      colour = change_direction
    ),
    size = 2.6,
    alpha = 0.90
  ) +

  stat_summary(
    aes(group = timepoint),
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.8,
    fill = "white",
    colour = "black"
  ) +

  scale_colour_manual(
    values = c(
      "Increase" = "#0072B2",
      "Decrease" = "#D55E00",
      "No meaningful change" =
        "grey50"
    ),
    name = "Post-exposure change"
  ) +

  labs(
    title = paste0(
      dataset_label,
      ": DNA concentration"
    ),

    subtitle = paste0(
      "Paired participants: n = ",
      length(
        paired_concentration_ids
      )
    ),

    x = NULL,

    y =
      "DNA concentration (ng/µL)",

    caption = paste(
      "Lines connect measurements from",
      "the same participant.",
      "White diamonds indicate medians."
    )
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      ),

    panel.grid.minor =
      element_blank(),

    legend.position =
      "right"
  )


# ============================================================
# 7. Plasma-normalised DNA yield
# ============================================================

paired_yield_ids <-
  plot_data %>%

  filter(
    !is.na(
      dna_yield_ng_per_ml_plasma
    )
  ) %>%

  distinct(
    participant_id,
    timepoint
  ) %>%

  count(
    participant_id,
    name = "n_timepoints"
  ) %>%

  filter(
    n_timepoints == 2
  ) %>%

  pull(participant_id)


paired_yield_data <-
  plot_data %>%

  filter(
    participant_id %in%
      paired_yield_ids
  )


yield_change <-
  paired_yield_data %>%

  select(
    participant_id,
    timepoint,
    dna_yield_ng_per_ml_plasma
  ) %>%

  pivot_wider(
    names_from = timepoint,
    values_from =
      dna_yield_ng_per_ml_plasma
  ) %>%

  mutate(

    yield_change_ng_per_ml =
      `Post-exposure` -
      `Pre-exposure`,

    change_direction =
      case_when(

        yield_change_ng_per_ml >
          yield_threshold ~
          "Increase",

        yield_change_ng_per_ml <
          -yield_threshold ~
          "Decrease",

        TRUE ~
          "No meaningful change"
      )
  )


paired_yield_data <-
  paired_yield_data %>%

  left_join(
    yield_change %>%
      select(
        participant_id,
        yield_change_ng_per_ml,
        change_direction
      ),
    by = "participant_id"
  )


yield_summary <-
  yield_change %>%

  count(change_direction) %>%

  mutate(
    percent =
      round(
        100 * n / sum(n),
        1
      )
  )


cat(
  "\nPaired yield participants: ",
  length(paired_yield_ids),
  "\n",
  sep = ""
)

print(yield_summary)


p2 <- ggplot(
  paired_yield_data,
  aes(
    x = timepoint,
    y = dna_yield_ng_per_ml_plasma,
    group = participant_id
  )
) +

  geom_boxplot(
    aes(group = timepoint),
    width = 0.45,
    outlier.shape = NA,
    alpha = 0.15,
    colour = "black"
  ) +

  geom_line(
    aes(
      colour = change_direction
    ),
    alpha = 0.60,
    linewidth = 0.7
  ) +

  geom_point(
    aes(
      colour = change_direction
    ),
    size = 2.6,
    alpha = 0.90
  ) +

  stat_summary(
    aes(group = timepoint),
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.8,
    fill = "white",
    colour = "black"
  ) +

  scale_colour_manual(
    values = c(
      "Increase" = "#0072B2",
      "Decrease" = "#D55E00",
      "No meaningful change" =
        "grey50"
    ),
    name = "Post-exposure change"
  ) +

  labs(
    title = paste0(
      dataset_label,
      ": plasma-normalised DNA yield"
    ),

    subtitle = paste0(
      "Paired participants: n = ",
      length(
        paired_yield_ids
      )
    ),

    x = NULL,

    y =
      "DNA yield (ng per mL starting plasma)",

    caption = paste(
      "Total DNA yield was divided by",
      "starting plasma volume.",
      "White diamonds indicate medians."
    )
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      ),

    panel.grid.minor =
      element_blank(),

    legend.position =
      "right"
  )


print(p1)
print(p2)


# ============================================================
# 8. Save summary tables
# ============================================================

write_tsv(
  concentration_change,
  file.path(
    output_dir,
    paste0(
      output_prefix,
      "_concentration_pre_post_change.tsv"
    )
  )
)

write_tsv(
  yield_change,
  file.path(
    output_dir,
    paste0(
      output_prefix,
      "_yield_pre_post_change.tsv"
    )
  )
)


# ============================================================
# 9. Save figures
# ============================================================

combined_pdf <- file.path(
  output_dir,
  paste0(
    output_prefix,
    "_pre_post_dna_isolation_plots.pdf"
  )
)


pdf(
  combined_pdf,
  width = plot_width,
  height = plot_height,
  onefile = TRUE
)

print(p1)
print(p2)

dev.off()


ggsave(
  filename = file.path(
    output_dir,
    paste0(
      output_prefix,
      "_pre_post_dna_concentration.pdf"
    )
  ),
  plot = p1,
  width = plot_width,
  height = plot_height,
  device = cairo_pdf
)


ggsave(
  filename = file.path(
    output_dir,
    paste0(
      output_prefix,
      "_pre_post_normalised_dna_yield.pdf"
    )
  ),
  plot = p2,
  width = plot_width,
  height = plot_height,
  device = cairo_pdf
)


cat("\n============================================\n")
cat("Plotting completed successfully\n")
cat("============================================\n")
cat("Results directory:\n")
cat(output_dir, "\n")
cat("============================================\n")