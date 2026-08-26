library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. File paths
# ============================================================

input_file <- paste0(
  "/quobyte/lasallegrp/projects/LACOFD/",
  "data/processed/wet_lab/dna_isolation/la_county/",
  "la_county_dna_isolation_analysis_ready.tsv"
)

output_dir <- paste0(
  "/quobyte/lasallegrp/projects/LACOFD/",
  "results/wet_lab_qc/dna_isolation/la_county"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. Read data
# ============================================================

la_extractions <- read_tsv(
  input_file,
  show_col_types = FALSE
)

cat("Input data dimensions:", dim(la_extractions), "\n")

# ============================================================
# 3. Prepare plotting data
# ============================================================

plot_data <- la_extractions %>%
  filter(
    !is.na(participant_id),
    !is.na(timepoint)
  ) %>%
  mutate(
    timepoint = case_when(
      timepoint == "pre" ~ "Pre-exposure",
      timepoint == "post" ~ "Post-exposure",
      TRUE ~ as.character(timepoint)
    ),
    timepoint = factor(
      timepoint,
      levels = c("Pre-exposure", "Post-exposure")
    )
  )

# ============================================================
# 4. Check duplicates
# ============================================================

duplicate_check <- plot_data %>%
  count(participant_id, timepoint) %>%
  filter(n > 1)

if (nrow(duplicate_check) > 0) {
  print(duplicate_check)
  stop("Some participants have more than one record per timepoint.")
}

# ============================================================
# 5. Plot 1: DNA concentration
# ============================================================

# Keep participants with concentration measured at both timepoints
paired_concentration_ids <- plot_data %>%
  filter(!is.na(dna_concentration_ng_ul)) %>%
  distinct(participant_id, timepoint) %>%
  count(participant_id, name = "n_timepoints") %>%
  filter(n_timepoints == 2) %>%
  pull(participant_id)

paired_concentration_data <- plot_data %>%
  filter(participant_id %in% paired_concentration_ids)

# Classify increase / decrease / no change
concentration_change <- paired_concentration_data %>%
  select(participant_id, timepoint, dna_concentration_ng_ul) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = dna_concentration_ng_ul
  ) %>%
  mutate(
    concentration_change_ng_ul =
      `Post-exposure` - `Pre-exposure`,
    change_direction = case_when(
      concentration_change_ng_ul > 0.01 ~ "Increase",
      concentration_change_ng_ul < -0.01 ~ "Decrease",
      TRUE ~ "No meaningful change"
    )
  )

paired_concentration_data <- paired_concentration_data %>%
  left_join(
    concentration_change %>%
      select(
        participant_id,
        concentration_change_ng_ul,
        change_direction
      ),
    by = "participant_id"
  )

cat(
  "Number of paired participants for concentration:",
  length(paired_concentration_ids), "\n"
)

print(
  concentration_change %>%
    count(change_direction) %>%
    mutate(percent = round(100 * n / sum(n), 1))
)

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
    aes(colour = change_direction),
    alpha = 0.60,
    linewidth = 0.7
  ) +
  geom_point(
    aes(colour = change_direction),
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
      "No meaningful change" = "grey50"
    ),
    name = "Post-exposure change"
  ) +
  labs(
    title = "LA County DNA concentration before and after exposure",
    subtitle = paste0(
      "Paired participants: n = ",
      length(paired_concentration_ids)
    ),
    x = NULL,
    y = "DNA concentration (ng/µL)",
    caption = paste(
      "Lines are from the same participant.",
      "Blue = increase, orange = decrease, grey = no meaningful change.",
      "White diamonds indicate medians."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ============================================================
# 6. Plot 2: Plasma-normalised DNA yield
# ============================================================

# Keep participants with normalised yield measured at both timepoints
paired_yield_ids <- plot_data %>%
  filter(!is.na(dna_yield_ng_per_ml_plasma)) %>%
  distinct(participant_id, timepoint) %>%
  count(participant_id, name = "n_timepoints") %>%
  filter(n_timepoints == 2) %>%
  pull(participant_id)

paired_yield_data <- plot_data %>%
  filter(participant_id %in% paired_yield_ids)

# Classify increase / decrease / no change
yield_change <- paired_yield_data %>%
  select(participant_id, timepoint, dna_yield_ng_per_ml_plasma) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = dna_yield_ng_per_ml_plasma
  ) %>%
  mutate(
    yield_change_ng_per_ml =
      `Post-exposure` - `Pre-exposure`,
    change_direction = case_when(
      yield_change_ng_per_ml > 0.1 ~ "Increase",
      yield_change_ng_per_ml < -0.1 ~ "Decrease",
      TRUE ~ "No meaningful change"
    )
  )

paired_yield_data <- paired_yield_data %>%
  left_join(
    yield_change %>%
      select(
        participant_id,
        yield_change_ng_per_ml,
        change_direction
      ),
    by = "participant_id"
  )

cat(
  "Number of paired participants for normalised yield:",
  length(paired_yield_ids), "\n"
)

print(
  yield_change %>%
    count(change_direction) %>%
    mutate(percent = round(100 * n / sum(n), 1))
)

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
    aes(colour = change_direction),
    alpha = 0.60,
    linewidth = 0.7
  ) +
  geom_point(
    aes(colour = change_direction),
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
      "No meaningful change" = "grey50"
    ),
    name = "Post-exposure change"
  ) +
  labs(
    title = "LA County plasma-normalised DNA yield",
    subtitle = paste0(
      "Paired participants: n = ",
      length(paired_yield_ids)
    ),
    x = NULL,
    y = "DNA yield (ng per mL starting plasma)",
    caption = paste(
      "Total DNA yield was divided by starting plasma volume.",
      "Blue = increase, orange = decrease, grey = no meaningful change.",
      "White diamonds indicate medians."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Show plots in session
print(p1)
print(p2)

# ============================================================
# 7. Save plots
# ============================================================

# Multi-page PDF
pdf_file <- file.path(
  output_dir,
  "la_county_pre_post_dna_isolation_plots_colored.pdf"
)

pdf(
  pdf_file,
  width = 14,
  height = 7,
  onefile = TRUE
)

print(p1)
print(p2)

dev.off()

# Individual PDFs
ggsave(
  filename = file.path(
    output_dir,
    "la_county_pre_post_dna_concentration_colored.pdf"
  ),
  plot = p1,
  width = 14,
  height = 7,
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "la_county_pre_post_normalised_dna_yield_colored.pdf"
  ),
  plot = p2,
  width = 14,
  height = 7,
  device = cairo_pdf
)

cat("Plots saved to:\n", output_dir, "\n")