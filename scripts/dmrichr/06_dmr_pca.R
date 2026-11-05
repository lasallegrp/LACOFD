#!/usr/bin/env Rscript
#=====================================================================================
#  06_dmr_pca.R
#  PCA of samples on the smoothed methylation of their DMRs, one plot per comparison.
#
#  Reads from each comparison folder made by 02_setup_comparisons.py:
#    <comparison>/DMRs/DMR_individual_smoothed_methylation.txt   (written by DMRichR)
#    <comparison>/sample_info.xlsx                                (Name + test covariate)
#  Sample columns are found by matching sample_info$Name, so the number of region
#  annotation columns does not matter.
#
#  Writes to <comparison>/DMRs/PCA/ (or --out_dir/<run>/<comparison>/):
#    <comparison>_DMR_PCA.pdf     PC1 vs PC2 (+ PC2 vs PC3 if --pcs 3), coloured by group
#    <comparison>_DMR_PCA.xlsx    sample scores, variance explained, and PC-group tests
#
#  Examples (from the repository root):
#    Rscript scripts/dmrichr/06_dmr_pca.R --run DMR_AZ_04
#    Rscript scripts/dmrichr/06_dmr_pca.R --run DMR_AZ_04,DMR_AZ_05 --shape Sex --label
#    Rscript scripts/dmrichr/06_dmr_pca.R --comparison_dir results/dmrichr/DMR_AZ_04/recruit_post_vs_pre
#    Rscript scripts/dmrichr/06_dmr_pca.R --run DMR_AZ_04 --regions background   # all tested regions
#=====================================================================================
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(openxlsx)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--run", type = "character", default = NA,
              help = "Comma list of run IDs, e.g. DMR_AZ_04,DMR_AZ_05"),
  make_option("--run_root", type = "character", default = "results/dmrichr",
              help = "Folder holding the runs [default %default]"),
  make_option("--comparisons", type = "character", default = NA,
              help = "Comma list of comparisons to plot [default: all finished in each run]"),
  make_option("--comparison_dir", type = "character", default = NA,
              help = "Instead of --run: one or more comparison folders, comma-separated"),
  make_option("--regions", type = "character", default = "dmr",
              help = "dmr = DMRs only; background = all tested regions [default %default]"),
  make_option("--colour", type = "character", default = NA,
              help = "sample_info column to colour by [default: the test covariate, i.e. 2nd column]"),
  make_option("--shape", type = "character", default = NA,
              help = "Optional sample_info column for point shape, e.g. Sex"),
  make_option("--label", action = "store_true", default = FALSE,
              help = "Label points with sample names"),
  make_option("--pcs", type = "integer", default = 2L,
              help = "2 = PC1 vs PC2; 3 = also PC2 vs PC3 [default %default]"),
  make_option("--scale", type = "logical", default = TRUE,
              help = "Scale regions to unit variance before PCA [default %default]"),
  make_option("--show_p", action = "store_true", default = FALSE,
              help = "Show the PC group-test p-values in the subtitle (always written to the .xlsx)"),
  make_option("--min_regions", type = "integer", default = 3L,
              help = "Skip comparisons with fewer regions than this [default %default]"),
  make_option("--out_dir", type = "character", default = NA,
              help = "Write here as <out_dir>/<run>/<comparison>/ [default: <comparison>/DMRs/PCA/]"),
  make_option("--width", type = "double", default = 7),
  make_option("--height", type = "double", default = 6)
)))

if (is.na(opt$run) && is.na(opt$comparison_dir)) stop("Give --run or --comparison_dir")
stopifnot(opt$regions %in% c("dmr", "background"), opt$pcs %in% 2:3)
split_opt <- function(x) if (is.na(x) || !nzchar(x)) NULL else trimws(strsplit(x, ",")[[1]])
meth_file <- if (opt$regions == "dmr") "DMRs/DMR_individual_smoothed_methylation.txt" else
  "DMRs/background_region_individual_smoothed_methylation.txt"

# ----------------------------- which comparison folders -----------------------------
dirs <- if (!is.na(opt$comparison_dir)) split_opt(opt$comparison_dir) else {
  keep <- split_opt(opt$comparisons)
  unlist(lapply(split_opt(opt$run), function(r) {
    rd <- file.path(opt$run_root, r)
    if (!dir.exists(rd)) { message("Run folder not found, skipped: ", rd); return(NULL) }
    d <- list.dirs(rd, recursive = FALSE)
    d <- d[file.exists(file.path(d, "sample_info.xlsx"))]          # comparison folders only
    if (!is.null(keep)) d <- d[basename(d) %in% keep]
    d
  }))
}
if (!length(dirs)) stop("No comparison folders found")

clean_level <- function(x) sub("^[0-9]+_", "", as.character(x))   # "1_Pre" -> "Pre"
pal <- c("#2166AC", "#B2182B", "#1B7837", "#762A83", "#E08214", "#4D4D4D")

# ----------------------------- one comparison ---------------------------------------
pca_one <- function(cdir) {
  comp <- basename(cdir); run <- basename(dirname(cdir))
  mf <- file.path(cdir, meth_file)
  if (!file.exists(mf)) { message(comp, ": no ", meth_file, " (DMRichR not finished?) -> skipped"); return(NULL) }
  si <- read.xlsx(file.path(cdir, "sample_info.xlsx"))
  if (names(si)[1] != "Name") stop(comp, ": first column of sample_info.xlsx must be Name")
  colour_by <- if (is.na(opt$colour)) names(si)[2] else opt$colour
  for (v in na.omit(c(colour_by, opt$shape)))
    if (!v %in% names(si)) stop(comp, ": column '", v, "' not in sample_info.xlsx (has: ",
                                paste(names(si), collapse = ", "), ")")

  x <- read.delim(mf, header = TRUE, check.names = FALSE)
  samp <- intersect(si$Name, names(x))
  miss <- setdiff(si$Name, names(x))
  if (length(miss)) message(comp, ": samples in sample_info but not in the methylation file: ",
                            paste(miss, collapse = ", "))
  if (length(samp) < 3) { message(comp, ": fewer than 3 samples -> skipped"); return(NULL) }
  m <- as.matrix(x[, samp, drop = FALSE])
  storage.mode(m) <- "double"
  n_regions_all <- nrow(m)
  m <- m[stats::complete.cases(m), , drop = FALSE]                        # drop regions with missing values
  if (opt$scale) m <- m[apply(m, 1, stats::var) > 0, , drop = FALSE]      # constant regions cannot be scaled
  if (nrow(m) < opt$min_regions) {
    message(comp, ": only ", nrow(m), " usable regions -> skipped"); return(NULL)
  }

  pca <- prcomp(t(m), center = TRUE, scale. = opt$scale)
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  npc <- min(ncol(pca$x), 5)
  sc <- data.frame(Name = rownames(pca$x), pca$x[, seq_len(npc), drop = FALSE], check.names = FALSE)
  sc <- merge(si, sc, by = "Name")
  sc$.group <- factor(clean_level(sc[[colour_by]]), levels = unique(clean_level(sort(unique(sc[[colour_by]])))))
  if (!is.na(opt$shape)) sc$.shape <- factor(sc[[opt$shape]])

  # does each PC separate the groups? (Wilcoxon for 2 groups, Kruskal-Wallis otherwise)
  tests <- do.call(rbind, lapply(seq_len(npc), function(k) {
    pc <- sc[[paste0("PC", k)]]
    p <- tryCatch(if (nlevels(sc$.group) == 2) wilcox.test(pc ~ sc$.group)$p.value
                  else kruskal.test(pc ~ sc$.group)$p.value, error = function(e) NA)
    data.frame(PC = paste0("PC", k), variance_pct = round(ve[k], 2), group_test_p = signif(p, 3))
  }))

  plot_pair <- function(a, b) {
    p <- ggplot(sc, aes(.data[[paste0("PC", a)]], .data[[paste0("PC", b)]], colour = .group)) +
      geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85")
    if (all(table(sc$.group) >= 4))
      p <- p + stat_ellipse(aes(fill = .group), geom = "polygon", type = "norm", alpha = 0.12,
                            colour = NA, show.legend = FALSE)
    p <- p + (if (is.na(opt$shape)) geom_point(size = 2.6) else geom_point(aes(shape = .shape), size = 2.6)) +
      scale_colour_manual(values = pal, name = colour_by) +
      scale_fill_manual(values = pal, guide = "none") +
      labs(x = sprintf("PC%d (%.1f%%)", a, ve[a]), y = sprintf("PC%d (%.1f%%)", b, ve[b]),
           shape = opt$shape,
           title = paste0(run, " / ", comp),
           subtitle = if (opt$show_p)
             sprintf("%s regions, %d samples; group difference on PC%d p = %s, PC%d p = %s",
                     format(nrow(m), big.mark = ","), nrow(sc), a,
                     format(tests$group_test_p[a]), b, format(tests$group_test_p[b]))
           else sprintf("%s regions, %d samples", format(nrow(m), big.mark = ","), nrow(sc)),
           caption = paste0(if (opt$regions == "dmr") "DMRs" else "All tested regions",
                            ", smoothed methylation", if (opt$scale) ", scaled" else "",
                            if (n_regions_all > nrow(m)) sprintf("; %d regions dropped (missing or constant)",
                                                                 n_regions_all - nrow(m)) else "")) +
      theme_classic(12) +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            axis.text = element_text(colour = "black"),
            axis.line = element_line(colour = "black", linewidth = 0.2),
            plot.subtitle = element_text(size = 9), plot.caption = element_text(size = 8, colour = "grey40"))
    if (opt$label) {
      p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_text_repel(aes(label = Name), size = 2.6, show.legend = FALSE, max.overlaps = Inf, seed = 1)
      else geom_text(aes(label = Name), size = 2.6, vjust = -0.8, show.legend = FALSE)
    }
    p
  }

  od <- if (is.na(opt$out_dir)) file.path(cdir, "DMRs", "PCA") else file.path(opt$out_dir, run, comp)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  tag <- if (opt$regions == "dmr") "DMR" else "background"
  pdf(file.path(od, paste0(comp, "_", tag, "_PCA.pdf")), width = opt$width, height = opt$height)
  print(plot_pair(1, 2))
  if (opt$pcs == 3 && npc >= 3) print(plot_pair(2, 3))
  invisible(dev.off())
  sc$.group <- NULL; sc$.shape <- NULL
  write.xlsx(list(Scores = sc, PC_tests = tests), file.path(od, paste0(comp, "_", tag, "_PCA.xlsx")),
             overwrite = TRUE)
  message(sprintf("%s / %s: %d regions, %d samples; PC1 %.1f%% (group p = %s) -> %s",
                  run, comp, nrow(m), nrow(sc), ve[1], format(tests$group_test_p[1]), od))
  data.frame(run = run, comparison = comp, regions = nrow(m), samples = nrow(sc),
             PC1_pct = round(ve[1], 1), PC1_group_p = tests$group_test_p[1],
             PC2_pct = round(ve[2], 1), PC2_group_p = tests$group_test_p[2])
}

res <- do.call(rbind, lapply(dirs, pca_one))
if (!is.null(res) && nrow(res)) {
  print(res, row.names = FALSE)
  if (!is.na(opt$out_dir)) {
    dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
    write.table(res, file.path(opt$out_dir, "PCA_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}