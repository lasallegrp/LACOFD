#!/usr/bin/env Rscript
#=====================================================================================
#  03_dmrichr_overlap.R
#  Overlap of DMRichR DMR sets, by genomic coordinates and/or gene symbol.
#
#  Two input layouts (--layout):
#
#  outcome_timepoint  (default; GROWELL layout)
#     <results_dir>/<outcome>/<timepoint>/[<subdir>/]DMRs/DMRs_annotated.xlsx
#     results_dir defaults to <project_dir>/results/dmrichr_output
#     e.g. results/apo/Baseline, results/apo/36-38wk, results/apo/Postpartum
#     Sets are discovered automatically; folders without a DMR file are skipped.
#     set_id = <outcome>_<timepoint> (non-alphanumerics -> "_", e.g. apo_36_38wk)
#     Comparisons:
#       within_outcome_across_timepoints   e.g. apo: Baseline / 36-38wk / Postpartum
#       within_timepoint_across_outcomes   e.g. Baseline: apo / egwg / ppwr / ...
#
#  config_sex  (original layout)
#     <project_dir>/config/dmrichr_sets.tsv  (set_prefix, contrast_dir, trait)
#     <dmrichr_root>/<contrast_dir>/cytosine_reports/<sex>/<cutoff>
#     Comparisons: within_sex_across_conditions, within_conditions_across_sex
#
#  custom    any grouping given in --groups_file (group_id <TAB> set_id), either layout
#  (--manifest with explicit paths: set_id, contrast, sex, path is still accepted.)
#
#  Every multi-set overlap pattern is also written to its own folder:
#    <out_dir>/<comparison>/<group>/<coord|gene>/overlaps/<pattern>/DMRs_annotated_<pattern>.xlsx
#
#  GWAS integration (optional; inputs = GWAS Catalog association export):
#    --gene_set <export.tsv>  Fisher test per set: GWAS genes among DMR genes, universe =
#                             genes of that set's DMRichR background; association-weighted
#                             Wilcoxon; outcomes + GWAS UpSet per group; evidence plots.
#    --snp_set  <export.tsv>  regioneR: DMRs within --snp_window bp of GWAS SNPs vs
#                             resampled background (candidate) regions.
#    Outputs in <out_dir>/gwas/ and in each group's gene/ folder.
#
#  Example (GROWELL):
#   Rscript 03_dmrichr_overlap.R --project_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/<project> \
#       --timepoints Baseline,36-38wk,Postpartum --overlap_mode both
#
#  With the regioneR permutation test (egwg vs ppwr at each timepoint):
#   Rscript 03_dmrichr_overlap.R --project_dir <project> --outcomes egwg,ppwr \
#       --timepoints Baseline,36-38wk,Postpartum --comparisons within_timepoint_across_outcomes \
#       --overlap_mode coord --perm_test                       # 100 x resampleRegions
#   add --perm_randomize randomize --perm_genome hg38   for genome-wide randomizeRegions
#=====================================================================================
suppressPackageStartupMessages({ library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--project_dir", type = "character", default = NA,
              help = "Project root containing results/ (and config/ for layout config_sex)"),
  make_option("--layout", type = "character", default = "outcome_timepoint",
              help = "outcome_timepoint | config_sex [default %default]"),
  # ---- outcome_timepoint layout ----
  make_option("--results_dir", type = "character", default = NA,
              help = "Folder holding <outcome>/<timepoint>/ [default <project_dir>/results/dmrichr_output]"),
  make_option("--outcomes", type = "character", default = NA,
              help = "Comma list of outcome folders, in display order [default: all discovered]"),
  make_option("--timepoints", type = "character", default = NA,
              help = "Comma list of timepoint folders, in display order, e.g. Baseline,36-38wk,Postpartum [default: all]"),
  make_option("--subdir", type = "character", default = "",
              help = "Extra path between <timepoint> and DMRs/, e.g. cutoff_005 [default none]"),
  make_option("--exclude_dirs", type = "character", default = "logs,metadata,DMRichR_overlap",
              help = "Folders in results_dir that are never treated as outcomes [default %default]"),
  # ---- config_sex layout ----
  make_option("--config_dir", type = "character", default = NA,
              help = "Config folder [default <project_dir>/config]"),
  make_option("--sets_file", type = "character", default = "dmrichr_sets.tsv",
              help = "Contrast table inside config_dir: set_prefix, contrast_dir, trait [default %default]"),
  make_option("--dmrichr_root", type = "character", default = NA,
              help = "DMRichR output root [default <project_dir>/results/dmrichr_output]"),
  make_option("--cutoff", type = "character", default = "cutoff_005",
              help = "Cutoff folder used for EVERY set (config_sex only) [default %default]"),
  make_option("--sexes", type = "character", default = "males,females,sex_combined",
              help = "Sex folders to include (config_sex only) [default %default]"),
  # ---- common ----
  make_option("--manifest", type = "character", default = NA,
              help = "Alternative to discovery: TSV with set_id, contrast, sex, path"),
  make_option("--out_dir", type = "character", default = NA,
              help = "[default <project_dir>/results/DMRichR_overlap (outcome_timepoint) or .../DMRichR_overlap_sex_stratified/<cutoff>]"),
  make_option("--comparisons", type = "character", default = NA,
              help = paste("Comma list. outcome_timepoint: within_outcome_across_timepoints,within_timepoint_across_outcomes;",
                           "config_sex: within_sex_across_conditions,within_conditions_across_sex [default: both for the layout]")),
  make_option("--groups_file", type = "character", default = NA,
              help = "Optional TSV group_id,set_id for custom groupings"),
  make_option("--per_overlap_folders", type = "logical", default = TRUE,
              help = "Write each overlap pattern's annotated DMRs to its own folder [default %default]"),
  make_option("--overlap_mode", type = "character", default = "both", help = "coord | gene | both"),
  make_option("--max_gap", type = "integer", default = 0L,
              help = "bp; DMRs this close count as overlapping (0 = must overlap)"),
  make_option("--min_overlap_bp", type = "integer", default = 1L,
              help = "bp of overlap required; if >1 overrides --max_gap"),
  make_option("--require_direction", action = "store_true", default = FALSE,
              help = "Count multi-set overlaps only when methylation direction agrees"),
  make_option("--exclude_sex_chr", action = "store_true", default = FALSE,
              help = "Drop chrX/chrY from DMRs and backgrounds"),
  make_option("--dmr_file", type = "character", default = "DMRs/DMRs_annotated.xlsx"),
  make_option("--background_file", type = "character", default = "DMRs/background_annotated.xlsx",
              help = "Tested candidate regions (coords; optional difference/pval). 'none' to skip"),
  make_option("--gene_universe", type = "character", default = "background",
              help = "Gene universe for gene-mode Fisher: 'background' (genes of all tested regions), a file, or 'none'"),
  make_option("--effect_col", type = "character", default = "auto",
              help = "Column giving direction: auto = difference > stat > beta"),
  make_option("--gene_col", type = "character", default = "geneSymbol"),
  make_option("--gene_sep", type = "character", default = NA,
              help = "Separator if a DMR lists several genes (e.g. ';')"),
  make_option("--max_pval", type = "double", default = NA, help = "Optional extra DMR pval filter"),
  make_option("--max_qval", type = "double", default = NA, help = "Optional DMR qval filter"),
  make_option("--sex_labels", type = "character", default = "males,females,sex_combined",
              help = "Sex values for male, female, combined (config_sex only)"),
  make_option("--upset_genes", type = "integer", default = 0L,
              help = "Gene names per bar on the coloured UpSet (*_upset_<mode>_colour.pdf; 0 = no names) [default %default]"),
  make_option("--upset_genes_scope", type = "character", default = "multi",
              help = "multi = label only overlap bars, all = label every bar [default %default]"),
  make_option("--upset_highlight_genes", type = "character", default = NA,
              help = "Optional gene list (one per line, or table with geneSymbol column); only these are labeled"),
  make_option("--scatter_genes", type = "integer", default = 10L,
              help = "Gene names to label on the cross-sex scatter (config_sex only; 0 = off) [default %default]"),
  make_option("--subthreshold_p", type = "double", default = 0.2,
              help = "Background p below which a non-significant other-sex effect counts as evidence"),
  # ---- regioneR permutation test (coord mode) ----
  make_option("--perm_test", action = "store_true", default = FALSE,
              help = "Run regioneR permutation test for every pair (both directions) in coord mode"),
  make_option("--perm_n", type = "integer", default = 100L,
              help = "Permutations per test; min p = 1/(n+1) [default %default]"),
  make_option("--perm_randomize", type = "character", default = "resample",
              help = paste("resample = regioneR::resampleRegions from A's tested candidate regions;",
                           "randomize = regioneR::randomizeRegions anywhere in --perm_genome [default %default]")),
  make_option("--perm_genome", type = "character", default = "hg38",
              help = "randomize only: genome name for regioneR (needs BSgenome) or a chrom.sizes file [default %default]"),
  make_option("--perm_mask", type = "character", default = NA,
              help = "randomize only: BED/TSV of regions to exclude (e.g. gaps, blacklist) [default none]"),
  make_option("--perm_alternative", type = "character", default = "greater",
              help = "greater (enrichment, matches the Fisher test) | less | auto [default %default]"),
  make_option("--perm_per_chr", type = "logical", default = TRUE,
              help = "Resample within chromosome (keeps A's chromosome distribution) [default %default]"),
  make_option("--perm_local_z", type = "logical", default = TRUE,
              help = "Compute regioneR local z-score (position specificity) [default %default]"),
  make_option("--perm_lz_window", type = "integer", default = NA,
              help = "Local-z window, bp [default 5 x mean DMR width, min 1000]"),
  make_option("--perm_lz_step", type = "integer", default = NA,
              help = "Local-z step, bp [default window/20]"),
  make_option("--perm_cores", type = "integer", default = 1L,
              help = "Cores for regioneR (parallel::mclapply) [default %default]"),
  make_option("--perm_seed", type = "integer", default = 1L, help = "Seed [default %default]"),
  # ---- GWAS gene set / SNP set ----
  make_option("--gene_set", type = "character", default = NA,
              help = "GWAS Catalog association export (or gene list with geneSymbol[, n_associations]) for the gene-set test"),
  make_option("--snp_set", type = "character", default = NA,
              help = "GWAS Catalog association export (or table chr, pos) for the SNP-proximity regioneR test"),
  make_option("--gwas_label", type = "character", default = "GWAS",
              help = "Short name for the GWAS set in plots/tables [default %default]"),
  make_option("--gwas_max_p", type = "double", default = NA,
              help = "Keep only associations with p <= this [default: all rows in the export]"),
  make_option("--snp_window", type = "integer", default = 50000L,
              help = "bp around each GWAS SNP for the SNP-proximity test [default %default]"),
  make_option("--harmonize_symbols", type = "logical", default = TRUE,
              help = "Map GWAS and DMRichR gene symbols to current symbols via org.Hs.eg.db [default %default]"),
  make_option("--gwas_expand_readthrough", type = "logical", default = TRUE,
              help = "Also credit read-through GWAS genes (A-B-C) to each component gene, when every part is a current symbol [default %default]"),
  make_option("--gwas_upset_genes", type = "integer", default = 10L,
              help = "Gene names per bar on the outcomes + GWAS UpSet, ranked by n associations [default %default]"),
  make_option("--gwas_plot_max", type = "integer", default = 60L,
              help = "Max genes on the GWAS evidence plots [default %default]"),
  make_option("--helper", type = "character", default = NA, help = "Path to helper_overlap.R")
)))

stopifnot(opt$layout %in% c("outcome_timepoint", "config_sex"))
is_ot <- opt$layout == "outcome_timepoint"
if (is.na(opt$project_dir) && is.na(opt$manifest) && !(is_ot && !is.na(opt$results_dir)))
  stop("Give --project_dir (or --results_dir for outcome_timepoint), or --manifest")
if (is_ot && is.na(opt$results_dir) && !is.na(opt$project_dir))
  opt$results_dir <- file.path(opt$project_dir, "results", "dmrichr_output")
if (is.na(opt$out_dir)) {
  if (is_ot) {
    opt$out_dir <- file.path(opt$results_dir, "DMRichR_overlap")
  } else {
    if (is.na(opt$project_dir)) stop("--out_dir is required without --project_dir")
    opt$out_dir <- file.path(opt$project_dir, "results", "DMRichR_overlap_sex_stratified", opt$cutoff)
  }
}
.a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
self <- if (length(.a)) normalizePath(sub("^--file=", "", .a[1])) else NA_character_
if (is.na(opt$helper)) opt$helper <- file.path(dirname(self), "helper_overlap.R")
source(opt$helper)
stopifnot(opt$overlap_mode %in% c("coord", "gene", "both"))
do_coord <- opt$overlap_mode %in% c("coord", "both")
do_gene  <- opt$overlap_mode %in% c("gene", "both")
sex_lab  <- setNames(strsplit(opt$sex_labels, ",")[[1]], c("M", "F", "C"))
if (opt$perm_test) {
  if (!do_coord) stop("--perm_test needs --overlap_mode coord or both")
  if (!requireNamespace("regioneR", quietly = TRUE))
    stop("--perm_test needs regioneR: BiocManager::install('regioneR') ",
         "(pixi/conda: bioconductor-regioner)")
  stopifnot(opt$perm_alternative %in% c("greater", "less", "auto"), opt$perm_n >= 10,
            opt$perm_randomize %in% c("resample", "randomize"))
  perm_genome <- NULL; perm_mask <- NA      # NA = no mask (NULL would make regioneR look for a masked BSgenome)
  if (opt$perm_randomize == "randomize") {
    perm_genome <- read_perm_genome(opt$perm_genome)
    if (!is.na(opt$perm_mask)) {
      mk <- read_table_any(opt$perm_mask)
      if (!"start" %in% names(mk)) names(mk)[1:3] <- c("chr", "start", "end")
      mc <- detect_col(mk, c("seqnames", "chr", "chrom"), "chromosome")
      perm_mask <- GRanges(mk[[mc]], IRanges(as.integer(mk$start), as.integer(mk$end)))
    }
    message("regioneR: randomizeRegions over ", if (is.character(perm_genome)) perm_genome else opt$perm_genome,
            if (is(perm_mask, "GRanges")) paste0(", mask = ", length(perm_mask), " regions") else ", no mask",
            " (genome-wide null; usually anti-conservative for DMRs)")
  }
  RNGkind("L'Ecuyer-CMRG")                    # reproducible streams under mclapply
  options(mc.cores = opt$perm_cores)
  if (tolower(opt$background_file) == "none")
    stop("--perm_test needs background files (the universe); remove --background_file none")
}
if (!is.na(opt$snp_set)) {
  if (!requireNamespace("regioneR", quietly = TRUE))
    stop("--snp_set needs regioneR (pixi/conda: bioconductor-regioner)")
  if (tolower(opt$background_file) == "none") stop("--snp_set needs background files (the universe)")
  RNGkind("L'Ecuyer-CMRG"); options(mc.cores = opt$perm_cores)
}
if (!is.na(opt$gene_set) && tolower(opt$background_file) == "none")
  stop("--gene_set needs background files (the gene universe)")
split_opt <- function(x) if (is.na(x) || !nzchar(x)) NULL else trimws(strsplit(x, ",")[[1]])
clean_id  <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

# Internally the two grouping dimensions keep the original column names so the helper
# functions work unchanged:  man$contrast = outcome,  man$sex = timepoint.
dim_names <- if (is_ot) c(contrast = "outcome", sex = "timepoint") else c(contrast = "contrast", sex = "sex")

# ----------------------------- discover sets (outcome_timepoint) ----------------------
discover_outcome_timepoint <- function(results_dir, dmr_file, subdir, outcomes, timepoints, exclude) {
  if (!dir.exists(results_dir)) stop("results_dir not found: ", results_dir)
  found <- setdiff(list.dirs(results_dir, recursive = FALSE, full.names = FALSE), exclude)
  if (!is.null(outcomes)) {
    miss <- setdiff(outcomes, found)
    if (length(miss)) message("Outcome folders not found, skipped: ", paste(miss, collapse = ", "))
    found <- outcomes[outcomes %in% found]                       # keep user order
  }
  rows <- list()
  for (o in found) {
    tp <- list.dirs(file.path(results_dir, o), recursive = FALSE, full.names = FALSE)
    if (!is.null(timepoints)) {
      extra <- setdiff(tp, timepoints)
      if (length(extra)) message(o, ": folders not in --timepoints, skipped: ", paste(extra, collapse = ", "))
      miss <- setdiff(timepoints, tp)
      if (length(miss)) message(o, ": no folder for timepoint(s): ", paste(miss, collapse = ", "))
      tp <- timepoints[timepoints %in% tp]
    }
    for (t in tp) {
      p <- if (nzchar(subdir)) file.path(results_dir, o, t, subdir) else file.path(results_dir, o, t)
      if (file.exists(file.path(p, dmr_file))) {
        rows[[length(rows) + 1]] <- data.frame(set_id = clean_id(paste(o, t, sep = "_")),
                                               contrast = o, sex = t, path = p,
                                               stringsAsFactors = FALSE)
      } else message("No ", dmr_file, " in ", p, " -> skipped")
    }
  }
  if (!length(rows)) stop("No DMR sets found under ", results_dir)
  message("Outcome folders used: ", paste(found, collapse = ", "),
          "  (excluded: ", paste(intersect(exclude, list.dirs(results_dir, FALSE, FALSE)), collapse = ", "), ")")
  man <- do.call(rbind, rows)
  if (anyDuplicated(man$set_id)) stop("set_id collision after cleaning names: ",
                                      paste(man$set_id[duplicated(man$set_id)], collapse = ", "))
  # factor order for plots/tables: user order if given, else discovery order
  man$sex <- as.character(man$sex); man$contrast <- as.character(man$contrast)
  man
}

# ----------------------------- load sets -----------------------------------------
if (!is.na(opt$manifest)) {
  man <- read_manifest(opt$manifest, c("set_id", "contrast", "sex", "path"))
  config_file <- opt$manifest
} else if (is_ot) {
  man <- discover_outcome_timepoint(opt$results_dir, opt$dmr_file, opt$subdir,
                                    split_opt(opt$outcomes), split_opt(opt$timepoints),
                                    split_opt(opt$exclude_dirs))
  config_file <- NULL
} else {
  if (is.na(opt$config_dir))   opt$config_dir   <- file.path(opt$project_dir, "config")
  if (is.na(opt$dmrichr_root)) opt$dmrichr_root <- file.path(opt$project_dir, "results", "dmrichr_output")
  config_file <- file.path(opt$config_dir, opt$sets_file)
  man <- dmr_sets_from_config(config_file, opt$dmrichr_root, opt$cutoff,
                              trimws(strsplit(opt$sexes, ",")[[1]]), sex_lab, opt$dmr_file)
}
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
{ mo <- man; names(mo)[match(names(dim_names), names(mo))] <- dim_names
  write.table(mo, file.path(opt$out_dir, "sets_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE) }
check_path_label(man, "sex"); check_path_label(man, "contrast")
man$dmr_path <- file.path(man$path, opt$dmr_file)
man$bg_path  <- if (tolower(opt$background_file) == "none") NA else file.path(man$path, opt$background_file)

dmrs <- setNames(lapply(seq_len(nrow(man)), function(i)
  read_dmr_set(man$dmr_path[i], man$set_id[i], opt$effect_col, opt$gene_col, opt$gene_sep,
               opt$exclude_sex_chr, opt$max_pval, opt$max_qval)), man$set_id)
bgs <- setNames(lapply(man$bg_path, read_background, exclude_sex_chr = opt$exclude_sex_chr), man$set_id)
if (any(vapply(bgs, is.null, TRUE)))
  message("No background for: ", paste(man$set_id[vapply(bgs, is.null, TRUE)], collapse = ", "),
          " -> region Fisher tests and other-sex effect lookup skipped for those sets")

stopifnot(opt$upset_genes_scope %in% c("multi", "all"), opt$upset_genes >= 0)
highlight <- NULL
if (!is.na(opt$upset_highlight_genes)) {
  hl <- read_table_any(opt$upset_highlight_genes)
  hc <- if (ncol(hl) == 1 && !any(c("geneSymbol", "gene_symbol") %in% names(hl))) {
    hl <- rbind(setNames(data.frame(names(hl)), names(hl)), hl); names(hl)[1]   # header-less list
  } else detect_col(hl, c("geneSymbol", "gene_symbol", "SYMBOL", "gene"), "gene")
  highlight <- unique(trimws(na.omit(as.character(hl[[hc]]))))
  message("UpSet labels restricted to ", length(highlight), " highlight genes")
}

gene_universe <- NULL
if (tolower(opt$gene_universe) == "background") {
  bg_genes <- unlist(lapply(bgs, function(b) if (!is.null(b) && !is.null(b$gene)) b$gene))
  if (length(bg_genes)) {
    gene_universe <- sort(unique(na.omit(bg_genes[bg_genes != ""])))
    message("Gene universe: ", length(gene_universe), " genes from background files")
  } else message("No gene symbols in background files; gene-mode Fisher tests skipped")
} else if (tolower(opt$gene_universe) != "none") {
  gu <- read_table_any(opt$gene_universe)
  gc <- detect_col(gu, c("geneSymbol", "gene_symbol", "SYMBOL", names(gu)[1]), "gene")
  gene_universe <- sort(unique(na.omit(as.character(gu[[gc]]))))
  message("Gene universe: ", length(gene_universe), " genes")
}

log_run_metadata(opt$out_dir, opt,
                 c(config_file, man$dmr_path, na.omit(man$bg_path), opt$groups_file,
                   if (file.exists(opt$gene_universe)) opt$gene_universe,
                   na.omit(c(opt$gene_set, opt$snp_set))),
                 self)
sizes <- data.frame(set_id = man$set_id, contrast = man$contrast, sex = man$sex,
                    n_dmrs = vapply(dmrs, nrow, 1L),
                    n_hyper = vapply(dmrs, function(d) sum(d$.dir > 0), 1L),
                    n_hypo = vapply(dmrs, function(d) sum(d$.dir < 0), 1L),
                    n_genes = vapply(dmrs, function(d) length(dmr_genes(d)), 1L),
                    n_background = vapply(bgs, function(b) if (is.null(b)) NA_integer_ else length(b), 1L),
                    effect_col = vapply(dmrs, function(d) attr(d, "effect_col"), ""))
names(sizes)[2:3] <- dim_names
print(sizes, row.names = FALSE)

# ----------------------------- GWAS gene set / SNP set -----------------------------
user_cols <- function(d) d[, !startsWith(names(d), "."), drop = FALSE]   # (also defined below)
gwas <- NULL; H <- identity; set_univ <- NULL; gw_label <- clean_id(opt$gwas_label)
gs_res <- NULL; snp_res <- NULL; ev_all <- NULL
gw_dir <- file.path(opt$out_dir, "gwas")
dimcols <- function(s) setNames(data.frame(man$contrast[man$set_id == s], man$sex[man$set_id == s],
                                           stringsAsFactors = FALSE), dim_names)
gin <- NULL
if (!is.na(opt$gene_set)) {
  dir.create(gw_dir, recursive = TRUE, showWarnings = FALSE)
  gin <- read_gwas_input(opt$gene_set, opt$gwas_max_p)
  gwas <- build_gwas_gene_set(gin, opt$harmonize_symbols, opt$gwas_expand_readthrough)
  message("\n", opt$gwas_label, " gene set: ", nrow(gwas$genes), " genes",
          if (!is.null(gwas$long)) paste0(" from ", length(unique(gwas$long$assoc_id)), " associations, ",
                                         length(unique(gwas$long$study)), " studies, ",
                                         length(unique(gwas$long$trait)), " traits"))
  # same symbol harmonization for the DMRichR genes (DMRs and backgrounds)
  dmr_syms <- unique(c(unlist(lapply(dmrs, dmr_genes)),
                       unlist(lapply(bgs, function(b) if (!is.null(b) && !is.null(b$gene))
                         split_genes(b$gene, opt$gene_sep)))))
  dmap <- symbol_map(dmr_syms, opt$harmonize_symbols)
  H <- function(x) apply_symbol_map(x, dmap)
  set_genes <- lapply(dmrs, function(d) unique(H(dmr_genes(d))))
  set_univ  <- lapply(bgs, function(b) if (is.null(b) || is.null(b$gene)) NULL
                      else unique(H(split_genes(b$gene, opt$gene_sep))))
  set_dirs  <- lapply(dmrs, gene_dir_table, H = H)

  # gene table of some DMR sets (harmonized symbols) + a GWAS column. GWAS genes are restricted
  # to the sets' background genes; the GWAS-only bar is never drawn. If file is given, draws the
  # UpSet: bars that include GWAS are labelled "GENE (n associations)", ranked by n associations.
  gwas_upset <- function(sids, file = NULL, title = "", cols = NULL) {
    sep <- if (is.na(opt$gene_sep)) NA else opt$gene_sep
    dlh <- lapply(dmrs[sids], function(d) { d$.gene <- harmonize_gene_col(d$.gene, H, sep); d })
    gth <- build_gene_table(dlh)
    gg <- intersect(gwas$genes$geneSymbol, unique(unlist(set_univ[sids])))
    if (is.null(gth) || !length(gg)) return(NULL)
    isg <- gth$geneSymbol %in% gg
    gth[[gw_label]] <- as.integer(isg)
    gth[[paste0("dir_", gw_label)]] <- NA_character_
    gth$pattern <- ifelse(isg, paste(gth$pattern, gw_label, sep = "&"), gth$pattern)
    gth$n_sets <- gth$n_sets + isg
    na <- gwas$genes$n_associations[match(gth$geneSymbol, gwas$genes$geneSymbol)]
    if (!is.null(file)) {
      alld <- do.call(rbind, lapply(dlh, function(d) d[!is.na(d$.gene), c(".gene", ".effect")]))
      eff <- tapply(abs(alld$.effect), alld$.gene, max)[gth$geneSymbol]
      gth$label_gene <- ifelse(isg, paste0(gth$geneSymbol, " (", na, ")"), gth$geneSymbol)
      if (is.null(cols)) cols <- setNames(rep("grey15", length(sids)), sids)
      plot_upset_labeled(gth, c(sids, gw_label), file,
                         paste0(title, " (", length(gg), " ", opt$gwas_label, " genes in universe); labels: gene (n associations)"),
                         gene_col = "label_gene", score = ifelse(isg, na, eff),
                         n_genes = opt$gwas_upset_genes, scope = "multi", highlight = NULL,
                         require_direction = FALSE,
                         set_sizes = c(vapply(dlh, function(d) length(dmr_genes(d)), 1L), setNames(length(gg), gw_label)),
                         set_cols = c(cols, setNames("#7F7F7F", gw_label)))
    }
    invisible(list(tab = gth, isg = isg, n_in_universe = length(gg)))
  }

  gs_res <- bind_rows_any(lapply(man$set_id, function(s) {
    if (is.null(set_univ[[s]])) { message("  ", s, ": no background genes -> gene-set test skipped"); return(NULL) }
    cbind(set_id = s, dimcols(s), gene_set_enrichment_one(set_genes[[s]], set_univ[[s]], gwas$genes))
  }))
  if (nrow(gs_res)) {
    gs_res$fisher_FDR_run    <- p.adjust(gs_res$fisher_p, "BH")
    gs_res$gw_fisher_FDR_run <- p.adjust(gs_res$gw_fisher_p, "BH")
    print(gs_res[, c("set_id", "dmr_genes", "gwas_genes_in_universe", "overlap", "expected",
                     "fisher_p", "fisher_FDR_run")], row.names = FALSE)
  }
  # every DMR gene (any set) that is a GWAS gene, with direction and GWAS evidence
  ev_all <- bind_rows_any(lapply(man$set_id, function(s) {
    x <- set_dirs[[s]]; x <- x[x$geneSymbol %in% gwas$genes$geneSymbol, ]
    if (!nrow(x)) return(NULL)
    cbind(set_id = s, dimcols(s), x)
  }))
  ov_genes <- if (nrow(ev_all)) {
    names(ev_all)[names(ev_all) == "dir"] <- "dmr_direction"
    merge(ev_all, gwas$genes, by = "geneSymbol", all.x = TRUE)
  } else data.frame()
  if (nrow(ov_genes)) ov_genes <- ov_genes[order(-ov_genes$n_associations, ov_genes$geneSymbol, ov_genes$set_id), ]
  if (nrow(ev_all)) {
    names(ev_all)[names(ev_all) == "dmr_direction"] <- "dir"; ev_all$set <- ev_all$set_id
    plot_gwas_evidence(ev_all, gwas$genes, man$set_id, file.path(gw_dir, paste0(gw_label, "_gene_evidence.pdf")),
                       paste0(opt$gwas_label, " genes among DMR genes"), opt$gwas_plot_max, opt$gwas_label,
                       set_groups = man$contrast)
  }
  # outcomes + GWAS UpSets, one per outcome (across timepoints) and one per timepoint
  # (across outcomes); written to gwas/upset/
  udir <- file.path(gw_dir, "upset"); dir.create(udir, showWarnings = FALSE)
  for (o in unique(man$contrast)) {
    mg <- man[man$contrast == o, ]
    if (nrow(mg) < 1) next
    gwas_upset(mg$set_id, file.path(udir, paste0(clean_id(o), "_across_", dim_names[["sex"]], "s_with_", gw_label, ".pdf")),
               paste0(o, ": DMR genes across ", dim_names[["sex"]], "s + ", opt$gwas_label, " genes"),
               set_colours(transform(mg, contrast = sex), "custom"))
  }
  for (t in unique(man$sex)) {
    mg <- man[man$sex == t, ]
    if (nrow(mg) < 1) next
    gwas_upset(mg$set_id, file.path(udir, paste0(clean_id(t), "_across_", dim_names[["contrast"]], "s_with_", gw_label, ".pdf")),
               paste0(t, ": DMR genes across ", dim_names[["contrast"]], "s + ", opt$gwas_label, " genes"),
               set_colours(mg, "within_sex_across_conditions"))
  }
  if (nrow(gs_res))
    plot_gene_set_enrichment(gs_res, file.path(gw_dir, paste0(gw_label, "_gene_set_enrichment.pdf")),
                             paste0(opt$gwas_label, " gene-set enrichment per DMR set"))
  smap <- rbind(cbind(source = opt$gwas_label, gwas$map), cbind(source = "DMRichR", dmap))
  write_xlsx_safe(list(Genes = gwas$genes, Trait_summary = gwas$traits,
                       Associations = if (!is.null(gwas$long)) gin else NULL,
                       Gene_by_association = gwas$long,
                       Symbol_changes = smap[smap$status != "current", ]),
                  file.path(gw_dir, paste0(gw_label, "_gene_set.xlsx")))
  write_xlsx_safe(list(Enrichment = gs_res, DMR_genes_in_GWAS = ov_genes),
                  file.path(gw_dir, paste0(gw_label, "_gene_set_enrichment.xlsx")))
}

if (!is.na(opt$snp_set)) {
  dir.create(file.path(gw_dir, "regioneR_snp"), recursive = TRUE, showWarnings = FALSE)
  sin <- if (!is.null(gin) && normalizePath(opt$snp_set) == normalizePath(opt$gene_set)) gin
         else read_gwas_input(opt$snp_set, opt$gwas_max_p)
  snp <- gwas_snp_gr(sin)
  if (opt$exclude_sex_chr) snp <- snp[!as.character(seqnames(snp)) %in% c("chrX", "chrY")]
  bg1 <- Filter(Negate(is.null), bgs)
  if (length(bg1) && !any(grepl("^chr", seqlevels(bg1[[1]]))))
    seqlevels(snp) <- sub("^chr", "", seqlevels(snp))
  W <- reduce(GRanges(seqnames(snp), IRanges(pmax(1L, start(snp) - opt$snp_window),
                                             end(snp) + opt$snp_window)))
  message("\n", opt$gwas_label, " SNP set: ", length(snp), " SNP positions -> ", length(W),
          " merged windows (+/- ", opt$snp_window, " bp, ", round(sum(as.numeric(width(W))) / 1e6, 1), " Mb)")
  sr <- list(); sp <- list(); near <- list()
  for (s in man$set_id) {
    message("  regioneR SNP proximity: ", s)
    r <- perm_snp_set(to_gr(dmrs[[s]]), bgs[[s]], W, opt$perm_n, opt$perm_per_chr,
                      opt$perm_alternative, opt$perm_seed)
    if (is.null(r)) { message("    skipped (no background or too few regions)"); next }
    sr[[s]] <- cbind(set_id = s, dimcols(s), r$summary)
    sp[[s]] <- data.frame(randomized_A = s, permuted = r$permuted)
    pdf(file.path(gw_dir, "regioneR_snp", paste0("snp_perm_", s, ".pdf")), width = 7, height = 5)
    plot(r$pt, main = paste0(s, " DMRs within ", opt$snp_window / 1000, " kb of ", opt$gwas_label, " SNPs"))
    dev.off()
    d <- dmrs[[s]]; A <- to_gr(d); hit <- which(overlapsAny(A, W))
    if (length(hit)) {
      dn <- distanceToNearest(A[hit], snp)
      j <- subjectHits(dn)[match(seq_along(hit), queryHits(dn))]
      near[[s]] <- cbind(set_id = s, dimcols(s),
                         nearest_snp = snp$snp[j], nearest_snp_pos = paste0(seqnames(snp)[j], ":", start(snp)[j]),
                         distance_bp = mcols(dn)$distance[match(seq_along(hit), queryHits(dn))],
                         gwas_trait = snp$trait[j], gwas_p = snp$gwas_p[j], gwas_study = snp$study[j],
                         n_snps_within_window = countOverlaps(A[hit] + opt$snp_window, snp),
                         user_cols(d[hit, ]))
    }
  }
  snp_res <- bind_rows_any(sr)
  if (nrow(snp_res)) {
    snp_res$perm_FDR_run <- p.adjust(snp_res$perm_p, "BH")
    print(snp_res[, c("set_id", "n_dmrs_tested", "observed", "perm_mean", "fold_enrichment", "z", "perm_p")],
          row.names = FALSE)
    pmx <- transform(snp_res, randomized_A = set_id,
                     fixed_B = paste0(opt$gwas_label, " SNPs +/-", opt$snp_window / 1000, "kb"),
                     group = paste0(opt$gwas_label, " SNP proximity"))
    prx <- transform(bind_rows_any(sp), fixed_B = pmx$fixed_B[1], group = pmx$group[1])
    plot_perm_summary(pmx, prx, file.path(gw_dir, paste0(gw_label, "_snp_proximity_summary.pdf")),
                      paste0("DMRs near ", opt$gwas_label, " SNPs vs resampled background regions"),
                      fdr_col = "perm_FDR_run",
                      xlab = paste0("DMRs within ", opt$snp_window / 1000, " kb of a SNP (grey = permuted null)"),
                      ylab = "DMR set")
  }
  snp_df <- data.frame(chr = as.character(seqnames(snp)), pos = start(snp), as.data.frame(mcols(snp)))
  write_xlsx_safe(list(Permutation = snp_res, DMRs_near_SNPs = bind_rows_any(near), SNPs = snp_df),
                  file.path(gw_dir, paste0(gw_label, "_snp_proximity.xlsx")))
}

# ----------------------------- define groups -------------------------------------
# Internal family keys stay as in the original script (helpers such as set_colours()
# may switch on them); `label` is the user-facing folder/summary name.
fam_label <- if (is_ot) {
  c(within_sex_across_conditions = "within_timepoint_across_outcomes",
    within_conditions_across_sex = "within_outcome_across_timepoints",
    custom = "custom")
} else {
  c(within_sex_across_conditions = "within_sex_across_conditions",
    within_conditions_across_sex = "within_conditions_across_sex",
    custom = "custom")
}
comps <- split_opt(opt$comparisons)
if (is.null(comps)) comps <- unname(fam_label[1:2])
# accept either the layout-specific or the original names
comps <- ifelse(comps %in% fam_label, names(fam_label)[match(comps, fam_label)], comps)
bad <- setdiff(comps, names(fam_label)[1:2])
if (length(bad)) stop("Unknown --comparisons: ", paste(bad, collapse = ", "),
                      " (use ", paste(fam_label[1:2], collapse = " and/or "), ")")
groups <- list()
if ("within_sex_across_conditions" %in% comps) for (s in unique(man$sex))
  groups[[clean_id(s)]] <- list(family = "within_sex_across_conditions", sets = man$set_id[man$sex == s])
if ("within_conditions_across_sex" %in% comps) for (ct in unique(man$contrast))
  groups[[clean_id(ct)]] <- list(family = "within_conditions_across_sex", sets = man$set_id[man$contrast == ct])
if (!is.na(opt$groups_file)) {
  gf <- read.delim(opt$groups_file, stringsAsFactors = FALSE)
  miss <- setdiff(gf$set_id, man$set_id)
  if (length(miss)) message("groups_file set_ids not found, dropped: ", paste(miss, collapse = ", "))
  gf <- gf[gf$set_id %in% man$set_id, ]
  for (g in unique(gf$group_id)) groups[[paste0("custom_", clean_id(g))]] <-
    list(family = "custom", sets = gf$set_id[gf$group_id == g])
}
groups <- groups[vapply(groups, function(g) length(g$sets) >= 2, TRUE)]
if (!length(groups)) stop("No group has >= 2 sets")

# strip helper columns for user-facing tables
user_cols <- function(d) d[, !startsWith(names(d), "."), drop = FALSE]

# one folder per overlap pattern with its annotated DMRs
save_per_overlap <- function(per_pat, base) {
  if (!isTRUE(opt$per_overlap_folders)) return(invisible())
  for (nm in names(per_pat)) {
    df <- per_pat[[nm]]
    if (is.null(df) || !nrow(df)) next
    pname <- sub("^ov_", "", nm)
    sub_dir <- file.path(base, "overlaps", pname)
    dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
    gcol <- intersect(c(opt$gene_col, "geneSymbol"), names(df))[1]
    genes <- if (!is.na(gcol)) data.frame(geneSymbol = sort(unique(na.omit(as.character(df[[gcol]])))))
             else data.frame(geneSymbol = character(0))
    write_xlsx_safe(list(DMRs = df, Genes = genes),
                    file.path(sub_dir, paste0("DMRs_annotated_", pname, ".xlsx")))
  }
}

# effect of a set at each locus: DMR effect if significant there, else background effect
effect_at_loci <- function(loci, set_id, dmr_map) {
  gr <- GRanges(loci$chr, IRanges(loci$start, loci$end))
  out <- data.frame(sig = loci[[set_id]] == 1, tested = NA, effect = NA_real_, bg_pval = NA_real_)
  d <- dmr_map[dmr_map$.set_id == set_id, ]
  e <- tapply(d$.effect, d$locus_id, function(v) v[which.max(abs(v))])
  out$effect[out$sig] <- e[loci$locus_id[out$sig]]
  bg <- bgs[[set_id]]
  if (!is.null(bg)) {
    h <- findOverlaps(gr, bg)
    out$tested <- seq_len(nrow(loci)) %in% queryHits(h)
    if (!is.null(bg$effect) && length(h)) {
      hp <- if (!is.null(bg$pval)) bg$pval[subjectHits(h)] else rep(1, length(h))
      best <- tapply(seq_along(h), queryHits(h), function(k) k[which.min(hp[k])])
      qi <- as.integer(names(best)); si <- subjectHits(h)[unlist(best)]
      fill <- !out$sig[qi]
      out$effect[qi[fill]] <- bg$effect[si[fill]]
      if (!is.null(bg$pval)) out$bg_pval[qi] <- bg$pval[si]
    }
  }
  out
}

# ----------------------------- within_conditions_across_sex: sex classification ------
# (config_sex layout only; skipped for outcome_timepoint)
classify_sex <- function(loci, sets, man_g) {
  sx <- setNames(names(sex_lab)[match(man_g$sex, sex_lab)], man_g$set_id)
  if (any(is.na(sx))) return(NULL)
  sid <- setNames(names(sx), sx)                       # M -> set_id
  has <- function(k) if (k %in% names(sid)) loci[[sid[[k]]]] == 1 else rep(FALSE, nrow(loci))
  M <- has("M"); F <- has("F"); C <- has("C")
  loci$sex_pattern <- paste0(ifelse(M, "M", ""), ifelse(F, "F", ""), ifelse(C, "C", ""))
  loci$sex_class <- ifelse(M & F, "Shared (M and F)",
                    ifelse(M, "Male-only", ifelse(F, "Female-only", "Combined-only")))
  if (all(c("M", "F") %in% names(sid))) {
    em <- effect_at_loci(loci, sid[["M"]], built$dmr_map)
    ef <- effect_at_loci(loci, sid[["F"]], built$dmr_map)
    loci$M_effect <- em$effect; loci$F_effect <- ef$effect
    loci$M_tested <- em$tested; loci$F_tested <- ef$tested
    loci$M_bg_pval <- em$bg_pval; loci$F_bg_pval <- ef$bg_pval
    other <- function(sig_eff, oth_eff, oth_tested, oth_p) {
      ifelse(is.na(oth_tested), "no background provided",
      ifelse(!oth_tested, "not tested in other sex",
      ifelse(is.na(oth_eff), "tested, no effect estimate",
      ifelse(!is.na(oth_p) & oth_p < opt$subthreshold_p & sign(oth_eff) == sign(sig_eff),
             "same direction, sub-threshold",
      ifelse(!is.na(oth_p) & oth_p < opt$subthreshold_p, "opposite direction, sub-threshold",
             "no evidence in other sex")))))
    }
    loci$other_sex_evidence <- NA_character_
    i <- loci$sex_class == "Male-only"
    loci$other_sex_evidence[i] <- other(em$effect[i], ef$effect[i], ef$tested[i], ef$bg_pval[i])
    i <- loci$sex_class == "Female-only"
    loci$other_sex_evidence[i] <- other(ef$effect[i], em$effect[i], em$tested[i], em$bg_pval[i])
    i <- loci$sex_class == "Shared (M and F)"
    loci$other_sex_evidence[i] <- ifelse(sign(em$effect[i]) == sign(ef$effect[i]),
                                         "both significant, same direction",
                                         "both significant, OPPOSITE direction")
  }
  loci
}

# ----------------------------- run each group ------------------------------------
all_patterns <- list(); all_pairwise <- list(); all_sexclass <- list(); all_perm <- list(); all_permuted <- list()

for (gname in names(groups)) {
  g <- groups[[gname]]; sets <- g$sets
  flab <- fam_label[[g$family]]
  gdir <- file.path(opt$out_dir, flab, gname)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  dl <- dmrs[sets]; man_g <- man[match(sets, man$set_id), ]
  # outcome_timepoint: sets within one outcome share a contrast, so colour them by timepoint
  set_cols <- if (is_ot && g$family == "within_conditions_across_sex")
    set_colours(transform(man_g, contrast = sex), "custom") else set_colours(man_g, g$family)
  message("\n== ", flab, " / ", gname, ": ", paste(sets, collapse = " | "))
  built <- NULL

  if (do_coord) {
    cdir <- file.path(gdir, "coord"); dir.create(cdir, showWarnings = FALSE)
    built <- build_loci(dl, opt$max_gap, opt$min_overlap_bp)
    loci <- built$loci
    pc <- pattern_counts(loci, sets, opt$require_direction)
    pw <- pairwise_coord(dl, bgs[sets], loci, opt$max_gap, opt$min_overlap_bp)
    plot_upset(loci, sets, file.path(cdir, paste0(gname, "_upset_coord.pdf")),
               paste0(gname, " loci", if (opt$require_direction) " (direction-concordant)"),
               opt$require_direction, set_cols)
    {
      score <- tapply(abs(built$dmr_map$.effect), built$dmr_map$locus_id, max)[loci$locus_id]
      plot_upset_labeled(loci, sets, file.path(cdir, paste0(gname, "_upset_coord_colour.pdf")),
                         paste0(gname, " loci", if (opt$require_direction) " (direction-concordant)"),
                         gene_col = "genes", score = score, n_genes = opt$upset_genes,
                         scope = opt$upset_genes_scope, highlight = highlight,
                         require_direction = opt$require_direction,
                         set_sizes = vapply(dl, nrow, 1L), set_cols = set_cols)
    }
    plot_pairwise_heatmap(pw, "jaccard_loci", file.path(cdir, paste0(gname, "_pairwise_jaccard_coord.pdf")),
                          paste(gname, "- locus Jaccard"), sets)
    if (any(!is.na(pw$fisher_OR)))
      plot_pairwise_heatmap(pw, "fisher_OR", file.path(cdir, paste0(gname, "_pairwise_OR_coord.pdf")),
                            paste(gname, "- enrichment OR (tested-in-both universe)"), sets)
    # per-intersection annotated DMR tables (all original DMRichR columns kept)
    multi <- pc$pattern[pc$n_sets >= 2]
    per_pat <- lapply(multi, function(p) {
      lid <- loci$locus_id[loci$pattern == p &
                           (!opt$require_direction | loci$direction_concordant %in% TRUE)]
      rows <- built$dmr_map[built$dmr_map$locus_id %in% lid, c(".dmr_id", "locus_id")]
      bind_rows_any(lapply(sets[sets %in% strsplit(p, "&")[[1]]], function(s) {
        d <- dl[[s]]; k <- match(rows$.dmr_id, d$.dmr_id); k <- k[!is.na(k)]
        if (!length(k)) return(NULL)
        cbind(Overlap = p, Source = s, Locus_ID = rows$locus_id[match(d$.dmr_id[k], rows$.dmr_id)],
              user_cols(d[k, ]))
      }))
    })
    names(per_pat) <- paste0("ov_", gsub("&", "+", multi))
    save_per_overlap(per_pat, cdir)
    sheets <- c(list(Patterns = pc, Pairwise = pw, Loci = loci), per_pat)
    if (opt$perm_test) {
      perm <- perm_pairwise(dl, bgs[sets], opt$perm_n, opt$max_gap, opt$min_overlap_bp,
                            opt$perm_randomize, perm_genome, perm_mask,
                            opt$perm_per_chr, opt$perm_alternative, opt$perm_local_z,
                            opt$perm_lz_window, opt$perm_lz_step, opt$perm_seed,
                            pdf_dir = file.path(cdir, "regioneR"))
      if (!is.null(perm)) {
        pr <- attr(perm, "permuted"); attr(perm, "permuted") <- NULL
        all_permuted[[gname]] <- cbind(group = gname, pr)
        plot_perm_summary(cbind(group = gname, perm), cbind(group = gname, pr),
                          file.path(cdir, "regioneR", paste0(gname, "_permutation_summary.pdf")),
                          paste0(gname, " - regioneR overlap test (", opt$perm_randomize, "Regions)"),
                          fdr_col = "perm_FDR_group")
        sheets <- c(list(Patterns = pc, Pairwise = pw, Permutation = perm, Loci = loci), per_pat)
        all_perm[[gname]] <- cbind(comparison = flab, group = gname, perm)
        print(perm[, c("randomized_A", "fixed_B", "observed", "perm_mean", "fold_enrichment",
                       "z", "perm_p")], row.names = FALSE)
      }
    }
    if (g$family == "within_conditions_across_sex" && !is_ot) {
      sc <- classify_sex(loci, sets, man_g)
      if (!is.null(sc)) {
        sheets$Sex_classification <- sc
        sc_sum <- as.data.frame(table(sex_class = sc$sex_class,
                                      other_sex_evidence = sc$other_sex_evidence, useNA = "ifany"))
        sc_sum <- sc_sum[sc_sum$Freq > 0, ]
        sheets$Sex_class_summary <- sc_sum
        all_sexclass[[gname]] <- cbind(group = gname, contrast = man_g$contrast[1], sc)
        # effect-size concordance between sexes at loci significant in either sex
        e <- sc[sc$sex_class != "Combined-only" & !is.na(sc$M_effect) & !is.na(sc$F_effect), ]
        # loci that cannot be plotted: significant in one sex, no overlapping region in the other
        drop_M <- sum(sc$sex_class == "Male-only" & is.na(sc$F_effect))
        drop_F <- sum(sc$sex_class == "Female-only" & is.na(sc$M_effect))
        n_M <- sum(sc$sex_class == "Male-only"); n_F <- sum(sc$sex_class == "Female-only")
        if (nrow(e) >= 5) {
          rho <- suppressWarnings(cor.test(e$M_effect, e$F_effect, method = "spearman"))
          e$marker <- ifelse(e$sex_class == "Shared (M and F)", "significant in both sexes",
                             "significant in one sex only (other sex: candidate region, not significant)")
          p <- ggplot(e, aes(M_effect, F_effect, colour = sex_class, shape = marker)) +
            geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
            geom_point(size = 2, stroke = .8, fill = "#1A1A1A") +
            scale_shape_manual(values = c("significant in both sexes" = 21,
                                          "significant in one sex only (other sex: candidate region, not significant)" = 1)) +
            scale_colour_manual(values = c("Male-only" = "#2C7BB6", "Female-only" = "#D7191C",
                                           "Shared (M and F)" = "#1A1A1A")) +
            labs(x = "Methylation difference, males (%)", y = "Methylation difference, females (%)",
                 colour = NULL, shape = NULL, title = paste(gname, "- cross-sex effect concordance"),
                 subtitle = sprintf(paste0("Spearman rho = %.2f (p = %.2g), n = %d loci plotted\n",
                                           "Not plotted (no overlapping region in other sex): ",
                                           "%d of %d male-only, %d of %d female-only"),
                                    rho$estimate, rho$p.value, nrow(e), drop_M, n_M, drop_F, n_F),
                 caption = paste0("Hollow = value on the non-significant axis comes from that sex's DMRichR ",
                                  "background (candidate) regions.\nIf candidates require the difference to pass ",
                                  "the cutoff, values near 0 cannot appear on that axis.")) +
            theme_bw(11) +
            theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
                  plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
                  plot.subtitle = element_text(size = 9.5))
          if (opt$scatter_genes > 0 || !is.null(highlight)) {
            # label priority: significant in both sexes first, then largest combined effect
            e$first_gene <- vapply(strsplit(ifelse(is.na(e$genes), "", e$genes), ";"),
                                   function(v) if (length(v)) paste(head(v, 2), collapse = "/") else "", "")
            e$rank_key <- (e$sex_class == "Shared (M and F)") * 1e6 + abs(e$M_effect) + abs(e$F_effect)
            cand <- e[e$first_gene != "", ]
            top <- head(cand[order(-cand$rank_key), ], opt$scatter_genes)
            if (!is.null(highlight)) {
              hl_rows <- cand[vapply(strsplit(cand$genes, ";"), function(v) any(v %in% highlight), TRUE), ]
              top <- unique(rbind(top, hl_rows))
            }
            if (nrow(top)) {
              p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
                ggrepel::geom_text_repel(data = top, inherit.aes = FALSE,
                                         mapping = aes(M_effect, F_effect, label = first_gene),
                                         size = 3, fontface = "italic", colour = "black",
                                         min.segment.length = 0, segment.colour = "grey50",
                                         box.padding = .4, max.overlaps = Inf, seed = 1)
              else geom_text(data = top, aes(M_effect, F_effect, label = first_gene), inherit.aes = FALSE,
                             size = 3, fontface = "italic", vjust = -0.8)
              sheets$Scatter_labeled_genes <- top[, intersect(c("locus_id", "chr", "start", "end", "genes",
                                                                "sex_class", "M_effect", "F_effect",
                                                                "other_sex_evidence"), names(top))]
            }
          }
          ggsave(file.path(cdir, paste0(gname, "_cross_sex_effect_concordance.pdf")), p, width = 7, height = 7)
        }
      }
    }
    write_xlsx_safe(sheets, file.path(cdir, paste0(gname, "_overlap_coord.xlsx")))
    all_patterns[[paste0(gname, "_coord")]] <- cbind(comparison = flab, group = gname, mode = "coord", pc)
    all_pairwise[[paste0(gname, "_coord")]] <- cbind(comparison = flab, group = gname, mode = "coord", pw)
    print(pc[, c("pattern", "n", "n_direction_concordant")], row.names = FALSE)
  }

  if (do_gene) {
    gdir2 <- file.path(gdir, "gene"); dir.create(gdir2, showWarnings = FALSE)
    gt <- build_gene_table(dl)
    if (!is.null(gt)) {
      pcg <- pattern_counts(gt, sets, opt$require_direction)
      pwg <- pairwise_gene(dl, gene_universe, if (do_coord) built$loci else NULL)
      plot_upset(gt, sets, file.path(gdir2, paste0(gname, "_upset_gene.pdf")),
                 paste0(gname, " genes", if (opt$require_direction) " (direction-concordant)"),
                 opt$require_direction, set_cols)
      {
        alld <- do.call(rbind, lapply(dl, function(d) d[!is.na(d$.gene), c(".gene", ".effect")]))
        score <- tapply(abs(alld$.effect), alld$.gene, max)[gt$geneSymbol]
        plot_upset_labeled(gt, sets, file.path(gdir2, paste0(gname, "_upset_gene_colour.pdf")),
                           paste0(gname, " genes", if (opt$require_direction) " (direction-concordant)"),
                           gene_col = "geneSymbol", score = score, n_genes = opt$upset_genes,
                           scope = opt$upset_genes_scope, highlight = highlight,
                           require_direction = opt$require_direction,
                           set_sizes = vapply(dl, function(d) length(dmr_genes(d)), 1L),
                           set_cols = set_cols)
      }
      plot_pairwise_heatmap(pwg, "jaccard", file.path(gdir2, paste0(gname, "_pairwise_jaccard_gene.pdf")),
                            paste(gname, "- gene Jaccard"), sets)
      multi <- pcg$pattern[pcg$n_sets >= 2]
      per_pat <- lapply(multi, function(p) {
        gg <- gt$geneSymbol[gt$pattern == p & (!opt$require_direction | gt$direction_concordant %in% TRUE)]
        bind_rows_any(lapply(strsplit(p, "&")[[1]], function(s) {
          d <- dl[[s]]; d <- d[!is.na(d$.gene) & d$.gene %in% gg, ]
          if (!nrow(d)) return(NULL)
          cbind(Overlap = p, Source = s, user_cols(d))
        }))
      })
      names(per_pat) <- paste0("ov_", gsub("&", "+", multi))
      save_per_overlap(per_pat, gdir2)
      gw_sheet <- NULL
      if (!is.null(gwas)) {          # GWAS genes among this group's DMR genes (plots are in gwas/)
        gu <- gwas_upset(sets)
        if (!is.null(gu) && any(gu$isg)) {
          gw_sheet <- merge(gu$tab[gu$isg, c("geneSymbol", "pattern", paste0("dir_", sets))], gwas$genes, by = "geneSymbol")
          gw_sheet <- gw_sheet[order(-gw_sheet$n_associations), ]
        }
      }
      write_xlsx_safe(c(list(Patterns = pcg, Pairwise = pwg, Genes = gt, GWAS_genes = gw_sheet), per_pat),
                      file.path(gdir2, paste0(gname, "_overlap_gene.xlsx")))
      all_patterns[[paste0(gname, "_gene")]] <- cbind(comparison = flab, group = gname, mode = "gene", pcg)
      all_pairwise[[paste0(gname, "_gene")]] <- cbind(comparison = flab, group = gname, mode = "gene", pwg)
    }
  }
}

# ----------------------------- cross-group summary -------------------------------
if (length(all_sexclass)) {
  sc <- bind_rows_any(all_sexclass)
  sc$other_sex_evidence[is.na(sc$other_sex_evidence)] <- "n/a (combined-only)"
  lv <- c("both significant, same direction", "both significant, OPPOSITE direction",
          "same direction, sub-threshold", "opposite direction, sub-threshold",
          "no evidence in other sex", "not tested in other sex", "tested, no effect estimate",
          "no background provided", "n/a (combined-only)")
  sc$other_sex_evidence <- factor(sc$other_sex_evidence, intersect(lv, unique(sc$other_sex_evidence)))
  sc$sex_class <- factor(sc$sex_class, c("Shared (M and F)", "Male-only", "Female-only", "Combined-only"))
  p <- ggplot(sc, aes(sex_class, fill = other_sex_evidence)) + geom_bar() +
    facet_wrap(~contrast, nrow = 1) +
    scale_fill_manual(values = c("#1A1A1A", "#E66101", "#5E9ED6", "#FDB863", "#D9D9D9",
                                 "#F7F7F7", "#BDBDBD", "#969696", "#7FBF7B")[seq_along(levels(sc$other_sex_evidence))],
                      drop = FALSE) +
    labs(x = NULL, y = "DMR loci", fill = "Evidence in the other sex",
         title = "Within condition, across sex - sex specificity of DMR loci",
         subtitle = "Male-/Female-only loci split by what the other sex shows at the same region") +
    theme_bw(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 2))
  ggsave(file.path(opt$out_dir, "within_conditions_across_sex", "sex_specificity_summary.pdf"), p, width = 11, height = 6)
}
pw_all <- bind_rows_any(all_pairwise)
if (nrow(pw_all)) pw_all$fisher_FDR_run <- p.adjust(pw_all$fisher_p, "BH")   # across every test in this run
summ <- list(Set_sizes = sizes, Patterns = bind_rows_any(all_patterns), Pairwise = pw_all,
             GWAS_gene_set = gs_res, GWAS_snp_proximity = snp_res)
if (length(all_perm)) {
  pm <- bind_rows_any(all_perm)
  pm$perm_FDR_run <- p.adjust(pm$perm_p, "BH")
  summ$Permutation <- pm
  # null distribution vs observed for every test in the run (FDR across the run)
  pr_all <- bind_rows_any(all_permuted)
  plot_perm_summary(transform(pm, group = paste(comparison, group, sep = " / ")),
                    transform(pr_all, group = paste(pm$comparison[match(pr_all$group, pm$group)],
                                                    group, sep = " / ")),
                    file.path(opt$out_dir, "regioneR_permutation_summary.pdf"),
                    paste0("regioneR overlap test (", opt$perm_randomize, "Regions), all groups"),
                    fdr_col = "perm_FDR_run")
  # overview: z-score per directed pair, faceted by group
  pm$label <- paste0(sprintf("%.1f", pm$z), ifelse(pm$perm_FDR_run < 0.05, "*", ""))
  p <- ggplot(pm, aes(fixed_B, randomized_A, fill = z)) + geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 3.2) +
    scale_fill_gradient2(low = "#5E3C99", mid = "#F7F7F7", high = "#E66101", midpoint = 0) +
    facet_wrap(~ paste(comparison, group, sep = " / "), scales = "free") +
    labs(x = "B (fixed)", y = "A (resampled from its candidate regions)", fill = "z",
         title = "regioneR permutation test of DMR overlap",
         caption = sprintf("%sRegions, %d permutations; * FDR < 0.05 across all tests in this run",
                           opt$perm_randomize, opt$perm_n)) +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(opt$out_dir, "regioneR_permutation_overview.pdf"), p,
         width = max(7, 3.5 * min(3, length(unique(pm$group)))),
         height = max(5, 3 * ceiling(length(unique(pm$group)) / 3)))
}
write_xlsx_safe(summ, file.path(opt$out_dir, "overlap_summary_all_groups.xlsx"))
message("\nDone. Outputs in ", normalizePath(opt$out_dir))