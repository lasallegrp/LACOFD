#!/usr/bin/env Rscript
# ============================================================================
# 02_enrichment_plot.R
#
# Enrichment of DMRs (e.g. from DMRichR), set up the same way as the comethyl
# 12a/12b scripts so DMR and co-methylation results are directly comparable.
#
# GENE-BASED (primary)
#   1. clusterProfiler ORA -- enricher() on the offline gene-set cache
#      (GO BP/MF/CC, KEGG, Reactome, Hallmark; the same RDS as comethyl 12b),
#      gene list = DMR genes, universe = BACKGROUND genes.
#   2. Enrichr (online) -- same gene list, background = BACKGROUND genes
#      (Speedrichr; enrichR >= 3.3). GWAS Catalog and any --extra_dbs here.
#
# REGION-BASED (secondary, optional)
#   3. GREAT (offline rGREAT rule) -- DMRs vs background regions, binomial
#      test over regulatory domains (same engine as comethyl 12b).
#
# BACKGROUND IS USED BY DEFAULT (--use_background TRUE). The right background
# for DMRs is the set of regions dmrseq actually TESTED, not the genome:
# testing against all genes inflates terms for genes that simply carry many
# testable CpGs. Provide it as --background_regions (BED/xlsx/csv/tsv of the
# tested regions) or --background_genes (a gene list made with the SAME
# annotation method as the DMR genes). Without it the script stops, unless
# you explicitly pass --use_background FALSE (genome-wide; logged as such).
#
# HOW GENES ARE ASSIGNED (--gene_assignment)
#   great      DMR genes and background genes both from the GREAT rule
#              (basalPlusExt 5 kb/1 kb, up to 1 Mb, up to 2 genes per region),
#              offline -- identical to comethyl 12a/15. Needs --dmrs and
#              --background_regions with coordinates.
#   annotation DMR genes from the geneSymbol column (DMRichR's own nearest-gene
#              annotation); background from --background_genes, or from a
#              geneSymbol column in --background_regions. Foreground and
#              background must come from the same annotation.
#   auto       (default) great when DMR and background coordinates are given,
#              otherwise annotation.
#
# PLOTS: --mode publication (one panel per database) or presentation (one
# combined overview). Plotted terms are filtered by --pval_cutoff on raw p,
# or on adjusted p with --use_adjusted. GO terms are slimmed with rrvgo.
# FULL RESULTS: one workbook per method with every tested term and its genes.
#
# enrichR is never library()'d (its .onAttach can fail offline); it is only
# called as enrichR::enrichr(), so the script runs on compute nodes and
# Enrichr is skipped with a logged reason.
#
# USAGE
#   # GREAT gene assignment, gene-based + GREAT, background = tested regions
#   Rscript 02_enrichment_plot.R \
#     --dmrs DMRs/DMRs_annotated.xlsx \
#     --background_regions DMRs/backgroundRegions.bed \
#     --gene_set_cache /path/reference/gene_sets_hg38.rds \
#     --outdir enrichment --prefix PPWR6m
#
#   # DMRichR's own annotation for genes; background gene list
#   Rscript 02_enrichment_plot.R --genes DMRs/DMRs_annotated.xlsx \
#     --background_genes background_genes.txt --gene_assignment annotation \
#     --gene_set_cache /path/reference/gene_sets_hg38.rds
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggsci)
  library(glue)
  library(purrr)
  library(Hmisc)
  library(openxlsx)
})

# ---- CLI --------------------------------------------------------------------
option_list <- list(
  # inputs
  make_option(c("-g", "--genes"), type = "character", default = NULL,
              help = "DMR gene list: .txt (one symbol per line) or .xlsx/.csv/.tsv with a geneSymbol column. Optional if --dmrs has a geneSymbol column or --gene_assignment great."),
  make_option("--dmrs", type = "character", default = NULL,
              help = "DMR regions: .bed or .xlsx/.csv/.tsv with seqnames|chr, start, end (DMRs_annotated.xlsx works). Needed for GREAT."),
  make_option("--background_regions", type = "character", default = NULL,
              help = "Background = regions tested by dmrseq (.bed or table with seqnames|chr, start, end; may also carry geneSymbol)."),
  make_option("--background_genes", type = "character", default = NULL,
              help = "Background gene list (.txt or table with geneSymbol), same annotation method as the DMR genes."),
  make_option("--use_background", type = "character", default = "TRUE",
              help = "Use the background for clusterProfiler, Enrichr and GREAT [default %default]. FALSE = genome-wide (not recommended)."),
  make_option("--gene_assignment", type = "character", default = "auto",
              help = "auto | great | annotation [default %default]. See header."),
  make_option("--bed_zero_based", type = "character", default = "TRUE",
              help = "TRUE = .bed starts are 0-based (BED standard; DMRichR writes BED this way) and get +1 [default %default]."),
  # outputs
  make_option(c("-o", "--outdir"), type = "character", default = "enrichment_results"),
  make_option("--prefix", type = "character", default = "enrichment"),
  # methods
  make_option("--gene_set_cache", type = "character", default = NULL,
              help = "Gene-set cache RDS (same file as comethyl 12b). Without it, GO falls back to enrichGO (offline) and KEGG to enrichKEGG (needs internet); GREAT is skipped."),
  make_option("--collections", type = "character", default = "all",
              help = "Collections from the cache: comma list of GO_BP,GO_MF,GO_CC,KEGG,Reactome,Hallmark, or all [default %default]."),
  make_option("--do_clusterprofiler", type = "character", default = "TRUE"),
  make_option("--do_enrichr", type = "character", default = "TRUE"),
  make_option("--do_great", type = "character", default = "TRUE",
              help = "Region-based GREAT; runs only when DMR + background coordinates and the cache are available [default %default]."),
  make_option("--enrichr_dbs", type = "character",
              default = "GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,Reactome_Pathways_2024,MSigDB_Hallmark_2020,GWAS_Catalog_2025",
              help = "Enrichr libraries; checked against Enrichr's current list [default %default]."),
  make_option("--extra_dbs", type = "character", default = NULL, help = "Additional Enrichr libraries, comma-separated."),
  make_option("--skip_gwas", action = "store_true", default = FALSE, help = "Drop GWAS Catalog libraries from Enrichr."),
  make_option("--skip_go", action = "store_true", default = FALSE, help = "Drop GO collections."),
  make_option("--skip_kegg", action = "store_true", default = FALSE, help = "Drop KEGG collections."),
  make_option("--min_gs_size", type = "integer", default = 10L, help = "clusterProfiler minimum gene-set size [default %default]."),
  make_option("--max_gs_size", type = "integer", default = 500L, help = "clusterProfiler maximum gene-set size [default %default]."),
  # GREAT settings (should match comethyl 12a/12b/15)
  make_option("--genome", type = "character", default = "hg38"),
  make_option("--tss_source", type = "character", default = "GREAT:hg38"),
  make_option("--great_mode", type = "character", default = "basalPlusExt"),
  make_option("--basal_upstream", type = "double", default = 5000),
  make_option("--basal_downstream", type = "double", default = 1000),
  make_option("--extension", type = "double", default = 1e6),
  make_option("--great_min_gene_set_size", type = "double", default = 5),
  make_option("--min_region_hits", type = "double", default = 5),
  make_option("--min_gene_hits", type = "double", default = 2),
  make_option("--fdr_cut", type = "double", default = 0.05),
  # Enrichr
  make_option("--enrichr_retries", type = "integer", default = 3L),
  make_option("--enrichr_sleep", type = "double", default = 2),
  make_option("--min_bg_genes", type = "integer", default = 200L),
  # plots
  make_option("--mode", type = "character", default = "publication"),
  make_option("--top_n", type = "integer", default = NULL),
  make_option("--pval_cutoff", type = "double", default = 0.05),
  make_option("--use_adjusted", action = "store_true", default = FALSE),
  make_option("--slim_threshold", type = "double", default = 0.7),
  make_option("--annoDb", type = "character", default = "org.Hs.eg.db"),
  make_option("--organism", type = "character", default = "hsa", help = "KEGG organism for the enrichKEGG fallback."),
  make_option("--fig_format", type = "character", default = "pdf"),
  make_option("--skip_full_results", action = "store_true", default = FALSE),
  make_option("--helper_file", type = "character", default = NULL,
              help = "comethyl helper.R (GREAT + gene-set cache functions) [default: next to this script]."),
  # accepted for old command lines; see messages below
  make_option("--go_method", type = "character", default = NULL),
  make_option("--kegg_method", type = "character", default = NULL),
  make_option("--go_dbs", type = "character", default = NULL),
  make_option("--kegg_db", type = "character", default = NULL),
  make_option("--gwas_db", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

tb <- function(x, nm) {
  v <- toupper(trimws(as.character(x)))
  if (v %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
  if (v %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
  stop("--", nm, " must be TRUE or FALSE", call. = FALSE)
}
for (nm in c("use_background", "bed_zero_based", "do_clusterprofiler", "do_enrichr", "do_great")) opt[[nm]] <- tb(opt[[nm]], nm)
opt$mode <- match.arg(opt$mode, c("publication", "presentation"))
opt$gene_assignment <- match.arg(opt$gene_assignment, c("auto", "great", "annotation"))
opt$fig_format <- match.arg(tolower(opt$fig_format), c("pdf", "png"))
if (is.null(opt$top_n)) opt$top_n <- if (opt$mode == "presentation") 3 else 7
if (!is.null(opt$go_method) || !is.null(opt$kegg_method) || !is.null(opt$go_dbs) || !is.null(opt$kegg_db)) {
  message("Note: --go_method/--kegg_method/--go_dbs/--kegg_db are no longer used. GO/KEGG run through clusterProfiler ",
          "on the gene-set cache with the background; Enrichr's own GO/KEGG libraries are set with --enrichr_dbs.")
}
if (!is.null(opt$gwas_db)) opt$enrichr_dbs <- paste(c(setdiff(trimws(strsplit(opt$enrichr_dbs, ",")[[1]]),
                                                              grep("GWAS", trimws(strsplit(opt$enrichr_dbs, ",")[[1]]), value = TRUE)),
                                                      opt$gwas_db), collapse = ",")

script_dir <- tryCatch({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else getwd()
}, error = function(e) getwd())
helper_file <- if (!is.null(opt$helper_file)) opt$helper_file else file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) stop("helper.R not found at ", helper_file, ". Pass --helper_file (the comethyl helper.R).", call. = FALSE)
source(helper_file)

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(opt$outdir, paste0(opt$prefix, "_run_log.txt"))
L <- function(...) append_log(log_file, ...)

# ---- Readers ----------------------------------------------------------------
read_table_any <- function(path) {
  validate_file_exists(path, path)
  ext <- tolower(tools::file_ext(sub("[.]gz$", "", path)))
  if (ext %in% c("xlsx", "xls")) return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  if (ext == "csv") return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  if (ext == "bed") {
    d <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
    d <- d[!grepl("^(track|browser)", d[[1]]), , drop = FALSE]
    out <- data.frame(chr = as.character(d[[1]]), start = as.integer(d[[2]]) + if (opt$bed_zero_based) 1L else 0L,
                      end = as.integer(d[[3]]), stringsAsFactors = FALSE)
    return(out)
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}
coords_of <- function(d, label) {
  cn <- tolower(colnames(d))
  ci <- match(c("chr", "seqnames", "chrom", "chromosome")[c("chr", "seqnames", "chrom", "chromosome") %in% cn][1], cn)
  si <- match("start", cn); ei <- match("end", cn)
  if (any(is.na(c(ci, si, ei)))) return(NULL)
  out <- data.frame(chr = as.character(d[[ci]]), start = as.integer(d[[si]]), end = as.integer(d[[ei]]), stringsAsFactors = FALSE)
  out <- out[!is.na(out$chr) & !is.na(out$start) & !is.na(out$end), , drop = FALSE]
  if (!all(grepl("^chr", out$chr))) out$chr <- ifelse(grepl("^chr", out$chr), out$chr, paste0("chr", out$chr))
  out <- unique(out)
  out$RegionID <- paste0(label, "_", seq_len(nrow(out)))
  out
}
symbols_of <- function(d) {
  if (is.character(d)) return(unique(trimws(d[nzchar(trimws(d))])))
  gc <- intersect(c("geneSymbol", "gene_symbol", "SYMBOL", "symbol", "gene"), colnames(d))[1]
  if (is.na(gc)) return(NULL)
  g <- unlist(strsplit(as.character(d[[gc]]), "[|;,]"))
  g <- trimws(g); unique(g[!is.na(g) & nzchar(g) & g != "NA"])
}
read_genes <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("txt", "")) return(symbols_of(readLines(path, warn = FALSE)))
  g <- symbols_of(read_table_any(path))
  if (is.null(g)) stop("No geneSymbol column in ", path, call. = FALSE)
  g
}

dmr_tab <- if (!is.null(opt$dmrs)) read_table_any(opt$dmrs) else NULL
bg_tab  <- if (!is.null(opt$background_regions)) read_table_any(opt$background_regions) else NULL
dmr_regions <- if (!is.null(dmr_tab)) coords_of(dmr_tab, "DMR") else NULL
bg_regions  <- if (!is.null(bg_tab)) coords_of(bg_tab, "BG") else NULL
if (!is.null(opt$dmrs) && is.null(dmr_regions)) L("NOTE: --dmrs has no seqnames/chr, start, end columns; GREAT unavailable.")
if (!is.null(opt$background_regions) && is.null(bg_regions)) stop("--background_regions needs seqnames/chr, start, end columns.", call. = FALSE)

assignment <- opt$gene_assignment
if (assignment == "auto") assignment <- if (!is.null(dmr_regions) && (!is.null(bg_regions) || !opt$use_background)) "great" else "annotation"
L("Gene assignment: ", assignment, if (opt$gene_assignment == "auto") " (auto)" else "")

# ---- Foreground / background genes -------------------------------------------
extended_tss <- NULL
need_ext <- assignment == "great" || (opt$do_great && !is.null(dmr_regions) && !is.null(bg_regions) && !is.null(opt$gene_set_cache))
if (need_ext) {
  extended_tss <- great_build_extended_tss(opt$tss_source, genome = opt$genome, mode = opt$great_mode,
                                           basal_upstream = opt$basal_upstream, basal_downstream = opt$basal_downstream,
                                           extension = opt$extension, verbose = TRUE)
}
fg_entrez <- bg_entrez <- NULL
assoc_out <- list()
if (assignment == "great") {
  if (is.null(dmr_regions)) stop("--gene_assignment great needs --dmrs with coordinates.", call. = FALSE)
  a_fg <- great_region_gene_associations(dmr_regions, extended_tss)
  fg_symbols <- unique(a_fg$gene_symbol); fg_entrez <- unique(a_fg$gene_entrezID)
  assoc_out$DMR_gene_associations <- a_fg
  if (opt$use_background) {
    if (is.null(bg_regions)) stop("--use_background TRUE with GREAT assignment needs --background_regions (the regions dmrseq tested).", call. = FALSE)
    a_bg <- great_region_gene_associations(bg_regions, extended_tss)
    bg_symbols <- unique(a_bg$gene_symbol); bg_entrez <- unique(a_bg$gene_entrezID)
  } else bg_symbols <- NULL
} else {
  fg_symbols <- if (!is.null(opt$genes)) read_genes(opt$genes) else if (!is.null(dmr_tab)) symbols_of(dmr_tab) else NULL
  if (is.null(fg_symbols) || !length(fg_symbols)) stop("No DMR genes: give --genes, or --dmrs with a geneSymbol column.", call. = FALSE)
  bg_symbols <- NULL
  if (opt$use_background) {
    bg_symbols <- if (!is.null(opt$background_genes)) read_genes(opt$background_genes) else if (!is.null(bg_tab)) symbols_of(bg_tab) else NULL
    if (is.null(bg_symbols) || !length(bg_symbols)) {
      stop("--use_background TRUE but no background genes. Give --background_genes (same annotation as the DMR genes), ",
           "or --dmrs + --background_regions to use --gene_assignment great, or set --use_background FALSE.", call. = FALSE)
    }
  }
}
L("DMR regions: ", if (is.null(dmr_regions)) "none" else nrow(dmr_regions),
  "; background regions: ", if (is.null(bg_regions)) "none" else nrow(bg_regions))
L("DMR genes: ", length(fg_symbols), "; background genes: ", if (is.null(bg_symbols)) "none (genome-wide)" else length(bg_symbols))
if (!is.null(bg_symbols)) {
  out_fg <- setdiff(fg_symbols, bg_symbols)
  if (length(out_fg) > 0) {
    L("WARNING: ", length(out_fg), " of ", length(fg_symbols), " DMR genes are not in the background and are dropped",
      if (length(out_fg) > 0.05 * length(fg_symbols)) " -- check that DMR and background genes use the same annotation" else "",
      ": ", paste(head(out_fg, 10), collapse = ", "), if (length(out_fg) > 10) ", ..." else "")
    fg_symbols <- intersect(fg_symbols, bg_symbols)
  }
}
if (!length(fg_symbols)) stop("No DMR genes left after intersecting with the background.", call. = FALSE)

full_wb_new <- function() if (opt$skip_full_results) NULL else openxlsx::createWorkbook()
add_sheet <- function(wb, df, name) {
  if (is.null(wb) || is.null(df) || !nrow(df)) return(invisible(NULL))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if ("geneID" %in% names(df) && !"Genes" %in% names(df)) { df$geneID <- gsub("/", ";", df$geneID); names(df)[names(df) == "geneID"] <- "Genes" }
  nm <- substr(gsub("[^A-Za-z0-9]", "_", name), 1, 31); base <- nm; i <- 1
  while (nm %in% openxlsx::sheets(wb)) { i <- i + 1; nm <- substr(paste0(base, "_", i), 1, 31) }
  openxlsx::addWorksheet(wb, nm); openxlsx::writeData(wb, nm, df)
}
save_wb <- function(wb, file) if (!is.null(wb) && length(openxlsx::sheets(wb))) {
  openxlsx::saveWorkbook(wb, file.path(opt$outdir, file), overwrite = TRUE); L("Saved ", file)
}
gl_wb <- openxlsx::createWorkbook()
add_sheet(gl_wb, data.frame(gene_symbol = sort(fg_symbols)), "DMR_genes")
if (!is.null(bg_symbols)) add_sheet(gl_wb, data.frame(gene_symbol = sort(bg_symbols)), "Background_genes")
for (nm in names(assoc_out)) add_sheet(gl_wb, assoc_out[[nm]], nm)
save_wb(gl_wb, paste0(opt$prefix, "_gene_lists.xlsx"))

# ---- Collections -------------------------------------------------------------
gs <- NULL; collections <- parse_collections(opt$collections)
if (opt$skip_go) collections <- setdiff(collections, c("GO_BP", "GO_MF", "GO_CC"))
if (opt$skip_kegg) collections <- setdiff(collections, "KEGG")
if (!is.null(opt$gene_set_cache)) gs <- load_gene_set_cache(opt$gene_set_cache, collections)

pcol_cp <- if (opt$use_adjusted) "p.adjust" else "pvalue"
flat <- function(df, label, term_col, p_col, n = opt$top_n, cutoff = opt$pval_cutoff) {
  if (is.null(df) || !nrow(df) || !p_col %in% names(df)) return(NULL)
  d <- df[!is.na(df[[p_col]]) & df[[p_col]] < cutoff, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(d[[p_col]]), , drop = FALSE][seq_len(min(n, nrow(d))), , drop = FALSE]
  dplyr::tibble(Database = label, Term = as.character(d[[term_col]]), `-log10.p-value` = -log10(d[[p_col]]))
}
overview <- list()
add_overview <- function(method, coll, df, padj_col) {
  overview[[length(overview) + 1]] <<- data.frame(method = method, collection = coll,
    n_tested = if (is.null(df)) 0 else nrow(df),
    n_adj_sig = if (is.null(df) || !padj_col %in% names(df)) 0 else sum(df[[padj_col]] < opt$fdr_cut, na.rm = TRUE))
}

# ---- 1. clusterProfiler (gene-based, universe = background) ------------------
cp_res <- list()
if (opt$do_clusterprofiler) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    L("clusterProfiler not installed; skipping gene-based ORA.")
  } else if (!is.null(gs)) {
    if (assignment == "great") {
      sym2e <- stats::setNames(assoc_out$DMR_gene_associations$gene_entrezID, assoc_out$DMR_gene_associations$gene_symbol)
      fg_e <- unique(unname(sym2e[fg_symbols])); uni <- bg_entrez
    } else {
      m1 <- symbols_to_entrez(fg_symbols, gs); fg_e <- m1$entrez
      L("clusterProfiler: ", m1$n_mapped, "/", m1$n_in, " DMR genes mapped to Entrez")
      uni <- if (!is.null(bg_symbols)) { m2 <- symbols_to_entrez(bg_symbols, gs)
        L("clusterProfiler universe: ", m2$n_mapped, "/", m2$n_in, " background genes mapped to Entrez"); m2$entrez } else NULL
    }
    if (!is.null(uni)) fg_e <- intersect(fg_e, uni)
    eg2sym <- stats::setNames(names(gs$symbol_map), gs$symbol_map)
    if (!is.null(extended_tss)) { sm <- attr(extended_tss, "symbol_map"); eg2sym[names(sm)] <- sm }
    cp_wb <- full_wb_new()
    for (coll in collections) {
      s <- gs$collections[[coll]]
      t2g <- data.frame(term = rep(names(s), lengths(s)), gene = unlist(s, use.names = FALSE), stringsAsFactors = FALSE)
      t2n <- data.frame(term = names(gs$term_names[[coll]]), name = unname(gs$term_names[[coll]]), stringsAsFactors = FALSE)
      res <- tryCatch(suppressMessages(clusterProfiler::enricher(fg_e, TERM2GENE = t2g, TERM2NAME = t2n, universe = uni,
                                                                 minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size,
                                                                 pvalueCutoff = 1, qvalueCutoff = 1)),
                      error = function(e) { L("ERROR clusterProfiler ", coll, ": ", conditionMessage(e)); NULL })
      df <- if (is.null(res)) NULL else as.data.frame(res)
      if (!is.null(df) && nrow(df)) {
        df$Genes <- vapply(strsplit(df$geneID, "/"), function(g) { x <- eg2sym[g]; x[is.na(x)] <- g[is.na(x)]; paste(unname(x), collapse = ";") }, "")
        df$collection <- coll; cp_res[[coll]] <- df; add_sheet(cp_wb, df, coll)
      }
      add_overview("clusterProfiler", coll, df, "p.adjust")
    }
    save_wb(cp_wb, paste0(opt$prefix, "_clusterProfiler_full_results.xlsx"))
  } else {
    L("No --gene_set_cache: GO via enrichGO (offline) and KEGG via enrichKEGG (needs internet), universe = background.")
    ids <- function(sym) { if (is.null(sym)) return(NULL)
      m <- suppressMessages(clusterProfiler::bitr(sym, "SYMBOL", "ENTREZID", OrgDb = opt$annoDb)); unique(m$ENTREZID) }
    fg_e <- if (!is.null(fg_entrez)) fg_entrez else ids(fg_symbols)
    uni <- if (!is.null(bg_entrez)) bg_entrez else ids(bg_symbols)
    cp_wb <- full_wb_new()
    for (ont in c("BP", "MF", "CC")) {
      coll <- paste0("GO_", ont); if (!coll %in% collections) next
      res <- tryCatch(clusterProfiler::enrichGO(fg_e, OrgDb = opt$annoDb, ont = ont, universe = uni, readable = TRUE,
                                                minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size, pvalueCutoff = 1, qvalueCutoff = 1),
                      error = function(e) { L("ERROR enrichGO ", ont, ": ", conditionMessage(e)); NULL })
      df <- if (is.null(res)) NULL else as.data.frame(res)
      if (!is.null(df) && nrow(df)) { df$Genes <- gsub("/", ";", df$geneID); df$collection <- coll; cp_res[[coll]] <- df; add_sheet(cp_wb, df, coll) }
      add_overview("clusterProfiler", coll, df, "p.adjust")
    }
    if ("KEGG" %in% collections) {
      res <- tryCatch(clusterProfiler::enrichKEGG(fg_e, organism = opt$organism, universe = uni, pvalueCutoff = 1, qvalueCutoff = 1,
                                                  minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size),
                      error = function(e) { L("ERROR enrichKEGG (needs internet): ", conditionMessage(e)); NULL })
      if (!is.null(res)) res <- tryCatch(clusterProfiler::setReadable(res, OrgDb = opt$annoDb, keyType = "ENTREZID"), error = function(e) res)
      df <- if (is.null(res)) NULL else as.data.frame(res)
      if (!is.null(df) && nrow(df)) { df$Genes <- gsub("/", ";", df$geneID); df$collection <- "KEGG"; cp_res$KEGG <- df; add_sheet(cp_wb, df, "KEGG") }
      add_overview("clusterProfiler", "KEGG", df, "p.adjust")
    }
    save_wb(cp_wb, paste0(opt$prefix, "_clusterProfiler_full_results.xlsx"))
  }
}

# ---- 2. Enrichr (gene-based, background = background genes) -----------------
enr_res <- list()
if (opt$do_enrichr) {
  dbs <- unique(trimws(c(strsplit(opt$enrichr_dbs, ",")[[1]], if (!is.null(opt$extra_dbs)) strsplit(opt$extra_dbs, ",")[[1]])))
  dbs <- dbs[nzchar(dbs)]
  if (opt$skip_gwas) dbs <- dbs[!grepl("GWAS", dbs)]
  if (opt$skip_go) dbs <- dbs[!grepl("^GO_", dbs)]
  if (opt$skip_kegg) dbs <- dbs[!grepl("^KEGG", dbs)]
  options(enrichR.sites.base.address = "https://maayanlab.cloud/", enrichR.base.address = "https://maayanlab.cloud/Enrichr/",
          speedrichr.base.address = "https://maayanlab.cloud/speedrichr/api/", enrichR.live = TRUE, enrichR.quiet = TRUE,
          modEnrichR.use = TRUE, enrichR.sites = c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr", "OxEnrichr"))
  reason <- NULL; bg_use <- NULL
  if (!requireNamespace("enrichR", quietly = TRUE)) reason <- "enrichR not installed"
  if (is.null(reason) && opt$use_background) {
    if (length(bg_symbols) < opt$min_bg_genes) reason <- paste0("background has ", length(bg_symbols), " genes (< --min_bg_genes)")
    else if (!(utils::packageVersion("enrichR") >= "3.3" && "background" %in% names(formals(enrichR::enrichr)))) {
      reason <- paste0("enrichR ", utils::packageVersion("enrichR"), " cannot use a background (need >= 3.3); not running without one")
    } else bg_use <- bg_symbols
  }
  if (is.null(reason)) {
    avail <- tryCatch(enrichR::listEnrichrDbs(), error = function(e) e)
    if (inherits(avail, "error") || is.null(avail) || !"libraryName" %in% colnames(avail)) {
      reason <- paste0("Enrichr unreachable (", if (inherits(avail, "error")) conditionMessage(avail) else "no library list", ")")
    } else {
      bad <- setdiff(dbs, avail$libraryName); if (length(bad)) L("Enrichr: unknown libraries dropped: ", paste(bad, collapse = ", "))
      dbs <- intersect(dbs, avail$libraryName); if (!length(dbs)) reason <- "none of the requested Enrichr libraries exist"
    }
  }
  if (!is.null(reason)) {
    L("Enrichr skipped: ", reason)
  } else {
    res <- NULL; last <- NULL
    for (k in seq_len(opt$enrichr_retries)) {
      res <- tryCatch(if (is.null(bg_use)) enrichR::enrichr(genes = fg_symbols, databases = dbs)
                      else enrichR::enrichr(genes = fg_symbols, databases = dbs, background = bg_use, include_overlap = TRUE),
                      error = function(e) { last <<- conditionMessage(e); NULL })
      if (!is.null(res) && length(res)) break
      if (k < opt$enrichr_retries) Sys.sleep(opt$enrichr_sleep)
    }
    if (is.null(res) || !length(res)) L("ERROR Enrichr failed: ", if (is.null(last)) "empty result" else last) else {
      enr_wb <- full_wb_new()
      for (nm in names(res)) {
        df <- res[[nm]]; if (is.null(df) || !nrow(df)) next
        enr_res[[nm]] <- df; add_sheet(enr_wb, df, nm); add_overview("Enrichr", nm, df, "Adjusted.P.value")
      }
      save_wb(enr_wb, paste0(opt$prefix, "_Enrichr_full_results.xlsx"))
      L("Enrichr: ", length(enr_res), " libraries; background = ", if (is.null(bg_use)) "none" else paste(length(bg_use), "genes"))
    }
  }
}

# ---- 3. GREAT (region-based, background = background regions) ----------------
great_res <- list()
if (opt$do_great) {
  why <- if (is.null(gs)) "no --gene_set_cache" else if (is.null(dmr_regions)) "no DMR coordinates (--dmrs)" else
         if (is.null(bg_regions)) "no --background_regions (GREAT always needs the tested regions as background)" else NULL
  if (!is.null(why)) {
    L("GREAT skipped: ", why)
  } else {
    sets <- list(); univ <- unique(extended_tss$gene_id)
    for (coll in collections) {
      s <- lapply(gs$collections[[coll]], intersect, univ)
      s <- s[lengths(s) >= opt$great_min_gene_set_size]
      if (length(s)) names(s) <- paste0(coll, "::", names(s)); sets <- c(sets, s)
    }
    bg_gr <- GenomicRanges::GRanges(bg_regions$chr, IRanges::IRanges(bg_regions$start, bg_regions$end))
    prep <- great_fast_prepare(extended_tss, bg_gr, sets, opt$great_min_gene_set_size, Inf)
    tb1 <- great_fast_test(prep, list(DMRs = dmr_regions))$DMRs
    hit <- attr(tb1, "hit_genes"); sym <- attr(extended_tss, "symbol_map")
    tb1$collection <- sub("::.*$", "", tb1$id)
    gr_wb <- full_wb_new()
    for (coll in collections) {
      d <- tb1[tb1$collection == coll, , drop = FALSE]; if (!nrow(d)) next
      d$term_id <- sub("^[^:]+::", "", d$id)
      nmv <- gs$term_names[[coll]][d$term_id]; d$Description <- ifelse(is.na(nmv), d$term_id, unname(nmv))
      d$Genes <- vapply(d$id, function(i) { g <- intersect(sets[[i]], hit); x <- sym[g]; x[is.na(x)] <- g[is.na(x)]; paste(sort(unname(x)), collapse = ";") }, "")
      tested <- d$observed_region_hits >= opt$min_region_hits
      d$fdr <- NA_real_; d$fdr[tested] <- stats::p.adjust(d$p_value[tested], "BH")
      d$significant <- tested & !is.na(d$fdr) & d$fdr < opt$fdr_cut & d$observed_gene_hits >= opt$min_gene_hits
      d <- d[order(is.na(d$fdr), d$fdr, d$p_value), c("collection", "term_id", "Description", "observed_region_hits",
              "fold_enrichment", "p_value", "fdr", "observed_gene_hits", "gene_set_size", "p_value_hyper", "significant", "Genes")]
      great_res[[coll]] <- d; add_sheet(gr_wb, d, paste0("GREAT_", coll)); add_overview("GREAT", coll, d, "fdr")
    }
    save_wb(gr_wb, paste0(opt$prefix, "_GREAT_full_results.xlsx"))
    L("GREAT: ", nrow(dmr_regions), " DMRs vs ", nrow(bg_regions), " background regions")
  }
}

if (length(overview)) utils::write.table(do.call(rbind, overview), file.path(opt$outdir, paste0(opt$prefix, "_enrichment_overview.tsv")),
                                         sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Plot data ----------------------------------------------------------------
slim_go <- function(df, ont) {
  if (is.null(df) || !nrow(df)) return(NULL)
  d <- df[!is.na(df[[pcol_cp]]) & df[[pcol_cp]] < opt$pval_cutoff, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  fallback <- function(msg) { if (!is.null(msg)) L("GO ", ont, ": ", msg, "; plotting top terms without slimming")
    flat(df, paste0("GO ", ont), "Description", pcol_cp) }
  if (!all(grepl("^GO:", d$ID))) return(fallback("term IDs are not GO IDs (MSigDB cache)"))
  if (!requireNamespace("rrvgo", quietly = TRUE)) return(fallback("rrvgo not installed"))
  if (nrow(d) < 2) return(fallback(NULL))
  red <- tryCatch({
    sm <- rrvgo::calculateSimMatrix(d$ID, orgdb = opt$annoDb, ont = ont, method = "Rel")
    rrvgo::reduceSimMatrix(sm, stats::setNames(-log10(d[[pcol_cp]]), d$ID), threshold = opt$slim_threshold, orgdb = opt$annoDb)
  }, error = function(e) { L("GO ", ont, " slimming failed: ", conditionMessage(e)); NULL })
  if (is.null(red) || !nrow(red)) return(fallback(NULL))
  par <- unique(as.character(red$parent)); par <- par[par %in% d$ID]
  d <- d[d$ID %in% par, , drop = FALSE]
  flat(d, paste0("GO ", ont), "Description", pcol_cp)
}
cp_plot <- dplyr::bind_rows(
  lapply(c("BP", "MF", "CC"), function(o) slim_go(cp_res[[paste0("GO_", o)]], o)),
  lapply(setdiff(names(cp_res), c("GO_BP", "GO_MF", "GO_CC")), function(cn) flat(cp_res[[cn]], cn, "Description", pcol_cp)))
pcol_en <- if (opt$use_adjusted) "Adjusted.P.value" else "P.value"
enr_plot <- lapply(names(enr_res), function(nm) {
  f <- flat(enr_res[[nm]], nm, "Term", pcol_en); if (!is.null(f)) f$Term <- sub(" \\(GO:[0-9]+\\)$", "", f$Term); f })
names(enr_plot) <- names(enr_res)
pcol_gr <- if (opt$use_adjusted) "fdr" else "p_value"
great_plot <- lapply(names(great_res), function(cn) flat(great_res[[cn]], paste("GREAT", cn), "Description", pcol_gr))
names(great_plot) <- names(great_res)

# ---- Plotting -------------------------------------------------------------------
font_cfg <- if (opt$mode == "presentation") list(axis_text = 20, axis_title = 20, legend_text = 18, legend_title = 18) else
  list(axis_text = 12, axis_title = 14, legend_text = 12, legend_title = 14)
xlab_p <- if (opt$use_adjusted) expression("-log"[10]("adj. p")) else expression("-log"[10](p))
make_barplot <- function(df, title = NULL, subtitle = NULL, fonts = font_cfg) {
  df <- df %>% dplyr::mutate(Term = stringr::str_wrap(Hmisc::capitalize(stringr::str_trim(Term)), 50),
                             Term = factor(Term, levels = unique(Term[order(forcats::fct_rev(Database), `-log10.p-value`)])))
  ggplot2::ggplot(df, ggplot2::aes(Term, `-log10.p-value`, fill = Database)) +
    ggplot2::geom_col(color = "black", width = 0.8) + ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = c(0, 0, 0.1, 0)) + ggsci::scale_fill_d3() +
    ggplot2::labs(y = xlab_p, title = title, subtitle = subtitle) + ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = fonts$axis_text, color = "black"),
                   axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_text(size = fonts$axis_title),
                   legend.text = ggplot2::element_text(size = fonts$legend_text), legend.title = ggplot2::element_text(size = fonts$legend_title),
                   plot.title = ggplot2::element_text(size = fonts$axis_title + 2, face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(size = fonts$legend_text - 2, hjust = 0.5))
}
save_plot <- function(p, name, w, h) {
  path <- file.path(opt$outdir, paste0(opt$prefix, "_", gsub("[^A-Za-z0-9_]+", "_", name), ".", opt$fig_format))
  ggplot2::ggsave(path, p, device = opt$fig_format, width = w, height = h, units = "in",
                  dpi = if (opt$mode == "presentation") 150 else 300, bg = "white")
  L("Saved ", basename(path))
}
bg_note <- if (is.null(bg_symbols)) "Background: genome (no background)" else
  paste0("Background: ", length(bg_symbols), " genes (", if (assignment == "great") "GREAT genes of tested regions" else "background gene list", ")")
plotted <- 0
if (opt$mode == "presentation") {
  comb <- dplyr::bind_rows(cp_plot, enr_plot[grepl("GWAS", names(enr_plot))])
  if (!nrow(comb)) comb <- dplyr::bind_rows(enr_plot)
  if (nrow(comb)) { save_plot(make_barplot(comb, subtitle = bg_note), "overview", 11, max(4, nrow(comb) * 0.5)); plotted <- 1 }
} else {
  plot_group <- function(df, name, title) if (!is.null(df) && nrow(df)) {
    save_plot(make_barplot(df, title, bg_note), name, 9, max(3, nrow(df) * 0.35 + 1)); plotted <<- plotted + 1 }
  plot_group(cp_plot[grepl("^GO ", cp_plot$Database), ], "GO", "Gene Ontology (clusterProfiler)")
  for (cn in setdiff(unique(cp_plot$Database), grep("^GO ", unique(cp_plot$Database), value = TRUE)))
    plot_group(cp_plot[cp_plot$Database == cn, ], cn, paste0(cn, " (clusterProfiler)"))
  for (nm in names(enr_plot)) plot_group(enr_plot[[nm]], paste0("Enrichr_", nm), paste0(nm, " (Enrichr)"))
  for (cn in names(great_plot)) plot_group(great_plot[[cn]], paste0("GREAT_", cn), paste0(cn, " (GREAT, region-based)"))
}
if (plotted == 0) L("No terms passed --pval_cutoff ", opt$pval_cutoff, " in any database -- no plots (full results still written).")

write_lines_safe(c(paste0(names(opt), "\t", vapply(opt, function(x) if (is.null(x)) "" else paste(x, collapse = ","), "")),
                   paste0("gene_assignment_used\t", assignment),
                   paste0("n_dmr_genes\t", length(fg_symbols)),
                   paste0("n_background_genes\t", if (is.null(bg_symbols)) "genome" else length(bg_symbols)),
                   if (!is.null(gs)) paste0("gene_sets_", names(gs$meta), "\t", vapply(gs$meta, function(x) paste(x, collapse = ","), "")) else character(0),
                   if (!is.null(extended_tss)) paste0("great_", names(attr(extended_tss, "great_params")), "\t",
                                                      unlist(attr(extended_tss, "great_params"))) else character(0)),
                 file.path(opt$outdir, paste0(opt$prefix, "_run_parameters.txt")))
L("Done.")
