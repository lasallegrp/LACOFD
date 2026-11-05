#=====================================================================================
#  helper_overlap.R
#  Shared functions for 03_dmrichr_overlap.R and 04_dmrichr_vs_comethyl.R
#  (merge into the project's helper.R; nothing here is script-specific)
#=====================================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(openxlsx)
  library(ggplot2)
})

# ------------------------------------------------------------------------------------
# Run metadata: parameters, session info, git commit, md5 of inputs
# ------------------------------------------------------------------------------------
log_run_metadata <- function(out_dir, opt, input_files, script_name) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  opt <- opt[setdiff(names(opt), "help")]
  par <- data.frame(parameter = names(opt),
                    value = vapply(opt, function(x) paste(x, collapse = ","), ""),
                    stringsAsFactors = FALSE)
  git <- tryCatch(suppressWarnings(system2("git", c("-C", dirname(script_name), "rev-parse", "HEAD"),
                          stdout = TRUE, stderr = FALSE)), error = function(e) NA_character_)
  if (!is.null(attr(git, "status"))) git <- NA_character_
  if (!length(git)) git <- NA_character_
  par <- rbind(par, data.frame(parameter = c("script", "git_commit", "run_time"),
                               value = c(script_name, git[1], format(Sys.time()))))
  write.table(par, file.path(out_dir, "run_parameters.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  input_files <- unique(input_files[file.exists(input_files)])
  md5 <- data.frame(file = input_files, md5 = unname(tools::md5sum(input_files)))
  write.table(md5, file.path(out_dir, "input_md5.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
  invisible(par)
}

script_path <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) normalizePath(sub("^--file=", "", a[1])) else NA_character_
}

# ------------------------------------------------------------------------------------
# Generic readers
# ------------------------------------------------------------------------------------
read_table_any <- function(path, sheet = 1) {
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    read.xlsx(path, sheet = sheet)
  } else if (ext %in% c("csv")) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, sep = "\t", check.names = FALSE, data.table = FALSE))
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# ------------------------------------------------------------------------------------
# DMRichR sets from config/dmrichr_sets.tsv (set_prefix, contrast_dir, trait)
#   path = <dmrichr_root>/<contrast_dir>/cytosine_reports/<sex>/<cutoff>
#   set_id = <set_prefix>_<M|F|C>; missing DMR files are reported and dropped
# ------------------------------------------------------------------------------------
dmr_sets_from_config <- function(sets_path, dmrichr_root, cutoff, sexes, sex_labels, dmr_file) {
  cfg <- read.delim(sets_path, stringsAsFactors = FALSE, comment.char = "#")
  miss <- setdiff(c("set_prefix", "contrast_dir"), names(cfg))
  if (length(miss)) stop(sets_path, " missing columns: ", paste(miss, collapse = ", "))
  abbr <- setNames(names(sex_labels), sex_labels)        # males -> M, ...
  if (any(!sexes %in% names(abbr))) stop("sexes must be among: ", paste(sex_labels, collapse = ","))
  man <- do.call(rbind, lapply(seq_len(nrow(cfg)), function(i) data.frame(
    set_id = paste0(cfg$set_prefix[i], "_", abbr[sexes]), method = "dmrichr",
    contrast = cfg$contrast_dir[i], sex = sexes, cutoff = cutoff,
    trait = if ("trait" %in% names(cfg)) cfg$trait[i] else NA,
    path = file.path(dmrichr_root, cfg$contrast_dir[i], "cytosine_reports", sexes, cutoff),
    stringsAsFactors = FALSE)))
  ok <- file.exists(file.path(man$path, dmr_file))
  if (any(!ok)) message("MISSING DMR files, skipped:\n  ",
                        paste(file.path(man$path, dmr_file)[!ok], collapse = "\n  "))
  man <- man[ok, ]
  if (!nrow(man)) stop("No DMR files found under ", dmrichr_root, " for ", cutoff)
  rownames(man) <- NULL
  man
}

# ------------------------------------------------------------------------------------
# Comethyl networks from config/comethyl_networks.tsv
#   network_id, network_type (combined|stratified|consensus), sex, annotated_regions, me_trait
#   Paths are relative to comethyl_root and may contain {variant}.
#   variant: a full folder name (v1_all_pcs) or shorthand V1 / V2 (matched as v1_* / v2_*).
# ------------------------------------------------------------------------------------
resolve_variant_path <- function(template, root, variant) {
  p <- file.path(root, template)
  if (!grepl("{variant}", p, fixed = TRUE)) return(p)
  if (grepl("^[vV][0-9]+$", variant)) {
    hits <- Sys.glob(sub("{variant}", paste0(tolower(variant), "_*"), p, fixed = TRUE))
    if (length(hits) == 1) return(hits)
    if (length(hits) > 1) stop("Variant ", variant, " is ambiguous for ", template, ": ",
                               paste(basename(dirname(hits)), collapse = ", "))
    return(sub("{variant}", paste0(tolower(variant), "_<not found>"), p, fixed = TRUE))
  }
  sub("{variant}", variant, p, fixed = TRUE)
}

networks_from_config <- function(net_path, comethyl_root, variant, types = NULL) {
  nc <- read.delim(net_path, stringsAsFactors = FALSE, comment.char = "#")
  miss <- setdiff(c("network_id", "network_type", "sex", "annotated_regions", "me_trait"), names(nc))
  if (length(miss)) stop(net_path, " missing columns: ", paste(miss, collapse = ", "))
  if (!is.null(types)) nc <- nc[nc$network_type %in% types, ]
  nc$annotated_regions <- vapply(nc$annotated_regions, resolve_variant_path, "", comethyl_root, variant)
  nc$me_trait <- vapply(nc$me_trait, resolve_variant_path, "", comethyl_root, variant)
  nc$variant <- basename(dirname(nc$annotated_regions))
  ok <- file.exists(nc$annotated_regions) & file.exists(nc$me_trait)
  if (any(!ok)) message("MISSING comethyl files, network skipped: ",
                        paste(nc$network_id[!ok], collapse = ", "), "\n  ",
                        paste(unique(c(nc$annotated_regions[!ok & !file.exists(nc$annotated_regions)],
                                       nc$me_trait[!ok & !file.exists(nc$me_trait)])), collapse = "\n  "))
  nc <- nc[ok, ]
  if (!nrow(nc)) stop("No comethyl networks found for variant ", variant)
  rownames(nc) <- NULL
  nc
}

read_manifest <- function(path, required) {
  m <- read.delim(path, stringsAsFactors = FALSE, comment.char = "#")
  miss <- setdiff(required, names(m))
  if (length(miss)) stop("Manifest ", path, " missing columns: ", paste(miss, collapse = ", "))
  if (anyDuplicated(m$set_id)) stop("Duplicate set_id in manifest: ",
                                    paste(m$set_id[duplicated(m$set_id)], collapse = ", "))
  m
}

# Check a manifest label against the path it came from (labels are path-derived
# in this project; a mismatch means the manifest row is wrong)
check_path_label <- function(m, col, path_col = "path") {
  bad <- !mapply(function(lbl, p) grepl(paste0("/", lbl, "(/|$)"), p), m[[col]], m[[path_col]])
  if (any(bad)) warning("Manifest '", col, "' not found as a path component for: ",
                        paste(m$set_id[bad], collapse = ", "), call. = FALSE)
  invisible(!bad)
}

detect_col <- function(df, candidates, what, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) {
    if (required) stop("Could not find a ", what, " column (tried: ",
                       paste(candidates, collapse = ", "), ")")
    return(NA_character_)
  }
  hit[1]
}

# ------------------------------------------------------------------------------------
# DMRichR DMR set -> standardized data.frame
#   adds: .chr .start .end .effect .dir (+1/-1) .gene .set_id .dmr_id
# ------------------------------------------------------------------------------------
read_dmr_set <- function(path, set_id, effect_col = "auto", gene_col = "geneSymbol",
                         gene_sep = NA, exclude_sex_chr = FALSE,
                         max_pval = NA, max_qval = NA) {
  df <- read_table_any(path)
  chr_c <- detect_col(df, c("seqnames", "chr", "chrom"), "chromosome")
  if (effect_col == "auto")
    effect_col <- detect_col(df, c("difference", "stat", "beta"), "effect/direction")
  if (!effect_col %in% names(df)) stop("effect_col '", effect_col, "' not in ", path)
  if (!gene_col %in% names(df)) {
    warning(set_id, ": gene column '", gene_col, "' missing; gene mode will be empty", call. = FALSE)
    df[[gene_col]] <- NA_character_
  }
  if (!is.na(max_pval) && "pval" %in% names(df)) df <- df[df$pval <= max_pval, , drop = FALSE]
  if (!is.na(max_qval) && "qval" %in% names(df)) df <- df[df$qval <= max_qval, , drop = FALSE]
  df$.chr    <- as.character(df[[chr_c]])
  df$.start  <- as.integer(df$start)
  df$.end    <- as.integer(df$end)
  df$.effect <- as.numeric(df[[effect_col]])
  df$.dir    <- sign(df$.effect)
  g <- as.character(df[[gene_col]])
  g[g %in% c("", "NA")] <- NA
  df$.gene   <- g
  if (exclude_sex_chr) df <- df[!df$.chr %in% c("chrX", "chrY", "X", "Y"), , drop = FALSE]
  df$.set_id <- set_id
  df$.dmr_id <- sprintf("%s_DMR%05d", set_id, seq_len(nrow(df)))
  attr(df, "effect_col") <- effect_col
  attr(df, "gene_sep") <- gene_sep
  df
}

dmr_genes <- function(df) {
  sep <- attr(df, "gene_sep")
  g <- df$.gene[!is.na(df$.gene)]
  if (!is.na(sep)) g <- trimws(unlist(strsplit(g, sep, fixed = TRUE)))
  sort(unique(g[g != ""]))
}

to_gr <- function(df, chr = ".chr", start = ".start", end = ".end") {
  GRanges(df[[chr]], IRanges(as.integer(df[[start]]), as.integer(df[[end]])))
}

# Background (all tested candidate regions). Coordinates required; an effect
# column (difference/stat/beta) and pval are used if present.
read_background <- function(path, exclude_sex_chr = FALSE) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  df <- read_table_any(path)
  if (!"start" %in% names(df) && ncol(df) >= 3) names(df)[1:3] <- c("chr", "start", "end")
  chr_c <- detect_col(df, c("seqnames", "chr", "chrom", "V1"), "chromosome")
  if (exclude_sex_chr) df <- df[!df[[chr_c]] %in% c("chrX", "chrY", "X", "Y"), , drop = FALSE]
  gr <- GRanges(df[[chr_c]], IRanges(as.integer(df$start), as.integer(df$end)))
  eff <- detect_col(df, c("difference", "stat", "beta"), "effect", required = FALSE)
  if (!is.na(eff)) gr$effect <- as.numeric(df[[eff]])
  if ("pval" %in% names(df)) gr$pval <- as.numeric(df$pval)
  gc <- detect_col(df, c("geneSymbol", "gene_symbol", "SYMBOL"), "gene", required = FALSE)
  if (!is.na(gc)) gr$gene <- as.character(df[[gc]])
  gr
}

# ------------------------------------------------------------------------------------
# Overlap criteria. GenomicRanges cannot combine maxgap and minoverlap, so:
#   min_overlap_bp > 1  -> require that many bp of overlap
#   else max_gap > 0    -> count DMRs within max_gap bp as overlapping
#   else                -> direct overlap (>= 1 bp)
# ------------------------------------------------------------------------------------
ov_hits <- function(q, s, max_gap = 0, min_overlap_bp = 1) {
  if (min_overlap_bp > 1) {
    findOverlaps(q, s, minoverlap = as.integer(min_overlap_bp))
  } else if (max_gap > 0) {
    findOverlaps(q, s, maxgap = as.integer(max_gap))
  } else {
    findOverlaps(q, s)
  }
}

# Connected components of an overlap graph (simple union-find, no igraph needed)
components_from_hits <- function(n, qh, sh) {
  parent <- seq_len(n)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (k in seq_along(qh)) {
    a <- find(qh[k]); b <- find(sh[k])
    if (a != b) parent[max(a, b)] <- min(a, b)
  }
  vapply(seq_len(n), find, 1L)
}

# ------------------------------------------------------------------------------------
# Build overlap loci across sets (coordinate mode)
#   Every DMR is a node; DMRs that overlap (under the criteria) are joined.
#   A locus is a connected component. Returns per-locus membership + direction.
# ------------------------------------------------------------------------------------
build_loci <- function(dmr_list, max_gap = 0, min_overlap_bp = 1) {
  all <- do.call(rbind, lapply(dmr_list, function(d)
    d[, c(".dmr_id", ".set_id", ".chr", ".start", ".end", ".effect", ".dir", ".gene")]))
  gr <- to_gr(all)
  h <- ov_hits(gr, gr, max_gap, min_overlap_bp)
  comp <- components_from_hits(length(gr), queryHits(h), subjectHits(h))
  ids <- match(comp, unique(comp[order(as.character(seqnames(gr)), start(gr))]))
  all$locus_id <- sprintf("L%06d", ids)

  sets <- names(dmr_list)
  loci <- do.call(rbind, lapply(split(all, all$locus_id), function(x) {
    data.frame(locus_id = x$locus_id[1], chr = x$.chr[1],
               start = min(x$.start), end = max(x$.end),
               n_dmrs = nrow(x),
               genes = paste(sort(unique(na.omit(x$.gene))), collapse = ";"),
               stringsAsFactors = FALSE)
  }))
  mem <- sapply(sets, function(s) loci$locus_id %in% all$locus_id[all$.set_id == s])
  if (is.null(dim(mem))) mem <- matrix(mem, nrow = 1, dimnames = list(NULL, sets))
  dirm <- sapply(sets, function(s) {
    x <- all[all$.set_id == s, ]
    d <- tapply(x$.dir, x$locus_id, function(v) if (all(v > 0)) 1 else if (all(v < 0)) -1 else 0)
    out <- rep(NA_real_, nrow(loci)); names(out) <- loci$locus_id
    out[names(d)] <- d
    out
  })
  if (is.null(dim(dirm))) dirm <- matrix(dirm, nrow = 1, dimnames = list(NULL, sets))
  colnames(dirm) <- sets
  loci$n_sets <- rowSums(mem)
  loci$pattern <- apply(mem, 1, function(r) paste(sets[r], collapse = "&"))
  loci$direction_concordant <- apply(dirm, 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < 2) return(NA)
    all(v != 0) && length(unique(v)) == 1
  })
  dir_lab <- matrix(c("hypo", "mixed", "hyper")[dirm + 2], nrow = nrow(dirm),
                    dimnames = list(NULL, paste0("dir_", sets)))
  loci <- cbind(loci, as.data.frame(mem) * 1L, dir_lab, stringsAsFactors = FALSE)
  rownames(loci) <- NULL
  list(loci = loci, dmr_map = all, sets = sets)
}

# Gene-mode equivalent of build_loci
build_gene_table <- function(dmr_list) {
  sets <- names(dmr_list)
  gl <- lapply(dmr_list, dmr_genes)
  genes <- sort(unique(unlist(gl)))
  if (!length(genes)) return(NULL)
  mem <- sapply(sets, function(s) genes %in% gl[[s]])
  if (is.null(dim(mem))) mem <- matrix(mem, nrow = 1, dimnames = list(NULL, sets))
  dirm <- sapply(sets, function(s) {
    d <- dmr_list[[s]]; d <- d[!is.na(d$.gene), ]
    v <- tapply(d$.dir, d$.gene, function(v) if (all(v > 0)) 1 else if (all(v < 0)) -1 else 0)
    out <- rep(NA_real_, length(genes)); names(out) <- genes
    out[intersect(names(v), genes)] <- v[intersect(names(v), genes)]
    out
  })
  if (is.null(dim(dirm))) dirm <- matrix(dirm, nrow = 1, dimnames = list(NULL, sets))
  tab <- data.frame(geneSymbol = genes, n_sets = rowSums(mem),
                    pattern = apply(mem, 1, function(r) paste(sets[r], collapse = "&")),
                    stringsAsFactors = FALSE)
  tab$direction_concordant <- apply(dirm, 1, function(v) {
    v <- v[!is.na(v)]; if (length(v) < 2) return(NA); all(v != 0) && length(unique(v)) == 1
  })
  dir_lab <- matrix(c("hypo", "mixed", "hyper")[dirm + 2], nrow = nrow(dirm),
                    dimnames = list(NULL, paste0("dir_", sets)))
  cbind(tab, as.data.frame(mem) * 1L, dir_lab, stringsAsFactors = FALSE)
}

# Intersection counts in UpSet style (exclusive patterns)
pattern_counts <- function(tab, sets, require_direction = FALSE) {
  t2 <- tab
  if (require_direction) t2 <- t2[t2$n_sets == 1 | t2$direction_concordant %in% TRUE, ]
  cnt <- as.data.frame(table(pattern = t2$pattern), stringsAsFactors = FALSE)
  names(cnt)[2] <- "n"
  conc <- tapply(tab$direction_concordant, tab$pattern, function(v) sum(v %in% TRUE))
  disc <- tapply(tab$direction_concordant, tab$pattern, function(v) sum(v %in% FALSE))
  cnt$n_direction_concordant <- ifelse(grepl("&", cnt$pattern), conc[cnt$pattern], NA)
  cnt$n_direction_discordant <- ifelse(grepl("&", cnt$pattern), disc[cnt$pattern], NA)
  cnt$n_sets <- lengths(strsplit(cnt$pattern, "&"))
  cnt[order(-cnt$n_sets, -cnt$n), ]
}

# ------------------------------------------------------------------------------------
# Pairwise statistics
# ------------------------------------------------------------------------------------
# Region-mode Fisher: universe = background regions tested in BOTH sets
region_fisher <- function(a, b, bg_a, bg_b, max_gap = 0, min_overlap_bp = 1) {
  if (is.null(bg_a) || is.null(bg_b)) return(c(OR = NA, p = NA, universe = NA))
  u <- reduce(c(granges(bg_a)[overlapsAny(bg_a, bg_b)], granges(bg_b)[overlapsAny(bg_b, bg_a)]))
  ina <- overlapsAny(u, a); inb <- overlapsAny(u, b)
  tb <- table(factor(ina, c(TRUE, FALSE)), factor(inb, c(TRUE, FALSE)))
  ft <- fisher.test(tb, alternative = "greater")
  c(OR = unname(ft$estimate), p = ft$p.value, universe = length(u))
}

gene_fisher <- function(ga, gb, universe) {
  if (is.null(universe)) return(c(OR = NA, p = NA, universe = NA))
  ga <- intersect(ga, universe); gb <- intersect(gb, universe)
  tb <- table(factor(universe %in% ga, c(TRUE, FALSE)), factor(universe %in% gb, c(TRUE, FALSE)))
  ft <- fisher.test(tb, alternative = "greater")
  c(OR = unname(ft$estimate), p = ft$p.value, universe = length(universe))
}

pairwise_coord <- function(dmr_list, bg_list, loci, max_gap, min_overlap_bp) {
  sets <- names(dmr_list); out <- list()
  for (i in seq_along(sets)) for (j in seq_along(sets)) if (i < j) {
    A <- dmr_list[[sets[i]]]; B <- dmr_list[[sets[j]]]
    ga <- to_gr(A); gb <- to_gr(B)
    h <- ov_hits(ga, gb, max_gap, min_overlap_bp)
    conc <- if (length(h)) mean(A$.dir[queryHits(h)] == B$.dir[subjectHits(h)]) else NA
    shared <- sum(loci[[sets[i]]] == 1 & loci[[sets[j]]] == 1)
    either <- sum(loci[[sets[i]]] == 1 | loci[[sets[j]]] == 1)
    f <- region_fisher(ga, gb, bg_list[[sets[i]]], bg_list[[sets[j]]], max_gap, min_overlap_bp)
    out[[length(out) + 1]] <- data.frame(
      set_A = sets[i], set_B = sets[j], n_A = nrow(A), n_B = nrow(B),
      n_A_overlapping_B = length(unique(queryHits(h))),
      n_B_overlapping_A = length(unique(subjectHits(h))),
      shared_loci = shared, jaccard_loci = ifelse(either > 0, shared / either, NA),
      pct_direction_concordant = round(100 * conc, 1),
      fisher_OR = f[["OR"]], fisher_p = f[["p"]], universe_regions = f[["universe"]])
  }
  res <- do.call(rbind, out)
  res$fisher_FDR <- p.adjust(res$fisher_p, "BH")
  res
}

pairwise_gene <- function(dmr_list, gene_universe, coord_loci = NULL) {
  sets <- names(dmr_list); gl <- lapply(dmr_list, dmr_genes); out <- list()
  for (i in seq_along(sets)) for (j in seq_along(sets)) if (i < j) {
    a <- gl[[sets[i]]]; b <- gl[[sets[j]]]; s <- intersect(a, b)
    f <- gene_fisher(a, b, gene_universe)
    # how many shared genes are NOT backed by a shared coordinate locus?
    gene_only <- NA
    if (!is.null(coord_loci) && length(s)) {
      both <- coord_loci[coord_loci[[sets[i]]] == 1 & coord_loci[[sets[j]]] == 1, ]
      lg <- unique(unlist(strsplit(both$genes, ";")))
      gene_only <- sum(!s %in% lg)
    }
    out[[length(out) + 1]] <- data.frame(
      set_A = sets[i], set_B = sets[j], n_genes_A = length(a), n_genes_B = length(b),
      shared_genes = length(s), jaccard = length(s) / length(union(a, b)),
      shared_genes_without_shared_locus = gene_only,
      fisher_OR = f[["OR"]], fisher_p = f[["p"]], universe_genes = f[["universe"]])
  }
  res <- do.call(rbind, out)
  res$fisher_FDR <- p.adjust(res$fisher_p, "BH")
  res
}

# ------------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------------
plot_upset <- function(tab, sets, file, title, require_direction = FALSE, set_cols = NULL) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    message("UpSetR not installed; skipping ", basename(file)); return(invisible(NULL))
  }
  t2 <- tab
  if (require_direction) t2 <- t2[t2$n_sets == 1 | t2$direction_concordant %in% TRUE, ]
  m <- as.data.frame(t2[, sets, drop = FALSE])
  if (nrow(m) == 0 || ncol(m) < 2) return(invisible(NULL))
  pdf(file, width = 8, height = 5, onefile = FALSE)
  print(UpSetR::upset(m, sets = rev(sets), keep.order = TRUE, order.by = "freq",
                      nintersects = NA, mainbar.y.label = title,
                      sets.x.label = "Total per set", text.scale = 1.2,
                      sets.bar.color = if (is.null(set_cols)) "gray23" else unname(set_cols[rev(sets)]),
                      main.bar.color = "grey25", matrix.color = "grey20"))
  dev.off()
}

plot_pairwise_heatmap <- function(pw, value, file, title, sets) {
  d <- rbind(pw[, c("set_A", "set_B", value)],
             setNames(pw[, c("set_B", "set_A", value)], c("set_A", "set_B", value)))
  d$set_A <- factor(d$set_A, sets); d$set_B <- factor(d$set_B, rev(sets))
  lab <- if (value == "jaccard_loci" || value == "jaccard") sprintf("%.3f", d[[value]]) else
    formatC(d[[value]], digits = 2, format = "g")
  if ("fisher_FDR" %in% names(pw)) {
    fdr <- c(pw$fisher_FDR, pw$fisher_FDR)
    lab <- paste0(lab, ifelse(!is.na(fdr) & fdr < 0.05, "*", ""))
  }
  p <- ggplot(d, aes(set_A, set_B, fill = .data[[value]])) + geom_tile(colour = "white") +
    geom_text(aes(label = lab), size = 3.5) +
    scale_fill_gradient(low = "#F2F2F2", high = "#2C7BB6", na.value = "white") +
    labs(x = NULL, y = NULL, fill = value, title = title,
         caption = "* Fisher FDR < 0.05 (one-sided, enrichment)") +
    theme_minimal(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file, p, width = 6, height = 5)
}

# ------------------------------------------------------------------------------------
# Excel writing with safe (<=31 char, unique) sheet names
# ------------------------------------------------------------------------------------
write_xlsx_safe <- function(lst, file) {
  lst <- lst[!vapply(lst, is.null, TRUE)]
  full <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", names(lst))
  # Excel: <= 31 chars, unique case-insensitively. Truncate, then on a clash replace the
  # tail with ~k (so the suffix is never cut off again).
  nm <- character(0)
  for (x in full) {
    cand <- substr(x, 1, 31); k <- 1
    while (tolower(cand) %in% tolower(nm)) {
      suf <- paste0("~", k)
      cand <- paste0(substr(x, 1, 31 - nchar(suf)), suf); k <- k + 1
    }
    nm <- c(nm, cand)
  }
  names(lst) <- nm
  lst <- lapply(lst, function(x) if (nrow(x) == 0) data.frame(note = "no rows") else x)
  # if any name was shortened, add a lookup sheet: sheet -> full name
  if (any(nm != full)) lst <- c(lst, list(Sheet_index = data.frame(sheet = nm, full_name = full)))
  write.xlsx(lst, file, asTable = FALSE, overwrite = TRUE, firstRow = TRUE)
}

# Bind data frames with differing columns (from the original 03 script)
bind_rows_any <- function(dflist) {
  dflist <- dflist[vapply(dflist, function(d) !is.null(d) && nrow(d) > 0, TRUE)]
  if (!length(dflist)) return(data.frame())
  all_cols <- Reduce(union, lapply(dflist, names))
  do.call(rbind, lapply(dflist, function(df) {
    miss <- setdiff(all_cols, names(df)); if (length(miss)) df[miss] <- NA; df[all_cols]
  }))
}

# ------------------------------------------------------------------------------------
# UpSet plot with gene names per intersection (ggplot2 + grid only)
#   tab:        loci table (coord) or gene table (gene mode), with set 0/1 columns + pattern
#   gene_col:   "genes" (";"-joined, coord) or "geneSymbol" (gene mode)
#   score:      numeric per row used to rank genes (e.g. max |difference|); higher = shown first
#   n_genes:    max names per bar; scope "multi" labels only overlaps, "all" labels every bar
#   highlight:  optional character vector; only these genes are eligible for labels
# ------------------------------------------------------------------------------------
upset_gene_labels <- function(tab, gene_col, score, n_genes, scope = "multi", highlight = NULL) {
  pats <- unique(tab$pattern)
  setNames(vapply(pats, function(p) {
    i <- which(tab$pattern == p)
    if (scope == "multi" && !grepl("&", p)) return("")
    g <- strsplit(ifelse(is.na(tab[[gene_col]][i]), "", tab[[gene_col]][i]), ";", fixed = TRUE)
    df <- data.frame(gene = unlist(g), score = rep(score[i], lengths(g)), stringsAsFactors = FALSE)
    df <- df[df$gene != "", ]
    if (!is.null(highlight)) df <- df[df$gene %in% highlight, ]
    if (!nrow(df)) return("")
    s <- tapply(df$score, df$gene, max, na.rm = TRUE)
    s <- sort(s, decreasing = TRUE)
    top <- names(s)[seq_len(min(n_genes, length(s)))]
    more <- length(s) - length(top)
    paste(c(top, if (more > 0) paste0("+", more, " more")), collapse = "\n")
  }, ""), pats)
}

# Consistent colours: within_conditions_across_sex by sex, within_sex_across_conditions (and custom) by contrast
set_colours <- function(man_g, family) {
  sex_pal <- c(males = "#2C7BB6", females = "#D7191C", sex_combined = "#4D4D4D")
  okabe <- c("#009E73", "#CC79A7", "#56B4E9", "#0072B2", "#F0E442", "#000000")
  if (family == "within_conditions_across_sex" && all(man_g$sex %in% names(sex_pal))) {
    cols <- sex_pal[man_g$sex]
  } else {
    lv <- unique(man_g$contrast)
    cols <- okabe[(match(man_g$contrast, lv) - 1) %% length(okabe) + 1]
  }
  setNames(unname(cols), man_g$set_id)
}

# methylation direction colours (match RdBu heatmaps: red = hyper, blue = hypo)
dir_pal <- c(Hyper = "#B2182B", Hypo = "#2166AC", `Discordant / mixed` = "#BDBDBD")

# direction class per row: same direction in every set it appears in -> Hyper/Hypo
direction_class <- function(tab, sets) {
  d <- as.matrix(tab[, paste0("dir_", sets), drop = FALSE])
  apply(d, 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) && all(v == "hyper")) "Hyper" else if (length(v) && all(v == "hypo")) "Hypo" else "Discordant / mixed"
  })
}

plot_upset_labeled <- function(tab, sets, file, title, gene_col, score, n_genes = 5,
                               scope = "multi", highlight = NULL, require_direction = FALSE,
                               set_sizes = NULL, set_cols = NULL) {
  t2 <- tab; sc <- score
  if (require_direction) {
    k <- t2$n_sets == 1 | t2$direction_concordant %in% TRUE
    t2 <- t2[k, ]; sc <- sc[k]
  }
  if (!nrow(t2)) return(invisible(NULL))
  if (is.null(set_cols)) set_cols <- setNames(rep("grey15", length(sets)), sets)
  t2$dir_class <- factor(direction_class(t2, sets), names(dir_pal))
  cnt <- as.data.frame(table(pattern = t2$pattern), stringsAsFactors = FALSE)
  cnt <- cnt[order(-cnt$Freq), ]
  lv <- cnt$pattern
  lab <- if (n_genes > 0) upset_gene_labels(t2, gene_col, sc, n_genes, scope, highlight) else
    setNames(rep("", length(lv)), lv)
  cnt$lab <- lab[cnt$pattern]
  cnt$pattern <- factor(cnt$pattern, lv)
  stk <- as.data.frame(table(pattern = t2$pattern, dir_class = t2$dir_class), stringsAsFactors = FALSE)
  stk <- stk[stk$Freq > 0, ]
  stk$pattern <- factor(stk$pattern, lv)
  stk$dir_class <- factor(stk$dir_class, rev(names(dir_pal)))
  n_lines <- max(1, vapply(strsplit(cnt$lab, "\n"), length, 1L))

  ylab_sets <- if (is.null(set_sizes)) sets else paste0(sets, " (", set_sizes[sets], ")")
  mat <- expand.grid(set = sets, pattern = lv, stringsAsFactors = FALSE)
  mat$on <- mapply(function(s, p) s %in% strsplit(p, "&")[[1]], mat$set, mat$pattern)
  mat$y <- match(mat$set, rev(sets))
  mat$pattern <- factor(mat$pattern, lv)
  mat$col <- ifelse(mat$on, set_cols[mat$set], "grey88")
  seg <- do.call(rbind, lapply(lv, function(p) {
    y <- mat$y[mat$pattern == p & mat$on]
    if (length(y) < 2) return(NULL)
    data.frame(pattern = factor(p, lv), ymin = min(y), ymax = max(y))
  }))

  th <- theme_minimal(11) + theme(panel.grid.major.x = element_blank(),
                                  panel.grid.minor = element_blank())
  sx <- scale_x_discrete(limits = lv, drop = FALSE)
  p_bar <- ggplot(stk, aes(pattern, Freq, fill = dir_class)) + geom_col(width = .7) +
    geom_text(data = cnt, aes(pattern, Freq, label = Freq), inherit.aes = FALSE, vjust = -0.4, size = 3.4) +
    sx + scale_fill_manual(values = dir_pal, breaks = names(dir_pal), name = "Methylation direction") +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    labs(x = NULL, y = "Intersection size", title = title) +
    th + theme(axis.text.x = element_blank(), legend.position = c(.98, .95), legend.justification = c(1, 1),
               legend.background = element_rect(fill = "white", colour = NA))
  p_lab <- ggplot(cnt, aes(pattern, 1)) +
    geom_text(aes(label = lab), vjust = 1, size = 2.7, fontface = "italic", lineheight = .9) +
    sx + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) + th +
    theme(axis.text = element_blank(), panel.grid = element_blank(), legend.position = "none")
  p_mat <- ggplot(mat, aes(pattern, y)) +
    { if (!is.null(seg)) geom_segment(data = seg, aes(x = pattern, xend = pattern, y = ymin, yend = ymax),
                                      inherit.aes = FALSE, linewidth = 1, colour = "grey40") } +
    geom_point(colour = mat$col, size = 3.8) + sx +
    scale_y_continuous(breaks = seq_along(sets), labels = rev(ylab_sets),
                       limits = c(.5, length(sets) + .5)) +
    labs(x = NULL, y = NULL) + th +
    theme(axis.text.x = element_blank(), panel.grid = element_blank(), legend.position = "none",
          axis.text.y = element_text(face = "bold"))

  g <- rbind(ggplotGrob(p_bar), ggplotGrob(p_lab), ggplotGrob(p_mat), size = "max")
  panels <- grep("^panel", g$layout$name)
  rows <- sort(unique(g$layout$t[panels]))
  h_lab <- if (n_genes > 0) max(.6, .16 * n_lines) else .1
  g$heights[rows] <- unit(c(3, h_lab, .35 * length(sets)), c("in", "in", "in"))
  w <- max(6.5, 1.0 * length(lv) + 2.5)
  pdf(file, width = w, height = 3 + h_lab + .35 * length(sets) + 1.4)
  grid::grid.newpage(); grid::grid.draw(g)
  dev.off()
}
# ------------------------------------------------------------------------------------
# regioneR permutation test of region overlap (used by 03_dmrichr_overlap.R --perm_test)
#
#   Statistic = number of A DMRs overlapping >= 1 B DMR, using the same overlap rule as
#   the rest of the script (ov_hits: max_gap / min_overlap_bp) = n_A_overlapping_B.
#   B is held fixed; both directions are run because set sizes differ.
#   A and B are restricted to the space tested in both sets (backgrounds overlap).
#
#   randomize = "resample"  (default) regioneR::resampleRegions
#       random A = |A| regions drawn from A's own DMRichR candidate regions that were
#       also tested in B. Like-for-like null (real widths, CpG-testable loci).
#   randomize = "randomize" regioneR::randomizeRegions
#       random A = A's regions placed at random positions in the genome (optionally
#       masked). Ignores where DMRs could be called -> usually anti-conservative.
# ------------------------------------------------------------------------------------
read_perm_genome <- function(genome) {
  # "hg38" etc. is passed to regioneR (needs the matching BSgenome package);
  # a file is read as chrom.sizes (chr <tab> length) -> GRanges
  if (file.exists(genome)) {
    cs <- read.delim(genome, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
    return(GRanges(cs[[1]], IRanges(1L, as.integer(cs[[2]]))))
  }
  genome
}

perm_overlap_pair <- function(A, B, bg_A, bg_B, ntimes = 100, max_gap = 0, min_overlap_bp = 1,
                              randomize = "resample", genome = "hg38", mask = NA,
                              per_chr = TRUE, alternative = "greater", local_z = TRUE,
                              lz_window = NA, lz_step = NA, seed = 1) {
  if (is.null(bg_A) || is.null(bg_B)) return(NULL)
  uA <- granges(bg_A)[overlapsAny(bg_A, bg_B)]            # A candidates tested in both
  A  <- granges(A)[overlapsAny(A, uA)]
  B  <- granges(B)[overlapsAny(B, granges(bg_A))]
  if (!length(A) || !length(B)) return(NULL)
  ev <- list(numOverlaps = function(A, B, ...)
    length(unique(queryHits(ov_hits(A, B, max_gap, min_overlap_bp)))))
  if (randomize == "resample") {
    if (length(uA) <= length(A)) return(NULL)
    if (per_chr) {                     # per-chromosome resampling needs room on every chr
      nA <- table(as.character(seqnames(A))); nU <- table(as.character(seqnames(uA)))
      if (any(nU[names(nA)] < nA)) per_chr <- FALSE
    }
    rargs <- list(randomize.function = regioneR::resampleRegions,
                  randomize.function.name = "resampleRegions", universe = uA,
                  per.chromosome = per_chr)
  } else {
    rargs <- list(randomize.function = regioneR::randomizeRegions,
                  randomize.function.name = "randomizeRegions", genome = genome,
                  mask = mask, per.chromosome = per_chr, allow.overlaps = TRUE)
  }
  set.seed(seed)
  pt <- do.call(regioneR::permTest,
                c(list(A = A, B = B, ntimes = ntimes, evaluate.function = ev,
                       alternative = alternative, verbose = FALSE), rargs))
  r <- pt$numOverlaps
  lz <- NULL; z0 <- NA; z_edge <- NA
  if (local_z && is.finite(r$zscore)) {
    w <- if (is.na(lz_window)) max(1000, round(5 * mean(width(A)))) else lz_window
    s <- if (is.na(lz_step)) max(1, round(w / 20)) else lz_step
    lz <- tryCatch(regioneR::localZScore(pt = pt, A = A, B = B, window = w, step = s),
                   error = function(e) { message("localZScore failed: ", conditionMessage(e)); NULL })
    if (!is.null(lz)) {
      sh <- lz$numOverlaps$shifts; zz <- lz$numOverlaps$shifted.z.scores
      z0 <- zz[which.min(abs(sh))]; z_edge <- mean(zz[abs(sh) == max(abs(sh))])
    }
  }
  perm <- r$permuted
  list(pt = pt, lz = lz, permuted = perm,
       summary = data.frame(randomization = randomize, n_A_tested = length(A), n_B_tested = length(B),
                            universe_A = if (randomize == "resample") length(uA) else NA,
                            per_chromosome = per_chr,
                            observed = r$observed, perm_mean = mean(perm), perm_sd = sd(perm),
                            fold_enrichment = ifelse(mean(perm) > 0, r$observed / mean(perm), NA),
                            z = r$zscore, perm_p = r$pval, alternative = r$alternative,
                            ntimes = ntimes, localz_at_0 = z0, localz_at_edges = z_edge,
                            stringsAsFactors = FALSE))
}

# All pairs of a group, both directions. Writes regioneR's own plots per direction to
# pdf_dir; the permuted counts are returned in attr(result, "permuted") for summary plots.
perm_pairwise <- function(dmr_list, bg_list, ntimes, max_gap, min_overlap_bp, randomize, genome,
                          mask, per_chr, alternative, local_z, lz_window, lz_step, seed,
                          pdf_dir = NULL) {
  sets <- names(dmr_list); out <- list(); perms <- list()
  if (!is.null(pdf_dir)) dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(sets)) for (j in seq_along(sets)) if (i != j) {
    a <- sets[i]; b <- sets[j]
    message("  regioneR (", randomize, "): ", a, " -> ", b)
    res <- perm_overlap_pair(to_gr(dmr_list[[a]]), to_gr(dmr_list[[b]]), bg_list[[a]], bg_list[[b]],
                             ntimes, max_gap, min_overlap_bp, randomize, genome, mask, per_chr,
                             alternative, local_z, lz_window, lz_step, seed)
    if (is.null(res)) {
      message("    skipped (no background, or too few regions in the shared universe)")
      next
    }
    out[[length(out) + 1]] <- cbind(randomized_A = a, fixed_B = b, res$summary)
    perms[[length(perms) + 1]] <- data.frame(randomized_A = a, fixed_B = b, permuted = res$permuted)
    if (!is.null(pdf_dir)) {
      pdf(file.path(pdf_dir, paste0("perm_", a, "__vs__", b, ".pdf")), width = 7, height = 5)
      plot(res$pt, main = paste0(a, " vs ", b, " (", randomize, "Regions, n = ", ntimes, ")"))
      if (!is.null(res$lz)) plot(res$lz, main = paste0("Local z-score: ", a, " vs ", b))
      dev.off()
    }
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  res$perm_FDR_group <- p.adjust(res$perm_p, "BH")
  attr(res, "permuted") <- do.call(rbind, perms)
  res
}

# Summary figure: null distribution (permuted overlaps) per directed pair, observed value
# as a red diamond, annotated with fold, z, p and FDR.
#   pm: permutation table (randomized_A, fixed_B, observed, fold_enrichment, z, perm_p, <fdr_col>, group)
#   pr: permuted counts (randomized_A, fixed_B, permuted, group)
plot_perm_summary <- function(pm, pr, file, title, fdr_col = "perm_FDR_group",
                              xlab = "A DMRs overlapping B DMRs (grey = permuted null)",
                              ylab = "A (randomized) -> B (fixed)") {
  key <- function(d) paste0(d$randomized_A, " -> ", d$fixed_B)
  pm$pair <- key(pm); pr$pair <- key(pr)
  lv <- rev(unique(pm$pair[order(pm$group, pm$randomized_A, pm$fixed_B)]))
  pm$pair <- factor(pm$pair, lv); pr$pair <- factor(pr$pair, lv)
  pm$sig <- ifelse(pm[[fdr_col]] < 0.05, "FDR < 0.05", "n.s.")
  pm$lab <- sprintf("fold %.2f | z = %.1f | p = %.2g | FDR = %.2g",
                    pm$fold_enrichment, pm$z, pm$perm_p, pm[[fdr_col]])
  xmax <- tapply(c(pr$permuted, pm$observed), c(as.character(pr$group), as.character(pm$group)), max)
  pm$xlab <- xmax[as.character(pm$group)] * 1.05
  p <- ggplot(pr, aes(permuted, pair)) +
    geom_violin(fill = "grey85", colour = "grey55", scale = "width", width = .8) +
    geom_boxplot(width = .15, outlier.shape = NA, fill = "white", colour = "grey35") +
    geom_point(data = pm, aes(observed, pair, colour = sig), shape = 18, size = 4.5) +
    geom_text(data = pm, aes(xlab, pair, label = lab), hjust = 0, size = 2.9, colour = "grey20") +
    scale_colour_manual(values = c("FDR < 0.05" = "#D7191C", "n.s." = "#2C7BB6"), drop = FALSE,
                        name = "Observed overlap") +
    scale_x_continuous(expand = expansion(mult = c(.02, .9))) +
    facet_wrap(~ group, ncol = 1, scales = "free") +
    labs(x = xlab, y = ylab,
         title = title,
         caption = sprintf("%s; %d permutations per test (min p = %.3g)",
                           unique(pm$randomization)[1], unique(pm$ntimes)[1], 1 / (unique(pm$ntimes)[1] + 1))) +
    theme_bw(10) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
                         plot.caption = element_text(hjust = 0, colour = "grey30"))
  n <- length(lv); ng <- length(unique(pm$group))
  ggsave(file, p, width = 11, height = max(3.5, 0.55 * n + 1.2 * ng + 1.5), limitsize = FALSE)
  invisible(p)
}

# ====================================================================================
# GWAS gene sets and SNP sets (used by 03_dmrichr_overlap.R --gene_set / --snp_set)
#   Input: a GWAS Catalog association export (riskAllele, pValue, mappedGenes,
#   efoTraits, accessionId, locations, ...) or a plain table (gene list for
#   --gene_set; chr + pos/start for --snp_set).
# ====================================================================================
read_gwas_input <- function(path, max_p = NA) {
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(tools::file_ext(path))
  df <- if (ext %in% c("xlsx", "xls", "csv")) read_table_any(path) else
    if (requireNamespace("data.table", quietly = TRUE))
      as.data.frame(data.table::fread(path, sep = "\t", quote = "\"", check.names = FALSE,
                                      data.table = FALSE, colClasses = "character"))
    else read.delim(path, quote = "\"", check.names = FALSE, stringsAsFactors = FALSE,
                    colClasses = "character")
  is_cat <- "mappedGenes" %in% names(df)
  if (is_cat) {
    df$pValue <- suppressWarnings(as.numeric(df$pValue))
    if (!is.na(max_p)) df <- df[!is.na(df$pValue) & df$pValue <= max_p, , drop = FALSE]
    df$assoc_id <- sprintf("A%05d", seq_len(nrow(df)))
  }
  attr(df, "type") <- if (is_cat) "gwas_catalog" else "table"
  df
}

split_genes <- function(g, sep = NA) {
  g <- as.character(g); g <- g[!is.na(g) & !g %in% c("", "NA", "-")]
  if (!is.na(sep)) g <- trimws(unlist(strsplit(g, sep, fixed = TRUE)))
  unique(g[g != ""])
}

# Current symbols via org.Hs.eg.db: current symbols kept; unique aliases mapped to the
# current symbol; ambiguous or unknown symbols kept as given (and flagged).
symbol_map <- function(sym, enable = TRUE) {
  sym <- unique(as.character(sym[!is.na(sym) & sym != ""]))
  ok <- enable && requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
        requireNamespace("AnnotationDbi", quietly = TRUE)
  out <- data.frame(original = sym, symbol = sym,
                    status = if (ok) "current" else "not checked", stringsAsFactors = FALSE)
  if (enable && !ok) message("org.Hs.eg.db/AnnotationDbi not installed: gene symbols NOT harmonized")
  if (!ok || !length(sym)) return(out)
  db <- org.Hs.eg.db::org.Hs.eg.db
  nf <- !sym %in% AnnotationDbi::keys(db, "SYMBOL")
  q <- sym[nf & sym %in% AnnotationDbi::keys(db, "ALIAS")]
  if (length(q)) {
    am <- suppressMessages(AnnotationDbi::select(db, keys = q, keytype = "ALIAS", columns = "SYMBOL"))
    am <- unique(am[!is.na(am$SYMBOL), c("ALIAS", "SYMBOL")])
    n <- table(am$ALIAS)
    uni <- am[am$ALIAS %in% names(n)[n == 1], ]
    i <- match(uni$ALIAS, out$original)
    out$symbol[i] <- uni$SYMBOL; out$status[i] <- "alias -> current symbol"
    out$status[out$original %in% names(n)[n > 1]] <- "ambiguous alias (kept as given)"
  }
  out$status[nf & out$status == "current"] <- "not in org.Hs.eg.db (kept as given)"
  out
}
apply_symbol_map <- function(x, map) {
  i <- match(x, map$original)
  ifelse(is.na(i), x, map$symbol[i])
}
# harmonize a gene column that may hold several genes joined by sep
harmonize_gene_col <- function(x, H, sep = NA) {
  if (is.na(sep)) return(ifelse(is.na(x), NA, H(x)))
  vapply(strsplit(ifelse(is.na(x), "", x), sep, fixed = TRUE), function(v) {
    v <- trimws(v); if (!length(v)) NA_character_ else paste(unique(H(v)), collapse = sep)
  }, "")
}

# GWAS Catalog / gene list -> per-gene evidence table
#   intergenic "A - B" -> both genes; interaction "A x B" -> both genes
# Read-through names (e.g. DTX2P1-UPK3BP1-PMS2P11, STX16-NPEPL1): split on "-" only when
# EVERY part is itself a current gene symbol in org.Hs.eg.db (so HLA-DRB1, NKX2-1,
# SLC9A3-AS1, LYPLAL1-DT stay whole). Returns a named list: symbol -> component genes.
readthrough_components <- function(sym) {
  sym <- unique(sym[grepl("-", sym)])
  if (!length(sym) || !requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) return(list())
  cur <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")
  parts <- strsplit(sym, "-", fixed = TRUE)
  ok <- vapply(parts, function(p) length(p) >= 2 && all(p %in% cur), TRUE)
  setNames(parts[ok], sym[ok])
}

build_gwas_gene_set <- function(df, harmonize = TRUE, expand_readthrough = TRUE) {
  if (identical(attr(df, "type"), "gwas_catalog")) {
    mg <- ifelse(is.na(df$mappedGenes), "", as.character(df$mappedGenes))
    sp <- strsplit(mg, ",|;| - | x ")
    long <- data.frame(assoc_id = rep(df$assoc_id, lengths(sp)),
                       gene_original = trimws(unlist(sp)),
                       intergenic = rep(grepl(" - ", mg), lengths(sp)), stringsAsFactors = FALSE)
    long <- unique(long[!long$gene_original %in% c("", "-", "NA"), ])
    long$from_readthrough <- NA_character_
    if (expand_readthrough) {
      rt <- readthrough_components(long$gene_original)
      if (length(rt)) {
        add <- do.call(rbind, lapply(names(rt), function(g) {
          x <- long[long$gene_original == g, c("assoc_id", "intergenic")]
          data.frame(assoc_id = rep(x$assoc_id, each = length(rt[[g]])),
                     gene_original = rep(rt[[g]], nrow(x)),
                     intergenic = rep(x$intergenic, each = length(rt[[g]])),
                     from_readthrough = g, stringsAsFactors = FALSE)
        }))
        long <- rbind(long, add)
        long <- long[!duplicated(long[, c("assoc_id", "gene_original")]), ]
        message("  read-through names expanded to their component genes: ",
                paste(names(rt), collapse = ", "))
      }
    }
    map <- symbol_map(long$gene_original, harmonize)
    long$geneSymbol <- apply_symbol_map(long$gene_original, map)
    k <- match(long$assoc_id, df$assoc_id)
    col <- function(nm) if (nm %in% names(df)) as.character(df[[nm]][k]) else rep(NA_character_, length(k))
    long$pValue <- df$pValue[k]; long$study <- col("accessionId"); long$pubmed <- col("pubmedId")
    long$trait <- col("efoTraits"); long$snp <- col("riskAllele")
    genes <- do.call(rbind, lapply(split(long, long$geneSymbol), function(x) {
      p <- suppressWarnings(min(x$pValue, na.rm = TRUE))
      data.frame(geneSymbol = x$geneSymbol[1],
                 original_symbols = paste(sort(unique(x$gene_original)), collapse = ";"),
                 n_associations = length(unique(x$assoc_id)),
                 n_studies = length(unique(x$study)),
                 n_publications = length(unique(x$pubmed)),
                 n_traits = length(unique(x$trait)),
                 min_p = ifelse(is.finite(p), p, NA),
                 n_gw_significant = length(unique(x$assoc_id[!is.na(x$pValue) & x$pValue < 5e-8])),
                 only_intergenic = all(x$intergenic),
                 via_readthrough = paste(sort(unique(na.omit(x$from_readthrough))), collapse = ";"),
                 traits = paste(sort(unique(x$trait)), collapse = " | "),
                 studies = paste(sort(unique(x$study)), collapse = ";"),
                 snps = paste(sort(unique(sub("-[^-]*$", "", x$snp))), collapse = ";"),
                 stringsAsFactors = FALSE)
    }))
    genes <- genes[order(-genes$n_associations, genes$geneSymbol), ]
    tr <- do.call(rbind, lapply(split(long, long$trait), function(x) data.frame(
      trait = x$trait[1], n_associations = length(unique(x$assoc_id)),
      n_studies = length(unique(x$study)), n_genes = length(unique(x$geneSymbol)),
      stringsAsFactors = FALSE)))
    tr <- tr[order(-tr$n_associations), ]
  } else {
    gc <- detect_col(df, c("geneSymbol", "gene_symbol", "SYMBOL", "gene", names(df)[1]), "gene")
    keep <- !trimws(as.character(df[[gc]])) %in% c("", "NA", "-")
    g0 <- trimws(as.character(df[[gc]]))[keep]
    map <- symbol_map(g0, harmonize)
    nac <- detect_col(df, c("n_associations"), "n", required = FALSE)
    na <- if (is.na(nac)) rep(1L, length(g0)) else as.integer(df[[nac]][keep])
    na[is.na(na)] <- 1L
    genes <- aggregate(list(n_associations = na), list(geneSymbol = apply_symbol_map(g0, map)), max)
    genes$n_studies <- NA; genes$n_traits <- NA; genes$min_p <- NA; genes$n_gw_significant <- NA
    genes$traits <- NA
    long <- NULL; tr <- NULL
  }
  rownames(genes) <- NULL
  list(genes = genes, long = long, map = map, traits = tr)
}

# SNP positions (GWAS Catalog 'locations' = GRCh38 chr:pos) or a table with chr + pos/start
gwas_snp_gr <- function(df) {
  if (identical(attr(df, "type"), "gwas_catalog")) {
    loc <- ifelse(is.na(df$locations), "", as.character(df$locations))
    sp <- strsplit(loc, ",|;| x ")
    d <- data.frame(assoc_id = rep(df$assoc_id, lengths(sp)), loc = trimws(unlist(sp)),
                    stringsAsFactors = FALSE)
    d <- d[grepl("^(chr)?[0-9XYMxy]+:[0-9]+$", d$loc), ]
    k <- match(d$assoc_id, df$assoc_id)
    chr <- paste0("chr", toupper(sub("^chr", "", sub(":.*", "", d$loc), ignore.case = TRUE)))
    gr <- GRanges(chr, IRanges(as.integer(sub(".*:", "", d$loc)), width = 1))
    gr$assoc_id <- d$assoc_id
    gr$snp <- if ("riskAllele" %in% names(df)) as.character(df$riskAllele[k]) else NA
    gr$trait <- if ("efoTraits" %in% names(df)) as.character(df$efoTraits[k]) else NA
    gr$gwas_p <- df$pValue[k]
    gr$study <- if ("accessionId" %in% names(df)) as.character(df$accessionId[k]) else NA
  } else {
    cc <- detect_col(df, c("chr", "seqnames", "chrom", "CHR"), "chromosome")
    pc <- detect_col(df, c("pos", "position", "start", "BP"), "position")
    chr <- as.character(df[[cc]]); chr <- ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
    gr <- GRanges(chr, IRanges(as.integer(df[[pc]]), width = 1))
    sc <- detect_col(df, c("snp", "rsid", "SNP", "riskAllele"), "snp", required = FALSE)
    gr$snp <- if (is.na(sc)) NA else as.character(df[[sc]])
    gr$trait <- NA; gr$gwas_p <- NA; gr$study <- NA
  }
  gr
}

# Gene-set test for one DMR set.
#   U = genes of the set's DMRichR background (universe); D = DMR genes in U;
#   G = GWAS genes in U. One-sided Fisher (all associations, and p < 5e-8 only), plus an
#   association-weighted check: among testable GWAS genes, do DMR genes carry more GWAS
#   associations than the rest (one-sided Wilcoxon)?
gene_set_enrichment_one <- function(D, U, gs) {
  U <- unique(U); D <- intersect(unique(D), U)
  G <- intersect(gs$geneSymbol, U)
  Ggw <- intersect(gs$geneSymbol[which(gs$n_gw_significant > 0)], U)
  f <- function(G) {
    if (!length(U) || !length(D) || !length(G)) return(c(n = length(G), ov = 0, exp = NA, OR = NA, p = NA))
    tb <- table(factor(U %in% D, c(TRUE, FALSE)), factor(U %in% G, c(TRUE, FALSE)))
    ft <- fisher.test(tb, alternative = "greater")
    c(n = length(G), ov = sum(D %in% G), exp = length(D) * length(G) / length(U),
      OR = unname(ft$estimate), p = ft$p.value)
  }
  a <- f(G); b <- f(Ggw)
  na <- gs$n_associations[match(G, gs$geneSymbol)]; inD <- G %in% D
  w <- if (sum(inD) >= 1 && sum(!inD) >= 1)
    suppressWarnings(wilcox.test(na[inD], na[!inD], alternative = "greater")$p.value) else NA
  ov <- G[inD][order(-na[inD])]
  data.frame(universe_genes = length(U), dmr_genes = length(D),
             gwas_genes_in_universe = a[["n"]], overlap = a[["ov"]], expected = a[["exp"]],
             fold = ifelse(isTRUE(a[["exp"]] > 0), a[["ov"]] / a[["exp"]], NA),
             fisher_OR = a[["OR"]], fisher_p = a[["p"]],
             gw_gwas_genes_in_universe = b[["n"]], gw_overlap = b[["ov"]], gw_expected = b[["exp"]],
             gw_fold = ifelse(isTRUE(b[["exp"]] > 0), b[["ov"]] / b[["exp"]], NA),
             gw_fisher_OR = b[["OR"]], gw_fisher_p = b[["p"]],
             mean_assoc_dmr_genes = if (any(inD)) mean(na[inD]) else NA,
             mean_assoc_other_gwas_genes = if (any(!inD)) mean(na[!inD]) else NA,
             assoc_wilcox_p = w,
             overlap_genes = if (length(ov)) paste0(ov, " (", na[match(ov, G)], ")", collapse = "; ") else "",
             stringsAsFactors = FALSE)
}

# hyper / hypo / mixed per (harmonized) gene for one DMR set
gene_dir_table <- function(d, H) {
  sep <- attr(d, "gene_sep"); if (is.null(sep)) sep <- NA
  d <- d[!is.na(d$.gene), c(".gene", ".dir")]
  if (!nrow(d)) return(data.frame(geneSymbol = character(0), dir = character(0)))
  g <- if (!is.na(sep)) strsplit(d$.gene, sep, fixed = TRUE) else as.list(d$.gene)
  x <- data.frame(geneSymbol = H(trimws(unlist(g))), dir = rep(d$.dir, lengths(g)),
                  stringsAsFactors = FALSE)
  x <- x[x$geneSymbol != "", ]
  a <- tapply(x$dir, x$geneSymbol, function(v) if (all(v > 0)) "hyper" else if (all(v < 0)) "hypo" else "mixed")
  data.frame(geneSymbol = names(a), dir = unname(a), stringsAsFactors = FALSE)
}

# Evidence plot: DMR genes that are GWAS genes. Left: number of GWAS associations
# (fill = strongest p, label = studies / traits). Right: UpSet-style dot matrix of the DMR sets
# (coloured = gene is a DMR gene in that set, red hyper / blue hypo; grey = not), rows joined by
# a line and shaded alternately across both panels; vertical lines separate set_groups.
plot_gwas_evidence <- function(ev, info, sets, file, title, max_genes = 60, label = "GWAS",
                               set_groups = NULL) {
  if (is.null(ev) || !nrow(ev)) return(invisible(NULL))
  info <- info[info$geneSymbol %in% ev$geneSymbol, ]
  nset <- table(ev$geneSymbol)
  info$n_sets <- as.integer(nset[info$geneSymbol])
  info <- info[order(-info$n_associations, -info$n_sets, info$geneSymbol), ]
  dropped <- max(0, nrow(info) - max_genes)
  info <- head(info, max_genes); ev <- ev[ev$geneSymbol %in% info$geneSymbol, ]
  lv <- rev(info$geneSymbol); ng <- length(lv)
  info$nlp <- -log10(pmax(info$min_p, 1e-300))
  info$lab <- ifelse(is.na(info$n_studies), "", paste0(info$n_studies, " st / ", info$n_traits, " tr"))
  stripes <- data.frame(pos = seq(1, ng, by = 2))
  stripe_layer <- geom_rect(data = stripes, aes(xmin = -Inf, xmax = Inf, ymin = pos - .5, ymax = pos + .5),
                            inherit.aes = FALSE, fill = "grey93", colour = NA)
  p1 <- ggplot(info, aes(n_associations, geneSymbol)) + stripe_layer +
    geom_col(aes(fill = nlp), width = .7) +
    geom_text(aes(label = lab), hjust = -0.1, size = 2.5, colour = "grey25") +
    scale_y_discrete(limits = lv) +
    scale_x_continuous(expand = expansion(mult = c(0, .35))) +
    scale_fill_gradient(low = "#C6DBEF", high = "#08306B", na.value = "grey70", name = "-log10(best p)") +
    labs(x = paste0(label, " associations (n)"), y = NULL, title = title,
         subtitle = paste0("DMR genes (>= 1 set) that are ", label, " genes",
                           if (dropped) paste0("; top ", max_genes, " shown (", dropped, " more in table)") else ""),
         caption = "st = studies, tr = traits") +
    theme_bw(10) + theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom",
                         panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
                         plot.caption = element_text(hjust = 0, colour = "grey30"))
  grid_df <- expand.grid(set = sets, geneSymbol = lv, stringsAsFactors = FALSE)
  k <- match(paste(grid_df$set, grid_df$geneSymbol), paste(ev$set, ev$geneSymbol))
  grid_df$dir <- ev$dir[k]; grid_df$on <- !is.na(k)
  grid_df$x <- match(grid_df$set, sets); grid_df$y <- match(grid_df$geneSymbol, lv)
  seg <- do.call(rbind, lapply(split(grid_df[grid_df$on, ], grid_df$geneSymbol[grid_df$on]), function(z)
    if (nrow(z) >= 2) data.frame(y = z$y[1], xmin = min(z$x), xmax = max(z$x)) else NULL))
  vl <- if (!is.null(set_groups)) which(head(set_groups, -1) != tail(set_groups, -1)) + .5 else numeric(0)
  p2 <- ggplot(grid_df, aes(x, y)) + stripe_layer +
    { if (length(vl)) geom_vline(xintercept = vl, colour = "grey55", linewidth = .4) } +
    { if (!is.null(seg) && nrow(seg)) geom_segment(data = seg, aes(x = xmin, xend = xmax, y = y, yend = y),
                                                    inherit.aes = FALSE, colour = "grey35", linewidth = .7) } +
    geom_point(data = grid_df[!grid_df$on, ], colour = "grey85", size = 2.4) +
    geom_point(data = grid_df[grid_df$on, ], aes(colour = dir), size = 3.6) +
    scale_colour_manual(values = c(hyper = "#B2182B", hypo = "#2166AC", mixed = "#7F7F7F"),
                        name = "DMR direction", drop = FALSE) +
    scale_x_continuous(breaks = seq_along(sets), labels = sets, limits = c(.5, length(sets) + .5),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_len(ng), limits = c(.5, ng + .5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = " ", subtitle = " ", caption = " ") + theme_bw(10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
          panel.grid = element_blank())
  p1 <- p1 + scale_y_discrete(limits = lv, expand = c(0, .5))
  g <- cbind(ggplotGrob(p1), ggplotGrob(p2), size = "max")
  panels <- grep("^panel", g$layout$name)
  cols <- sort(unique(g$layout$l[panels]))
  wr <- max(1.2, .3 * length(sets))
  g$widths[cols] <- unit(c(3, wr), "in")
  pdf(file, width = 3 + wr + 3.4, height = 3 + .19 * ng)
  grid::grid.newpage(); grid::grid.draw(g); dev.off()
  invisible(info)
}

# Gene-set enrichment summary: fold enrichment per set, all associations vs p < 5e-8
plot_gene_set_enrichment <- function(res, file, title) {
  mk <- function(test, ov, ex, fdr) data.frame(set = res$set_id, test = test, overlap = ov,
                                               expected = ex, fdr = fdr, stringsAsFactors = FALSE)
  d <- rbind(mk("All associations", res$overlap, res$expected, res$fisher_FDR_run),
             mk("p < 5e-8 associations only", res$gw_overlap, res$gw_expected, res$gw_fisher_FDR_run))
  d <- d[!is.na(d$expected) & d$expected > 0, ]
  if (!nrow(d)) return(invisible(NULL))
  d$fold <- d$overlap / d$expected
  d$sig <- ifelse(!is.na(d$fdr) & d$fdr < 0.05, "FDR < 0.05", "n.s.")
  d$set <- factor(d$set, rev(res$set_id))
  p <- ggplot(d, aes(fold, set)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_segment(aes(x = 1, xend = fold, yend = set), colour = "grey65") +
    geom_point(aes(size = overlap, colour = sig)) +
    geom_text(aes(label = sprintf("%d / %.1f", overlap, expected)), hjust = -0.35, size = 2.8) +
    facet_wrap(~ test) +
    scale_colour_manual(values = c("FDR < 0.05" = "#D7191C", "n.s." = "#2C7BB6"), name = NULL, drop = FALSE) +
    scale_size_area(max_size = 6, name = "Overlapping genes") +
    scale_x_continuous(expand = expansion(mult = c(.05, .25))) +
    labs(x = "Fold enrichment (observed / expected GWAS genes among DMR genes)", y = NULL, title = title,
         caption = paste("Label: observed / expected. Universe = genes annotated to each set's DMRichR",
                         "background regions; one-sided Fisher, BH across all sets.")) +
    theme_bw(10) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, colour = "grey30"))
  ggsave(file, p, width = 10, height = max(3.5, .35 * length(unique(d$set)) + 2))
  invisible(p)
}

# regioneR: are DMRs closer to GWAS SNPs than resampled background (candidate) regions?
#   B = SNP +/- window (merged); statistic = number of DMRs overlapping a window.
perm_snp_set <- function(A, bg, W, ntimes = 100, per_chr = TRUE, alternative = "greater", seed = 1) {
  if (is.null(bg)) return(NULL)
  U <- granges(bg); A <- granges(A)[overlapsAny(A, U)]
  if (!length(A) || length(U) <= length(A)) return(NULL)
  if (per_chr) {
    nA <- table(as.character(seqnames(A))); nU <- table(as.character(seqnames(U)))
    if (any(nU[names(nA)] < nA)) per_chr <- FALSE
  }
  ev <- list(numOverlaps = function(A, B, ...) sum(overlapsAny(A, B)))
  set.seed(seed)
  pt <- regioneR::permTest(A = A, B = W, universe = U, per.chromosome = per_chr, ntimes = ntimes,
                           randomize.function = regioneR::resampleRegions,
                           randomize.function.name = "resampleRegions",
                           evaluate.function = ev, alternative = alternative, verbose = FALSE)
  r <- pt$numOverlaps; perm <- r$permuted
  list(pt = pt, permuted = perm,
       summary = data.frame(randomization = "resample", n_dmrs_tested = length(A), universe = length(U),
                            per_chromosome = per_chr, observed = r$observed,
                            perm_mean = mean(perm), perm_sd = sd(perm),
                            fold_enrichment = ifelse(mean(perm) > 0, r$observed / mean(perm), NA),
                            z = r$zscore, perm_p = r$pval, alternative = r$alternative,
                            ntimes = ntimes, stringsAsFactors = FALSE))
}