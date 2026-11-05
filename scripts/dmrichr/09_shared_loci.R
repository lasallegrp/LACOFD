#!/usr/bin/env Rscript
#=====================================================================================
#  09_shared_loci.R
#  Region sets for GREAT enrichment of loci SHARED between DMRichR comparisons, made
#  from the coordinate output of 05_dmrichr_overlap.R (<group>_overlap_coord.xlsx,
#  sheet Loci). The region-based counterpart of 09_shared_genes.R.
#
#  Foreground: loci where DMRs from every comparison of the group overlap (--shared all)
#              or from at least two (--shared any2); optionally direction-concordant only,
#              optionally split into shared-hyper and shared-hypo loci.
#  Background: genome territory tested in EVERY comparison of the group = base-pair
#              intersection of each comparison's DMRs/backgroundRegions.bed. Shared loci
#              that extend past it are added to it (counted in selection.txt), so no
#              foreground region is silently dropped by GREAT.
#
#  Writes to --out_dir (BED, 0-based, as 07_enrichment_plot.R expects by default):
#    <list>_loci.bed            for 07_enrichment_plot.R --dmrs
#    background_regions.bed     for --background_regions
#    lists.tsv                  list, n_loci, run (TRUE if >= --min_items)
#=====================================================================================
suppressPackageStartupMessages({ library(optparse); library(readxl); library(GenomicRanges) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--overlap_xlsx", type = "character", help = "<group>_overlap_coord.xlsx from 05_dmrichr_overlap.R"),
  make_option("--manifest", type = "character", help = "manifest.tsv from 04_make_overlap_sets.py (set_id -> path)"),
  make_option("--out_dir", type = "character"),
  make_option("--background_bed", type = "character", default = "DMRs/backgroundRegions.bed"),
  make_option("--shared", type = "character", default = "all", help = "all | any2 [default %default]"),
  make_option("--concordant", action = "store_true", default = FALSE),
  make_option("--split_direction", action = "store_true", default = FALSE),
  make_option("--min_items", type = "integer", default = 20L)
)))
for (a in c("overlap_xlsx", "manifest", "out_dir")) if (is.null(opt[[a]])) stop("--", a, " is required")
stopifnot(opt$shared %in% c("all", "any2"))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- shared loci (1-based, from DMRs_annotated.xlsx coordinates) --------------------
lo <- as.data.frame(read_excel(opt$overlap_xlsx, sheet = "Loci"))
dcols <- grep("^dir_", names(lo), value = TRUE)
sets <- sub("^dir_", "", dcols)
if (length(sets) < 2) stop("Expected >= 2 dir_<set> columns in the Loci sheet")
need <- if (opt$shared == "all") length(sets) else 2L
lo <- lo[lo$n_sets >= need, ]
if (opt$concordant) lo <- lo[lo$direction_concordant %in% c(TRUE, "TRUE"), ]
conc <- lo$direction_concordant %in% c(TRUE, "TRUE")
ldir <- apply(lo[, dcols, drop = FALSE], 1, function(v) { v <- unique(na.omit(v)); if (length(v) == 1) v else "mixed" })
lists <- list(shared = lo)
if (opt$split_direction) {
  lists$shared_hyper <- lo[conc & ldir == "hyper", ]
  lists$shared_hypo  <- lo[conc & ldir == "hypo", ]
}

# ---- background: territory tested in every comparison -------------------------------
man <- read.delim(opt$manifest, stringsAsFactors = FALSE)
miss <- setdiff(sets, man$set_id)
if (length(miss)) stop("Sets not in ", opt$manifest, ": ", paste(miss, collapse = ", "))
read_bed <- function(f) {
  b <- read.delim(f, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  b <- b[!grepl("^(track|browser)", b$V1), ]
  reduce(GRanges(b$V1, IRanges(as.integer(b$V2) + 1L, as.integer(b$V3))))   # BED 0-based -> 1-based
}
bgs <- lapply(sets, function(s) {
  f <- file.path(man$path[man$set_id == s], opt$background_bed)
  if (!file.exists(f)) stop("Missing background for ", s, ": ", f)
  g <- read_bed(f)
  message(sprintf("  %s: %s background regions, %.1f Mb", s, format(length(g), big.mark = ","), sum(width(g)) / 1e6))
  g
})
bg <- Reduce(function(a, b) GenomicRanges::intersect(a, b, ignore.strand = TRUE), bgs)
fg <- GRanges(lists$shared$chr, IRanges(lists$shared$start, lists$shared$end))
mid <- GRanges(seqnames(fg), IRanges(floor((start(fg) + end(fg)) / 2), width = 1))
n_out <- sum(!overlapsAny(mid, bg))
if (n_out) message("  ", n_out, " shared loci have a midpoint outside the common background; added to it")
bg <- reduce(c(bg, fg))
message(sprintf("  Common background: %s regions, %.1f Mb", format(length(bg), big.mark = ","), sum(width(bg)) / 1e6))

write_bed <- function(chr, start1, end, f)
  write.table(data.frame(chr, as.integer(start1) - 1L, as.integer(end)), f,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write_bed(as.character(seqnames(bg)), start(bg), end(bg), file.path(opt$out_dir, "background_regions.bed"))
for (n in names(lists)) {
  x <- lists[[n]]
  write_bed(x$chr, x$start, x$end, file.path(opt$out_dir, paste0(n, "_loci.bed")))
  write.table(x, file.path(opt$out_dir, paste0(n, "_loci_table.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
summ <- data.frame(list = names(lists), n = vapply(lists, nrow, 1L),
                   run = vapply(lists, nrow, 1L) >= opt$min_items, stringsAsFactors = FALSE)
write.table(summ, file.path(opt$out_dir, "lists.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("overlap_xlsx:", opt$overlap_xlsx), paste("sets:", paste(sets, collapse = ", ")),
             paste("shared:", opt$shared), paste("concordant:", opt$concordant),
             sprintf("common background: %d regions, %.1f Mb", length(bg), sum(width(bg)) / 1e6),
             paste("shared loci added to background (midpoint outside):", n_out)),
           file.path(opt$out_dir, "selection.txt"))
print(summ, row.names = FALSE)
