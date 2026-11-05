#!/bin/bash
#SBATCH --job-name=dmr_enrich
#SBATCH --account=lasallegrp
#SBATCH --partition=high
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/quobyte/lasallegrp/projects/LACOFD/results/dmrichr/logs/enrichment_%j.log
#=====================================================================================
#  07_run_enrichment.sbatch
#  Gene-set / GREAT enrichment (07_enrichment_plot.R + helper.R) for every finished
#  comparison in one or more DMRichR runs.
#
#  Usage (from anywhere):
#    mkdir -p /quobyte/lasallegrp/projects/LACOFD/results/dmrichr/logs
#    sbatch scripts/dmrichr/07_run_enrichment.sbatch DMR_AZ_04
#    sbatch scripts/dmrichr/07_run_enrichment.sbatch DMR_AZ_04,DMR_AZ_05
#    sbatch scripts/dmrichr/07_run_enrichment.sbatch DMR_AZ_04 pre_recruit_vs_instructor,post_recruit_vs_instructor
#
#  Optional environment variables (put before sbatch, e.g. MIN_DMRS=20 sbatch ...):
#    MIN_DMRS=50      skip comparisons with fewer DMRs than this
#    FORCE=1          rerun comparisons that already have DMRs/enrichment/
#    EXTRA_ARGS="..." extra options passed to 07_enrichment_plot.R (e.g. "--mode presentation")
#
#  Output per comparison: <run>/<comparison>/DMRs/enrichment/  (+ enrichment.log there)
#  Summary for the job:   results/dmrichr/logs/enrichment_<jobid>.log
#=====================================================================================
set -uo pipefail

RUNS="${1:?Give run ID(s), e.g. DMR_AZ_04 or DMR_AZ_04,DMR_AZ_05}"
COMPS="${2:-}"                                   # optional comma list of comparisons
MIN_DMRS="${MIN_DMRS:-5}"
FORCE="${FORCE:-0}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

PROJECT=/quobyte/lasallegrp/projects/LACOFD
RUN_ROOT=$PROJECT/results/dmrichr
PIXI_MANIFEST=$PROJECT/reproducibility/dmrichr/pixi/pixi.toml
RSCRIPT=$PROJECT/scripts/dmrichr/07_enrichment_plot.R
HELPER=$PROJECT/scripts/dmrichr/helper_overlap.R
GENE_SETS=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_dmrich_George/reference/gene_sets_hg38.rds
# Local GREAT gene table, needed if rGREAT cannot download it on compute nodes (leave empty to skip)
export GREAT_GENE_TABLE="${GREAT_GENE_TABLE:-}"

for f in "$PIXI_MANIFEST" "$RSCRIPT" "$HELPER" "$GENE_SETS"; do
  [ -f "$f" ] || { echo "Missing: $f"; exit 1; }
done

eval "$(pixi shell-hook --manifest-path "$PIXI_MANIFEST")"
echo "Job ${SLURM_JOB_ID:-local} on $(hostname) at $(date)"
echo "R: $(which Rscript)   pixi env: $PIXI_MANIFEST"
for p in clusterProfiler rGREAT rrvgo org.Hs.eg.db ggsci Hmisc; do
  Rscript -e "q(status = !requireNamespace('$p', quietly = TRUE))" >/dev/null 2>&1 \
    || echo "WARNING: R package $p not in this env (that part of the analysis will be skipped or fail)"
done
echo

declare -a done_list=() skip_list=() fail_list=()

IFS=',' read -ra RUN_ARR <<< "$RUNS"
for RUN in "${RUN_ARR[@]}"; do
  RD=$RUN_ROOT/$RUN
  [ -d "$RD" ] || { echo "Run folder not found: $RD"; skip_list+=("$RUN (no folder)"); continue; }

  for C in "$RD"/*/; do
    C=${C%/}; COMP=$(basename "$C")
    [ -f "$C/sample_info.xlsx" ] || continue                       # not a comparison folder
    if [ -n "$COMPS" ] && [[ ",$COMPS," != *",$COMP,"* ]]; then continue; fi

    DMRS=$C/DMRs/DMRs_annotated.xlsx
    BG=$C/DMRs/backgroundRegions.bed
    if [ ! -f "$DMRS" ] || [ ! -f "$BG" ]; then
      echo "[$RUN/$COMP] DMRichR not finished -> skipped"; skip_list+=("$RUN/$COMP (not finished)"); continue
    fi
    N=$(grep -vc '^track\|^#' "$C/DMRs/DMRs.bed" 2>/dev/null || true); N=${N:-0}
    if [ "$N" -lt "$MIN_DMRS" ]; then
      echo "[$RUN/$COMP] $N DMRs (< $MIN_DMRS) -> skipped"; skip_list+=("$RUN/$COMP ($N DMRs)"); continue
    fi
    OUT=$C/DMRs/enrichment
    if [ -d "$OUT" ] && [ "$FORCE" != "1" ]; then
      echo "[$RUN/$COMP] $OUT exists -> skipped (FORCE=1 to rerun)"; skip_list+=("$RUN/$COMP (done before)"); continue
    fi
    mkdir -p "$OUT"

    echo "[$RUN/$COMP] $N DMRs -> $OUT  ($(date +%H:%M))"
    # shellcheck disable=SC2086
    if Rscript "$RSCRIPT" \
         --dmrs "$DMRS" \
         --background_regions "$BG" \
         --gene_set_cache "$GENE_SETS" \
         --helper_file "$HELPER" \
         --outdir "$OUT" \
         --prefix "$COMP" \
         $EXTRA_ARGS > "$OUT/enrichment.log" 2>&1; then
      done_list+=("$RUN/$COMP")
    else
      echo "[$RUN/$COMP] FAILED, see $OUT/enrichment.log (last lines below)"; tail -n 5 "$OUT/enrichment.log"
      fail_list+=("$RUN/$COMP")
    fi
  done
done

echo
echo "===== Summary $(date) ====="
echo "Done    (${#done_list[@]}): ${done_list[*]:-none}"
echo "Skipped (${#skip_list[@]}):"; for s in ${skip_list[@]+"${skip_list[@]}"}; do echo "  $s"; done
echo "Failed  (${#fail_list[@]}): ${fail_list[*]:-none}"
[ ${#fail_list[@]} -eq 0 ]
