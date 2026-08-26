#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# LACOFD repository reorganization
#
# Run from:
# /quobyte/lasallegrp/projects/LACOFD
# ============================================================

PROJECT_ROOT="/quobyte/lasallegrp/projects/LACOFD"

cd "${PROJECT_ROOT}"

echo
echo "============================================"
echo "LACOFD repository reorganization"
echo "============================================"
echo "Project root:"
pwd
echo


# ============================================================
# Helper functions
# ============================================================

move_file() {

    src="$1"
    dest="$2"

    if [ -f "${src}" ]; then
        echo "MOVE:"
        echo "  ${src}"
        echo "  -> ${dest}"

        mkdir -p "$(dirname "${dest}")"
        mv "${src}" "${dest}"
    else
        echo "SKIP: file not found: ${src}"
    fi
}


move_dir() {

    src="$1"
    dest="$2"

    if [ -d "${src}" ]; then
        echo "MOVE DIRECTORY:"
        echo "  ${src}"
        echo "  -> ${dest}"

        mkdir -p "$(dirname "${dest}")"
        mv "${src}" "${dest}"
    else
        echo "SKIP: directory not found: ${src}"
    fi
}


# ============================================================
# 1. Create target directory structure
# ============================================================

echo
echo "Creating directory structure..."
echo

mkdir -p config/schema

mkdir -p \
    data/raw/wet_lab/dna_isolation/arizona/DNAEXT_AZ_01 \
    data/raw/wet_lab/dna_isolation/la_county/DNAEXT_LA_01 \
    data/raw/wet_lab/dna_isolation/la_county/DNAEXT_LA_02

mkdir -p \
    data/interim/wet_lab/dna_isolation/arizona \
    data/interim/wet_lab/dna_isolation/la_county

mkdir -p \
    data/processed/wet_lab/dna_isolation/arizona \
    data/processed/wet_lab/dna_isolation/la_county \
    data/processed/wet_lab/dna_isolation/combined

mkdir -p \
    scripts/wet_lab/dna_isolation

mkdir -p \
    results/wet_lab_qc/dna_isolation/arizona \
    results/wet_lab_qc/dna_isolation/la_county \
    results/wet_lab_qc/dna_isolation/combined

mkdir -p \
    results/sequencing_qc/EPI_PILOT_01/library_complexity


# ============================================================
# 3. Move project config out of analysis/
# ============================================================

move_file \
    "analysis/config/config.yaml" \
    "config/config.yaml"

move_file \
    "analysis/config/samples.tsv" \
    "config/samples.tsv"

move_file \
    "analysis/config/units.tsv" \
    "config/units.tsv"

move_file \
    "analysis/config/schema/config.schema.yaml" \
    "config/schema/config.schema.yaml"


# ============================================================
# 4. Move renv.lock to project root
# ============================================================

move_file \
    "analysis/renv.lock" \
    "renv.lock"


# ============================================================
# 5. Reorganize DNA isolation scripts
# ============================================================

move_file \
    "scripts/dna_isolation/01_clean_dna_isolation.R" \
    "scripts/wet_lab/dna_isolation/01_clean_dna_isolation.R"

move_file \
    "scripts/dna_isolation/02_plot_dna_isolation_pre_post.R" \
    "scripts/wet_lab/dna_isolation/02_plot_dna_isolation_pre_post.R"


# ============================================================
# 6. Reorganize cleaned DNA isolation data
# ============================================================

CURRENT_DNA_DIR="data/processed/wet_lab/dna_isolation/la_county"

# Cleaned dataset belongs in interim/
move_file \
    "${CURRENT_DNA_DIR}/la_county_dna_isolation_cleaned.tsv" \
    "data/interim/wet_lab/dna_isolation/la_county/la_county_dna_isolation_cleaned.tsv"

# Analysis-ready dataset stays under processed/
if [ -f "${CURRENT_DNA_DIR}/la_county_dna_isolation_analysis_ready.tsv" ]; then
    echo "KEEP:"
    echo "  ${CURRENT_DNA_DIR}/la_county_dna_isolation_analysis_ready.tsv"
fi


# ============================================================
# 7. Preserve redundant older datasets temporarily
# ============================================================
#
# Rather than deleting these now, move them into an archive.
# Once you confirm they contain no unique information,
# you can remove the archive.
# ============================================================

mkdir -p \
    data/interim/wet_lab/dna_isolation/la_county/archive_old_outputs

move_file \
    "${CURRENT_DNA_DIR}/la_county_dna_isolation_aliquot_level_cleaned.tsv" \
    "data/interim/wet_lab/dna_isolation/la_county/archive_old_outputs/la_county_dna_isolation_aliquot_level_cleaned.tsv"

move_file \
    "${CURRENT_DNA_DIR}/la_county_dna_isolation_extraction_level.tsv" \
    "data/interim/wet_lab/dna_isolation/la_county/archive_old_outputs/la_county_dna_isolation_extraction_level.tsv"


# ============================================================
# 8. Reorganize Bioanalyzer
# ============================================================

BIOA_ROOT="data/raw/wet_lab/bioanalyzer/la_county"

for n in 1 2 3
do

    RUN_ID=$(printf "BIOA_LA_%02d" "${n}")

    mkdir -p \
        "${BIOA_ROOT}/${RUN_ID}/submission_form" \
        "${BIOA_ROOT}/${RUN_ID}/reports"

    # Move submission forms
    if [ -d "${BIOA_ROOT}/submission_forms/submission_${n}" ]; then

        find \
            "${BIOA_ROOT}/submission_forms/submission_${n}" \
            -maxdepth 1 \
            -type f \
            -exec mv {} "${BIOA_ROOT}/${RUN_ID}/submission_form/" \;

    fi

    # Move Bioanalyzer reports/electropherograms
    if [ -d "${BIOA_ROOT}/electropherograms/submission_${n}" ]; then

        find \
            "${BIOA_ROOT}/electropherograms/submission_${n}" \
            -maxdepth 1 \
            -type f \
            -exec mv {} "${BIOA_ROOT}/${RUN_ID}/reports/" \;

    fi

done


# ============================================================
# 9. Remove empty old Bioanalyzer directories
# ============================================================

find \
    "${BIOA_ROOT}/submission_forms" \
    -type d \
    -empty \
    -delete \
    2>/dev/null || true

find \
    "${BIOA_ROOT}/electropherograms" \
    -type d \
    -empty \
    -delete \
    2>/dev/null || true


# ============================================================
# 10. Move pilot library-complexity metrics
# ============================================================

if [ -d "analysis/SmallSequencingBatchForQC" ]; then

    find \
        "analysis/SmallSequencingBatchForQC" \
        -maxdepth 1 \
        -type f \
        -exec mv {} \
        "results/sequencing_qc/EPI_PILOT_01/library_complexity/" \;

fi


# ============================================================
# 11. Rename Epigenerator pilot run
# ============================================================

if [ -d "pipeline_runs/epigenerator/SmallSequencingBatchForQC" ] && \
   [ ! -d "pipeline_runs/epigenerator/EPI_PILOT_01" ]; then

    mv \
        "pipeline_runs/epigenerator/SmallSequencingBatchForQC" \
        "pipeline_runs/epigenerator/EPI_PILOT_01"

    echo "Renamed Epigenerator pilot to EPI_PILOT_01"

elif [ -d "pipeline_runs/epigenerator/EPI_PILOT_01" ]; then

    echo "SKIP: EPI_PILOT_01 already exists."

fi


# ============================================================
# 12. Clean empty legacy directories
# ============================================================

rmdir \
    analysis/config/schema \
    2>/dev/null || true

rmdir \
    analysis/config \
    2>/dev/null || true

rmdir \
    analysis/SmallSequencingBatchForQC \
    2>/dev/null || true

rmdir \
    scripts/dna_isolation \
    2>/dev/null || true


# ============================================================
# 13. Remove accidental Rplots.pdf
# ============================================================

if [ -f "Rplots.pdf" ]; then

    echo "Removing accidental Rplots.pdf"
    rm "Rplots.pdf"

fi


# ============================================================
# 14. Add README placeholders to important batch directories
# ============================================================

for batch_dir in \
    data/raw/wet_lab/dna_isolation/arizona/DNAEXT_AZ_01 \
    data/raw/wet_lab/dna_isolation/la_county/DNAEXT_LA_01 \
    data/raw/wet_lab/dna_isolation/la_county/DNAEXT_LA_02
do

    if [ ! -f "${batch_dir}/README.md" ]; then

        touch "${batch_dir}/README.md"

        echo "Created placeholder:"
        echo "  ${batch_dir}/README.md"

    fi

done


for bioa_dir in \
    "${BIOA_ROOT}/BIOA_LA_01" \
    "${BIOA_ROOT}/BIOA_LA_02" \
    "${BIOA_ROOT}/BIOA_LA_03"
do

    if [ ! -f "${bioa_dir}/README.md" ]; then

        touch "${bioa_dir}/README.md"

        echo "Created placeholder:"
        echo "  ${bioa_dir}/README.md"

    fi

done


# ============================================================
# 15. Final report
# ============================================================

echo
echo "============================================"
echo "Reorganization completed"
echo "============================================"
echo
echo "Recommended checks:"
echo
echo "  tree -L 6"
echo
echo "  git status"
echo
echo "Check extraction files:"
echo
echo "  tree data/raw/wet_lab/dna_isolation"
echo
echo "Check Bioanalyzer:"
echo
echo "  tree data/raw/wet_lab/bioanalyzer/la_county"
echo
echo "Check Epigenerator:"
echo
echo "  ls pipeline_runs/epigenerator"
echo
echo "============================================"