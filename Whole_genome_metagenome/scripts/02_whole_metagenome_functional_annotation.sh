#!/bin/bash
set -euo pipefail

# =========================================================
# WHOLE METAGENOME FUNCTIONAL ANALYSIS
# Includes:
#   - BlastKOALA
#   - KEGG Decoder
#   - dbCAN
#   - antiSMASH
# =========================================================
cd Whole_genome_metagenome/results
echo "Starting whole metagenome functional analysis..."

THREADS=16

ASSEMBLY_DIR="results/filtered_contigs_WGM"
OUTDIR="results/02_whole_metagenome_functional_analysis"

mkdir -p ${OUTDIR}

# =========================================================
# KEGG ANALYSIS
# =========================================================

echo "BlastKOALA annotation performed manually."

mkdir -p ${OUTDIR}/kegg/blastkoala

echo "Place KO annotation files in:"
echo "  ${OUTDIR}/kegg/blastkoala/"

KEGG-decoder \
  -i ${OUTDIR}/kegg/blastkoala/combined_ko.tsv \
  -o ${OUTDIR}/kegg/kegg_decoder_output.tsv \
  -v static \
  > ${OUTDIR}/kegg/kegg_decoder.log 2>&1

# =========================================================
# dbCAN ANALYSIS
# =========================================================

echo "Running dbCAN..."

run_dbcan \
  assembly_proteins.faa \
  protein \
  --out_dir ${OUTDIR}/dbcan \
  > ${OUTDIR}/dbcan.log 2>&1


# =========================================================
# antiSMASH
# =========================================================

echo "antiSMASH analysis performed manually via web server."

mkdir -p ${OUTDIR}/antismash
