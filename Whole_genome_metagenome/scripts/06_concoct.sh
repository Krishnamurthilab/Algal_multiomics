#!/bin/bash
set -euo pipefail

# =========================================================
# CONCOCT BINNING
# Input : contigs + BAM file (from mapping step)
# Output: CONCOCT bins
# =========================================================
cd Whole_genome_metagenome/results
echo "Starting CONCOCT binning..."

# -------------------------
# INPUTS
# -------------------------
CONTIGS="filtered_contigs_WGM/S1.fasta"
BAM="03_mapping/S1_mapped.sorted.bam"

OUTDIR="06_concoct"
PREFIX="(basename)"

THREADS=16
CHUNK_SIZE=10000

mkdir -p ${OUTDIR}

# -------------------------
# 1. CUT CONTIGS INTO FRAGMENTS
# -------------------------
echo "Cutting contigs..."

cut_up_fasta.py \
  ${CONTIGS} \
  -c ${CHUNK_SIZE} \
  -o 0 \
  --merge_last \
  -b ${OUTDIR}/cuts.bed \
  > ${OUTDIR}/cuts.fa

# -------------------------
# 2. COVERAGE TABLE
# -------------------------
echo "Generating coverage table..."

concoct_coverage_table.py \
  ${OUTDIR}/cuts.bed \
  ${BAM} \
  > ${OUTDIR}/coverage.tsv

# -------------------------
# 3. RUN CONCOCT
# -------------------------
echo "Running CONCOCT..."

mkdir -p ${OUTDIR}/concoct_output

concoct \
  --composition_file ${OUTDIR}/cuts.fa \
  --coverage_file ${OUTDIR}/coverage.tsv \
  -b ${OUTDIR}/concoct_output \
  -t ${THREADS} \
  > ${OUTDIR}/concoct.log 2>&1

# -------------------------
# 4. MERGE CUTUP CLUSTERS
# -------------------------
echo "Merging clusters..."

merge_cutup_clustering.py \
  ${OUTDIR}/concoct_output/clustering_gt1000.csv \
  > ${OUTDIR}/concoct_output/clustering_merged.csv

# -------------------------
# 5. EXTRACT BIN FASTA FILES
# -------------------------
echo "Extracting bins..."

mkdir -p ${OUTDIR}/bins

extract_fasta_bins.py \
  ${CONTIGS} \
  ${OUTDIR}/concoct_output/clustering_merged.csv \
  --output_path ${OUTDIR}/bins

# -------------------------
# DONE
# -------------------------
echo "CONCOCT binning completed successfully!"
echo "Outputs:"
echo "  - ${OUTDIR}/bins/"
echo "  - ${OUTDIR}/concoct.log"

