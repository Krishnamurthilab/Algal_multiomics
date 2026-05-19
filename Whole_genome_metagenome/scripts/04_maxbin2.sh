#!/bin/bash
set -euo pipefail

# =========================================================
# MaxBin2 Binning
# Input: contigs + depth file from Bowtie2 mapping step
# Output: genome bins (MaxBin2)
# =========================================================
cd Whole_genome_metagenome/results
echo "Starting MaxBin2 binning..."

# -------------------------
# INPUTS (edit via config later if needed)
# -------------------------
CONTIGS="filtered_contigs_WGM/*.fasta"
DEPTH="03_mapping/*_depth.txt"

OUTDIR="04_maxbin2"
PREFIX="(basename)"

THREADS=16

mkdir -p ${OUTDIR}

# -------------------------
# RUN MAXBIN2
# -------------------------
echo "Running MaxBin2..."

run_MaxBin.pl \
  -contig ${CONTIGS} \
  -abund ${DEPTH} \
  -out ${OUTDIR}/${PREFIX}_maxbin2 \
  -thread ${THREADS} \
  -min_contig_length 1500 \
  > ${OUTDIR}/maxbin2.log 2>&1

# -------------------------
# OPTIONAL CLEANING / ORGANIZATION
# -------------------------
mkdir -p ${OUTDIR}/bins

mv ${OUTDIR}/${PREFIX}_maxbin2.*.fasta ${OUTDIR}/bins/ 2>/dev/null || true

# -------------------------
# DONE
# -------------------------
echo "MaxBin2 binning completed successfully!"
echo "Outputs:"
echo "  - ${OUTDIR}/bins/"
echo "  - ${OUTDIR}/maxbin2.log"
