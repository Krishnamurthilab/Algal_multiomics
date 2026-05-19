#!/bin/bash
set -euo pipefail

# =========================================================
# MetaBAT2 Binning
# Input : contigs + depth file (from mapping step)
# Output: MAG bins
# =========================================================
cd Whole_genome_metagenome/results
echo "Starting MetaBAT2 binning..."

# -------------------------
# INPUTS
# -------------------------
CONTIGS="filtered_contigs_WGM/*.fasta"
DEPTH="03_mapping/*_depth.txt"

OUTDIR="05_metabat2"
PREFIX="(basename)"

THREADS=16
MIN_CONTIG=1500

mkdir -p ${OUTDIR}

# -------------------------
# RUN METABAT2
# -------------------------
echo "Running MetaBAT2..."

metabat2 \
  -i ${CONTIGS} \
  -a ${DEPTH} \
  -o ${OUTDIR}/${PREFIX}_bin \
  -t ${THREADS} \
  --minContig ${MIN_CONTIG} \
  --unbinned \
  > ${OUTDIR}/metabat2.log 2>&1

# -------------------------
# ORGANIZE OUTPUT
# -------------------------
mkdir ${OUTDIR}/bins
mv ${OUTDIR}/${PREFIX}_bin.*.fa ${OUTDIR}/bins/ 2>/dev/null || true

# -------------------------
# DONE
# -------------------------
echo "MetaBAT2 binning completed successfully!"
echo "Outputs:"
echo "  - ${OUTDIR}/bins/"
echo "  - ${OUTDIR}/metabat2.log"
