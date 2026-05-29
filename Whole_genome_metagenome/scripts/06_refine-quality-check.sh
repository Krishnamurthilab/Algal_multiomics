#!/bin/bash
set -euo pipefail

# =========================================================
#
# STEP 06: BIN REFINEMENT (DAS Tool) + QUALITY CHECK (CheckM)
# Input : bins from MetaBAT2 + MaxBin2 + CONCOCT
# Output: refined MAGs + quality statistics
# =========================================================

echo "Starting DAS Tool + CheckM pipeline..."

THREADS=16
SAMPLE="S1"

# =========================================================
# INPUT DIRECTORIES
# =========================================================

METABAT_DIR="results/04_metabat2/bins"
MAXBIN_DIR="results/03_maxbin2/bins"
CONCOCT_DIR="results/05_concoct/bins"

CONTIGS="results/Assembly/S1.fasta"

OUTDIR="results/06_refinement_quality"
mkdir -p ${OUTDIR}

# =========================================================
# STEP 1: CREATE DAS TOOL INPUT TABLES
# =========================================================

echo "Preparing DAS Tool inputs..."

mkdir -p ${OUTDIR}/dastool_inputs

# -------------------------
# MetaBAT2 scaffolds2bin
# -------------------------
for f in ${METABAT_DIR}/*.fa; do
  bin=$(basename "$f" .fa)
  grep ">" "$f" | sed 's/>//' | awk -v b="$bin" '{print $1"\t"b}'
done > ${OUTDIR}/dastool_inputs/metabat.scaffolds2bin.tsv

# -------------------------
# MaxBin2 scaffolds2bin
# -------------------------
for f in ${MAXBIN_DIR}/*.fasta; do
  bin=$(basename "$f" .fasta)
  grep ">" "$f" | sed 's/>//' | awk -v b="$bin" '{print $1"\t"b}'
done > ${OUTDIR}/dastool_inputs/maxbin.scaffolds2bin.tsv

# -------------------------
# CONCOCT scaffolds2bin
# -------------------------
awk -F',' 'NR>1 {print $1"\tbin_"$2}' \
  ${CONCOCT_DIR}/../concoct_output/clustering_merged.csv \
  > ${OUTDIR}/dastool_inputs/concoct.scaffolds2bin.tsv

# =========================================================
# STEP 2: RUN DAS TOOL
# =========================================================

echo "Running DAS Tool..."

DAS_Tool \
  -i ${OUTDIR}/dastool_inputs/metabat.scaffolds2bin.tsv,${OUTDIR}/dastool_inputs/maxbin.scaffolds2bin.tsv,${OUTDIR}/dastool_inputs/concoct.scaffolds2bin.tsv \
  -l metabat2,maxbin2,concoct \
  -c ${CONTIGS} \
  -o ${OUTDIR}/${SAMPLE}_DAS \
  --search_engine diamond \
  --threads ${THREADS} \
  > ${OUTDIR}/dastool.log 2>&1

# =========================================================
# STEP 3: CHECKM QUALITY ASSESSMENT
# =========================================================

echo "Running CheckM..."

CHECKM_OUT=${OUTDIR}/checkm_output
mkdir -p ${CHECKM_OUT}

checkm lineage_wf \
  -x fa \
  -t ${THREADS} \
  ${OUTDIR}/${SAMPLE}_DAS_DASTool_bins \
  ${CHECKM_OUT} \
  > ${OUTDIR}/checkm.log 2>&1

# =========================================================
# STEP 4: SUMMARY TABLE
# =========================================================

echo "Generating summary..."

checkm qa \
  ${CHECKM_OUT}/lineage.ms \
  ${CHECKM_OUT} \
  -o 2 \
  > ${OUTDIR}/checkm_summary.tsv

# =========================================================
# DONE
# =========================================================

echo "DAS Tool + CheckM completed successfully!"
echo "Outputs:"
echo "  - ${OUTDIR}/${SAMPLE}_DAS_DASTool_bins/"
echo "  - ${CHECKM_OUT}/"
echo "  - ${OUTDIR}/checkm_summary.tsv"
