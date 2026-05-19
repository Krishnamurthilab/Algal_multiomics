#!/bin/bash
set -euo pipefail

# =========================================================
# MAG TAXONOMY (GTDB-Tk + CAT PACK)
# Input : MAGs from DAS Tool
# Output: taxonomy assignments (GTDB + CAT)
# =========================================================
cd Whole_genome_metagenome/results
echo "Starting MAG taxonomy assignment..."

THREADS=16

MAG_DIR="07_mag_refinement_quality/NCBI_FCS_GX_cleaned_high_quality_MAGs/cleaned_MAGs"
OUTDIR="09_mag_taxonomy"
mkdir -p ${OUTDIR}

# =========================================================
# STEP 1: GTDB-Tk classification
# =========================================================
export GTDBTK_DATA_PATH=/path/to/gtdbtk_db

echo "Running GTDB-Tk..."

gtdbtk classify_wf \
  --genome_dir ${MAG_DIR} \
  --out_dir ${OUTDIR}/gtdbtk_out \
  --cpus ${THREADS} \
  --mash_db ${GTDB_DATA_PATH} \
  > ${OUTDIR}/gtdbtk_out/gtdbtk.log 2>&1

# =========================================================
# STEP 2: CAT PACK TAXONOMY
# =========================================================
CAT_gtdb_db="/path/to/cat_db"
CAT_gtdb_tax="/path/to/cat_gtdb_tax"

echo "Running CAT pack..."

CAT_pack bins \
  -b ${MAG_DIR}/ \
  -d ${CAT_gtdb_db} \
  -t ${CAT_gtdb_tax} \
  -o ${OUTDIR}/cat_bins_output \
  -n ${THREADS}

# =========================================================
# STEP 3: SUMMARY PREPARATION
# =========================================================

echo "Generating taxonomy summary..."

mkdir -p ${OUTDIR}/summary

cp ${OUTDIR}/gtdbtk_out/classify/*.summary.tsv \
   ${OUTDIR}/summary/ 2>/dev/null || true

cp ${OUTDIR}/cat_bins_output/*.tsv \
   ${OUTDIR}/summary/ 2>/dev/null || true

# =========================================================
# STEP 4: UBCG PHYLOGENOMIC ANALYSIS
# =========================================================

echo "Running UBCG phylogenomic analysis..."

UBCG_OUT="${OUTDIR}/ubcg"
THREADS="16"
mkdir -p ${UBCG_OUT}
mkdir -p ${UBCG_OUT}/ucg_profiles

# ---------------------------------------------------------
# STEP 4.1: Generate UBCG profiles (.ucg)
# ---------------------------------------------------------

echo "Generating UBCG profiles..."

GENOME_DIR="09_mag_taxonomy/MAG_ref_genomes"


for genome in ${GENOME_DIR}/*.fasta; do

  base=$(basename ${genome} .fasta)

  ubcg.py \
    -i ${genome} \
    -o ${UBCG_OUT}/ucg_profiles/${base} \
    -t ${THREADS}

done

# ---------------------------------------------------------
# STEP 4.2: Build phylogenomic tree
# ---------------------------------------------------------

echo "Constructing UBCG phylogenomic tree..."

build_ucg_tree.py \
  -i ${UBCG_OUT}/ucg_profiles \
  -o ${UBCG_OUT}/phylogeny \
  -t ${THREADS}

echo "UBCG phylogeny completed."

echo "Tree output:"
echo "  ${UBCG_OUT}/phylogeny/"

# =========================================================
# DONE
# =========================================================

echo "MAG taxonomy completed!"
echo "Outputs:"
echo "  - ${OUTDIR}/gtdbtk_out/"
echo "  - ${OUTDIR}/cat_bins_output/"
echo "  - ${OUTDIR}/summary/"
echo "  - ${OUTDIR}/ubcg/phylogeny"
