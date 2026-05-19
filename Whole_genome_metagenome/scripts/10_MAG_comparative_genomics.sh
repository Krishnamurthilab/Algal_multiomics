#!/bin/bash
set -euo pipefail

# =========================================================
# COMPARATIVE GENOMICS + FUNCTIONAL PROFILING
# Includes:
#   1. OrthoANIu v1.2 (ANI)
#   2. EzAAI v1.2.4 (AAI)
#   3. KEGG Decoder (metabolic inference)
#   4. dbCAN (CAZyme annotation)
#   5. antiSMASH (BGC prediction)
#   6. UBCG (phylogenomic tree)
# =========================================================

echo "Starting comparative genomics analysis..."

THREADS=16

MAG_DIR="07_mag_refinement_quality/NCBI_FCS_GX_cleaned_high_quality_MAGs/cleaned_MAGs"
OUTDIR="10_mag_comparative_genomics"

mkdir -p ${OUTDIR}

# ==========================================================================
# STEP 1: OrthoANIu [EZAAI orthoaniu stand alone pipeline (ANI calculation)]
# ==========================================================================

echo "Running OrthoANIu..."

OUTDIR="10_mag_comparative_genomics/orthoani"

mkdir -p ${OUTDIR}

java -jar /path/to/bi-orthoaniu.jar \
  -fd ${MAG_DIR} \
  -n 16 \
  -o ${OUTDIR} \
  -t ${OUTDIR}/tmp \
  -u /path/to/usearch

echo "OrthoANIu completed."
echo "Output directory: ${OUTDIR}"


# =========================================================
# STEP 2: EzAAI (Extract → Calculate → Cluster)
# =========================================================

echo "Running EzAAI pipeline..."

MAG_DIR="07_mag_refinement_quality/NCBI_FCS_GX_cleaned_high_quality_MAGs/cleaned_MAGs"
OUTDIR="10_mag_comparative_genomics/ezaai"

mkdir -p ${OUTDIR}/db
mkdir -p ${OUTDIR}/out

# =========================================================
# STEP 2.1: EXTRACT CDS PROFILE DATABASES
# =========================================================

echo "EzAAI extract step..."

ezaai extract \
  -i ${MAG_DIR} \
  -o ${OUTDIR}/db \
  -l ${OUTDIR}/labels.tsv

# Output: *.db files per genome

# =========================================================
# STEP 2.2: CALCULATE AAI MATRIX
# =========================================================

echo "EzAAI calculate step..."

ezaai calculate \
  -i ${OUTDIR}/db \
  -j ${OUTDIR}/db \
  -o ${OUTDIR}/out/aai.tsv


# =========================================================
# STEP 3: KEGG FUNCTIONAL PROFILING
# BlastKOALA → KEGG Decoder
# =========================================================

echo "KEGG analysis..."

KEGG_OUT="${OUTDIR}/kegg_decoder"
mkdir -p ${KEGG_OUT}

# ---------------------------------------------------------
# STEP 3.1: BlastKOALA annotation
# ---------------------------------------------------------

echo "BlastKOALA was performed manually via KEGG web server."
echo "Place KO annotation files in:"
echo "  ${KEGG_OUT}/blastkoala/"

mkdir -p ${KEGG_OUT}/blastkoala

# Expected input:
# *.txt or KO assignment files from BlastKOALA

# ---------------------------------------------------------
# STEP 3.2: Combine KO annotations
# ---------------------------------------------------------

echo "Combining KO annotations from all MAGs..."

cat ${KEGG_OUT}/blastkoala/*.tsv \
  > ${KEGG_OUT}/combined_ko.tsv

# ---------------------------------------------------------
# STEP 3.3: KEGG Decoder
# ---------------------------------------------------------

echo "Running KEGG Decoder..."

KEGG-decoder \
  -i ${KEGG_OUT}/combined_ko.tsv \
  -o ${KEGG_OUT}/kegg_decoder_output.tsv \
  -v static \
  > ${KEGG_OUT}/kegg_decoder.log 2>&1

echo "KEGG Decoder completed."


# =========================================================
# STEP 4: dbCAN (CAZyme annotation)
# =========================================================

echo "Running dbCAN CAZyme annotation..."

DBCAN_OUT="${OUTDIR}/dbcan"

mkdir -p ${DBCAN_OUT}

FAA_DIR="10_mag_comparative_genomics/dbcan"

for faa in ${FAA_DIR}/*.faa; do

  base=$(basename ${faa} .faa)

  echo "Processing ${base}..."

  run_dbcan \
    ${faa} \
    protein \
    --out_dir ${DBCAN_OUT}/${base} \
    > ${DBCAN_OUT}/${base}.log 2>&1

done

echo "dbCAN completed."


#==========================================================
# STEP 5: antiSMASH (NOTE: webserver-based)
# =========================================================

echo "antiSMASH is web-based (manual upload required)."
echo "Place results manually in:"
echo "  ${OUTDIR}/antismash/"

mkdir -p ${OUTDIR}/antismash

#==========================================================
# DONE
# =========================================================

echo "Comparative genomics completed!"
echo "Outputs:"
echo "  - OrthoANIu"
echo "  - EzAAI"
echo "  - KEGG Decoder"
echo "  - dbCAN"
echo "  - antiSMASH (manual)"

