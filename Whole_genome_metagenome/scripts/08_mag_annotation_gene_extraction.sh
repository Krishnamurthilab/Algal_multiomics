set -euo pipefail

# =========================================================
# FUNCTIONAL ANNOTATION (BAKTA + BARRNAP)
# Input : high-quality MAGs (from DAS Tool)
# Output: annotated genomes + rRNA prediction
# =========================================================
cd Whole_genome_metagenome/results

echo "Starting functional annotation..."

THREADS=16

MAG_DIR="07_mag_refinement_quality/NCBI_FCS_GX_cleaned_high_quality_MAGs/cleaned_MAGs"
OUTDIR="08_mag_annotation"
mkdir -p ${OUTDIR}

# =========================================================
# STEP 1: BAKTA ANNOTATION
# =========================================================

echo "Running Bakta..."

mkdir -p ${OUTDIR}/bakta

for genome in ${MAG_DIR}/*.fa; do
  name=$(basename ${genome} .fa)

  bakta \
    --db /path/to/bakta_db \
    --threads ${THREADS} \
    --output ${OUTDIR}/bakta/${name} \
    ${genome} \
    > ${OUTDIR}/bakta/${name}.log 2>&1
done

# =========================================================
# STEP 2: BARRNAP (rRNA prediction)
# =========================================================

echo "Running Barrnap..."

mkdir -p ${OUTDIR}/barrnap

for genome in ${MAG_DIR}/*.fa; do
  name=$(basename ${genome} .fa)

  barrnap \
    --threads ${THREADS} \
    ${genome} \
    > ${OUTDIR}/barrnap/${name}.rRNA.gff
done
# =========================================================
# DONE
# =========================================================

echo "Functional annotation completed!"
echo "Outputs:"
echo "  - ${OUTDIR}/bakta/"
echo "  - ${OUTDIR}/barrnap/"
