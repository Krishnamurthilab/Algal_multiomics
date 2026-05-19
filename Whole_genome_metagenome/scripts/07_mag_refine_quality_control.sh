
#!/bin/bash
set -euo pipefail

# =========================================================
# BIN REFINEMENT (DAS_Tool-1.1.7) + QUALITY CHECK (CheckM v1.2.2)
# Covergage estimation (coverm 0.7.0)
# Assembly statistics estimation (quast-5.3.0)
# Input : bins from MetaBAT2 + MaxBin2 + CONCOCT
# Output: refined MAGs + quality statistics
# =========================================================
cd Whole_genome_metagenome_results
echo "Starting DAS Tool + CheckM pipeline + coverm + QUAST..."

THREADS=16
SAMPLE="(basename)"

# =========================================================
# INPUT DIRECTORIES
# =========================================================

METABAT_DIR="05_metabat2/bins"
MAXBIN_DIR="04_maxbin2/bins"
CONCOCT_DIR="06_concoct/bins"

CONTIGS="filtered_contigs_WGM/*.fasta"

OUTDIR="07_mag_refinement_quality"
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
# STEP 5: HIGH-QUALITY MAGs SELECTION + FCS-GX CLEANING
# =========================================================

echo "Selecting high-quality MAGs (CheckM filtering)..."

CHECKM_TABLE="07_mag_refinement_quality/checkm_summary.tsv"

HQ_DIR="07_refinement_quality/high_quality_MAGs"
FCS_DIR="NCBI_FCS_GX_cleaned_high_quality_bins"
TRIMMED_DIR="${FCS_DIR}/trimmed_headers"

mkdir -p ${HQ_DIR}
mkdir -p ${FCS_DIR}
mkdir -p ${TRIMMED_DIR}

# ---------------------------------------------------------
# Step 5.1: Select HQ MAGs (Completeness >90%, Contamination <5%)
# ---------------------------------------------------------
## Bin having completeness >90%, and contamination <5% were manually selected and placed in HQ_DIR

# ---------------------------------------------------------
# Step 5.2: Remove FCS contamination + standardize headers
# ---------------------------------------------------------

echo "Running NCBI FCS-GX contamination screening..."

fcs_gx run \
	--in-dir ${HQ_DIR} \
	--out_dir ${FCS_DIR} \
	--threads 16

fcs_gx clean \	
	--in-dir ${HQ_DIR} \	
	--report ${FCS_DIR}/fcs_gx_report.tsv \
	--out-dir ${FCS_DIR}
echo "FCS-GX screening completed."

echo "Trimming headers after FCS-GX cleaning..."

for fasta in ${HQ_DIR}/*.fasta; do
  name=$(basename ${fasta} .fasta)

  # remove long/complex headers → simplify to bin name
  awk '/^>/ {print ">"FILENAME"_"NR; next} {print}' ${fasta} \
    > ${TRIMMED_DIR}/${name}.fasta

done

echo "HQ MAG selection and FCS-GX cleanup completed."


# =========================================================
# STEP 6: COVERM (FINAL COVERAGE PROFILE for high_quality_bins)
# =========================================================
CONTIGS="07_mag_refinement_quality/NCBI_FCS_GX_cleaned_high_quality_MAGS/cleaned_MAGs"
READ1="../data/reads/*_R1.fastq"
READ2="../data/reads/*_R2.fastq"
OUTDIR="07_mag_refinement_quality/coverm"
echo "Running CoverM..."

mkdir -p ${OUTDIR}/coverm

coverm contig \
  -1 ${READ1} \
  -2 ${READ2} \
  -r ${CONTIGS} \
  --methods mean \
  --output-format dense \
  -t ${THREADS} \
  > ${OUTDIR}/coverm/coverage.tsv


# =========================================================
# STEP 7: QUAST
# =========================================================
QUAST_INPUTS="07_mag_refinement_quality/NCBI_FCS_GX_cleaned_high_quality_MAGs/cleaned_MAGs"
OUTDIR="07_refinement_quality/
echo "Running QUAST..."

QUAST_OUT=${OUTDIR}/quast_output
mkdir -p ${QUAST_OUT}

quast.py \
  ${OUTDIR}/${QUAST_INPUTS}/*.fa \
  -o ${QUAST_OUT} \
  -t ${THREADS} \
  > ${OUTDIR}/quast_output/quast.log 2>&1

# =========================================================
# DONE
# =========================================================

echo "DAS Tool + CheckM + QUAST + NCBI_FCS_GX + coverm completed successfully!"
echo "Outputs:"
