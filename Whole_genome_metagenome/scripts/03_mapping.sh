#!/bin/bash
set -euo pipefail

# =========================================================
# READ MAPPING + COVERAGE GENERATION
# Required for: MetaBAT2 / CONCOCT binning
# =========================================================
cd Whole_genome_metagenome/results
echo "Starting read mapping step..."

# -------------------------
# INPUTS (edit or use config file later)
# -------------------------
CONTIGS="filtered_contigs_WGM/*.fa"
READ1="../data/reads/*_R1.fastq.gz"
READ2="../data/reads/*_R2.fastq.gz"

OUTDIR="03_mapping"
THREADS=16

mkdir -p ${OUTDIR}

# -------------------------
# 1. INDEX ASSEMBLY
# -------------------------
echo "Building Bowtie2 index..."
bowtie2-build ${CONTIGS} ${OUTDIR}/contigs_index \
  > ${OUTDIR}/bowtie2_build.log 2>&1

# -------------------------
# 2. MAP READS
# -------------------------
echo "Mapping reads..."
bowtie2 \
  -x ${OUTDIR}/contigs_index \
  -1 ${READ1} \
  -2 ${READ2} \
  -p ${THREADS} \
  | samtools view -bS - \
  > ${OUTDIR}/mapped.bam

# -------------------------
# 3. SORT BAM
# -------------------------
echo "Sorting BAM..."
samtools sort -@ ${THREADS} \
  -o ${OUTDIR}/mapped.sorted.bam \
  ${OUTDIR}/mapped.bam

rm ${OUTDIR}/mapped.bam

# -------------------------
# 4. INDEX BAM
# -------------------------
echo "Indexing BAM..."
samtools index ${OUTDIR}/mapped.sorted.bam

# -------------------------
# 5. DEPTH FILE (MetaBAT2 input)
# -------------------------
echo "Generating depth file..."
jgi_summarize_bam_contig_depths \
  --outputDepth ${OUTDIR}/depth.txt \
  ${OUTDIR}/mapped.sorted.bam

echo "Mapping step completed successfully"
echo "Outputs:"
echo "  - ${OUTDIR}/mapped.sorted.bam"
echo "  - ${OUTDIR}/depth.txt"
