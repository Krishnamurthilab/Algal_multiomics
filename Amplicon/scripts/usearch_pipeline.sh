#!/bin/bash
set -euo pipefail 
# =====================================================
# AMPLICON ANALYSIS PIPELINE
# MetaBiomelab
# =====================================================

# Software: USEARCH v11.0.667
#Update USEARCH path according to local installation
USEARCH=/path/to/usearch11.0.667_i86linux32

# RDP Classifier 2.13
# MEGAN CE version 6.24.22

# =====================================================
# A1. Paired-end read merging
# =====================================================

# Create output directory
mkdir -p merged

# Loop through all forward reads
for r1 in *_R1.fastq
do
    # Extract sample name
    sample=$(basename ${r1} _R1.fastq)
    
    echo "Merging sample: ${sample}"

    # Merge paired-end reads
    usearch -fastq_mergepairs ${sample}_R1.fastq \
             -reverse ${sample}_R2.fastq \
             -fastqout merged/${sample}_merged.fq

done

# =====================================================
# A2. Quality filtering
# =====================================================

mkdir -p filtered

for merged_file in merged/*_merged.fq
do
    sample=$(basename ${merged_file} _merged.fq)

    echo "Filtering sample: ${sample}"

    $USEARCH \
        -fastq_filter ${merged_file} \
        -fastq_trunclen 150 \
        -fastq_maxee 0.5 \
        -fastaout filtered/${sample}_filtered.fasta

done

# =====================================================
# A3. Chimera removal
# =====================================================

cat filtered/*_filtered.fasta > all_filtered.fasta

$USEARCH \
    -fastx_uniques all_filtered.fasta \
    -fastaout uniques.fasta \
    -sizeout \
    -relabel Uniq


# =====================================================
# A4. OTU clustering
# =====================================================

$USEARCH \
    -cluster_otus uniques.fasta \
    -otus otus.fasta \
    -relabel Otu

echo "Amplicon processing completed successfully!"

# =====================================================
# A5. OTU taxonomic classification using RDP
# =====================================================

# Update RDP classifier path according to local installation
RDP=/path/to/rdp_classifier.jar

echo "Assigning taxonomy to OTUs"

java -Xmx4g -jar $RDP \
    -c 0.5 \
    -o otu_taxonomy.txt \
    -h otu_taxonomy_hierarchy.txt \
    otus.fasta

# =====================================================
# A6. Generation of MEGAN-compatible file
# =====================================================

# Example:
# Commands/tools used for MEGAN import
