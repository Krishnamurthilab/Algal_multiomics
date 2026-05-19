
#!/bin/bash

set -euo pipefail

# ==========================================
# Software paths
# ==========================================

FASTQC=/path/to/fastqc
TRIMMOMATIC=/path/to/trimmomatic.jar
SPADES=/path/to/spades.py
DIAMOND=/path/to/diamond

# ==========================================
# B1. Quality assessment
# ==========================================
## move to directory Whole_genome_metagenome/data 
cd Whole_genome_metagenome/data
mkdir -p fastqc_output

fastqc reads/*fastq.gz -o fastqc_output


# ==========================================
# B2. Quality filtering using Trimmomatic
# ==========================================

mkdir -p trimmed_reads

for r1 in reads/*_R1.fastq
do
    sample=$(basename ${r1} _R1.fastq)

    echo "Running Trimmomatic for sample: ${sample}"

    java -jar $TRIMMOMATIC PE \
        ${sample}_R1.fastq.gz \
        ${sample}_R2.fastq.gz \
        trimmed_reads/${sample}_F_paired.fastq \
        trimmed_reads/${sample}_F_unpaired.fastq \
        trimmed_reads/${sample}_R_paired.fastq \
        trimmed_reads/${sample}_R_unpaired.fastq \
        ILLUMINACLIP:all_adapter.fa:2:30:10 \
        LEADING:33 \
        TRAILING:33 \
        MINLEN:50 \
	CROP:240 \
	HEADCROP:20 \

done
# ==========================================
# B3. Assembly using MetaSPAdes
# ==========================================
## move to directory Whole_genome_metagenome
cd Whole_genome_metagenome/results
mkdir -p whole_genome_metagenome_assembly

for r1 in ../data/trimmed_reads/*_F_paired.fastq 
do
    sample=$(basename ${r1} _F_paired.fastq)

    echo "Running MetaSPAdes assembly for sample: ${sample}"

    $METASPADES \
        -1 ../data/trimmed_reads/${sample}_F_paired.fastq \
        -2 ../data/trimmed_reads/${sample}_R_paired.fastq \
        -s ../data/trimmed_reads/${sample}_F_unpaired.fastq \
        -s ../data/trimmed_reads/${sample}_R_unpaired.fastq \
        -o whole_genome_metagenome_assembly/${sample}

done

# ==========================================
# B4. Removal of short contigs using Bioawk
# ==========================================

mkdir -p filtered_contigs_WGM

for assembly_dir in whole_genome_metagenome_assembly/*
do
    sample=$(basename ${assembly_dir})

    echo "Filtering contigs for sample: ${sample}"

    $BIOAWK -c fastx \
    'length($seq) > 1000 { print ">"$name; print $seq }' \
    ${assembly_dir}/final.contigs.fa \
    > filtered_contigs/${sample}_filtered_contigs.fasta

done

echo "Shotgun metagenome preprocessing completed successfully!"

# ====================================
# MEGAN + DIAMOND tool based analysis
# ====================================
echo "read alignemnt using blastx tool"

diamond blastx -d nr -q results/filtered_contigs_WGM/*fasta -o results/megan/*{basename}.daa -f 100

echo "read meganization"

daa-meganizer -i /results/filtered_contigs_WGM/*{basename}.daa -mdb megan-map-feb2022.db

echo "manual import of meganized data into megan CE tool"
echo "OTUs and functional gene table export"


# ================================================
# OTUs table import into microbiome analyst server
# =================================================

echo " Alpha and Beta diversity analysis"
echo " Co-occurance network plot analysis"
echo " statistical analysis"