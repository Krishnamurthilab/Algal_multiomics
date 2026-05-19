# Amplicon Analysis Workflow

This includes complete workflow information for 16S rRNA gene amplicon analysis using USEARCH, RDP_classifeir and MEGAN CE GUI tool.

## Pipeline overview
Workflow includes
### A1. Paired-end read merging
#### For samples with replicate barcode sequencing files, paired-end reads from individual replicates were merged independently using USEARCH prior to pooling and downstream processing.
usearch -fastq_mergepairs ${r1} \
         -reverse ${r2} \
         -fastqout merged.fq

### A2. Quality filtering
usearch -fastq_filter merged.fq \
         -fastq_trunclen 150 \
         -fastq_maxee 0.5 \
         -fastaout filtered.fasta

### A3. Chimera removal
usearch -fastx_uniques filtered.fasta \
         -fastaout uniques.fasta \
         -sizeout \
         -relabel Uniq

### A4. OTU clustering
usearch -cluster_otus uniques.fasta \
         -otus otus.fasta \
         -relabel Otu

### A5. Taxonomic classification
java -Xmx4g -jar rdp_classifier.jar \
     -c 0.5 \
     -o otu_taxonomy.txt \
     -h otu_taxonomy_hierarchy.txt \
     otus.fasta

### A6. Generation of MEGAN-compatible file
Taxonomic classification results were imported into MEGAN CE (version 6.24.22) for downstream visualization and analysis.

File → Import → Text (CSV) Format
### A7. Downstream taxonomic visualization and statistical analysis
Taxonomic abundance profiles generated in MEGAN Community Edition were used for downstream visualization, including:

stacked bar plots
bubble plots

Taxon abundance tables exported from MEGAN CE were subsequently analyzed using MicrobiomeAnalyst for:

alpha and beta diversity analysis
taxonomic profiling
co-occurrence network analysis
multivariate statistical analyses
graphical visualization
Others

## Software used

- USEARCH v11.0.667
-RDP Classifier v2.13
-MEGAN CE 6.24.22 (GUI)

##Example Repository Structure

amplicon/
amplicon/
│
├── scripts/
│   └── usearch_pipeline.sh
│
├── metadata/
│
├── results/
│   ├── merged/
│   ├── filtered/
│   ├── otus/
│   └── taxonomy/
│
└── README.md
# Default parameters were used unless otherwise specified in the associated manuscript.
