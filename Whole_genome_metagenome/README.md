# Whole Genome Metagenomics and Comparative Genomics Pipeline
This repository complete workflow for whole genome metagenome analysis including community profiling, systematic analysis of metagenome assembled genomes, and functional ecology of microbial community.
## Overview

The workflow included:

* Whole metagenome assembly
* Community-level functional annotation
* Genome binning using multiple algorithms
* MAG refinement and quality assessment
* Contamination screening and genome curation
* Taxonomic classification of MAGs
* Functional annotation and gene extraction
* Comparative genomics and metabolic reconstruction

The workflow was designed for genome-resolved metagenomic analysis and integrates multiple complementary approaches to recover, refine, classify, and functionally characterize metagenome-assembled genomes (MAGs).

---

# Workflow Summary

```text
01_assembly.sh
02_whole_metagenome_functional_annotation.sh
03_mapping.sh
04_maxbin2.sh
05_metabat2.sh
06_concoct.sh
07_refine-quality-check.sh
08_mag_annotation_gene_extraction.sh
09_MAGs_taxonomy.sh
10_MAG_comparative_genomics.sh
```

---

# Workflow Description

# 01_assembly.sh
## Quality controll of whole genome metagenome raw sequence
### QC check using fastqc v0.11.9

fastqc reads/*fastq.gz -o fastqc_output

### Quality filtering using Trimmomatic

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

## Denovo assembly
Performs de novo assembly of quality-filtered metagenomic reads to generate assembled contigs for downstream metagenomic analyses.

### Input

* Paired-end Illumina reads
* Quality-filtered FASTQ files

### Major Steps

* Metagenome assembly
* Contig filtering

### Typical Tools

* metaSPAdes v3.13.1

### script used
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

#### Output

* Assembled contigs (`*.fasta`)

### Filtering short contigs using Bioawk v20110810
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

#### Output
* Filtered contig assemblies

## Megan compatible data generation
### read alignment
diamond blastx -d nr -q results/filtered_contigs_WGM/*fasta -o results/megan/*{basename}.daa -f 100

### Data Meganization
daa-meganizer -i /results/filtered_contigs_WGM/*{basename}.daa -mdb megan-map-feb2022.db

### Data import in MEGAN 
File → Import → Text (CSV) Format

### Downstream taxonomic visualization and statistical analysis
*Taxonomic abundance profiles generated in MEGAN Community Edition were used for downstream visualization, including:

* stacked bar plots
* bubble plots

*Taxon abundance and tables exported from MEGAN CE were subsequently analyzed using MicrobiomeAnalyst for:

* alpha and beta diversity analysis
* taxonomic profiling
* co-occurrence network analysis
* multivariate statistical analyses
* graphical visualization
* Others

 
------

## 02_whole_metagenome_functional_annotation.sh

### Purpose

Performs community-level functional profiling directly on whole metagenome assemblies independent of MAG recovery.

This analysis evaluates the overall metabolic and biosynthetic potential of the microbial community.

### Analyses Included

### DIAMOND + MEGAN

* Functional annotation of metagenomic reads and contigs
* community-level taxonomic and functional profiling

### BlastKOALA

* Manual KEGG web server annotation
* KO (KEGG Orthology) assignment generation

### KEGG Decoder

* Metabolic pathway reconstruction
* Comparative pathway completeness analysis
* Static heatmap generation using combined KO profiles

KEGG-decoder \
  -i ${OUTDIR}/kegg/blastkoala/combined_ko.tsv \
  -o ${OUTDIR}/kegg/kegg_decoder_output.tsv \
  -v static \
  > ${OUTDIR}/kegg/kegg_decoder.log 2>&1

### dbCAN based cazyme profiling

* Carbohydrate-active enzyme (CAZyme) annotation

run_dbcan \
  assembly_proteins.faa \
  protein \
  --out_dir ${OUTDIR}/dbcan \
  > ${OUTDIR}/dbcan.log 2>&1

### antiSMASH server based analysis of biosynthetic gene cluster

* Biosynthetic gene cluster prediction
* Secondary metabolite biosynthetic potential assessment

## Workflow strategy 


Whole metagenome assembly
        ↓
DIAMOND alignment
        ↓
MEGAN taxonomic and functional annotation
        ↓
Community-level profiling
        ↓
BLASTKOALA annotation
        ↓
Combined KO matrix generation
        ↓
KEGG Decoder analysis
        ↓
Metabolic heatmap reconstruction


### Important Notes

#### BlastKOALA

BlastKOALA annotation was performed manually using the KEGG web server.

Expected output files from BlastKOALA should be placed in:


results/02_whole_metagenome_functional_annotation/kegg/blastkoala/


#### antiSMASH

antiSMASH analyses were performed using the antiSMASH web server.

### Output

* DIAMOND alignment files
* MEGAN taxonomic profiles
* Functional annotation tables
* Community composition summaries
* KO annotation tables
* KEGG metabolic pathway matrices
* Metabolic heatmaps
* CAZyme annotations
* Biosynthetic gene cluster predictions

---

# 03_mapping.sh

## Purpose

Maps metagenomic reads back to assembled contigs and calculates coverage depth required for metagenomic binning.

## Tools

* Bowtie2
* SAMtools
* CoverM
* jgi_summarize_bam_contig_depths

## Workflow

### Step 1: Bowtie2 index generation

* Contig indexing for read mapping

bowtie2-build ${CONTIGS} ${OUTDIR}/contigs_index \
  > ${OUTDIR}/bowtie2_build.log 2>&1


### Step 2: Read mapping

* Paired-end read alignment to assembled contigs

bowtie2 \
  -x ${OUTDIR}/contigs_index \
  -1 ${READ1} \
  -2 ${READ2} \
  -p ${THREADS} \
  | samtools view -bS - \
  > ${OUTDIR}/mapped.bam

### Step 3: SAM to BAM conversion

* Conversion of alignment files into compressed BAM format

samtools sort -@ ${THREADS} \
  -o ${OUTDIR}/mapped.sorted.bam \
  ${OUTDIR}/mapped.bam

rm ${OUTDIR}/mapped.bam

samtools index ${OUTDIR}/mapped.sorted.bam


### Step 4: Coverage estimation

* Contig depth estimation
* Coverage calculation

jgi_summarize_bam_contig_depths \
  --outputDepth ${OUTDIR}/depth.txt \
  ${OUTDIR}/mapped.sorted.bam


## Output

* SAM files
* BAM files
* Sorted BAM files
* Coverage profiles
* Depth tables

---

# 04_maxbin2.sh

## Purpose

Performs metagenome binning using MaxBin2.

## Tool

* MaxBin2

## Workflow

* Uses contig composition
* Uses abundance/depth information
* Generates genome bins from assembled contigs

## Input

* Filtered contigs
* Coverage depth files

## Used scipt
run_MaxBin.pl \
  -contig ${CONTIGS} \
  -abund ${DEPTH} \
  -out ${OUTDIR}/${PREFIX}_maxbin2 \
  -thread ${THREADS} \
  -min_contig_length 1500 \
  > ${OUTDIR}/maxbin2.log 2>&1

## Output

* Genome bins (`*.fasta`)
* MaxBin2 output files

---

# 05_metabat2.sh

## Purpose

Performs metagenome binning using MetaBAT2.

## Tool

* MetaBAT2

## Workflow

* Coverage-based binning
* Tetranucleotide frequency-based clustering
* Genome reconstruction from metagenomic contigs

## Input

* Filtered contigs
* Contig depth profiles

## Used scripts

metabat2 \
  -i ${CONTIGS} \
  -a ${DEPTH} \
  -o ${OUTDIR}/${PREFIX}_bin \
  -t ${THREADS} \
  --minContig ${MIN_CONTIG} \
  --unbinned \
  > ${OUTDIR}/metabat2.log 2>&1

## Output

* MetaBAT2 bins
* Unbinned contigs

---

# 06_concoct.sh

## Purpose

Performs metagenome binning using CONCOCT.

## Tools

* CONCOCT
* cut_up_fasta.py
* concoct_coverage_table.py
* merge_cutup_clustering.py
* extract_fasta_bins.py

## Workflow

### Step 1: Contig fragmentation

* Long contigs fragmented into smaller windows
cut_up_fasta.py \
  ${CONTIGS} \
  -c ${CHUNK_SIZE} \
  -o 0 \
  --merge_last \
  -b ${OUTDIR}/cuts.bed \
  > ${OUTDIR}/cuts.fa


### Step 2: Coverage table generation

* Coverage profile calculation from BAM files

concoct_coverage_table.py \
  ${OUTDIR}/cuts.bed \
  ${BAM} \
  > ${OUTDIR}/coverage.tsv

### Step 3: CONCOCT clustering

* Clustering based on sequence composition and abundance

concoct \
  --composition_file ${OUTDIR}/cuts.fa \
  --coverage_file ${OUTDIR}/coverage.tsv \
  -b ${OUTDIR}/concoct_output \
  -t ${THREADS} \
  > ${OUTDIR}/concoct.log 2>&1

### Step 4: Cluster merging

* Reconstruction of original contig clustering
merge_cutup_clustering.py \
  ${OUTDIR}/concoct_output/clustering_gt1000.csv \
  > ${OUTDIR}/concoct_output/clustering_merged.csv

### Step 5: Bin extraction

* Extraction of genome bins

extract_fasta_bins.py \
  ${CONTIGS} \
  ${OUTDIR}/concoct_output/clustering_merged.csv \
  --output_path ${OUTDIR}/bins

## Output

* CONCOCT genome bins
* Clustering files
* Coverage matrices

---

# 07_refine-quality-check.sh

## Purpose

Performs MAG refinement, quality assessment, contamination screening, and genome curation.

This step integrates outputs from multiple binning tools to generate refined high-quality MAGs.

## Analyses Included

### DAS Tool

* Consensus bin refinement
* Integration of MetaBAT2, MaxBin2, and CONCOCT bins

DAS_Tool \
  -i ${OUTDIR}/dastool_inputs/metabat.scaffolds2bin.tsv,${OUTDIR}/dastool_inputs/maxbin.scaffolds2bin.tsv,${OUTDIR}/dastool_inputs/concoct.scaffolds2bin.tsv \
  -l metabat2,maxbin2,concoct \
  -c ${CONTIGS} \
  -o ${OUTDIR}/${SAMPLE}_DAS \
  --search_engine diamond \
  --threads ${THREADS} \
  > ${OUTDIR}/dastool.log 2>&1

### CoverM

* Coverage estimation
* MAG abundance profiling
coverm contig \
  -1 ${READ1} \
  -2 ${READ2} \
  -r ${CONTIGS} \
  --methods mean \
  --output-format dense \
  -t ${THREADS} \
  > ${OUTDIR}/coverm/coverage.tsv

### CheckM

* MAG completeness estimation
* Contamination estimation
* Lineage-specific marker analysis

checkm lineage_wf \
  -x fa \
  -t ${THREADS} \
  ${OUTDIR}/${SAMPLE}_DAS_DASTool_bins \
  ${CHECKM_OUT} \
  > ${OUTDIR}/checkm.log 2>&1

### QUAST

* Assembly quality statistics
* Genome assembly metrics

quast.py \
  ${OUTDIR}/${QUAST_INPUTS}/*.fa \
  -o ${QUAST_OUT} \
  -t ${THREADS} \
  > ${OUTDIR}/quast_output/quast.log 2>&1

### NCBI FCS-GX

* Foreign contamination screening
* Removal of contaminant sequences
* Genome cleaning and curation

fcs_gx run \
	--in-dir ${HQ_DIR} \
	--out_dir ${FCS_DIR} \
	--threads 16

fcs_gx clean \	
	--in-dir ${HQ_DIR} \	
	--report ${FCS_DIR}/fcs_gx_report.tsv \
	--out-dir ${FCS_DIR}

## MAG Selection Criteria

High-quality MAGs were selected based on:

| Metric        | Threshold |
| ------------- | --------- |
| Completeness  | > 90%     |
| Contamination | < 5%      |

## Workflow Logic

```text
MetaBAT2 bins
MaxBin2 bins
CONCOCT bins
        ↓
DAS Tool refinement
        ↓
CheckM quality assessment
        ↓
High-quality MAG selection
        ↓
NCBI FCS-GX contamination screening
        ↓
Cleaned high-quality MAGs
```

## Output

* Refined MAGs
* High-quality MAGs
* CheckM reports
* QUAST statistics
* Cleaned genomes
* Coverage profiles

---

# 08_mag_annotation_gene_extraction.sh

## Purpose

Performs MAG annotation, protein prediction, and gene extraction.

## Tools

* Bakta
* Barrnap

## Workflow

### Bakta annotation

* Protein coding gene prediction
* Functional annotation
* Protein FASTA generation

for genome in ${MAG_DIR}/*.fa; do
  name=$(basename ${genome} .fa)

  bakta \
    --db /path/to/bakta_db \
    --threads ${THREADS} \
    --output ${OUTDIR}/bakta/${name} \
    ${genome} \
    > ${OUTDIR}/bakta/${name}.log 2>&1
done

### Barrnap

* rRNA gene identification
* Ribosomal RNA annotation

for genome in ${MAG_DIR}/*.fa; do
  name=$(basename ${genome} .fa)

  barrnap \
    --threads ${THREADS} \
    ${genome} \
    > ${OUTDIR}/barrnap/${name}.rRNA.gff
done

## Output

* Annotated MAGs
* Protein FASTA files (`*.faa`)
* Nucleotide gene files
* rRNA annotations
* Functional annotation tables

---

# 09_MAGs_taxonomy.sh

## Purpose

Performs taxonomic classification and phylogenomic analysis of recovered MAGs.

## Tools

* GTDB-Tk
* CAT pack
* UBCG

## GTDB-Tk

Performs genome-based taxonomic classification using the Genome Taxonomy Database.

### Important Configuration

GTDB-Tk database path must be configured before running:

```bash
export GTDBTK_DATA_PATH=/path/to/gtdbtk_database
```
gtdbtk classify_wf \
  --genome_dir ${MAG_DIR} \
  --out_dir ${OUTDIR}/gtdbtk_out \
  --cpus ${THREADS} \
  --mash_db ${GTDB_DATA_PATH} \
  > ${OUTDIR}/gtdbtk_out/gtdbtk.log 2>&1

## CAT pack

Performs contig and bin-level taxonomic classification.

CAT_pack bins \
  -b ${MAG_DIR}/ \
  -d ${CAT_gtdb_db} \
  -t ${CAT_gtdb_tax} \
  -o ${OUTDIR}/cat_bins_output \
  -n ${THREADS}

## UBCG Phylogenomics

Phylogenomic analysis was performed in two stages:

### Stage 1: UBCG profile generation

* Conserved bacterial core gene extraction
* `.ucg` profile generation

for genome in ${GENOME_DIR}/*.fasta; do

  base=$(basename ${genome} .fasta)

  ubcg.py \
    -i ${genome} \
    -o ${UBCG_OUT}/ucg_profiles/${base} \
    -t ${THREADS}

done

### Stage 2: Phylogenomic tree construction

* Concatenated marker alignment
* Phylogenetic tree inference

build_ucg_tree.py \
  -i ${UBCG_OUT}/ucg_profiles \
  -o ${UBCG_OUT}/phylogeny \
  -t ${THREADS}

## Workflow Logic

MAGs
  ↓
GTDB-Tk classification
  ↓
CAT pack taxonomy
  ↓
UBCG marker extraction (.ucg)
  ↓
Phylogenomic tree construction


## Output

* GTDB taxonomy assignments
* CAT taxonomy outputs
* UBCG profiles
* Phylogenomic trees

---

# 10_MAG_comparative_genomics.sh

## Purpose

Performs comparative genomic and functional analyses across recovered MAGs.

## Analyses Included

### OrthoANIu

* Average Nucleotide Identity (ANI)
* Genome similarity assessment
* Species-level genomic comparison

java -jar /path/to/bi-orthoaniu.jar \
  -fd ${MAG_DIR} \
  -n 16 \
  -o ${OUTDIR} \
  -t ${OUTDIR}/tmp \
  -u /path/to/usearch

### EzAAI

* Average Amino Acid Identity (AAI)
* Genus-level evolutionary relationships
* Hierarchical clustering

ezaai extract \
  -i ${MAG_DIR} \
  -o ${OUTDIR}/db \
  -l ${OUTDIR}/labels.tsv

ezaai calculate \
  -i ${OUTDIR}/db \
  -j ${OUTDIR}/db \
  -o ${OUTDIR}/out/aai.tsv

### BlastKOALA

* KO assignment generation
* KEGG functional annotation

### KEGG Decoder

* Metabolic reconstruction
* Pathway completeness estimation
* Comparative metabolic heatmap generation

cat ${KEGG_OUT}/blastkoala/*.tsv \
  > ${KEGG_OUT}/combined_ko.tsv

KEGG-decoder \
  -i ${KEGG_OUT}/combined_ko.tsv \
  -o ${KEGG_OUT}/kegg_decoder_output.tsv \
  -v static \
  > ${KEGG_OUT}/kegg_decoder.log 2>&1

### dbCAN

* CAZyme annotation
* Carbohydrate metabolism profiling
* Individual MAG-level CAZyme analysis

for faa in ${FAA_DIR}/*.faa; do

  base=$(basename ${faa} .faa)

  echo "Processing ${base}..."

  run_dbcan \
    ${faa} \
    protein \
    --out_dir ${DBCAN_OUT}/${base} \
    > ${DBCAN_OUT}/${base}.log 2>&1

done

### antiSMASH

* Biosynthetic gene cluster prediction
* Secondary metabolite biosynthetic potential analysis

---

# Detailed Comparative Genomics Workflow

## OrthoANIu Workflow

ANI analysis was performed using OrthoANIu v1.2 in directory mode.

### Workflow

```text
MAG directory
      ↓
OrthoANIu (-fd mode)
      ↓
Pairwise ANI matrix
```

### Important Notes

* USEARCH dependency required
* ANI matrix generated using directory mode (`-fd`)
* Matrix output used for comparative genomic analyses

---

## EzAAI Workflow

AAI analysis was performed in three stages.

### Stage 1: CDS profile extraction

* Prodigal-based CDS prediction
* `.db` profile generation

### Stage 2: Pairwise AAI calculation

* MMSeqs2-based AAI estimation
* All-vs-all genome comparisons

### Stage 3: Hierarchical clustering

* UPGMA clustering
* Newick tree generation

## Workflow Logic

```text
MAG genomes
    ↓
EzAAI extract
    ↓
CDS profile databases (.db)
    ↓
EzAAI calculate
    ↓
AAI matrix
    ↓
EzAAI cluster
    ↓
Hierarchical tree
```

---

## KEGG Functional Reconstruction Workflow

### Step 1: BlastKOALA annotation

* Manual KEGG web server annotation
* KO assignment generation

### Step 2: Combined KO matrix generation

* KO annotations merged across all MAGs

### Step 3: KEGG Decoder analysis

* Pathway reconstruction
* Pathway completeness estimation
* Heatmap generation using `-v static`

## Workflow Logic

```text
MAG proteins
      ↓
BlastKOALA annotation
      ↓
Combined KO matrix
      ↓
KEGG Decoder
      ↓
Metabolic heatmap
```

---

## dbCAN Workflow

CAZyme annotation was performed individually for each MAG.

### Workflow

```text
Individual MAG protein FASTA
          ↓
run_dbcan
          ↓
CAZyme annotation
```

### Important Notes

* Protein mode was used
* Individual MAG profiling performed in loop-based workflow
* Output generated separately for each MAG

---

## antiSMASH Workflow

antiSMASH analysis was performed manually using the antiSMASH web server.

### Output

* Biosynthetic gene clusters
* Secondary metabolite predictions

---

# Software Requirements

# Core Tools

## Assembly
* metaSPAdes

## Mapping and Coverage

* Bowtie2
* SAMtools
* CoverM

## Binning

* MetaBAT2
* MaxBin2
* CONCOCT
* DAS Tool

## Quality Assessment

* CheckM
* QUAST
* NCBI FCS-GX

## Annotation

* Bakta
* Barrnap

## Taxonomy and Phylogeny

* GTDB-Tk
* CAT pack
* UBCG

## Comparative Genomics

* OrthoANIu
* EzAAI
* KEGG Decoder
* dbCAN
* antiSMASH

## Whole Metagenome Taxonomic and Functional Annotation

* DIAMOND
* MEGAN

---

# Repository Structure

```text
Whole_genome_metagenome/
├── README.md
├── data/
│   ├── reads/
│   ├── trimmed_reads/
│   ├── fastqc_output/
│   └── filtered_contigs/
│
├── scripts/
│   ├── 01_assembly.sh
│   ├── 02_whole_metagenome_functional_annotation.sh
│   ├── 03_mapping.sh
│   ├── 04_maxbin2.sh
│   ├── 05_metabat2.sh
│   ├── 06_concoct.sh
│   ├── 07_mag_refine_quality_control.sh
│   ├── 08_mag_annotation_gene_extraction.sh
│   ├── 09_MAGs_taxonomy.sh
│   └── 10_MAG_comparative_genomics.sh
│
├── results/
│   ├── 02_whole_metagenome_functional_analysis/
│   │   ├── antismash/
│   │   ├── dbcan/
│   │   └── kegg/
│   │
│   ├── 03_mapping/
│   ├── 04_maxbin2/
│   ├── 05_metabat2/
│   ├── 06_concoct/
│   ├── 07_mag_refinement_quality/
│   │
│   ├── 08_mag_annotation/
│   │   ├── bakta/
│   │   └── barrnap/
│   │
│   ├── 09_mag_taxonomy/
│   │   ├── MAG_ref_genomes/
│   │   ├── cat_bins_output/
│   │   ├── gtdbtk_out/
│   │   └── ubcg/
│   │
│   ├── 10_mag_comparative_genomics/
│   ├── megan/
│   └── microbiome_analyst/
│
└── .gitignore
---

# General Notes

## Manual Web Server Analyses

The following analyses were performed manually using web servers:

| Tool            | Purpose                                              |
| --------------- | ---------------------------------------------------- |
| DIAMOND + MEGAN | Whole metagenome taxonomic and functional annotation |
| BlastKOALA      | KEGG orthology annotation                            |
| antiSMASH       | Biosynthetic gene cluster prediction                 |

---

## Recommended Citation

Please cite the original publications of the following tools if using this workflow:

* MetaBAT2
* MaxBin2
* CONCOCT
* DAS Tool
* CheckM
* GTDB-Tk
* CAT pack
* UBCG
* OrthoANIu
* EzAAI
* KEGG Decoder
* dbCAN
* antiSMASH
* Bakta
* CoverM
* Bowtie2

---

# Contact

For workflow-related issues, suggestions, or troubleshooting, please open an issue in this repository.
