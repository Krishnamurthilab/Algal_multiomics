** NOTE: This repository only provides pipeline and workflows implementated in the current research to analyse the multiomics datasets. The tools and softawre belong to the original developers cited in the manuscript.

# Software and Tool Usage Notice

This repository does not distribute or claim ownership of the third-party software and bioinformatics tools used in the analyses.

All analyses were performed using freely available open-source tools, packages, and publicly available software developed by their respective authors and communities.

This repository primarily provides:

- workflow organization
- analysis scripts
- parameter settings
- reproducible analytical pipelines
- downstream statistical analyses

Users should cite the original software/tools/packages according to their respective publications and licenses.

# Major Software and Packages Used
- USEARCH v11.0.667
- RDP classifier v2.13
- MEGAN 6.24.22
- MetaboAnalyst (webserver)
- MicrobiomeAnalyst (webserver)
- FASTQC v0.11.9
- Trimmomatic v0.39
- Metaspades v3.13.1
- DIAMOND + MEGAN tool v3.0
- MaxBin2 v2.2
- MetaBAT2 v1.7
- GTDB-Tk v2.4.1
- Concoct v1.1
- das v1.1.2
- checkM v1.0.18
- QUAST v5.0.2
- mixOmics
- coverM
- barrnap
- Bakta
- vegan
- ggplot2
- pheatmap


# Algal Multi-omics Analysis

## Overview

This repository contains the workflows, scripts, and analytical approaches used for multi-omics characterization of marine algal-associated microbial communities and metabolomes.

The study integrates:

- Amplicon sequencing analysis
- Metagenomic functional profiling
- LC-MS based metabolome analysis
- Multi-omics integration analysis

to investigate host-associated ecological patterns, microbial functional potential, and metabolite diversity across marine algal systems.

---

# Study Objectives

The primary objectives of this study were:

- Characterize algal-associated microbial communities
- Explore metagenomic functional diversity
- Identify ecologically relevant metabolites
- Investigate relationships between microbial functions and metabolite profiles
- Understand species- and location-specific ecological signatures
- Integrate metagenomic and metabolomic datasets using multivariate statistical approaches

---

# Workflow Overview

```text
Marine algal samples
          ↓
Amplicon sequencing
          ↓
Microbial diversity analysis
          ↓
Shotgun metagenome sequencing
          ↓
Functional annotation and profiling
          ↓
LC-MS metabolome profiling
          ↓
Manual metabolite curation and annotation
          ↓
Statistical analyses
    ├── Diversity analysis
    ├── NMDS ordination
    ├── envfit analysis
    ├── PLS-DA
    └── PERMANOVA
          ↓
Multi-omics integration
    ├── DIABLO
    ├── Procrustes analysis
    ├── Co-inertia / CCA
    └── Correlation analysis
```

---

# Repository Structure

```text
Algal_multiomics/
├── Amplicon/
├── Metagenome/
├── Metabolome/
├── Integration/
├── README.md
└── LICENSE
```

---

# 1. Amplicon Analysis

Amplicon sequencing analysis was performed to investigate microbial community composition associated with marine algal samples.

## Analyses Included

- Quality filtering
- Denoising
- Taxonomic assignment
- Alpha diversity
- Beta diversity
- Ordination analysis

---

# 2. Metagenome Analysis

Shotgun metagenomic analysis was performed to explore microbial functional potential and ecological adaptation.

## Analyses Included

- Quality control
- Assembly
- Functional annotation
- Taxonomic profiling
- COG analysis
- Metabolic pathway analysis
- Comparative functional profiling

---

# 3. Metabolome Analysis

LC-MS based metabolome profiling was performed to investigate algal-associated metabolite diversity.

## Metabolite Curation

The metabolite dataset was manually curated to remove:

- Duplicate metabolites
- Unmatched compounds
- Synthetic drug-related metabolites
- Human and animal-associated metabolites

## Metabolite Annotation Strategy

Each metabolite was manually investigated using:

- PubChem
- HMDB
- KEGG
- ChEBI
- MetaboLights
- LipidMaps
- Published literature

to determine:

- Biological source
- Chemical nature
- Functional and ecological relevance

Approximately 300 ecologically relevant metabolites associated with marine, algal, plant, and microbial systems were retained for downstream analyses.

## Statistical Analyses

The metabolome analysis included:

- Data normalization
- Bray-Curtis dissimilarity
- NMDS ordination
- Environmental fitting (`envfit`)
- Convex hull analysis
- PCA biplot analysis
- PLS-DA
- PERMANOVA

---

# 4. Multi-omics Integration Analysis

Integrated analysis of metagenomic and metabolomic datasets was performed to identify relationships between microbial functions and metabolite profiles.

## Integration Approaches

### DIABLO Analysis

Performed using `mixOmics` to:

- identify correlated features
- detect discriminative biomarkers
- visualize multi-omics relationships

### Procrustes Analysis

Used to compare ordination similarity between metagenomic and metabolomic datasets.

### Co-inertia / CCA Analysis

Performed to identify shared ecological gradients and associations between datasets.

### Additional Analyses

- Heatmap analysis
- Correlation analysis
- PLS-DA
- PERMANOVA

---

# Software and Tools

## Programming Languages

- R
- Bash

## Major R Packages

```r
vegan
ggplot2
mixOmics
ade4
reshape2
ggrepel
pheatmap
dplyr
```

---

# Key Statistical Methods

The study utilized several multivariate and ecological statistical approaches, including:

- Bray-Curtis dissimilarity
- NMDS
- PCA
- PLS-DA
- PERMANOVA
- envfit
- DIABLO
- Procrustes analysis
- CCA / Co-inertia analysis

---

# Biological Significance

The integrated analyses provide insights into:

- Algal-associated microbial ecology
- Functional adaptation of microbiomes
- Host-associated metabolite diversity
- Species-specific ecological signatures
- Microbe-metabolite interactions
- Marine ecological metabolomics

---

# Data Availability

Input datasets include:

- Metagenome abundance matrices
- Metabolome abundance matrices
- Metadata tables

All analyses were performed using curated datasets and reproducible R workflows provided in this repository.

---

Algal_multiomics/
├── README.md
├── LICENSE
│
├── Amplicon/
│   ├── README.md
│   ├── data/
│   ├── scripts/
│   │   ├── usearch_pipeline.sh
│   │
│   └── results/
│       ├── filtered/
│       ├── MEGAN/
│       ├── merged/
│       ├── microbiome_analyst/
│       ├── taxonomy
│       ├── Unique/
│       └── otus/
│
├── Metagenome/
│   ├── README.md
│   ├── data/
│   │   ├── reads/
│   │   ├── fastqc_output
│   │   ├──filtered_contigs
│   │   └── trimmed_reads
│   │
│   ├── scripts/
│   │   ├── 01_assembly.sh
│   │   ├── 02_whole_metagenome_functional_annotation.sh
│   │   ├── 03_mapping.sh
│   │   ├── 04_maxbin2.sh
│   │   ├── 05_metabat2.sh
│   │   ├── 06_concoct.sh
│   │   ├── 07_mag_refine_quality_control.sh
│   │   ├── 08_mag_annotation_gene_extraction.sh
│   │   ├── 09_MAGs_taxonomy.sh
│   │   └── 10_MAG_comparative_genomics.sh
│   │
│   └── results/
│       ├── filtered_contigs_WGM/
│       ├── megan/
│       ├── microbiome_analyst/
│       ├── whole_genome_metagenome_assembly/
│       ├── 02_whole_metagenome_functional_analysis/
│       │    ├── antismash/
│       │    ├── dbcan/
│       │    └── kegg/
│       │        └── blastkoala/
│       ├── whole_genome_metagenome_assembly/
│       ├── 03_mapping/
│       ├── 04_maxbin2/
│       │     └── bins/
│       │  
│       ├── 05_metabat2/
│       │     └── bins/
│       │  
│       ├── 06_concoct/
│       │     ├── bins/
│       │     └── concoct_output/
│       │  
│       ├── 07_mag_refinement_quality/
│       │   ├── checkm_output
│       │   ├── coverm
│       │   ├── dastool_inputs
│       │   ├── high_quality_MAGs
│       │   ├── NCBI_FCS_GX_cleaned_high_quality_MAGs
│       │   │    └── cleaned_MAGs
│       │   │  
│       │   └── quast_output
│       │  
│       ├── 08_mag_annotation/
│       │    ├── bakta/
│       │    └── barnap/
│       │  
│       ├── 09_mag_taxonomy/
│       │    ├── cat_bins_output/
│       │    ├── gtdbtk_out/
│       │    ├── MAG_ref_genomes/
│       │    └── ubcg/
│       │  
│       └── 10_mag_comparative_genomics/
│            ├── dbcan/
│            ├── ezaai/
│            ├── kegg_decoder/
│            └── orthoani/
│                 └── temp/
│
├── Metabolome/
│   ├── metabolome_README.md
│   ├── data/
│   │   ├── filtered_metabolome.csv
│   │   └── metadata.csv
│   │
│   ├── scripts/
│   │    └── metabolome_analysis.R
│   │ 
│   ├── R_analysis/
│   │   ├── NMDS_output_final/
│   │  
│   └── metaboanalyst
│   
└── Integration/
    ├── multiomics_integration_README.md
    ├── multiomics_integration.R
    │
    ├── data/
    │   ├── metagenome.csv
    │   ├── metabolome.csv
    │   └── metadata.csv
    │
    └── Integration_Output/
        ├── CCA/
        ├── DIABLO/
        ├── Heatmap/
        ├── PERMANOVA/
        ├── PLSDA/
        ├── Procrustes/


# Citation

If using this repository or workflow, please cite:

Kumar et al. — Algal multi-omics analysis workflow

---

# Authors

Pravin Kumar1, Shiva Sundharam S1,2, Gaurav Tripathi3, Manisha Yadav3, Jaswinder Singh Maras3, Kalyan De2,4, Sambhaji Mote4, Srinivasan Krishnamurthi1,2


1Microbial Type Culture Collection & Gene Bank (MTCC), CSIR-Institute of Microbial Technology, Sector-39A, Chandigarh-160036.

2Academy of Scientific and Innovative Research (AcSIR), Ghaziabad, Uttar Pradesh, India

3Department of Molecular and Cellular Medicine, Institute of Liver and Biliary Sciences (ILBS), New Delhi

4Biological Oceanography Division, CSIR National Institute of Oceanography (NIO), Goa, India

Correpondance: srinivasan.kmurthi@csir.res.in
