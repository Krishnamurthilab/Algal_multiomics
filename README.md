# Algal multiomics analysis
The repositpory contains complete and reproducible workflow used for multiomics analysis of macroalgal samples from Central West Coast of Goa, India.

## Overview
The study integrates v3-v4 amplicon sequencing, whole genome metagenome sequencing and untargeted LC-MS based metabolome analysis to investigate structural and functional community composition and factors influencing the community composition of the algae.

## 1. Amplicon Analysis
Amplicon sequencing analysis was performed to investigate microbial community composition associated with marine algal samples.

Analyses Included
Quality filtering
Denoising
Taxonomic assignment
Alpha diversity
Beta diversity
Ordination analysis

## 2. Metagenome Analysis
Shotgun metagenomic analysis was performed to explore microbial functional potential and ecological adaptation.

Analyses Included
Quality control
Assembly
Functional annotation
Taxonomic profiling
COG analysis
Metabolic pathway analysis
Comparative functional profiling

## 3. Metabolome Analysis
LC-MS based metabolome profiling was performed to investigate algal-associated metabolite diversity.

### Metabolite Curation
The metabolite dataset was manually curated to remove: 1. Duplicate metabolites, 2. Unmatched compounds, 3. Synthetic drug-related metabolites


### Metabolite annotation strategies
Each metabolites were manually search against PubChem, HMDB, KEGG, ChEBI, LipidMaps, and published literatures to determine 1. Biological source, 2. Chemical nature, 3. Functional and ecological relevance

Additional filteration steps were perforemd which excluded terrestrail, human, animal related origin source, resulting in a final set of 300 metabolites associated with marine, algal, plant, and microbial systems.

The metabolome analysis included:

Data normalization
Bray-Curtis dissimilarity
NMDS ordination
Environmental fitting (envfit)
Convex hull analysis
PCA biplot analysis
PLS-DA
PERMANOVA

## 4. Multi-omics Integration Analysis

Integrated analysis of metagenomic and metabolomic datasets was performed to identify relationships between microbial functions and metabolite profiles.

The workflow included:

DIABLO Analysis
Procrustes Analysis
Co-inertia / CCA Analysis


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
