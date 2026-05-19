# LC-MS Metabolome Analysis Workflow

## Overview

This repository contains the workflow used for LC-MS based metabolome profiling and statistical analysis of marine algal samples.

The study focused on identifying biologically relevant metabolites associated with marine, algal, plant, and microbial origin followed by diversity, ordination, environmental fitting, and multivariate statistical analysis to know the influence of algal surface metabolites in colonization of microbial community.

The workflow integrates:

- Manual metabolite curation
- Metabolite annotation and validation
- MetaboAnalyst-based diversity analysis
- R-based metabolome ordination and statistical analysis
- NMDS analysis
- Environmental fitting (`envfit`)
- PLS-DA analysis

# Experimental Workflow

Raw LC-MS metabolite abundance table
                    ↓
Preliminary filteration of abundance table
                    ↓
Removal of:
  • duplicates
  • unmatched compounds
                    ↓
Manual metabolite annotation
using PubChem, HMDB, KEGG, literature
                    ↓
removal of terrestrial, human, and animal origin compounds
                    ↓
Selection of biologically relevant metabolites (~300 metabolites retained)
                    ↓
R-based statistical analyses
    ├── Statistical analysis
    ├── ordination analysis
    ├── envfit metabolite vectors
    └── PLS-DA
                    ↓
MetaboAnalyst based analysis

# Step 1: Initial filteration 
* Manual removal of duplicate and no matched compounds

# Steps 2: Metabolite Source, Nature, and Functional Annotation

Following initial filtering, each metabolite was manually investigated to determine its:

- Biological source
- Chemical nature
- Ecological and functional relevance

## Annotation Strategy

Each metabolite was individually searched using Google and cross-validated through multiple metabolomics and chemical databases, including:

- PubChem
- HMDB (Human Metabolome Database)
- KEGG
- ChEBI
- MetaboLights
- LipidMaps
- FooDB
- Scientific literature and published reports

The annotation process focused on identifying:

### A. Biological Source

Metabolites were categorized based on reported biological origin, including:

- Human/Animal derived metabolites
- Terrestrial originate metabolites
- Marine-associated metabolites
- Algal-derived metabolites
- Plant-derived metabolites
- Microbial metabolites
- Secondary metabolites from symbiotic organisms

### B. Chemical Nature

Metabolites were further classified according to their biochemical characteristics, such as:

- Lipids and fatty acids
- Amino acid derivatives
- Phenolics
- Terpenoids
- Alkaloids
- Pigments
- Organic acids
- Osmolytes
- Sugars and sugar alcohols

### C. Functional and Ecological Relevance

Functional annotation was performed to identify ecological or biological roles of metabolites, including:

- Antioxidant activity
- Antimicrobial compounds
- Stress-response metabolites
- Osmoregulation-related metabolites
- Signaling molecules
- Nutrient-associated metabolites
- Host-microbe interaction compounds
- Secondary metabolite biosynthesis products

---

# Final Metabolite Selection

Based on the manual annotation and ecological relevance assessment, approximately 300 metabolites associated with marine, algal, plant, and microbial systems were retained for downstream analyses.

# R-based Statistical Analysis

## R Packages Used

```r
library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(mixOmics)
```

---

# NMDS Analysis Workflow

## 1. Data Preprocessing

The following preprocessing steps were applied:

- Matrix transposition
- Numeric conversion
- Removal of zero variance metabolites
- Removal of zero abundance samples
- Replacement of missing values
- Selection of top 200 most variable metabolites
- Log transformation (`log1p`)
- Wisconsin double standardization

---

## 2. Bray-Curtis Distance

Bray-Curtis dissimilarity matrix was calculated using:

```r
vegdist(method = "bray")
```

---

## 3. NMDS Ordination

Three-dimensional NMDS analysis was performed using:

```r
metaMDS(k = 3, trymax = 200)
```

---

## 4. Environmental Fitting

Environmental fitting analysis (`envfit`) was performed to identify metabolites significantly associated with sample clustering patterns.

Metabolites were ranked using:

- Correlation coefficient (`r²`)
- Statistical significance (`p-value`)

Top 50 significant metabolites were visualized as vectors.

---

## 5. Convex Hull Visualization

Species-level clustering patterns were visualized using convex hull polygons generated with:

```r
chull()
```

---

# Output Files

```text
NMDS_output_final/
├── matrix_filtered.csv
├── matrix_log1p.csv
├── matrix_wisconsin.csv
├── bray_curtis_distance_matrix.csv
├── NMDS_scores_k3.csv
├── envfit_all_metabolites.csv
├── envfit_top50_metabolites.csv
├── species_convex_hull.csv
├── NMDS_summary.txt
└── NMDS_plot.svg
```

---

# Repository Structure

```text
.
├── data/
│   ├── raw_metabolome.xlsx
│   ├── filtered_metabolome.csv
│   └── metadata.csv
│
├── scripts/
│   ├── 01_metabolite_filtering_notes.md
│   ├── 02_nmds_envfit_analysis.R
│   ├── 03_pca_biplot_analysis.R
│   └── 04_plsda_analysis.R
│
├── results/
│   ├── NMDS_output_final/
│   ├── biplot/
│   └── plsda/
│
└── README.md
```

---

# Biological Interpretation

The metabolome analysis revealed metabolite compositional differences associated with:

- Host species
- Geographic sampling location
- Ecological niche adaptation

Significant metabolites identified through `envfit` likely contribute to:

- Host-microbe interactions
- Marine ecological adaptation
- Stress response
- Osmoregulation
- Secondary metabolite production
- Chemical signaling
----------

# MetaboAnalyst Analysis

The filtered metabolite abundance matrix was uploaded to the MetaboAnalyst server for metabolome diversity analysis.

Analyses included:

- Data normalization
- Data transformation
- Ordination analysis
- Diversity analysis
- Clustering visualization

MetaboAnalyst server:

https://www.metaboanalyst.ca/

---
