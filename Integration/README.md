# Multi-omics Integration Analysis

To investigate relationships between metagenomic functional profiles and metabolomic compositions, multiple multivariate integration approaches were performed in R.

The integration analyses aimed to identify:

- Shared variation between metagenome and metabolome datasets
- Correlated functional and metabolic signatures
- Sample clustering consistency across omics layers
- Key metabolites and metagenomic features contributing to ecological differentiation

---

# R Packages Used

```r
library(vegan)
library(ggplot2)
library(ade4)
library(reshape2)
library(mixOmics)
library(pheatmap)
```

---

# Data Preparation

Two datasets were used:

1. Metagenome functional abundance matrix
2. Metabolome abundance matrix

Both datasets:

- Used samples as rows
- Features as columns
- Included only matched/common samples across datasets

Prior to integration analysis:

- Common samples were retained
- Features were variance filtered
- Top variable features were selected for visualization clarity
- Data scaling and normalization were applied where required

---

# 1. DIABLO Analysis (mixOmics)

DIABLO (`block.splsda`) analysis was performed using the `mixOmics` package to integrate metagenomic and metabolomic datasets.

## Objectives

- Identify correlated features across datasets
- Detect discriminative biomarkers
- Explore relationships between metabolites and metagenomic functions
- Visualize multi-omics associations

## Workflow

### Feature Selection

Top variable features were selected from:

- Metagenome dataset
- Metabolome dataset

using variance-based ranking.

### DIABLO Model

A supervised multi-block sparse PLS-DA model was constructed using:


block.splsda()


with:

- 2 components
- feature selection (`keepX`)
- scaling enabled

### Outputs Generated

- Sample clustering plots
- Variable plots
- Circos correlation plots
- Network plots
- Selected biomarker features
- Cross-omics correlation matrices

### Exported Files

DIABLO_SamplePlot_HQ.png
DIABLO_VariablePlot_HQ.png
DIABLO_CircosPlot_HQ.png
DIABLO_NetworkPlot_HQ.png
Selected_Metagenome_Features.csv
Selected_Metabolome_Features.csv
Feature_Correlations.csv

---

# 2. Procrustes Analysis

Procrustes analysis was performed using PCA ordinations from metagenome and metabolome datasets.

## Objectives

- Compare ordination similarity between datasets
- Evaluate concordance between metagenomic and metabolomic structures

## Workflow

### PCA Generation

PCA ordinations were generated separately for:

- Metagenome dataset
- Metabolome dataset

using:


rda()


### Procrustes Transformation

The ordinations were aligned using:


procrustes()


### Statistical Testing

Permutation-based significance testing was performed using:


protest()


with 999 permutations.

### Outputs


Procrustes_Plot.png
Procrustes_Plot.pdf


### Interpretation

- Lower residual distances indicate stronger agreement between datasets
- Significant Procrustes correlation suggests coupled metagenomic and metabolomic patterns

---

# 3. Co-inertia / Canonical Correspondence Analysis (CCA)

CCA analysis was performed to investigate relationships between metagenomic functions and metabolomic variables.

## Objectives

- Identify shared ecological gradients
- Determine metabolite-driven functional variation
- Explore co-structure between datasets

## Workflow

CCA was performed using:


cca()

where:

- Metagenome features were treated as response variables
- Metabolome features were treated as explanatory variables

### Statistical Testing

Permutation tests were performed for:

- Overall model significance
- Axis significance
- Individual variable significance

using:

```r
anova.cca()
```

### Outputs

```text
CCA_Plot.png
CCA_Plot.pdf
```

---

# 4. Heatmap Analysis

Heatmaps were generated using the `pheatmap` package to visualize abundance patterns of highly variable features.

## Workflow

- Top 50 variable features selected
- Feature scaling applied
- Hierarchical clustering performed on:
  - samples
  - features

### Outputs

```text
Heatmap_Metagenome_Top50.png
Heatmap_Metagenome_Top50.pdf
```

---

# Additional Statistical Analyses

The following additional analyses were also performed:

## PLS-DA Analysis

Partial Least Squares Discriminant Analysis (PLS-DA) was used to:

- investigate sample discrimination
- identify important metabolites/features
- evaluate group separation patterns

## PERMANOVA Analysis

PERMANOVA was performed using:

```r
adonis2()
```

to statistically evaluate differences in community/metabolome composition across:

- species
- locations
- experimental groups

based on Bray-Curtis dissimilarity matrices.

---

# Analyses Included in Manuscript

The manuscript primarily reports the following integration analyses:

- Procrustes analysis
- Co-inertia / CCA analysis
- DIABLO Circos plot

Additional analyses (PLS-DA, PERMANOVA, heatmaps, network plots) were performed as supporting exploratory analyses.

-----



# Repository structure

Integration/
├── README.md
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
    └── Procrustes/
