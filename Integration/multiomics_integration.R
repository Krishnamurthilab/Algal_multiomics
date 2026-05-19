# =========================================================
# Multi-omics Integration Analysis
# Metagenome + Metabolome Integration Workflow
# =========================================================

# =========================================================
# REQUIRED LIBRARIES
# =========================================================

# Install if needed
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!requireNamespace("mixOmics", quietly = TRUE))
  BiocManager::install("mixOmics")

packages <- c(
  "vegan",
  "ggplot2",
  "ade4",
  "reshape2",
  "mixOmics",
  "pheatmap",
  "dplyr"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
  
  library(pkg, character.only = TRUE)
}

# =========================================================
# CREATE OUTPUT DIRECTORIES
# =========================================================

dir.create("Integration_Output", showWarnings = FALSE)
dir.create("Integration_Output/DIABLO", showWarnings = FALSE)
dir.create("Integration_Output/Procrustes", showWarnings = FALSE)
dir.create("Integration_Output/CCA", showWarnings = FALSE)
dir.create("Integration_Output/Heatmap", showWarnings = FALSE)
dir.create("Integration_Output/PERMANOVA", showWarnings = FALSE)
dir.create("Integration_Output/PLSDA", showWarnings = FALSE)

# =========================================================
# LOAD DATA
# =========================================================

metaG <- read.csv(
  "metagenome.csv",
  row.names = 1,
  check.names = FALSE
)

metaB <- read.csv(
  "metabolome.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "metadata.csv",
  check.names = FALSE
)

# =========================================================
# FIX SAMPLE IDS
# =========================================================

rownames(metadata) <- metadata$SampleID

# =========================================================
# MATCH COMMON SAMPLES
# =========================================================

common_samples <- Reduce(
  intersect,
  list(
    rownames(metaG),
    rownames(metaB),
    metadata$SampleID
  )
)

metaG <- metaG[common_samples, , drop = FALSE]
metaB <- metaB[common_samples, , drop = FALSE]
metadata <- metadata[common_samples, ]

# =========================================================
# REMOVE NA
# =========================================================

metaG[is.na(metaG)] <- 0
metaB[is.na(metaB)] <- 0

# =========================================================
# REMOVE ZERO VARIANCE FEATURES
# =========================================================

metaG <- metaG[, apply(metaG, 2, var) > 0]
metaB <- metaB[, apply(metaB, 2, var) > 0]

# =========================================================
# SELECT TOP VARIABLE FEATURES
# =========================================================

metaG_var <- apply(metaG, 2, var)
metaB_var <- apply(metaB, 2, var)

metaG_top20 <- metaG[
  ,
  order(metaG_var,
        decreasing = TRUE)[1:min(20, ncol(metaG))],
  drop = FALSE
]

metaB_top20 <- metaB[
  ,
  order(metaB_var,
        decreasing = TRUE)[1:min(20, ncol(metaB))],
  drop = FALSE
]

# =========================================================
# DEFINE GROUPS
# Replace with your metadata column
# =========================================================

Y <- as.factor(metadata$Species)

# =========================================================
# SAVE FILTERED MATRICES
# =========================================================

write.csv(
  metaG_top20,
  "Integration_Output/metagenome_top20.csv"
)

write.csv(
  metaB_top20,
  "Integration_Output/metabolome_top20.csv"
)

# =========================================================
# =========================================================
# 1. DIABLO ANALYSIS
# =========================================================
# =========================================================

cat("\nRunning DIABLO analysis...\n")

# =========================================================
# PREPARE DATA LIST
# =========================================================

data_list <- list(
  metaG = metaG_top20,
  metaB = metaB_top20
)

# =========================================================
# KEEPX
# =========================================================

keepX_list <- list(
  metaG = rep(10, 2),
  metaB = rep(10, 2)
)

# =========================================================
# RUN DIABLO
# =========================================================

diablo_res <- block.splsda(
  X = data_list,
  Y = Y,
  ncomp = 2,
  keepX = keepX_list,
  scale = TRUE
)

# =========================================================
# SELECTED FEATURES
# =========================================================

selected_metaG <- selectVar(
  diablo_res,
  block = "metaG",
  comp = 1
)$name

selected_metaB <- selectVar(
  diablo_res,
  block = "metaB",
  comp = 1
)$name

# =========================================================
# EXPORT FEATURES
# =========================================================

write.csv(
  selected_metaG,
  "Integration_Output/DIABLO/Selected_Metagenome_Features.csv",
  row.names = FALSE
)

write.csv(
  selected_metaB,
  "Integration_Output/DIABLO/Selected_Metabolome_Features.csv",
  row.names = FALSE
)

# =========================================================
# SAMPLE PLOT
# =========================================================

png(
  "Integration_Output/DIABLO/DIABLO_SamplePlot.png",
  width = 3000,
  height = 2500,
  res = 600
)

plotIndiv(
  diablo_res,
  legend = TRUE,
  title = "DIABLO Sample Plot"
)

dev.off()

# =========================================================
# VARIABLE PLOT
# =========================================================

png(
  "Integration_Output/DIABLO/DIABLO_VariablePlot.png",
  width = 3000,
  height = 2500,
  res = 600
)

plotVar(
  diablo_res,
  var.names = TRUE,
  title = "DIABLO Variable Plot"
)

dev.off()

# =========================================================
# CIRCOS PLOT
# =========================================================

png(
  "Integration_Output/DIABLO/DIABLO_CircosPlot.png",
  width = 3000,
  height = 3000,
  res = 600
)

print(
  circosPlot(
    diablo_res,
    cutoff = 0.6,
    title = "DIABLO Circos Plot"
  )
)

dev.off()

# =========================================================
# NETWORK PLOT
# =========================================================

png(
  "Integration_Output/DIABLO/DIABLO_NetworkPlot.png",
  width = 3000,
  height = 2500,
  res = 600
)

print(
  network(
    diablo_res,
    cutoff = 0.6
  )
)

dev.off()

# =========================================================
# CORRELATION MATRIX
# =========================================================

if(length(selected_metaG) > 0 &
   length(selected_metaB) > 0){
  
  cor_matrix <- cor(
    metaG_top20[, selected_metaG, drop = FALSE],
    metaB_top20[, selected_metaB, drop = FALSE],
    use = "pairwise.complete.obs"
  )
  
  write.csv(
    round(cor_matrix, 3),
    "Integration_Output/DIABLO/Feature_Correlations.csv"
  )
}

# =========================================================
# =========================================================
# 2. PROCRUSTES ANALYSIS
# =========================================================
# =========================================================

cat("\nRunning Procrustes analysis...\n")

# =========================================================
# PCA
# =========================================================

pca_metaG <- rda(
  metaG_top20,
  scale = TRUE
)

pca_metaB <- rda(
  metaB_top20,
  scale = TRUE
)

# =========================================================
# PROCRUSTES
# =========================================================

proc_res <- procrustes(
  pca_metaG,
  pca_metaB,
  symmetric = TRUE
)

# =========================================================
# PLOT
# =========================================================

png(
  "Integration_Output/Procrustes/Procrustes_Plot.png",
  width = 2500,
  height = 2000,
  res = 600
)

plot(
  proc_res,
  kind = 1,
  main = "Procrustes Plot"
)

dev.off()

# =========================================================
# PROTEST
# =========================================================

proc_test <- protest(
  pca_metaG,
  pca_metaB,
  permutations = 999
)

capture.output(
  proc_test,
  file = "Integration_Output/Procrustes/Procrustes_statistics.txt"
)

# =========================================================
# =========================================================
# 3. CCA ANALYSIS
# =========================================================
# =========================================================

cat("\nRunning CCA analysis...\n")

cca_res <- cca(
  metaG_top20 ~ .,
  data = metaB_top20
)

# =========================================================
# SAVE SUMMARY
# =========================================================

capture.output(
  summary(cca_res),
  file = "Integration_Output/CCA/CCA_summary.txt"
)

# =========================================================
# PLOT
# =========================================================

png(
  "Integration_Output/CCA/CCA_Plot.png",
  width = 3000,
  height = 2500,
  res = 600
)

plot(
  cca_res,
  display = c("species", "sites", "bp"),
  main = "CCA Analysis"
)

dev.off()

# =========================================================
# SIGNIFICANCE TESTS
# =========================================================

cca_overall <- anova(
  cca_res,
  permutations = 999
)

cca_axis <- anova(
  cca_res,
  by = "axis",
  permutations = 999
)

cca_term <- anova(
  cca_res,
  by = "term",
  permutations = 999
)

capture.output(
  cca_overall,
  file = "Integration_Output/CCA/CCA_overall_test.txt"
)

capture.output(
  cca_axis,
  file = "Integration_Output/CCA/CCA_axis_test.txt"
)

capture.output(
  cca_term,
  file = "Integration_Output/CCA/CCA_term_test.txt"
)

# =========================================================
# =========================================================
# 4. HEATMAP
# =========================================================
# =========================================================

cat("\nGenerating heatmap...\n")

metaG_var50 <- apply(metaG, 2, var)

top50_features <- names(
  sort(metaG_var50,
       decreasing = TRUE)[1:min(50, ncol(metaG))]
)

metaG_top50 <- metaG[, top50_features]

metaG_scaled <- t(
  scale(
    t(metaG_top50)
  )
)

png(
  "Integration_Output/Heatmap/Heatmap_Metagenome_Top50.png",
  width = 3000,
  height = 2500,
  res = 600
)

pheatmap(
  metaG_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Top 50 Metagenome Features"
)

dev.off()

# =========================================================
# =========================================================
# 5. PLS-DA
# =========================================================
# =========================================================

cat("\nRunning PLS-DA...\n")

plsda_res <- plsda(
  metaB_top20,
  Y,
  ncomp = 2
)

png(
  "Integration_Output/PLSDA/PLSDA_plot.png",
  width = 3000,
  height = 2500,
  res = 600
)

plotIndiv(
  plsda_res,
  comp = c(1,2),
  group = Y,
  ellipse = TRUE,
  legend = TRUE,
  title = "PLS-DA"
)

dev.off()

# =========================================================
# =========================================================
# 6. PERMANOVA
# =========================================================
# =========================================================

cat("\nRunning PERMANOVA...\n")

dist_metaB <- vegdist(
  metaB_top20,
  method = "bray"
)

permanova_res <- adonis2(
  dist_metaB ~ Species + Location,
  data = metadata,
  permutations = 999
)

capture.output(
  permanova_res,
  file = "Integration_Output/PERMANOVA/PERMANOVA_results.txt"
)

# =========================================================
# FINAL SESSION INFO
# =========================================================

writeLines(
  capture.output(sessionInfo()),
  "Integration_Output/sessionInfo.txt"
)

cat("\n=====================================\n")
cat("Integration analysis completed.\n")
cat("=====================================\n")
