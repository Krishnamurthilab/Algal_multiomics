# =========================================================
# LC-MS METABOLOME NMDS + ENVFIT ANALYSIS
# =========================================================

# =========================================================
# LOAD LIBRARIES
# =========================================================

library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)

# =========================================================
# IMPORT DATA
# =========================================================

df <- read.csv(
  file.choose(),
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(file.choose())

# =========================================================
# FIX METADATA
# =========================================================

metadata <- metadata %>%
  rename(SampleID = Label)

rownames(metadata) <- metadata$SampleID

# =========================================================
# TRANSPOSE MATRIX
# Samples as rows
# =========================================================

mat <- t(df)

mat <- as.data.frame(mat)

# =========================================================
# CONVERT TO NUMERIC
# =========================================================

mat[] <- lapply(mat, as.numeric)

# =========================================================
# MATCH METADATA ORDER
# =========================================================

mat <- mat[metadata$SampleID, ]

# =========================================================
# REMOVE ZERO VARIANCE METABOLITES
# =========================================================

mat <- mat[, apply(mat, 2, var) > 0]

# =========================================================
# REMOVE ZERO SUM SAMPLES
# =========================================================

mat <- mat[rowSums(mat) > 0, ]

# =========================================================
# REPLACE NA VALUES
# =========================================================

mat[is.na(mat)] <- 0

# =========================================================
# SELECT TOP 200 VARIABLE METABOLITES
# =========================================================

var_rank <- apply(mat, 2, var)

mat <- mat[, order(var_rank,
                   decreasing = TRUE)[1:200]]

# =========================================================
# LOG TRANSFORMATION
# =========================================================

mat_log <- log1p(mat)

# =========================================================
# WISCONSIN NORMALIZATION
# =========================================================

mat_nmds <- wisconsin(mat_log)

# =========================================================
# BRAY-CURTIS DISTANCE
# =========================================================

dist_mat <- vegdist(
  mat_nmds,
  method = "bray"
)

# =========================================================
# NMDS ANALYSIS
# =========================================================

set.seed(123)

nmds3 <- metaMDS(
  dist_mat,
  k = 3,
  trymax = 200
)

# =========================================================
# STRESS VALUE
# =========================================================

print(nmds3$stress)

# =========================================================
# NMDS SCORES
# =========================================================

scores_nmds <- as.data.frame(
  scores(nmds3)
)

scores_nmds$SampleID <- rownames(scores_nmds)

scores_nmds <- merge(
  scores_nmds,
  metadata,
  by = "SampleID"
)

# =========================================================
# ENVFIT ANALYSIS
# =========================================================

fit <- envfit(
  nmds3,
  mat_nmds,
  permutations = 999
)

vec <- as.data.frame(
  scores(fit, "vectors")
)

vec$r2 <- fit$vectors$r

vec$pval <- fit$vectors$pvals

# =========================================================
# TOP 50 SIGNIFICANT METABOLITES
# =========================================================

vec_top50 <- vec %>%
  filter(pval < 0.05) %>%
  arrange(desc(r2)) %>%
  head(50)

# =========================================================
# CONVEX HULL
# =========================================================

hull_species <- scores_nmds %>%
  group_by(Species) %>%
  slice(chull(NMDS1, NMDS2))

# =========================================================
# PLOT
# =========================================================

arrow_scale <- 0.1

p <- ggplot(scores_nmds,
            aes(NMDS1, NMDS2)) +

  # Species hull
  geom_polygon(
    data = hull_species,
    aes(fill = Species,
        group = Species),
    alpha = 0.2,
    color = NA
  ) +

  # Sample points
  geom_point(
    aes(color = Species,
        shape = Location),
    size = 4
  ) +

  # Metabolite arrows
  geom_segment(
    data = vec_top50,
    aes(
      x = 0,
      y = 0,
      xend = NMDS1 * arrow_scale,
      yend = NMDS2 * arrow_scale
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "darkred",
    inherit.aes = FALSE
  ) +

  # Labels
  geom_text_repel(
    data = vec_top50,
    aes(
      x = NMDS1 * arrow_scale,
      y = NMDS2 * arrow_scale,
      label = rownames(vec_top50)
    ),
    size = 2.5,
    max.overlaps = 100
  ) +

  theme_bw(base_size = 14) +

  theme(
    panel.grid = element_blank()
  ) +

  labs(
    title = paste0(
      "NMDS (k = 3, Stress = ",
      round(nmds3$stress, 3),
      ")"
    )
  )

print(p)

# =========================================================
# CREATE OUTPUT DIRECTORY
# =========================================================

dir.create(
  "NMDS_output_final",
  showWarnings = FALSE
)

# =========================================================
# EXPORT MATRICES
# =========================================================

write.csv(
  mat,
  "NMDS_output_final/matrix_filtered.csv"
)

write.csv(
  mat_log,
  "NMDS_output_final/matrix_log1p.csv"
)

write.csv(
  mat_nmds,
  "NMDS_output_final/matrix_wisconsin.csv"
)

# =========================================================
# EXPORT DISTANCE MATRIX
# =========================================================

dist_df <- as.data.frame(
  as.matrix(dist_mat)
)

write.csv(
  dist_df,
  "NMDS_output_final/bray_curtis_distance_matrix.csv"
)

# =========================================================
# EXPORT NMDS SCORES
# =========================================================

write.csv(
  scores_nmds,
  "NMDS_output_final/NMDS_scores_k3.csv"
)

# =========================================================
# EXPORT ENVFIT RESULTS
# =========================================================

write.csv(
  vec,
  "NMDS_output_final/envfit_all_metabolites.csv"
)

write.csv(
  vec_top50,
  "NMDS_output_final/envfit_top50_metabolites.csv"
)

# =========================================================
# EXPORT HULL
# =========================================================

write.csv(
  hull_species,
  "NMDS_output_final/species_convex_hull.csv"
)

# =========================================================
# SAVE SUMMARY
# =========================================================

sink("NMDS_output_final/NMDS_summary.txt")

cat("NMDS analysis summary\n\n")

cat("Stress (k = 3):\n")

print(nmds3$stress)

cat("\nFull NMDS object:\n")

print(nmds3)

sink()

# =========================================================
# SAVE FIGURE
# =========================================================

ggsave(
  "NMDS_output_final/NMDS_plot.svg",
  plot = p,
  width = 10,
  height = 8
)

# =========================================================
# COMPLETED
# =========================================================

cat("\nNMDS analysis completed successfully.\n")
