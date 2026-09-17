
# ============================================================
# Spatial Transcriptomic Analysis of Human Breast Cancer
# 05 - Spatial variability analysis using Moran's I
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Moran's I dependency
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The 'ape' package is required for Moran's I analysis.")
}

# Load clustered Visium object
breast <- readRDS(
  "data/processed/breast_visium_clustered.rds"
)

# Use the SCTransform-normalised assay
DefaultAssay(breast) <- "SCT"

breast

# ============================================================
# Confirm spatial image and variable features
# ============================================================

Images(breast)

length(VariableFeatures(breast))

head(VariableFeatures(breast))

# ============================================================
# Select highly variable genes for spatial testing
# ============================================================

genes_for_morans <- VariableFeatures(breast)[1:1000]

length(genes_for_morans)

head(genes_for_morans)

# ============================================================
# Identify spatially variable genes using Moran's I
# ============================================================

breast <- FindSpatiallyVariableFeatures(
  breast,
  assay = "SCT",
  layer = "scale.data",
  features = genes_for_morans,
  image = "breast_cancer",
  selection.method = "moransi",
  verbose = TRUE
)

# ============================================================
# Retrieve genes ranked by spatial variability
# ============================================================

spatial_genes <- SpatiallyVariableFeatures(
  breast,
  method = "moransi"
)

head(spatial_genes, 20)

# ============================================================
# Extract Moran's I statistics
# ============================================================

morans_results <- SVFInfo(
  breast[["SCT"]],
  method = "moransi",
  status = TRUE
)

# Add gene names as a normal column
morans_results$gene <- rownames(morans_results)

# Create clearly named copies of the Seurat output columns
morans_results$variable_flag <- as.logical(
  unlist(morans_results[["variable"]])
)

morans_results$spatial_rank <- as.numeric(
  unlist(morans_results[["rank"]])
)

# Sort by Moran's I spatial-variable rank
morans_results <- morans_results[
  order(morans_results$spatial_rank, na.last = TRUE),
  ,
  drop = FALSE
]

head(morans_results, 20)

# ============================================================
# Save Moran's I results
# ============================================================

# Create compact Moran's I results table for export
morans_export <- data.frame(
  gene = morans_results$gene,
  MoransI_observed = as.numeric(morans_results$MoransI_observed),
  MoransI_p_value = as.numeric(morans_results$MoransI_p.value),
  spatially_variable = morans_results$variable_flag,
  spatial_rank = morans_results$spatial_rank
)

write.csv(
  morans_export,
  "results/morans_i_results_compact.csv",
  row.names = FALSE
)

# Keep spatially variable genes with a valid rank
top50_morans <- morans_results[
  morans_results$variable_flag == TRUE &
    !is.na(morans_results$spatial_rank),
  ,
  drop = FALSE
]

# Order and keep the top 50
top50_morans <- top50_morans[
  order(top50_morans$spatial_rank),
  ,
  drop = FALSE
]

top50_morans <- head(top50_morans, 50)

write.csv(
  top50_morans,
  "results/top50_spatially_variable_genes.csv",
  row.names = FALSE
)

top50_morans

# ============================================================
# Top 20 spatially variable genes
# ============================================================

top20_spatial <- morans_results[
  morans_results$variable_flag == TRUE &
    !is.na(morans_results$spatial_rank),
  ,
  drop = FALSE
]

top20_spatial <- top20_spatial[
  order(top20_spatial$spatial_rank),
  ,
  drop = FALSE
]

top20_spatial <- head(top20_spatial, 20)

top20_spatial

# ============================================================
# Select top six spatially variable genes
# ============================================================

top6_spatial_genes <- head(
  spatial_genes,
  6
)

top6_spatial_genes

# ============================================================
# Visualise top spatially variable genes
# ============================================================

top_spatial_plot <- SpatialFeaturePlot(
  breast,
  features = top6_spatial_genes,
  ncol = 3,
  alpha = c(0.1, 1),
  max.cutoff = "q95"
)

top_spatial_plot

ggsave(
  filename = "figures/16_top_spatially_variable_genes.png",
  plot = top_spatial_plot,
  width = 14,
  height = 9,
  dpi = 300
)

# ============================================================
# Examine biologically relevant spatial genes
# ============================================================

candidate_genes <- c(
  "CCL19",
  "CCL21",
  "COL11A1",
  "GREM1",
  "VEGFA",
  "MKI67",
  "EPCAM",
  "KRT19"
)

candidate_genes <- intersect(
  candidate_genes,
  rownames(breast)
)

candidate_genes

candidate_spatial_plot <- SpatialFeaturePlot(
  breast,
  features = candidate_genes,
  ncol = 4,
  max.cutoff = "q95"
)

candidate_spatial_plot

ggsave(
  filename = "figures/17_candidate_spatial_programmes.png",
  plot = candidate_spatial_plot,
  width = 16,
  height = 8,
  dpi = 300
)

# ============================================================
# Save object containing spatial variability analysis
# ============================================================

saveRDS(
  breast,
  file = "data/processed/breast_visium_spatial_variability.rds"
)

# ============================================================
# Record computational environment for reproducibility
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  "results/sessionInfo.txt"
)