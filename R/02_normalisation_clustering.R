# ============================================================
# Spatial Transcriptomic Analysis of Human Breast Cancer
# 02 - Normalisation, dimensionality reduction and clustering
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load QC-completed Visium object
breast <- readRDS(
  "data/processed/breast_visium_qc.rds"
)

# Confirm object loaded correctly
breast

# ============================================================
# SCTransform normalisation
# ============================================================

breast <- SCTransform(
  breast,
  assay = "Spatial",
  verbose = FALSE
)

DefaultAssay(breast) <- "SCT"

# ============================================================
# Principal component analysis
# ============================================================

breast <- RunPCA(
  breast,
  assay = "SCT",
  verbose = FALSE
)

# Examine variance explained across principal components
elbow_plot <- ElbowPlot(
  breast,
  ndims = 50
)

elbow_plot

# Save PCA elbow plot
ggsave(
  filename = "figures/04_pca_elbow.png",
  plot = elbow_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# ============================================================
# Nearest-neighbour graph
# ============================================================

set.seed(1234)

breast <- FindNeighbors(
  breast,
  reduction = "pca",
  dims = 1:30
)

# ============================================================
# Graph-based clustering
# ============================================================

breast <- FindClusters(
  breast,
  resolution = 0.5,
  random.seed = 1234
)

# ============================================================
# UMAP dimensionality reduction
# ============================================================

breast <- RunUMAP(
  breast,
  reduction = "pca",
  dims = 1:30,
  seed.use = 1234
)

# Visualise transcriptomic clusters in UMAP space
umap_clusters <- DimPlot(
  breast,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

umap_clusters

ggsave(
  filename = "figures/05_umap_clusters.png",
  plot = umap_clusters,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# Map transcriptomic clusters back onto the tissue
# ============================================================

spatial_clusters <- SpatialDimPlot(
  breast,
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 3
)

spatial_clusters

ggsave(
  filename = "figures/06_spatial_clusters.png",
  plot = spatial_clusters,
  width = 8,
  height = 8,
  dpi = 300
)

# Save clustered spatial transcriptomic object
saveRDS(
  breast,
  file = "data/processed/breast_visium_clustered.rds"
)