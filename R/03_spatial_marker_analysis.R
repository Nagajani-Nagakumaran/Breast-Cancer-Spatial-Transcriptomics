# ============================================================
# Spatial Transcriptomic Analysis of Human Breast Cancer
# 03 - Marker analysis and biological interpretation
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load clustered Visium object
breast <- readRDS(
  "data/processed/breast_visium_clustered.rds"
)

breast

# ============================================================
# Find cluster marker genes
# ============================================================

markers <- FindAllMarkers(
  breast,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# View top markers for each cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

top_markers

write.csv(
  markers,
  "results/all_cluster_markers.csv",
  row.names = FALSE
)

write.csv(
  top_markers,
  "results/top10_markers_per_cluster.csv",
  row.names = FALSE
)

# ============================================================
# Known marker genes of tumour, immune and stromal regions
# ============================================================

tumour_markers <- c("EPCAM", "KRT8", "KRT18", "KRT19")
immune_markers <- c("PTPRC", "CD3D", "CD3E", "CD68")
stroma_markers <- c("COL1A1", "COL1A2", "DCN", "VIM")

tumour_umap <- FeaturePlot(
  breast,
  features = tumour_markers,
  reduction = "umap",
  ncol = 2
)

immune_umap <- FeaturePlot(
  breast,
  features = immune_markers,
  reduction = "umap",
  ncol = 2
)

stroma_umap <- FeaturePlot(
  breast,
  features = stroma_markers,
  reduction = "umap",
  ncol = 2
)

tumour_umap
immune_umap
stroma_umap

ggsave("figures/07_tumour_markers_umap.png", tumour_umap, width = 10, height = 8, dpi = 300)
ggsave("figures/08_immune_markers_umap.png", immune_umap, width = 10, height = 8, dpi = 300)
ggsave("figures/09_stroma_markers_umap.png", stroma_umap, width = 10, height = 8, dpi = 300)

tumour_spatial <- SpatialFeaturePlot(
  breast,
  features = tumour_markers,
  ncol = 2
)

immune_spatial <- SpatialFeaturePlot(
  breast,
  features = immune_markers,
  ncol = 2
)

stroma_spatial <- SpatialFeaturePlot(
  breast,
  features = stroma_markers,
  ncol = 2
)

tumour_spatial
immune_spatial
stroma_spatial

ggsave("figures/10_tumour_markers_spatial.png", tumour_spatial, width = 10, height = 8, dpi = 300)
ggsave("figures/11_immune_markers_spatial.png", immune_spatial, width = 10, height = 8, dpi = 300)
ggsave("figures/12_stroma_markers_spatial.png", stroma_spatial, width = 10, height = 8, dpi = 300)

selected_markers <- c(
  "EPCAM", "KRT8", "KRT18", "KRT19",
  "PTPRC", "CD3D", "CD3E", "CD68",
  "COL1A1", "COL1A2", "DCN", "VIM"
)

marker_heatmap <- DoHeatmap(
  breast,
  features = selected_markers,
  group.by = "seurat_clusters"
) + NoLegend()

marker_heatmap

ggsave(
  "figures/13_marker_heatmap.png",
  marker_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)