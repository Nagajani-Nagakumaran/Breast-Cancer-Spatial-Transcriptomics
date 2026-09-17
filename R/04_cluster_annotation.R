# ============================================================
# Spatial Transcriptomic Analysis of Human Breast Cancer
# 04 - Cluster annotation and tumour microenvironment profiling
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

breast <- readRDS(
  "data/processed/breast_visium_clustered.rds"
)

DefaultAssay(breast) <- "SCT"

# Marker panels for broad tissue compartments

epithelial_markers <- c(
  "EPCAM", "KRT8", "KRT18", "KRT19",
  "GATA3", "FOXA1", "ESR1"
)

tcell_markers <- c(
  "PTPRC", "CD3D", "CD3E",
  "CD4", "CD8A"
)

myeloid_markers <- c(
  "LST1", "TYROBP", "FCER1G",
  "CD68", "C1QA", "C1QB"
)

fibroblast_markers <- c(
  "COL1A1", "COL1A2",
  "COL3A1", "DCN", "LUM"
)

endothelial_markers <- c(
  "PECAM1", "VWF", "EMCN"
)

proliferation_markers <- c(
  "MKI67", "TOP2A"
)

all_requested_markers <- unique(c(
  epithelial_markers,
  tcell_markers,
  myeloid_markers,
  fibroblast_markers,
  endothelial_markers,
  proliferation_markers
))

available_markers <- intersect(
  all_requested_markers,
  rownames(breast)
)

available_markers

cluster_dotplot <- DotPlot(
  breast,
  features = available_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  labs(
    x = "Marker gene",
    y = "Transcriptomic cluster"
  )

cluster_dotplot

ggsave(
  "figures/14_cluster_marker_dotplot.png",
  cluster_dotplot,
  width = 1000,
  height = 400,
  dpi = 300
)

breast <- AddModuleScore(
  breast,
  features = list(epithelial_markers),
  name = "Epithelial"
)

breast <- AddModuleScore(
  breast,
  features = list(tcell_markers),
  name = "Tcell"
)

breast <- AddModuleScore(
  breast,
  features = list(myeloid_markers),
  name = "Myeloid"
)

breast <- AddModuleScore(
  breast,
  features = list(fibroblast_markers),
  name = "Fibroblast"
)

programme_spatial <- SpatialFeaturePlot(
  breast,
  features = c(
    "Epithelial1",
    "Tcell1",
    "Myeloid1",
    "Fibroblast1"
  ),
  ncol = 2
)

programme_spatial

ggsave(
  "figures/15_spatial_tme_programmes.png",
  programme_spatial,
  width = 11,
  height = 9,
  dpi = 300
)