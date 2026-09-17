# ============================================================
# Spatial Transcriptomic Analysis of Human Breast Cancer
# 01 - Data loading and initial quality control
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Path to the 10x Visium breast cancer dataset
data_dir <- file.path(
  "data",
  "raw",
  "Visium_Human_Breast_Cancer"
)

# Confirm that R can find the files
list.files(data_dir)
list.files(file.path(data_dir, "spatial"))

# Load the 10x Visium dataset
breast <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "Visium_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "breast_cancer"
)

breast
breast@images
SpatialDimPlot(
  breast,
  images = "breast_cancer"
)

# Create processed-data folder if it does not already exist
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save imported Visium object
saveRDS(
  breast,
  file = "data/processed/breast_visium_raw.rds"
)

# Check that the file was saved
list.files("data/processed")

# Create processed-data folder if it does not already exist
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save imported Visium object
saveRDS(
  breast,
  file = "data/processed/breast_visium_raw.rds"
)

# Check that the file was saved
list.files("data/processed")

# Create processed data folder
dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save the Seurat object
saveRDS(
  breast,
  file = "data/processed/breast_visium_raw.rds"
)

# ============================================================
# Quality control
# ============================================================

# Calculate mitochondrial transcript percentage per Visium spot
breast[["percent.mt"]] <- PercentageFeatureSet(
  breast,
  pattern = "^MT-"
)

# Inspect QC metrics
summary(breast$nCount_Spatial)
summary(breast$nFeature_Spatial)
summary(breast$percent.mt)

# Create figures folder if needed
dir.create(
  "figures",
  showWarnings = FALSE
)

# QC distributions
qc_plot <- VlnPlot(
  breast,
  features = c(
    "nCount_Spatial",
    "nFeature_Spatial",
    "percent.mt"
  ),
  ncol = 3,
  pt.size = 0.1
)

qc_plot

ggsave(
  filename = "figures/01_qc_violin.png",
  plot = qc_plot,
  width = 12,
  height = 5,
  dpi = 300
)

qc_spatial_counts <- SpatialFeaturePlot(
  breast,
  features = "nCount_Spatial"
)

qc_spatial_features <- SpatialFeaturePlot(
  breast,
  features = "nFeature_Spatial"
)

qc_spatial_counts
qc_spatial_features

ggsave(
  "figures/02_spatial_counts.png",
  qc_spatial_counts,
  width = 7,
  height = 7,
  dpi = 300
)

ggsave(
  "figures/03_spatial_features.png",
  qc_spatial_features,
  width = 7,
  height = 7,
  dpi = 300
)

saveRDS(
  breast,
  file = "data/processed/breast_visium_qc.rds"
)

list.files("figures")
list.files("data/processed")

# ============================================================
# QC summary statistics
# ============================================================

summary(breast$nCount_Spatial)
summary(breast$nFeature_Spatial)
summary(breast$percent.mt)

quantile(
  breast$nCount_Spatial,
  probs = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)
)

quantile(
  breast$nFeature_Spatial,
  probs = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)
)

quantile(
  breast$percent.mt,
  probs = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)
)
saveRDS(
  breast,
  file = "data/processed/breast_visium_qc.rds"
)