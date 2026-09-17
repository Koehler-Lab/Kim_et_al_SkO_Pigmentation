## CellChatV2
# Perform cell-cell communication analysis and plot results

library(scDown)
library(Seurat)
library(SeuratObject)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
# Loading RDS, Dark Vs. Light
# lvl4 annotation transferred from Gopee 2024, renamed Melaoblast to Melanocyte
P200 <- readRDS(file = "P200_NewFetalSkinAnnotated2_Melanocyte.rds")

seurat_obj <- P200
output_dir <- "."
annotation_column <- "transferedCellType"
species <- "human"

sample_column <- NULL
annotation_selected <- NULL
group_column <- "orig.ident"
group_cmp <- list(unique(P200$orig.ident))
top_n <- 10

# Fix non-ASCII values in place across all columns
seurat_obj <- cleanSeuratMeta(seurat_obj)

# ------------------------------------------------------------------------------
# Figure S11d and Figure S13a–b (CellChat)
# ------------------------------------------------------------------------------
# Pairwise condition comparisons.
# Run CellChat V2
run_cellchatV2(
  seurat_obj = seurat_obj,
  output_dir = output_dir,
  annotation_column = annotation_column,
  species = species,
  group_column = group_column,
  group_cmp = group_cmp,
  cores = 1
)
