# Differential expression (Dark vs Light)
# - Cell types transferred from Gopee 2024 lvl4 annotation
# - Analyze cell types with ≥0.1% representation in both samples
# - Output: per-cell-type DEG tables + volcano plots; plus combined analysis across all eligible cell types

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(EnhancedVolcano)

# ------------------------------------------------------------------------------
# Load inputs (DEG Light vs Dark)
# ------------------------------------------------------------------------------
# Cell types with at least 0.1% representation in both L'Oreal P200 Dark vs Light samples
celltype <- read.table("P200_NewFetalSkinlvl4_.001celltype.txt", sep = "\t", stringsAsFactors = FALSE)[, 1]
length(celltype) # 25

# Loading RDS, Dark vs Light
# lvl4 annotation transferred from Gopee 2024; renamed Melaoblast -> Melanocyte
Y <- readRDS("P200_NewFetalSkinAnnotated2_Melanocyte.rds")

dge <- Y
Idents(dge) <- dge$orig.ident
DefaultAssay(dge) <- "RNA"
dgeall <- dge

# ------------------------------------------------------------------------------
# Figure 4e and Figure S11b
# Differential expression in each cell type (≥0.1% in both samples)
# ------------------------------------------------------------------------------
markerslist <- list()

for (i in seq_along(celltype)) {
  label <- celltype[i]
  print(c(i, label))

  dge <- subset(dgeall, transferedCellType == label)
  DefaultAssay(dge) <- "RNA"
  Idents(dge) <- dge$orig.ident

  print(table(Idents(dge)))

  markers1 <- FindMarkers(
    object = dge,
    ident.1 = "4-L-GEX",
    only.pos = FALSE,
    min.pct = -Inf,
    min.diff.pct = -Inf,
    logfc.threshold = -Inf,
    min.cells.feature = -Inf,
    min.cells.group = -Inf
  )

  write.table(
    x = markers1,
    file = paste0("markers/P200_LvsD_lvl4_C", i, "_allgenes_de_12.2023.txt"),
    col.names = TRUE,
    row.names = TRUE,
    quote = FALSE,
    sep = "\t"
  )

  markerslist[[i]] <- markers1
}

# ------------------------------------------------------------------------------
# Volcano plots per cell type (from saved DEG tables)
# ------------------------------------------------------------------------------
markerslist <- list()

for (i in seq_along(celltype)) {
  label <- celltype[i]
  print(c(i, label))

  markers1 <- read.table(
    file = paste0("markers/P200_LvsD_lvl4_C", i, "_allgenes_de_12.2023.txt"),
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE,
    sep = "\t"
  )

  a <- markers1
  markerslist[[i]] <- markers1

  png(filename = paste0("markers/P200_LvsD_lvl4_C", i, "_VolcanoEnhanced.png"), res = 300, height = 2000, width = 1600)
  par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))

  plot <- EnhancedVolcano(
    a,
    lab = rownames(a),
    x = "avg_log2FC",
    y = "p_val",
    FCcutoff = log2(1.5)
  )

  print(plot)
  dev.off()
}

# ------------------------------------------------------------------------------
# Figure S11a
# Differential expression across all eligible cell types combined (all 25)
# ------------------------------------------------------------------------------
labels <- celltype

dgeall$CellType <- as.character(dgeall$transferedCellType)
names(dgeall$CellType) <- names(dgeall$transferedCellType)

dge <- subset(dgeall, CellType %in% labels)
length(table(dge$CellType)) # 25

DefaultAssay(dge) <- "RNA"
Idents(dge) <- dge$orig.ident
print(table(Idents(dge)))

markers1 <- FindMarkers(
  object = dge,
  ident.1 = "4-L-GEX",
  only.pos = FALSE,
  min.pct = -Inf,
  min.diff.pct = -Inf,
  logfc.threshold = -Inf,
  min.cells.feature = -Inf,
  min.cells.group = -Inf
)

write.table(
  x = markers1,
  file = "markers/P200_LvsD_lvl4_.001celltypes_allgenes_de_12.2023.txt",
  col.names = TRUE,
  row.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

markerslist[["All25"]] <- markers1

markers1 <- read.table(
  file = "markers/P200_LvsD_lvl4_.001celltypes_allgenes_de_12.2023.txt",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE,
  sep = "\t"
)

a <- markers1
markerslist[["All25"]] <- markers1

png(filename = "markers/P200_LvsD_lvl4_.001celltypes_VolcanoEnhanced.png", res = 300, height = 2000, width = 1600)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))

plot <- EnhancedVolcano(
  a,
  lab = rownames(a),
  x = "avg_log2FC",
  y = "p_val",
  FCcutoff = log2(1.5)
)

print(plot)
dev.off()
