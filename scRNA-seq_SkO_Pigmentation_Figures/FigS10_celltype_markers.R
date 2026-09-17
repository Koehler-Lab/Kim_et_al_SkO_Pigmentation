### Dot Plots of top unique markers (fine and merged cell types)

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(pbapply)
library(Nebulosa)
library(ComplexHeatmap)

# Loading RDS, Dark Vs. Light
# lvl4 annotation transferred from Gopee 2024, renamed Melaoblast to Melanocyte
P200 <- readRDS("P200_NewFetalSkinAnnotated2_Melanocyte.rds")

dge <- P200
levels(Idents(dge))

levelsall <- levels(Idents(dge))[c(
  8:23, 40, 1:5, 38, 41, 6:7, 35, 25:34, 42:43, 36:37, 24, 39, 44:45
)]
length(levelsall) # [1] 45

Idents(dge) <- factor(Idents(dge), levels = levelsall)

### Figure S10b – Top 5 unique markers for each cell type
X <- dge
n <- 100

M1 <- pblapply(levels(Idents(X)), function(cID) {
  FC <- FoldChange(X, ident.1 = cID) # this does not consider p-value
  FC <- FC[order(FC$avg_log2FC, decreasing = TRUE), ]
  rownames(FC[seq_len(n), ])
})

M1

dup <- unlist(M1)[duplicated(unlist(M1))]
M2 <- M1

for (ii in seq_along(M1)) {
  M2[[ii]] <- M2[[ii]][which(!(M2[[ii]] %in% dup))]
  M2[[ii]] <- M2[[ii]][1:5]
}

M <- unique(unlist(M2))

options(warn = -1) # Temporarily suppress warnings
O <- DotPlot(X, features = M, col.min = 0, scale.by = "size") +
  scale_color_gradient2(low = rgb(1, 1, 1, 0), mid = rgb(1, 1, 1, 0)) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = 3)
  )
options(warn = 0) # Reset warning options to default

png(
  paste0("marker_files/Dotplot_MP5_unique.png"),
  width = 8000,
  height = 2200,
  res = 300
)

print(
  O +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    scale_size(range = c(1, 3))
)

dev.off()

### Figure S10a – Top 5 unique markers for each broad cell type
# merge similar cell types
dge <- P200
levels(Idents(dge))

levelsall <- levels(Idents(dge))[c(
  8:23, 40, 1:5, 38, 41, 6:7, 35, 25:34, 42:43, 36:37, 24, 39, 44:45
)]

Idents(dge) <- factor(Idents(dge), levels = levelsall)
levels(Idents(dge))

levels(Idents(dge)) <- gsub(":.*", "", levels(Idents(dge)))
levels(Idents(dge)) <- gsub("Early KC.*", "Early KC", levels(Idents(dge)))
levels(Idents(dge)) <- gsub("Placode.*", "Placode", levels(Idents(dge)))
levels(Idents(dge)) <- gsub("Satellite muscle", "Myocyte", levels(Idents(dge)))
levels(Idents(dge)) <- gsub("SkM myocyte/myoblast", "Myocyte", levels(Idents(dge)))

levels(Idents(dge))
#  [1] "Fb"            "Pericyte"      "Chondrocytes"  "EC"
#  [5] "NCProgenitors" "Periderm"      "Early KC"      "KC"
#  [9] "Placode"       "Melanocyte"    "Merkel cell"   "Glial"
# [13] "Neuronal"      "Myocyte"

# Select top 5 unique markers not repeated in any other cell type
X <- dge
n <- 100

M1 <- pblapply(levels(Idents(X)), function(cID) {
  FC <- FoldChange(X, ident.1 = cID) # this does not consider p-value
  FC <- FC[order(FC$avg_log2FC, decreasing = TRUE), ]
  rownames(FC[seq_len(n), ])
})

M1

dup <- unlist(M1)[duplicated(unlist(M1))]
M2 <- M1

for (ii in seq_along(M1)) {
  M2[[ii]] <- M2[[ii]][which(!(M2[[ii]] %in% dup))]
  M2[[ii]] <- M2[[ii]][1:5]
}

M <- unique(unlist(M2))

options(warn = -1) # Temporarily suppress warnings
O <- DotPlot(X, features = M, col.min = 0, scale.by = "size") +
  scale_color_gradient2(low = rgb(1, 1, 1, 0), mid = rgb(1, 1, 1, 0)) +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = 3)
  )
options(warn = 0) # Reset warning options to default

png(
  paste0("marker_files/Dotplot_MP5_merged_unique.png"),
  width = 3000,
  height = 1000,
  res = 300
)

print(
  O +
    theme(panel.border = element_blank(), panel.grid = element_blank()) +
    scale_size(range = c(1, 3))
)

dev.off()
