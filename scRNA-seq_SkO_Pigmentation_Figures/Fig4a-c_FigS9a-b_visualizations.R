# Generate UMAPs, sample-split UMAPs, cell-type proportion barplots,
# and melanocyte marker FeaturePlots

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(pbapply)
library(Nebulosa)
library(ComplexHeatmap)
library(EnhancedVolcano)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
# Loading RDS, Dark vs. Light
# lvl4 annotation transferred from Gopee 2024, renamed "Melaoblast" to "Melanocyte"
Y <- readRDS("P200_NewFetalSkinAnnotated2_Melanocyte.rds")

# ------------------------------------------------------------------------------
# Figure 4a – UMAP
# ------------------------------------------------------------------------------
png(filename = "Annotation/P200_NewFetalSkinAnnotated1_noaxes_biggerdots.png", height = 2000, width = 2000, res = 300)

UMAPPlot(Y, label = TRUE, repel = TRUE, pt.size = 1) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 0),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("")

dev.off()

png(filename = "Annotation/P200_NewFetalSkinAnnotated1_noaxes_SampleSplit.png", width = 1500 * length(unique(Y$orig.ident)), height = 1550, res = 300)

DimPlot(
  object = Y,
  label = FALSE,
  pt.size = 1,
  reduction = "umap",
  split.by = "orig.ident",
  ncol = length(unique(Y$orig.ident))
) +
  theme_light() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0),
    axis.text.y = element_text(size = 0),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(colour = NA, fill = NA),
    strip.text = element_text(color = "black", size = 16),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("")

dev.off()

# ------------------------------------------------------------------------------
# Figure 4b – Bar plot of cell type proportion
# ------------------------------------------------------------------------------
CP <- table(Y$orig.ident, Idents(Y))
CP <- (CP / rowSums(CP)) * 100
CP <- reshape2::melt(CP)
colnames(CP) <- c("Sample", "CT", "Proportion")

png(filename = "P200_NewFetalSkinAnnotatedCP.png", width = 3500, height = 1200, res = 300)

p <- ggplot(CP, aes(Proportion, Sample, fill = CT)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank()
  ) +
  labs(fill = "NewFetalSkin\nCell Type")

print(p)
dev.off()

# ------------------------------------------------------------------------------
# Figure 4c – FeaturePlot (melanocyte markers)
# ------------------------------------------------------------------------------
genes <- c("TYRP1", "DCT", "MLANA")

png(filename = "Annotation/P200-Feature-melanocyte_red.png", res = 300, height = 800 * 3, width = 800 * 2)

FeaturePlot(Y, features = genes, split.by = "orig.ident") & NoLegend() + NoAxes()

dev.off()

# ------------------------------------------------------------------------------
# Figure S9a – UMAP transferred from Gopee 2024
# ------------------------------------------------------------------------------
Y2 <- readRDS(file = "P200_NewFetalSkinAnnotated.rds")

png(filename = "P200_NewFetalSkinAnnotated.png", width = 2000, height = 2000, res = 300)

p <- UMAPPlot(Y2, label = TRUE, repel = TRUE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

print(p)
dev.off()

# ------------------------------------------------------------------------------
# Figure S9b – UMAP split by Light and Dark
# ------------------------------------------------------------------------------
png(filename = "P200_NewFetalSkinAnnotated-UMAP.png", width = 3000, height = 1500, res = 300)

DimPlot(
  object = Y2,
  label = TRUE,
  repel = TRUE,
  reduction = "umap",
  split.by = "orig.ident"
) +
  theme_light() +
  theme(
    legend.position = "none",
    strip.text = element_text(color = "black")
  ) +
  xlab("UMAP 1") +
  ylab("UMAP 2")

dev.off()
