### InferCNV

library(Seurat)
library(infercnv)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
# Loading RDS, Dark Vs. Light
# lvl4 annotation transferred from Gopee 2024, renamed Melaoblast to Melanocyte
group.combined <- readRDS("P200_NewFetalSkinAnnotated2_Melanocyte.rds")

# Define 'normal' reference cells and 'case' cells
# Light: normal reference
# Dark:  case
case_samples <- "6"
reference_samples <- "4"

##### Dark Vs Light
condition <- "DarkVsLight"
cat("Processing ", condition, "...\n")

# Create InferCNV object
DarkVsLight <- group.combined
counts_matrix <- DarkVsLight@assays$RNA@counts

# Create Cell_id annotation file
cell_id_anno <- cbind(rownames(DarkVsLight@meta.data), DarkVsLight@meta.data$group)
fname_cellanno <- paste0("InferCNV/cell_anno_", condition, ".txt")
write.table(cell_id_anno, file = fname_cellanno, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Creating input object for inferCNV
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_matrix,
  annotations_file = fname_cellanno,
  delim = "\t",
  gene_order_file = "InferCNV/GRCh38-2020-A_gene_order_file.txt",
  ref_group_names = reference_samples
)

# ------------------------------------------------------------------------------
# Figure S12a InferCNV
# ------------------------------------------------------------------------------
# Set options to avoid scientific notation, especially for subclusters
options(scipen = 100)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = paste0("InferCNV/samples_", condition),
  cluster_by_groups = TRUE,
  analysis_mode = "samples",
  denoise = TRUE,
  HMM = TRUE,
  HMM_type = "i6",
  num_threads = 60,
  no_plot = FALSE
)

saveRDS(infercnv_obj, file = paste0("InferCNV/samples_", condition, ".rds"))

# ------------------------------------------------------------------------------
# Figure S12b overlap between differential expression and InferCNV results
# ------------------------------------------------------------------------------
# load cell types that have at least 0.1% representation in both L'Oreal P200 Light Vs. Dark samples with cell types transferred from lvl4_annotation of new organoids-fetal skin reference data
celltype <- read.table("P200_NewFetalSkinlvl4_.001celltype.txt", sep = "\t", stringsAsFactors = FALSE)[, 1]
length(celltype) # 25

# load genes in Chr 1 gain region
gain <- read.table("../InferCNV/samples_DarkVsLight/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_genes.dat", header = TRUE)
table(gain$state)

gaingenes <- gain$gene

### Differential expression in each cell type with at least 0.1% representation in both samples
markerslist <- list()
geneslist <- list()
genes2list <- list()
for (i in seq_along(celltype)) {
  label <- celltype[i]

  gene <- read.table(
    paste0("markers/P200_LvsD_lvl4_C", i, "_de_minpct0.1_fc1.5_12.2023.txt"),
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE,
    sep = "\t"
  )

  geneslist[[i]] <- gene
  gene1 <- rownames(gene[which(gene$avg_log2FC > 0), ]) # higher in Light
  gene2 <- rownames(gene[which(gene$avg_log2FC < 0), ]) # higher in Dark

  gene2gain <- gene2[gene2 %in% gaingenes]
  print(c(i, label, length(gene2), length(gene2gain)))

  write.table(
    t(c(i, label, length(gene2), length(gene2gain), gene2gain)),
    "../InferCNV/samples_DarkVsLight/DEG1.5fc_dark_chr1gain.txt",
    quote = FALSE,
    sep = "\t",
    append = TRUE,
    row.names = F,
    col.names = F
  )

  gene1gain <- gene1[gene1 %in% gaingenes]
  print(c(i, label, length(gene1), length(gene1gain)))

  write.table(
    t(c(i, label, length(gene1), length(gene1gain), gene1gain)),
    "../InferCNV/samples_DarkVsLight/DEG1.5fc_light_doublecheck_chr1gain.txt",
    quote = FALSE,
    sep = "\t",
    append = TRUE,
    row.names = F,
    col.names = F
  )

  gene <- read.table(
    paste0("markers/P200_LvsD_lvl4_C", i, "_de_minpct0.1_log2fc0.25_12.2023.txt"),
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE,
    sep = "\t"
  )

  print(c(length(which(gene$avg_log2FC > 0)), length(which(gene$avg_log2FC < 0))))
  genes2list[[i]] <- gene

  gene1 <- rownames(gene[which(gene$avg_log2FC > 0), ])
  gene2 <- rownames(gene[which(gene$avg_log2FC < 0), ])

  gene2gain <- gene2[gene2 %in% gaingenes]
  print(c(i, label, length(gene2), length(gene2gain)))

  write.table(
    t(c(i, label, length(gene2), length(gene2gain), gene2gain)),
    "../InferCNV/samples_DarkVsLight/DEG0.25logfc_dark_chr1gain.txt",
    quote = FALSE,
    sep = "\t",
    append = TRUE,
    row.names = F,
    col.names = F
  )

  gene1gain <- gene1[gene1 %in% gaingenes]
  print(c(i, label, length(gene1), length(gene1gain)))

  write.table(
    t(c(i, label, length(gene1), length(gene1gain), gene1gain)),
    "../InferCNV/samples_DarkVsLight/DEG0.25logfc_light_doublecheck_chr1gain.txt",
    quote = FALSE,
    sep = "\t",
    append = TRUE,
    row.names = F,
    col.names = F
  )
}
