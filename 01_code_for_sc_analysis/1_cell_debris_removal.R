# --------------------------------------------#
# Quality Control and Cellular Debris Removal #  # Perpetrated by Manu Moradiellos
# --------------------------------------------#

# This scripts takes one of the outputs of Luis preliminary's analysis 
# and performs some further quality control steps as well as a reannotation
# just in case it was needed

# |- Packages -|  ----
set.seed(666)
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('Seurat'))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('qusage'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('patchwork'))
suppressPackageStartupMessages(library('SeuratWrappers'))
suppressPackageStartupMessages(library('viridis'))
suppressPackageStartupMessages(library('scales'))
suppressPackageStartupMessages(library('RColorBrewer'))
suppressPackageStartupMessages(library('ggpubr'))
suppressPackageStartupMessages(library('ggridges'))

# Load previous filtered object
seurat <- readRDS('./input/seurat_new_annotation.RDS')

# 1. Applying more filtering to previous object ----
## Tirosh et al.'s House-keeping genes 
hk_genes <- scan(file = '/resources/housekeeping.gene_cell_groups.txt',
                 what = 'character', skip = 2) 
hk_genes.found_macrophages <- which(toupper(rownames(seurat@assays$RNA)) %in% hk_genes) # Remove HKgenes that were not found
n.expressed.hk_genes <- Matrix::colSums(seurat@assays$RNA@data[hk_genes.found_macrophages, ] > 0) # Number of HK genes being expressed
seurat <- AddMetaData(object = seurat, metadata = n.expressed.hk_genes, col.name = 'n.exp.hk_genes') 
VlnPlot(seurat, features = c('nFeature_RNA', 'nCount_RNA', 'n.exp.hk_genes')) # Plot all qc features

## New filtering, I decided on 2500 < UMIS < 27500 & HK Genes > 50
seurat.2500_27500umis.50hk <- subset(x = seurat,
                                     subset = nCount_RNA < 27500 &
                                       nCount_RNA > 2500 & n.exp.hk_genes > 50)
VlnPlot(seurat.2500_27500umis.50hk, features = c('nFeature_RNA', 'nCount_RNA', 'n.exp.hk_genes')) # Compare changes

# 2. Renormalization and annotation of the newly filtered object

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(limma))

set.seed(1)

# 2.1. Normalization ----
## Global-scaling normalization method that normalizes the feature expression
## for each cell by the total expression, multiplies this by a scale factor of 10000
## the log-transforms the result
seurat.2500_27500umis.50hk <- NormalizeData(seurat.2500_27500umis.50hk, normalization.method = 'LogNormalize', scale.factor = 10000)

# This is not mandatory but maybe can be used to get an idea of what's going on
## It is also recommended by seurat's authors to check these subsets of features
## that exhibit high cell-to-cell variation in the dataset as focusing on these
## genes in downstream analysis helps to highlight possible biological signals 
seurat.2500_27500umis.50hk <- FindVariableFeatures(seurat.2500_27500umis.50hk, selection.method = 'vst', nfeatures = 2500) 
top10 <- head(VariableFeatures(seurat.2500_27500umis.50hk), 10)
p1 <- VariableFeaturePlot(seurat.2500_27500umis.50hk) + theme(legend.position = 'bottom')
p1 <- p1 + LabelPoints(plot = p1, points = top10, repel = T) + theme(legend.position = 'bottom')
p1

# 2.2. Scaling data to perform PCA ----
## Applying a linear transformation/scaling is a standard pre-processing step 
## prior to dimensional reduction techniques such as PCA
seurat.2500_27500umis.50hk <- ScaleData(seurat.2500_27500umis.50hk, features = rownames(seurat.2500_27500umis.50hk), vars.to.regress = NULL)

# 2.3. Perform PCA ----
## PCA calculated only using highly variable genes, with a total of 50 PCs to compute
seurat.2500_27500umis.50hk <- RunPCA(seurat.2500_27500umis.50hk, features = VariableFeatures(seurat.2500_27500umis.50hk) , npcs = 50)
p4 <- ElbowPlot(seurat.2500_27500umis.50hk, ndims = 50) + theme(legend.position="bottom")
p4 # Elbow Plot to see PCs, around ~16 we see reduction in p.c.e
seurat.2500_27500umis.50hk <- RunUMAP(seurat.2500_27500umis.50hk, reduction = 'pca', dims = 1:16, n.components = 2) 
DimPlot(seurat.2500_27500umis.50hk)

# 2.4. Cell annotation with SingleR ----
seurat_SCE <- as.SingleCellExperiment(seurat.2500_27500umis.50hk)
immgen <- ImmGenData(ensembl = F) # Get Immune cells information

cell_pred <- SingleR(test = seurat_SCE, ref = immgen, 
                     labels = immgen$label.fine, assay.type.test = 1) # Obtain annotation at cell level using immune data reference
seurat.2500_27500umis.50hk[['SingleR.labels.final']] <- cell_pred$labels # Add annotation
seurat.2500_27500umis.50hk$SingleR.labels.final.simplified <- sapply(seurat.2500_27500umis.50hk$SingleR.labels.final, function(label){
  gsub(strsplit(label, split = "\\(")[[1]][1], pattern = " $", replacement = "")
} ) # Annotation has very specific types of immune cells, so we just keep the general names (i.e., 'T cell (T.4FP3+25+)' to just 'T cell')

# Those annotation groups with really low number of cells are stored under the same umbrella group
## After filtering, they seem to be only stromal
mixed_group <-  names(table(seurat.2500_27500umis.50hk$SingleR.labels.final.simplified)[table(seurat.2500_27500umis.50hk$SingleR.labels.final.simplified) < 10]) 
seurat.2500_27500umis.50hk$SingleR.labels.final.simplified = sapply(seurat.2500_27500umis.50hk$SingleR.labels.final.simplified, function(x){
  ifelse(x %in% mixed_group, "Mixed group", x)
})

# 2.5. Further cell annotation with SingleR, correcting previous labels ----
## This section is to make some first changes to SingleR labelling as it
## is a bit too descriptive, so we simplify it to just keep general namings

# We create a new vector to store the fine-tuned labels
seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.corrected <- seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.simplified

# We got to differentiate better the different kinds of T-cells available,
# so we go back to the original labels proposed by SingleR and find those
# entries matching any T-cell and separate them according to basic markers/indicators
# appearing in the name. Those that cannot be identified from the get-go
# where decided to be grouped under other immature/early T-cells
t_cell_labels <- grep(pattern = "^T cell", x = unique(seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final), value = T) ## Value = TRUE to obtain the element itself and not its position
t_cells_labels_CD4 = grep(pattern = "T\\.CD4|T\\.4", t_cell_labels, value = T)
t_cells_labels_CD8 = grep(pattern = "T\\.CD8|T\\.8", t_cell_labels, value = T)
t_cells_labels_regs = grep(pattern = "Tregs|T\\.4FP3\\+25\\+)", t_cell_labels, value = T)
t_cells_labels_extra = t_cell_labels[!(t_cell_labels %in% c(t_cells_labels_CD4, t_cells_labels_CD8, t_cells_labels_regs))]

# We add this cells to the new meta data level, changing the group name accordingly
seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.corrected[seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final %in% t_cells_labels_CD4] <- 'T cells CD4'
seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.corrected[seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final %in% t_cells_labels_CD8] <- 'T cells CD8'
seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.corrected[seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final %in% t_cells_labels_regs] <- 'Tregs'
seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.corrected[seurat.2500_27500umis.50hks@meta.data$SingleR.labels.final %in% t_cells_labels_extra] <- 'T cells immature-early'

# Also, Jose Luis asked to merge the microglia cells to the macrophages group
microglia_cells <- grep(pattern = 'Microglia', unique(seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final), value = T)
seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final.corrected[seurat.2500_27500umis.50hk@meta.data$SingleR.labels.final %in% microglia_cells] <- 'Macrophages'

# 2.6. Compare UMAPs between the two labelings to check for possible differences in the 'clusters'
asa_palette <- c("B cells" = "#426600", "DC" = "#FF0010","Endothelial cells" = "#4C005C", "Epithelial cells" = "#FFCC99",
                 "Fibroblasts" = "#8F7C00", "ILC" = "#808080", "Macrophages" = "#9DCC00", "Mixed group" = "#F0A0FF", "Monocytes" = "#005C31", 
                 "Neutrophils" = "#191919", "NK cells" = "#C20088", "NKT" = "#FFA405","Stem cells" = "#003380",
                 "T cells CD4" = "#FFA8BB", "T cells CD8" = "#993F00", "T cells immature-early" = "#94FFB5", "Tgd" = "#2BCE48", "Tregs" = "#0075DC") # Pls use this for the UMAPs

DimPlot(seurat.2500_27500umis.50hk, cols = asa_palette, group.by = 'SingleR.labels.final.corrected', pt.size = 0.3) 

# 2.8. Obtain heatmap of principal markers of each group obtained by SingleR and the simplifying method ----
## Find markers for every cluster compared to all remaining cells, reporting just the positive ones
seurat.cell_types.markers <- FindAllMarkers(seurat.2500_27500umis.50hk, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top10markers_cluster <- seurat.cell_types.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC) 
top5markers_cluster <- seurat.cell_types.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC) # To make it a more readable plot
DoHeatmap(seurat.2500_27500umis.50hk, features = top5markers_cluster$gene) + NoLegend()

write.table(top10markers_cluster, sep = '\t', quote = F, file = '/output/top10genes_per_cluster.tsv')

# 2.9. Change Idents to new ones and save the object
Idents(seurat.2500_27500umis.50hk) <- seurat.2500_27500umis.50hk$SingleR.labels.final.corrected
saveRDS(seurat.2500_27500umis.50hk, file = '/output/seurat.2000_27500umis_50hk.idents.RDS')
