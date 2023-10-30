# -------------------------------------------------------------- #
# Macrophages Subclustering Analysis after final Cell Filtering  #  # Perpetrated by Manu Moradiellos
# -------------------------------------------------------------- #

# Main code of the macrophages subclustering and marker genes expression
# analysis. This was done to finish the single collaboration and these
# were the last tasks proposed by Jose Luis.

# |- Packages -|  ----
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('Seurat'))
suppressPackageStartupMessages(library('dplyr'))
suppressPackageStartupMessages(library('qusage'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('slingshot'))
suppressPackageStartupMessages(library('monocle3'))
suppressPackageStartupMessages(library('SeuratWrappers'))

# |- Import complete or macrophages Seurat Object -| ----
## Read single-cell experiment analyzed with Seurat by Luis, but with some further Cell Filtering to remove low quality samples and some high-readings
seurat <- readRDS('./seurat_obj.hypoxia_breast_cancer.2500_27500umis.50hk.RDS')

# |- 1. Macrophages Subclustering Analysis -| --------------
macrophages <- subset(seurat, ident = 'Macrophages')  # Focus on cell cluster of interest
marophages <- NormalizeData(macrophages, normalization.method = 'LogNormalize', scale.factor = 10000)
macrophages <- FindVariableFeatures(macrophages)      # Identify features high cell-to-cell variation, focus on those genes
macrophages <- RunPCA(macrophages, verbose = F)       # Allows for removal of noise in features
if (FALSE) {                                           # Used to find the appropriate number of PCs (only do it one time)
  macrophages <- JackStraw(macrophages, num.replicate = 100) 
  macrophages <- ScoreJackStraw(macrophages, dims = 1:20)
  JackStrawPlot(macrophages, dims = 1:20)   # Checking PCs and their p-values
}
macrophages <- FindNeighbors(macrophages, dims = 1:12) 
macrophages <- FindClusters(macrophages, resolution = 0.65) # Lower resolution to find less clusters, return similar number of clusters as expected
macrophages_oldumap <- macrophages  # Save copy before running again the UMAP on subcluster, just to compare to the main experiment UMAP
macrophages <- RunUMAP(macrophages, dims = 1:12) 
if (FALSE) { # Plot results
  DimPlot(macrophages, label = T, label.size = 6)
  DimPlot(macrophages_oldumap, label = T, label.size = 6) # Check distribution on all cells' UMAP
}

# Add macrophages subclusters info to original seurat object initial step, no refining on clusters
seurat$sub_cluster <- as.character(Idents(seurat))
seurat$sub_cluster[Cells(macrophages)] <- paste("Macro_",Idents(macrophages))
  
# Save subsetted and "reannotated" object 
saveRDS(macrophages, paste0(working_directory,'output/1.initial_subclustering_no.manual.corrections/', 'macrophages_initialclustering_nomanualrefining.sub_seurat_obj.RDS'))
saveRDS(macrophages_oldumap, paste0(working_directory, 'output/1.initial_subclustering_no.manual.corrections/','macrophages_oldumap_initialclustering_nomanualrefining.sub_seurat_obj.RDS'))
saveRDS(seurat, paste0(output_dir, 'seurat_new_annotation_macrophages.sub.RDS'))


# |- 2. Number of macrophages per condition and sub cluster |----
count_table <- as.data.frame.matrix(addmargins(table(macrophages_oldumap@active.ident, macrophages_oldumap@meta.data$hypoxia_label))) 

# |- 3. Use Known Macrophages Markers -| ----
## We can use previous knowledge to plot our markers 
## for the 'classical' macrophages types (a bit outdated)
# General Macrophages Markers (not useful for subclustering)
general_markers <- c('Csf1', 'Itgam', 'Adgre1', 'Cd68')
# M1 Markers (Pro-inflammatory)
m1_markers <- c('Cd86', 'Cd80', 'Stat1', 'Irf7', 'Isg15')
# M2 Markers (Anti-inflammatory, promote tumor cell proliferation)
m2_markers <- c('Vegfa', 'Tnf', 'Msr1', 'Mrc1', 'Tgfb1') ## Il10 not expressed, Msr1 is Cd204
# TAMs Markers (Tumor Associated Macrophages)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6472943/ Breast Human but some mentioned for mouse
tams_markers <- c('Ccl8', 'Siglec1', 'Ccl2', 'Cd163') #Ccl2 TAM recruitment, Cd163 TAM marker; Ccl8 Siglec1 in human breast

## Markers coming from the TAMs diversity review https://doi.org/10.1016/j.it.2022.04.008 
ifn_tams <- c('Stat1', 'Cxcl10', 'Ccl2', 'Ccl7', 'Ifit1', 'Ifit2', 'Ifit3', 'Ifitm3') # Interferon-primed
reg_tams <- c('Cd86', 'Cx3xr1', 'Itga4', 'Lgals9', 'Tgfbr1', 'Tgfbr2', 'H2-Aa', 'Apoe') # Immune Regulatory 
inflam_tams <- c('Cxcl1', 'Ccl3', 'Il1b') # Inflammatory cytokine-enriched
la_tams <- c('Apoe', 'Cd63', 'C1qa', 'Ctsb') # Lipid-associated 
angio_tams <- c('Adam8', 'Bnip3') # Pro-angiogenic
prolif_tams <- c('Cdk1', 'Mki67', 'Stmn1', 'Top2a', 'Tubb2a', 'Tubb5') # Proliferating TAMs

tims <- c('Ccr2', 'Ccl9', 'Fn1') # Tumor-infiltrating monocytes

FeaturePlot(macrophages, features = tims, order = T) # Use to plot list of markers

# |- 4. Unsupervised Macrophages Markers -| ----
## Obtain table of unsupervised markers
macro_markers <- FindAllMarkers(macrophages, min.pct = 0.25, logfc.threshold = 0.25) 
macro_markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 2, order_by = avg_log2FC) # Find two top positive markers per clusters
write_tsv(x = macro_markers,
          file = paste0(working_directory,'output/1.initial_subclustering_no.manual.corrections/', 'macro.subcluster_resolution_0.5_unsupervised_markers.tsv'))

## Heatmap of top 10 markers per subcluster
macro_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10
top.markers_heat <- DoHeatmap(macrophages, features = top10$gene) + NoLegend()
ggsave(filename = paste0(working_directory,'output/1.initial_subclustering_no.manual.corrections/', 'top.markers_macro.sub_resolution_0.5_unsupervised_markers_heatmap.png'),
       plot = top.markers_heat, dpi = 300)