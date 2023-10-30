# ---------------------------------------#
# Main Figures for Single Cell Analysis  #  # Perpetrated by Manu Moradiellos
# ---------------------------------------#

# Collection of all the single-cell analysis' plots made for the publication,
# most of these were later slightly edited with image editing software

# |- SETUP -| ----
## 1. Load Packages ----
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
suppressPackageStartupMessages(library('magick'))
suppressPackageStartupMessages(library('SCpubr')) 
# devtools::install_github('enblacar/SCpubr', ref = 'v1.1.2-dev-stable')

# For CellChat
suppressPackageStartupMessages(library('Seurat'))
suppressPackageStartupMessages(library('CellChat'))
suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('patchwork'))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library('svglite'))
suppressPackageStartupMessages(library('ComplexHeatmap'))
suppressPackageStartupMessages(library('NMF'))
suppressPackageStartupMessages(library('ggalluvial')) 

## A set of slightly modified functions from CellChat
source('/resources/netVisual_bubble_mod.R') 
source('/resources/ComplexHeatmap_heatmap_mod.R')
source('/resources/netVisual_heatmap_mod_v2.R')
source('/resources/cellchat_rankNet_mod.R')
source('/resources/netAnalysis_signalingRole_heatmap_mod.R')

## 2. Seurat objets and output dir setup -----
# Directory to save all plots and share them with JL
outdir <- '/output_plots/'

# Load object containing whole experiment
seurat <- readRDS('/input/seurat.2500_27500umis.50hk.RDS')

# Load object containing the macrophages subclustering analysis result (new subgroups and UMAP)
macro_obj <- readRDS('/macrophages_seurat_object.various_resolutions.RDS')
macro_obj$macro_sub_0.5 <- macro_obj$RNA_snn_res.0.5 # Resolution of choice for our experiment
levels(macro_obj$macro_sub_0.5) <- c('TAM01', 'TAM02', 'TAM03', 'TAM04', 'TAM05', 'TAM06')

# Load CellChat object, long name to indicate various processing steps 
# performed when converting the Seurat object into a CellChat one
cellchat_breast <- readRDS(file = '/cellchat_object.seurat.2500_27500umis.50hk_final.idents.cellclustersize_filtered.pop_size_T.RDS') 
cellchat_merged <- mergeCellChat(cellchat_breast, add.names = c('control', 'low', 'high'), merge.data = T) 
cellchat_merged_fixed <- liftCellChat(cellchat_merged, group.new = levels(cellchat_merged@idents$joint))

## 3. Colors used ----
# Color palette for UMAP of all cell types
asa_palette <- c('B cells' = '#426600',
                 'DC' = '#FF0010',
                 'Endothelial cells' = '#4C005C',
                 'Epithelial cells' = '#FFCC99',
                 'Fibroblasts' = '#8F7C00',
                 'ILC' = '#808080',
                 'Macrophages' = '#9DCC00',
                 'Mixed group' = '#F0A0FF',
                 'Monocytes' = '#005C31',
                 'Neutrophils' = '#191919',
                 'NK cells' = '#C20088',
                 'NKT' = '#FFA405',
                 'Stem cells' = '#003380',
                 'T cells CD4' = '#FFA8BB',
                 'T cells CD8' = '#993F00',
                 'T cells immature-early' = '#65c8db',
                 'Tgd' = '#2BCE48',
                 'Tregs' = '#0075DC'
                   )

# Color palette for Hypoxia Levels 
hypoxia_colors <- c('control' = '#1C1C1C',
                    'low' = '#579CFF',
                    'high' = '#FF5757'
                      )

# Color palette for UMAP of macrophages
macrofalco <- c('TAM01' = '#5BC0BE',
                'TAM02' = '#DC4141',
                'TAM03' = '#820263',
                'TAM04' = '#0AAE59',
                'TAM05' = '#F77F00',
                'TAM06' = '#004266'
                  )

#|- PLOTS -| ----

## 1. Whole Experiment UMAP (Fig. 4E and Supp. 11B Right) ----
p1 <- SCpubr::do_DimPlot(seurat,
                         colors.use = asa_palette,
                         pt.size = 0.3,
                         plot_cell_borders = F,
                         font.size = 12)

SCpubr::save_Plot(plot = p1,
                  figure_path = paste0(outdir, '1_umap_general/'),
                  file_name = '1_General_UMAP',
                  dpi = 300,
                  output_format = 'all')

p2 <- SCpubr::do_DimPlot(seurat,
                         split.by = 'SingleR.labels.new.corrected',
                         ncol = 5,
                         colors.use = asa_palette,
                         pt.size = 0.5,
                         plot_cell_borders = T,
                         border.size = 0.6,
                         na.value = 'grey90',
                         legend.position = 'none',
                         font.size = 9)
SCpubr::save_Plot(plot = p2,
                  figure_path = paste0(outdir, '2_umap_separate_clusters/'),
                  file_name = '2_Separate_Clusters_UMAP',
                  dpi = 300,
                  output_format = 'all',
                  width = 10)

# 3. Hypoxia Level UMAP (Supp. 11C) ----
## Have to choose one of the resulting plots...
p3_nonshuffled <- SCpubr::do_DimPlot(seurat,
                                  colors.use = hypoxia_colors,
                                  group.by = 'hypoxia_label',
                                  pt.size = 0.3,
                                  plot_cell_borders = F,
                                  font.size = 12,
                                  shuffle = F,
                                  order = c('high', 'low', 'control'))
SCpubr::save_Plot(plot = p3_nonshuffled,
                  figure_path = paste0(outdir, '3_umap_hypoxia_levels/'),
                  file_name = '3_Hypoxia_Levels_UMAP_notshuffled',
                  dpi = 300,
                  output_format = 'all')

## 4. Top 5 Markers Whole Experiment (Supp. 11B Left) ----
# Obtain differentially expressed genes per clusters to highlight as possible 'markers'
seurat.cell_types.markers <- FindAllMarkers(seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top5markers_cluster <- seurat.cell_types.markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC) 

# For visualization purposes we downsample the cell to show on the heatmap
downsampled_orderd_heatmap <- DoHeatmap(
  subset(seurat, downsample = 100),
  features = top5markers_cluster$gene,
  group.colors = asa_palette,
  size = 4
  ) +
  scale_fill_viridis(na.value = 'white') +
  NoLegend()

SCpubr::save_Plot(plot = downsampled_orderd_heatmap,
                  figure_path = paste0(outdir, '4_topmarkers_heatmap/'),
                  create_path = T,
                  file_name = '4_Top5markers_whole_experiment_downsamples100cells_withlegend',
                  dpi = 300,
                  output_format = 'all',
                  width = 15,
                  height = 15)

# 5. Macrophages Subclustering UMAP (Supp. 12B) ----
macro_umap1 <- SCpubr::do_DimPlot(macro_obj,
                   colors.use = macrofalco,
                   group.by = 'macro_sub_0.5',
                   pt.size = 1, 
                   plot_cell_borders = F,
                   font.size = 12,
                   shuffle = F)
SCpubr::save_Plot(macro_umap1,
                  figure_path = paste0(outdir, '5_macro_general'),
                  file_name = '5_macro_general_umap',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 10)

# 6. Macrophages Subclustering UMAP Separate Clusters (Supp. 12B) ----
macro_umap2 <- SCpubr::do_DimPlot(macro_obj,
                   split.by = 'macro_sub_0.5',
                   ncol = 3,
                   colors.use = macrofalco,
                   pt.size = 0.75,
                   plot_cell_borders = T,
                   border.size = 0.6,
                   na.value = 'grey90',
                   legend.position = 'none',
                   font.size = 12)
SCpubr::save_Plot(macro_umap2,
                  figure_path = paste0(outdir, '6_macro_separate_clusters'),
                  file_name = '6_macro_separate_clusters_umap',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 15)

# 7. Top 5 Markers Macrophages Subclusters (Supp. 12A) ----
Idents(macro_obj) <- macro_obj$macro_sub_0.5 
macro_markers <- FindAllMarkers(macro_obj, min.pct = 0.25, logfc.threshold = 0.25) 
macro_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) -> top5
top.markers_heat <- DoHeatmap(macro_obj, features = top5$gene) +
  scale_fill_viridis(na.value = 'white') +
  NoLegend()

SCpubr::save_Plot(top.markers_heat,
                  figure_path = paste0(outdir, '7_macro_subclusters_heatmap'),
                  file_name = '7_macro_subclusters_heatmap',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 21
                  ) 

# 8. Cell-cell Communication Heatmap (Fig. 5G and Fig. 5H) ----
# Just comparing Low vs High

### 8.1. Heatmaps ----
  ## Colors are the same as for networks. The top colored bar plot represents the
  ## sum of column of values displayed in the heatmap (incoming signal),
  ## the right colored bar plot represents the sum of row of values (outgoing)

heat_count <- netVisual_heatmap(cellchat_merged_fixed, comparison = c(2,3), remove.isolate = F,
                                    cluster.rows = F, title.name = paste0('Differential Number of Interactions ', ' Low vs. High'),
                                    font.size.title = 14, font.size = 12,
                                color.use = asa_palette) 
SCpubr::save_Plot(ComplexHeatmap::draw(heat_count),
                  figure_path = paste0(outdir, '8_low_vs_high_comms_diffs/heatmaps/'),
                  file_name = '8_low_vs_high_comms_diffs_count_ht_nonclustered',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 8,
                  height = 8
) 

heat_strength <- netVisual_heatmap(cellchat_merged_fixed, comparison = c(2,3), remove.isolate = F, measure = 'weight',
                                       cluster.rows = F, title.name = paste0('Differential Interaction Strength ', ' Low vs. High'),
                                       font.size.title = 14, font.size = 12, color.heatmap = c('#47115d', '#fdc330'),
                                   color.use = asa_palette) 

SCpubr::save_Plot(ComplexHeatmap::draw(heat_strength),
                  figure_path = paste0(outdir, '8_low_vs_high_comms_diffs/heatmaps/'),
                  file_name = '8_low_vs_high_comms_diffs_strength_ht_nonclustered',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 8,
                  height = 8
) 

# 9. Information Flow Plots (Fig. 5I)----
# Creates and save in correct format
to_test <- list('Low vs. High' = c(2,3), 'Control vs. Low' = c(1, 2), 'Control vs. High' = c(1,3))
colors_test <- list('Low vs. High' = c(hypoxia_colors['low'], hypoxia_colors['high']),
                    'Control vs. Low' = c(hypoxia_colors['control'], hypoxia_colors['low']),
                    'Control vs. High' = c(hypoxia_colors['control'], hypoxia_colors['high']))
comparison <- 1 # For the different comparisons to be made, 

## With measure = 'weight' or 'count' we can compare either the total interaction strength or the total number
## Remember that for 'weight' you need to obtain scaled and not-scaled results using show.raw = F or T
## All of this should be made into a better function but I don't have the time now

### 9.1.1 Low vs High. Weight/Interaction Strength not scaled ----
info_flow_group1 <- rankNet_mod(cellchat_merged_fixed, mode = "comparison",
                                measure = 'weight', stacked = T, do.stat = T,
                                comparison = to_test[[1]], font.size = 10,
                                show.raw = T, return.data = T, sources.use = c(1:18),
                                targets.use = c(1:18), color.use = colors_test[[1]])

info_flow_group2 <- rankNet_mod(cellchat_merged_fixed, mode = "comparison",
                                measure = 'weight', stacked = F, do.stat = T,
                                comparison = to_test[[1]], font.size = 10,
                                show.raw = T, return.data = T, sources.use = c(1:18),
                                targets.use = c(1:18), color.use = colors_test[[1]])


scaled_comparison <- info_flow_group1$gg.obj + info_flow_group2$gg.obj + 
  plot_annotation(title = paste0('Information Flow regarding Interaction Strength of ',
                                 names(to_test[1]), ' (not scaled)'),
                  theme = theme(plot.title = element_text(hjust = 0.5)) ) 

SCpubr::save_Plot(scaled_comparison,
                  figure_path = paste0(outdir, '9_information_plots/weight/not_scaled/'),
                  file_name = '9_information_plots_notscaled_low_high',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 17,
                  height = 11
) 


### 9.2.1 Low vs High. Weight/Interaction Strength scaled ----
info_flow_group1 <- rankNet_mod(cellchat_merged_fixed, mode = "comparison",
                                measure = 'weight', stacked = T, do.stat = T,
                                comparison = to_test[[1]], font.size = 10,
                                show.raw = F, return.data = T, sources.use = c(1:18),
                                targets.use = c(1:18), color.use = colors_test[[1]])

info_flow_group2 <- rankNet_mod(cellchat_merged_fixed, mode = "comparison",
                                measure = 'weight', stacked = F, do.stat = T,
                                comparison = to_test[[1]], font.size = 10,
                                show.raw = F, return.data = T, sources.use = c(1:18),
                                targets.use = c(1:18), color.use = colors_test[[1]])


scaled_comparison <- info_flow_group1$gg.obj + info_flow_group2$gg.obj + 
  plot_annotation(title = paste0('Information Flow regarding Interaction Strength of ',
                                 names(to_test[1])),
                  theme = theme(plot.title = element_text(hjust = 0.5)) ) 

SCpubr::save_Plot(scaled_comparison,
                  figure_path = paste0(outdir, '9_information_plots/weight/scaled/'),
                  file_name = '9_information_plots_notscaled_low_high',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 17,
                  height = 11
) 

# 10. Pathway Specific Plots (Fig. 5A)----
### 10.1 SPP1 Heatmap ----
pathways.show <- c('SPP1') # Which pathway to focus on
# Obtain maximum communication probability among all conditions
weight.max <- getMaxWeight(cellchat_breast, slot.name = c('netP'), attribute = pathways.show)

par(mfrow = c(1,3), xpd = T)
object.list <- list('control' = cellchat_breast$control,
                    'low' = cellchat_breast$low_hypoxia,
                    'high' = cellchat_breast$high_hypoxia)

ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap_mod_v2(object.list[[i]],
                               signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]),
                               max_value_dataset = weight.max
                               )
}

pathway_ht <- ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, 'cm'))
SCpubr::save_Plot(pathway_ht,
                  figure_path = paste0(outdir, '10_pathway_heatmaps/SPP1/'),
                  file_name = '10_SPP1_comm.prob_heatmap',
                  create_path = T,
                  dpi = 300,
                  output_format = 'all',
                  width = 17,
                  height = 7
)