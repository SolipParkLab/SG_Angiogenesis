# ------------------------------------------------------------------ #
# Communication Analysis Using CellChat after further Cell Filtering #  # Perpetrated by Manu Moradiellos
# ------------------------------------------------------------------ #

# Main code of the cell communication analysis requested.
# Using CellChat we'll study the main cell-cell interactions and intercommunication
# between the immune cells found in three conditions: control, low hypoxia and 
# high hypoxia (suppressed immune response) in tissues treated with antiVEGFA.

# |- Packages -|  -------------------------------------------------------------
library('Seurat')
library('CellChat')
library('ggplot2')
library('patchwork')
library('tidyverse')
library('svglite')
library('ComplexHeatmap')
source('/resources/netVisual_bubble_mod.R') # Slightly modified functions from CellChat
source('/resources/netVisual_heatmap_mod.R')
source('/resources/cellchat_rankNet_mod.R')
set.seed(1)

# |- Seurat Preprocess -| --------------
preprocess <- T # Whether to start from the seurat object or load an already
# created cellchat object to save some memory and time
if (preprocess) {
  # Read single-cell experiment prev. analyzed with Seurat by Luis, but with some further Cell Filtering to remove low quality samples and some high-readings
  seurat <- readRDS('./seurat_obj.hypoxia_breast_cancer.2500_27500umis.50hk.RDS')
  
  experiment_conditions <- SplitObject(seurat, split.by = 'hypoxia_label') # Split by control, low-hypoxia or high-hypoxia
  names(experiment_conditions) <- c('control', 'high_hypoxia', 'low_hypoxia')
  experiment_conditions <- experiment_conditions[c(1,3,2)] # Reorder to have the timeline right
  
  # Data initialization and CellChat object creation separately per experiment condition
  cellchat_breast <- setNames(lapply(names(experiment_conditions), function(condition) { 
    data.input <- GetAssayData(experiment_conditions[[condition]], assay = 'RNA', slot = 'data') # Normalized data matrix
    labels <- Idents(experiment_conditions[[condition]]) # Cell types defined from Seurat's resulting clusters
    cell_type_df <- data.frame(group = labels, row.names = names(labels)) # To identify the different cell clusters
    
    cellchat_cond <- createCellChat(object = data.input, meta = cell_type_df, group.by = 'group') 
    cellchat_cond <- addMeta(cellchat_cond, meta = cell_type_df, meta.name = 'group') # Add cell information
    cellchat_cond <- setIdent(cellchat_cond, ident.use = 'group') # Set "labels" as default cell identity
    cellchat_cond@DB <- CellChatDB.mouse # Import Ligand-Receptor Database on Mouse
    cellchat_cond <- subsetData(cellchat_cond) # Uses whole database to extract all signaling genes
    future::plan('multisession', workers = 4)
    cellchat_cond <- identifyOverExpressedGenes(cellchat_cond)        # Identify Genes and Ligand-Receptor interactions
    cellchat_cond <- identifyOverExpressedInteractions(cellchat_cond) # from the database and adds it to the object
    cellchat_cond <- computeCommunProb(cellchat_cond, population.size = T) # scRNAseq experiment was unsorted, uses Chromium 10x Standard Procedure
    cellchat_cond <- filterCommunication(cellchat_cond, min.cells = 10) # Filter out groups with low cell count in groups, maybe not representative 
    cellchat_cond <- computeCommunProbPathway(cellchat_cond)  # Infers probability from expression value of genes and reference database
    cellchat_cond <- aggregateNet(cellchat_cond) 
    print( data.frame('cell_type' = levels(cellchat_cond@idents) , 'group' = as.numeric(table(cellchat_cond@idents))) ) # To check "composition" differences
    return(cellchat_cond) 
  }), names(experiment_conditions) )
  saveRDS(cellchat_breast, file = paste0(working_directory, '1.cellchat_object.seurat.2000_27500umis.50hk_final.idents.cellclustersize_filtered.pop_size_T.RDS') )
} 

# |- Load CellChat Objects Prev. Obtained -|. ----------------

## From prior experience, I will only work with the object that 
## filters out those groups with less than 10 cell on them
# The groups affected are T-Cells Immature and mixed_group
                                                           
cellchat_breast <- readRDS(file = paste0(working_directory, '1.cellchat_object.seurat.2500_27500umis.50hk_final.idents.cellclustersize_filtered.pop_size_T.RDS') )
cellchat_merged <- mergeCellChat(cellchat_breast, add.names = c('control', 'low', 'high'), merge.data = T) # Object merging three conditions, the one above
                                                                                                           # is a list with separate objects in a list

# |- Dataframes of Cell Communications Found -| ---------
# Dataframe of ligand-receptor level cellular communication network
df.interaction_lr_net <- subsetCommunication(cellchat_merged)
dir.create(file.path(working_directory, '2.cellchat_cellcluster_filtered/tables'), recursive = T)
sapply( c('control', 'low', 'high') , function(test_cond){
  write_tsv(x = data.frame( df.interaction_lr_net[[test_cond]] ),
            file =  paste0(working_directory, '/2.cellchat_cellcluster_filtered/tables/cell_comm_ligand_receptor_table.', test_cond, '.csv') )
})

# Dataframe of cellular communication network at signaling pathway level
df.interaction_pathway_net <- subsetCommunication(cellchat_merged, slot.name = 'netP')
sapply( c('control', 'low', 'high') , function(test_cond){
  write_tsv(x = data.frame( df.interaction_pathway_net[[test_cond]] ),
            file =  paste0(working_directory, '/2.cellchat_cellcluster_filtered/tables/cell_comm_pathway_table.', test_cond, '.csv') )
})

# |- Barplots Total Num. and Strength Interactions -| ----
# Compare number of interactions and their strength per Test Condition 
barplot1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1:3),
                                title.name = 'Number of Interactions per Condition', size.text = 15, color.use = c('black', '#579CFF', '#FF5757') ) 
barplot2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1:3), measure = 'weight', color.use = c('black', '#579CFF', '#FF5757'),
                                title.name = 'Interaction Strength per Condition', size.text = 15) 
dir.create(file.path(working_directory, '2.cellchat_cellcluster_filtered/plots'), recursive = T)
svglite(filename = paste0(working_directory, '/2.cellchat_cellcluster_filtered/plots/general_comparison_barplots.svg') )
barplot1 + barplot2 + plot_layout(guides = 'collect') # Messes up x-axis order without guides collect
dev.off()

## Relist the cell clusters of interest, some aren't required anymore
## Number according to their level within the seurat object -checked by hand-
relist <- T
if(relist){
  dendritic_cells <- setNames(2, 'dendritic_cells')
  macrophages <- setNames(7, 'macrophages')
  monocytes <- setNames(9, 'monocytes')
  neutrophils <- setNames(10, 'neutrophils')
  # t-cells
  t_cells_cd8 <- setNames(15, 't_cells_cd8')
  t_cells_cd4 <- setNames(14, 't_cells_cd4')
  t_rex <- setNames(18, 't_cells_regs')
  # tumor cells
  epithelial_cells <- setNames(4, 'epithelial_cells')
  endothelial_cells <- setNames(3, 'endothelial_cells')
  stem_cells <- setNames(13, "stem_cells")
  ilc_cells <- setNames(6, "ilc_cells")
  b_cells <- setNames(1, "B cells")
  fibroblasts <- setNames(5, "fibroblasts")
  mixed_group <- setNames(8, "mixed_group")
}

cells_interest <- c(dendritic_cells, macrophages, monocytes, neutrophils,
                    t_cells_cd8, t_cells_cd4, t_rex, epithelial_cells,
                    endothelial_cells, stem_cells, ilc_cells, b_cells,
                    fibroblasts) # Number of ident for loops
cells_interest_subset <- c('Epithelial cells',  'Endothelial cells', 'Macrophages', 'DC',
                           'T cells CD8', 'Neutrophils', 'T cells CD4', 'Monocytes', 'ILC',
                           'Stem cells', 'NK cells', 'Fibroblasts', 'Tregs', 'B cells',
                           'Mixed Group')

#- 2.5 Information Flow comparison #----

# Red enriched in one condition, blue in the other and blank in none after statistical test 
# For now using show.raw = F so both graphs have the same results but 
control_color <- 'black'
low_color <- '#579CFF'
high_color <- '#FF5757'

to_test <- list('Low vs. High' = c(2,3), 'Control vs. Low' = c(1, 2), 'Control vs. High' = c(1,3))
colors_test <- list('Low vs. High' = c(low_color, high_color), 'Control vs. Low' = c(control_color, low_color), 'Control vs. High' = c(control_color, high_color))

# To solve problem with merged groups, as some are missing cell groups due to low number of cells
cellchat_merged_fixed <- liftCellChat(cellchat_merged, group.new = levels(cellchat_merged@idents$joint)) 

comparison <- 1 # For the different comparisons to be made
## With measure = 'weight' or 'count' we can compare either the total interaction strength or the total number
## Remember that for 'weight' you need to obtain scaled and not-scaled results using show.raw = F or T
info_flow <- rankNet_mod(cellchat_merged_fixed, mode = "comparison", measure = 'weight', stacked = T, do.stat = T,
                         comparison = to_test[[comparison]], font.size = 10, show.raw = T,
                         return.data = T, sources.use = cells_interest, targets.use = cells_interest, color.use = colors_test[[comparison]])

info_flow_2 <- rankNet_mod(cellchat_merged_fixed, mode = "comparison", measure = 'weight', stacked = F, do.stat = T,
                           comparison = to_test[[comparison]], font.size = 10, show.raw = T,
                           return.data = T, sources.use = cells_interest, targets.use = cells_interest, color.use = colors_test[[comparison]], ylim = 0.04)
info_flow$gg.obj + info_flow_2$gg.obj + plot_annotation(title = paste0('Information Flow regarding Interaction Strength of ', names(to_test[comparison]), ' (not scaled)'),
                                                        theme = theme(plot.title = element_text(hjust = 0.5)) ) 

info_flow$gg.obj + info_flow_2$gg.obj + plot_annotation(title = paste0('Information Flow regarding Number of Interactions of ', names(to_test[comparison])),
                                                        theme = theme(plot.title = element_text(hjust = 0.5)) )


#- 2. Pathway-level Plots #----
## Current pathways of interest: SPP1, MIF, CCL, TNF, TGF-Beta
pathways_interest <- c('CCL', 'MIF', 'TGFb', 'SPP1', 'TNF')

#- 2.3 Heatmaps #----
## Probably the best way to showcase general trends within conditions 

pathways.show <- c('SPP1')
weight.max <- getMaxWeight(cellchat_breast, slot.name = c('netP'), attribute = pathways.show)

## Both of these functions were found elsewhere, needed for correct visualization
zeros_after_period <- function(x) {
  if (isTRUE(all.equal(round(x),x))) return (0) # y would be -Inf for integer values
  y <- log10(abs(x)-floor(abs(x)))   
  ifelse(isTRUE(all.equal(round(y),y)), -y-1, -ceiling(y))}

check_min_nonzero <- function(pathway = pathways.show){
  c <- cellchat_breast$control@netP$prob[, , pathway]
  c[c == 0] <- NA
  l <- cellchat_breast$low_hypoxia@netP$prob[, , pathway]
  l[l == 0] <- NA
  h <- cellchat_breast$high_hypoxia@netP$prob[, , pathway]
  h[h == 0] <- NA
  return( min(c(c,l,h), na.rm = T) )
}

cutre_jitmap <- function(cellchat.object.list = cellchat_breast, cond, pathway = pathways.show,
                         max_prob = weight.max, min_prob = weight.min){
  whatevs <- data.frame( cellchat.object.list[[cond]]@netP$prob[, , pathway] )
  #whatevs <- whatevs[!rownames(whatevs) %in% c('NKT', 'Tgd', 'T cells immature-early', 'mixed_group'), ] # To be removed from analysis
  df <- whatevs %>% rownames_to_column() %>% 
    gather(colname, value, -rowname)  %>% 
    filter( !(rowname %in% c('NKT', 'Tgd', 'T cells immature-early', 'mixed_group'))) %>% 
    filter( !(colname %in% c('NKT', 'Tgd', 'T.cells.immature.early', 'mixed_group'))) # To be ignored from analysis
  df[df == 0] <- NA
  ht <- ggplot(df, aes(x = colname, y = rowname, fill = value)) +
    geom_tile(color = 'white') + coord_equal() +
    #scale_fill_gradient2(na.value = 'white', low = 'white' , mid = '#f7e7e9', high = '#b2182b', limits = c(0, round( (max_prob + as.numeric(paste0('0.',paste(c(rep('0', times = zeros_after_period(max_prob) + 1 )), collapse = ''), '1'))), zeros_after_period(max_prob) + 1) )) +
    scale_fill_gradient2(na.value = 'white', low = 'white' , mid = '#f7e7e9', high = '#b2182b', limits = c(0, max_prob + as.numeric(paste0('0.',paste(c(rep('0', times = zeros_after_period(max_prob) + 1 )), collapse = ''), '1')) )) +
    labs(title = paste0(pathway, ' signaling ' , cond), x = 'Target (Receiver)', y = 'Source (Sender)') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5)) 
  return(ht)
}

# Heatmap of various pathways of interest
pathways.show <- c('SPP1')
if ( pathways.show %in% c('SPP1', 'TNF', 'FGF', 'PDGF', 'ANGPTL') ) {
  weight.max <- max(cellchat_breast$high_hypoxia@netP$prob[, , pathways.show])
  weight.min <- min(cellchat_breast$high_hypoxia@netP$prob[, , pathways.show])
} else {
  weight.max <- getMaxWeight(cellchat_breast, slot.name = c('netP'), attribute = pathways.show)
  weight.min <- check_min_nonzero(pathways.show)
}

(   (cutre_jitmap(cond = 'control') + theme(axis.title.x = element_blank() ) ) +
    (cutre_jitmap(cond = 'low_hypoxia') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() ) ) + 
    (cutre_jitmap(cond = 'high_hypoxia') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank() ) ) ) +
  plot_layout(guides = 'collect')


cutre_jitmap(cond = 'low_hypoxia')  + 
  (cutre_jitmap(cond = 'high_hypoxia') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank() ) )  +
  plot_layout(guides = 'collect')


#- 2.4 Violin Plots Signaling Genes Expression #----
cellchat_merged@meta$datasets <- factor(cellchat_merged@meta$datasets, levels = c('control', 'low', 'high'))
plotGeneExpression(cellchat_merged, signaling = pathways.show, y.max = 6,
                   group.by = 'group', split.by = 'datasets',
                   color.use = c(control_color, low_color, high_color))
