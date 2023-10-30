# ---------------------------------------------------------------- #
# Differential Expression Analysis between Experimental Conditions #  # Perpetrated by Manu Moradiellos
# ---------------------------------------------------------------- #

# DEA performed by request of José Luis 

# |- Packages -|  -------------------------------------------------------------
library('Seurat')
library('DESeq2')
library('openxlsx')
library('stringr')

set.seed(666)

# Read single-cell experiment analyzed with Seurat by Luis, but with some further
# cell filtering done to remove low quality samples and some high-readings
seurat <- readRDS('./input/seurat_obj.hypoxia_breast_cancer.2500_27500umis.50hk.RDS')

# Create specific levels for each cell group per condition
seurat$Celltype.Hypoxia <- paste(Idents(seurat), seurat$hypoxia_label, sep = '_')
seurat$final_celltype <- Idents(seurat) # Copy again just in case
Idents(seurat) <- seurat$Celltype.Hypoxia

#Set styles for xlsx files
redStyle <- createStyle(fontColour = "#B60A1C", bgFill = "#FFF06A", textDecoration = c("BOLD"))
greenStyle <- createStyle(fontColour = "#309143", bgFill = "#FFF06A", textDecoration = c("BOLD"))
blackStyle_ita <- createStyle(fontColour = "black", bgFill = "#FFF06A", textDecoration = c("italic"))
blackStyle_bold <- createStyle(fontColour = "black", bgFill = "#FFF06A", textDecoration = c("bold"))


## Should create a function but this will be faster right now 
## Some failed and need to be rerun
#already_made <- c('T cells immature-early', 'Mixed group') # Neutrophils not available to use in any Low vs. comparison
#already_made <- c('B cells', 'DC', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'ILC', 'Macrophages', 'Mixed group', 'Monocytes', 
already_made <- c('NK cells', 'NKT', 'Stem cells', 'T cells CD4', 'T cells immature-early', 'Tgd', 'Tregs', 'Neutrophils') #'T cells CD8'

## Some parts are based on Bollito's code as I liked the way they did it: https://github.com/cnio-bu/bollito/blob/main/scripts/step5_degs.R
# I don't test cell groups that have less than three cells for that condition

# Low vs. High #----
if (TRUE){
  for (cell_cluster in levels(seurat$final_celltype)[ !(levels(seurat$final_celltype) %in% already_made )] ){ # Mixed group doesn't have enough cells 
    ident_low <- paste0(cell_cluster, '_low')
    ident_high <- paste0(cell_cluster, '_high')
    print(sprintf("Current DEA %s vs %s", ident_low, ident_high))
  
    # To show positive Log2FC as over-expressed genes in High this is the right order (use as "reference")
    clusterX.markers <- FindMarkers(seurat, ident.1 = ident_high, ident.2 = ident_low, min.pct = 0, logfc.threshold = 0, test.use = 'DESeq2') #min expressed
  
    wb <- createWorkbook() # Create an Excel file for JL to check 
    addWorksheet(wb, "DE analysis")
    writeData(wb, "DE analysis", clusterX.markers, rowNames = TRUE)
    conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                          rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2<0, $F2<0.05)",
                          style = greenStyle)
    conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                          rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2>0, $F2<0.05)",
                          style = redStyle)
    legend <- createComment(comment = c("","Red means a positive LogFold\n",
                                        "(Overexpressed in High Hypoxia)\n\n",
                                        "Green means a negative LogFold\n",
                                        "(Overexpressed in Low Hypoxia)"),
                            style = c(redStyle, redStyle, blackStyle_ita, greenStyle, blackStyle_ita))
    extra_legend <- createComment(comment = c("", "p_val: ", "unadjusted p-value\n",
                                              "avg_logFC: ", "log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group (High)\n",
                                              "pct.1: ", "The percentage of cells where the feature is detected in the first group (High Hypoxia)\n",
                                              "pct.2: ", "The percentage of cells where the feature is detected in the second group (Low Hypoxia)\n",
                                              "p_val_adj: ", "Adjusted p-value, based on bonferroni correction using all features in the dataset"),
                                  style = c(blackStyle_bold, blackStyle_bold, blackStyle_ita,
                                            blackStyle_bold, blackStyle_ita,
                                            blackStyle_bold, blackStyle_ita,
                                            blackStyle_bold, blackStyle_ita,
                                            blackStyle_bold, blackStyle_ita))
    writeComment(wb, "DE analysis", col = 8, row = 2, comment = legend)
    writeComment(wb, "DE analysis", col = 11, row = 2, comment = extra_legend)
    saveWorkbook(wb, paste0("./output/low_vs_high.dea/", str_replace_all(cell_cluster, ' ', '_'),".DE.Low_vs_High.xlsx"), overwrite = TRUE)    
    }
}

## Control vs. Low #----
for (cell_cluster in levels(seurat$final_celltype)[ !(levels(seurat$final_celltype) %in% already_made )] ){ # Mixed group doesn't have enough cells 
  ident_low <- paste0(cell_cluster, '_low')
  ident_control <- paste0(cell_cluster, '_control')
  print(sprintf("Current DEA %s vs %s", ident_low, ident_control))
  
  # To show positive Log2FC as over-expressed genes in High this is the right order (use as "reference")
  clusterX.markers <- FindMarkers(seurat, ident.1 = ident_control, ident.2 = ident_low, min.pct = 0, logfc.threshold = 0, test.use = 'DESeq2') #min expressed
  
  wb <- createWorkbook() # Create an Excel file for JL to check 
  addWorksheet(wb, "DE analysis")
  writeData(wb, "DE analysis", clusterX.markers, rowNames = TRUE)
  conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                        rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2<0, $F2<0.05)",
                        style = greenStyle)
  conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                        rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2>0, $F2<0.05)",
                        style = redStyle)
  legend <- createComment(comment = c("","Red means a positive LogFold\n",
                                      "(Overexpressed in High Hypoxia)\n\n",
                                      "Green means a negative LogFold\n",
                                      "(Overexpressed in Low Hypoxia)"),
                          style = c(redStyle, redStyle, blackStyle_ita, greenStyle, blackStyle_ita))
  extra_legend <- createComment(comment = c("", "p_val: ", "unadjusted p-value\n",
                                            "avg_logFC: ", "log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group (High)\n",
                                            "pct.1: ", "The percentage of cells where the feature is detected in the first group (High Hypoxia)\n",
                                            "pct.2: ", "The percentage of cells where the feature is detected in the second group (Low Hypoxia)\n",
                                            "p_val_adj: ", "Adjusted p-value, based on bonferroni correction using all features in the dataset"),
                                style = c(blackStyle_bold, blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita))
  writeComment(wb, "DE analysis", col = 8, row = 2, comment = legend)
  writeComment(wb, "DE analysis", col = 11, row = 2, comment = extra_legend)
  saveWorkbook(wb, paste0('./output/control_vs_low.dea/', str_replace_all(cell_cluster, ' ', '_'),".DE.Control_vs_Low.xlsx"), overwrite = TRUE)    
}

## Control vs. High #----
for (cell_cluster in levels(seurat$final_celltype)[ !(levels(seurat$final_celltype) %in% c('Mixed group', 'T cells immature-early'))] ){ # Mixed group doesn't have enough cells 
  ident_high <- paste0(cell_cluster, '_high')
  ident_control <- paste0(cell_cluster, '_control')
  print(sprintf("Current DEA %s vs %s", ident_high, ident_control))
  
  # To show positive Log2FC as over-expressed genes in High this is the right order (use as "reference")
  clusterX.markers <- FindMarkers(seurat, ident.1 = ident_control, ident.2 = ident_high, min.pct = 0, logfc.threshold = 0, test.use = 'DESeq2') #min expressed
  
  wb <- createWorkbook() # Create an Excel file for JL to check 
  addWorksheet(wb, "DE analysis")
  writeData(wb, "DE analysis", clusterX.markers, rowNames = TRUE)
  conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                        rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2<0, $F2<0.05)",
                        style = greenStyle)
  conditionalFormatting(wb, "DE analysis", cols = 1:(ncol(clusterX.markers)+1),
                        rows = 2:(nrow(clusterX.markers) + 1), rule = "AND($C2>0, $F2<0.05)",
                        style = redStyle)
  legend <- createComment(comment = c("","Red means a positive LogFold\n",
                                      "(Overexpressed in High Hypoxia)\n\n",
                                      "Green means a negative LogFold\n",
                                      "(Overexpressed in Low Hypoxia)"),
                          style = c(redStyle, redStyle, blackStyle_ita, greenStyle, blackStyle_ita))
  extra_legend <- createComment(comment = c("", "p_val: ", "unadjusted p-value\n",
                                            "avg_logFC: ", "log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group (High)\n",
                                            "pct.1: ", "The percentage of cells where the feature is detected in the first group (High Hypoxia)\n",
                                            "pct.2: ", "The percentage of cells where the feature is detected in the second group (Low Hypoxia)\n",
                                            "p_val_adj: ", "Adjusted p-value, based on bonferroni correction using all features in the dataset"),
                                style = c(blackStyle_bold, blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita,
                                          blackStyle_bold, blackStyle_ita))
  writeComment(wb, "DE analysis", col = 8, row = 2, comment = legend)
  writeComment(wb, "DE analysis", col = 11, row = 2, comment = extra_legend)
  saveWorkbook(wb, paste0("./output/control_vs_high.dea/", str_replace_all(cell_cluster, ' ', '_'),".DE.Control_vs_High.xlsx"), overwrite = TRUE)    
}