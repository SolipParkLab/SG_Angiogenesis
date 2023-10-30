
#' Visualization of network using heatmap
#'
#' This heatmap can be used to show differential number of interactions or interaction strength in the cell-cell communication network between two datasets;
#' the number of interactions or interaction strength in a single dataset
#' the inferred cell-cell communication network in single dataset, defined by `signaling`
#'
#' When show differential number of interactions or interaction strength in the cell-cell communication network between two datasets, the width of edges represent the relative number of interactions or interaction strength.
#' Red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
#'
#' The top colored bar plot represents the sum of column of values displayed in the heatmap. The right colored bar plot represents the sum of row of values.
#'
#'
#' @param object A merged CellChat object or a single CellChat object
#' @param comparison a numerical vector giving the datasets for comparison in object.list; e.g., comparison = c(1,2)
#' @param measure "count" or "weight". "count": comparing the number of interactions; "weight": comparing the total interaction weights (strength)
#' @param signaling a character vector giving the name of signaling networks in a single CellChat object
#' @param slot.name the slot name of object. Set is to be "netP" if input signaling is a pathway name; Set is to be "net" if input signaling is a ligand-receptor pair
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap A vector of two colors corresponding to max/min values, or a color name in brewer.pal only when the data in the heatmap do not contain negative values
#' @param title.name the name of the title
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @param font.size.title font size of the title
#' @param cluster.rows whether cluster rows
#' @param cluster.cols whether cluster columns
#' @param sources.use a vector giving the index or the name of source cell groups
#' @param targets.use a vector giving the index or the name of target cell groups.
#' @param remove.isolate whether remove the isolate nodes in the communication network
#' @param row.show,col.show a vector giving the index or the name of row or columns to show in the heatmap
#' @importFrom methods slot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_barplot rowAnnotation
#' @return  an object of ComplexHeatmap
#' @export

import::from(ComplexHeatmap, 'Heatmap', 'HeatmapAnnotation', 'anno_barplot', 'rowAnnotation')
import::from(grDevices, 'colorRampPalette')
import::from(RColorBrewer, 'brewer.pal')
import::from(methods, 'slot')
import::from(grid, 'gpar')
source('/local/mmoradiellos/work_stuff/post-TFM/breast.cancer_single_cell_collaboration/14_Thesis_figures/ComplexHeatmap_heatmap_mod.R')

netVisual_heatmap_mod_v2 <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL,
                              legend.title = 'Relative values', max_value_dataset = 1, min_value_dataset = 0){
  # obj1 <- object.list[[comparison[1]]]
  # obj2 <- object.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = legend.title
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  
  net <- net.diff
  
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  
  mat <- net
  og_max <- max(mat)
  og_min <- min(mat)
  
  # EDIT ----
  if (! ('Mixed group' %in% rownames(mat)) ){
    mat <- rbind(mat, 'Mixed group' = rep(0, ncol(mat)))
    mat <- cbind(mat, 'Mixed group' = rep(0, nrow(mat)))
    mat <- mat[rownames(mat)[order(rownames(mat))] , colnames(mat)[order(colnames(mat))] ,drop=FALSE]
  }
  
  # Trying to color the scale correctly, but the annotation needs to remain as before
  mat <- rbind(mat, 'dummy_value' = rep(max_value_dataset, ncol(mat)))
  mat <- cbind(mat, 'dummy_value' = rep(max_value_dataset, nrow(mat)))
  
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  
  
  if (min_value_dataset < 0) {
    color.heatmap.use = colorRamp3(c(min_value_dataset, 0, max_value_dataset), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    colorbar.break <- c(min_value_dataset, max_value_dataset)
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min_value_dataset, max_value_dataset), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min_value_dataset, max_value_dataset), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)), )(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- ComplexHeatmap::HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  # EDIT ANNOTATIONS ----
  ## Calculate true total values and add dummy value to not break scale
  incoming_signal <- rowSums(abs(mat[-c(ncol(mat)), -c(ncol(mat))]))
  incoming_signal['dummy_value'] <- 0
  outgoing_signal <- colSums(abs(mat[-c(ncol(mat)), -c(ncol(mat))]))
  outgoing_signal['dummy_value'] <- 0
  
  ha1 = ComplexHeatmap::rowAnnotation(Strength = ComplexHeatmap::anno_barplot(incoming_signal,
                                                                              border = FALSE,
                                                                              gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = ComplexHeatmap::HeatmapAnnotation(Strength = ComplexHeatmap::anno_barplot(outgoing_signal,
                                                                                  border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, 
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}
