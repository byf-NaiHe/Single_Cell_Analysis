# adjust based on sscClust::ssc.plot.headmap
suppressPackageStartupMessages({
library("ComplexHeatmap")
library("circlize")
library("gridBase")
library("grid")
library("RColorBrewer")
library("zoo")
})

heatmap_sm.2 = function (obj, assay.name = "exprs", out.prefix = NULL, ncell.downsample = NULL, 
  ave.by = NULL, columns = NULL, columns.order = NULL, 
  rows = NULL,
  gene.desc = NULL, 
  colSet = list(), pdf.width = 16, pdf.height = 15, do.scale = TRUE, 
  z.lo = -2.5, z.hi = 2.5, z.step = 1, exp.title = "Exp", 
  do.clustering.row = T, do.clustering.col = T, dend.col = FALSE, 
  dend.row = FALSE, clustering.distance = "spearman", clustering.method = "complete", 
  k.row = 1, k.col = 1, palette.name = NULL, annotation_legend_param = list(), 
  ann.bar.height = 1.5, mytitle = "", ...) 
{
  if (!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && 
    ("geneID" %in% colnames(gene.desc))) {
    obj <- obj[gene.desc$geneID, ]
  }
  if (!is.null(ncell.downsample) && ncell.downsample < ncol(obj)) {
    obj <- obj[, sample(seq_len(ncol(obj)), ncell.downsample)]
  }
  n <- nrow(obj)
  m <- ncol(obj)
  if (n < 3) {
    loginfo(sprintf("Too few genes: n=%s", n))
    return(NULL)
  }
  if (m < 3) {
    loginfo(sprintf("Too few samples: m=%s", m))
    return(NULL)
  }
  if (is.null(ave.by)) {
    obj <- ssc.assay.hclust(obj, assay.name = assay.name, 
      order.col = if (is.logical(dend.col) && FALSE == 
        dend.col) 
        do.clustering.col
      else FALSE, order.row = if (is.logical(dend.row) && 
        FALSE == dend.row) 
        do.clustering.row
      else FALSE, clustering.distance = "spearman", clustering.method = "complete", 
      k.row = 1, k.col = 1)
  }
  else {
    obj <- ssc.average.cell(obj, assay.name = assay.name, 
      column = ave.by, ret.type = "sce")
    columns <- intersect(ave.by, columns)
    columns.order <- intersect(ave.by, columns.order)
  }
  #
  ha.col <- NULL
  annDF <- data.frame()
  if (!is.null(columns)) {
    if (!is.null(columns.order)) {
      obj <- ssc.order(obj, columns.order = columns.order)
    }
    annDF <- as.data.frame(colData(obj)[columns])
    if (length(colSet) == 0) {
      for (i in seq_along(columns)) {
        x <- columns[i]
        if (class(colData(obj)[, x]) == "numeric") {
          if (all(colData(obj)[, x] <= 1) && all(colData(obj)[, 
            x] >= 0)) {
            Y.level <- c(0, 1)
          }
          else {
            Y.level <- pretty(colData(obj)[, x], n = 8)
          }
          colSet[[x]] <- colorRamp2(seq(Y.level[1], 
            Y.level[length(Y.level)], length = 7), rev(brewer.pal(n = 7, 
            name = "RdYlBu")), space = "LAB")
          annotation_legend_param[[x]] <- list(color_bar = "continuous", 
            legend_direction = "horizontal", legend_width = unit(4, 
              "cm"), legend_height = unit(2, "cm"))
        }
        else {
          group.value <- sort(unique(colData(obj)[, 
            x]))
          colSet[[x]] <- structure(auto.colSet(length(group.value), 
            name = "Accent"), names = group.value)
        }
      }
    }
    g.show.legend <- T
    ha.col <- ComplexHeatmap::HeatmapAnnotation(df = annDF, 
      col = colSet, show_legend = g.show.legend, simple_anno_size = unit(ann.bar.height, 
        "cm"), annotation_legend_param = annotation_legend_param)
  }
  # patch: add row annotation on right side
  ha.row <- NULL
  annDF <- data.frame()
  if (!is.null(rows)) {
    annDF <- as.data.frame(rowData(obj)[rows])
    if (length(colSet) == 0) {
      for (i in seq_along(rows)) {
        x <- rows[i]
        if (class(rowData(obj)[, x]) == "numeric") {
          if (all(rowData(obj)[, x] <= 1) && all(rowData(obj)[, 
            x] >= 0)) {
            Y.level <- c(0, 1)
          }
          else {
            Y.level <- pretty(rowData(obj)[, x], n = 8)
          }
          colSet[[x]] <- colorRamp2(seq(Y.level[1], 
            Y.level[length(Y.level)], length = 7), rev(brewer.pal(n = 7, 
            name = "RdYlBu")), space = "LAB")
          annotation_legend_param[[x]] <- list(color_bar = "continuous", 
            legend_direction = "horizontal", legend_width = unit(4, 
              "cm"), legend_height = unit(2, "cm"))
        }
        else {
          group.value <- sort(unique(rowData(obj)[, 
            x]))
          colSet[[x]] <- structure(auto.colSet(length(group.value), 
            name = "Accent"), names = group.value)
        }
      }
    }
    g.show.legend <- T
    ha.row <- ComplexHeatmap::rowAnnotation(df = annDF, 
      col = colSet, show_legend = g.show.legend, simple_anno_size = unit(ann.bar.height, 
        "cm"), annotation_legend_param = annotation_legend_param)
  }
  #
  obj <- ssc.order(obj, columns.order = NULL, gene.desc = gene.desc)
  dat.plot <- as.matrix(assay(obj, assay.name))
  rownames(dat.plot) <- unname(rowData(obj)$display.name)
  if (do.scale) {
    rowM <- rowMeans(dat.plot, na.rm = T)
    rowSD <- apply(dat.plot, 1, sd, na.rm = T)
    dat.plot <- sweep(dat.plot, 1, rowM)
    dat.plot <- sweep(dat.plot, 1, rowSD, "/")
    if (!is.null(z.lo)) {
      dat.plot[dat.plot < z.lo] <- z.lo
    }
    if (!is.null(z.hi)) {
      dat.plot[dat.plot > z.hi] <- z.hi
    }
  }else {
    tmp.var <- pretty((dat.plot), n = 8)
    if (is.null(z.lo)) {
      z.lo <- tmp.var[1]
    }
    if (is.null(z.hi)) {
      z.hi <- tmp.var[length(tmp.var)]
    }
    if (is.null(z.step)) {
      z.step <- tmp.var[2] - tmp.var[1]
    }
  }
  if (!is.null(out.prefix)) {
    pdf(sprintf("%s.pdf", out.prefix), width = pdf.width, 
      height = pdf.height)
  }
  par(mar = c(4, 12, 4, 4))
  plot.new()
  title(main = mytitle, cex.main = 2)
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  if (is.null(palette.name)) {
    exp.palette <- rev(brewer.pal(n = 7, name = ifelse(do.scale, 
      "RdBu", "RdYlBu")))
  } else {
      if(palette.name=="YlOrRd"){
	  exp.palette <- sscVis:::getColorPaletteFromNameContinuous("YlOrRd")
      }else{
	exp.palette <- rev(brewer.pal(n = 7, name = palette.name))
      }
  }
  # smooth value
  if (ncol(dat.plot)>50){
    dat.plot = t(apply(dat.plot, 1, function(x){rollmean(x, 50, fill="extend")}))
  }
  #
  ht <- ComplexHeatmap::Heatmap(dat.plot, name = exp.title, 
    col = colorRamp2(seq(z.lo, z.hi, length = 100), colorRampPalette(exp.palette)(100), 
      space = "LAB"), column_names_gp = grid::gpar(fontsize = 12 * 
      28/max(m, 32)), row_names_gp = grid::gpar(fontsize = 18 * 
      28/max(n, 32)), show_heatmap_legend = T, row_names_max_width = unit(10, 
      "cm"), cluster_columns = dend.col, cluster_rows = dend.row, 
    row_dend_reorder = FALSE, column_dend_reorder = FALSE, 
    heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
      grid_height = unit(0.8, "cm"), at = seq(z.lo, z.hi, 
        z.step), title_gp = grid::gpar(fontsize = 14, 
        fontface = "bold"), label_gp = grid::gpar(fontsize = 12), 
      color_bar = "continuous"), top_annotation = ha.col, 
     right_annotation =ha.row,
    ...)
  ComplexHeatmap::draw(ht, newpage = FALSE)
  if (!is.null(ha.col)) {
    for (i in seq_along(names(ha.col@anno_list))) {
      ComplexHeatmap::decorate_annotation(names(ha.col@anno_list)[i], 
        {
          grid.text(names(ha.col@anno_list)[i], unit(-4, 
            "mm"), gp = grid::gpar(fontsize = 14), just = "right")
        })
    }
  }
  if (!is.null(out.prefix)) {
    dev.off()
  }
}



