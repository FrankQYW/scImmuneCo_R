#' GSVA for one cell type
#'
#' @param seurat_object anotated seurat object
#' @param cell_column Cell type column
#' @param sample_column Sample column
#' @param condition_column Condition column
#' @param cell Cell type to subset
#' @importFrom Seurat AverageExpression
#' @return A Seurat object with GSVA results
#' @export
gsva_cell_type <- function(seurat_object, cell_column = 'predicted.celltype.l1',
                           sample_column = 'GEO_ident',
                           condition_column = 'Groups', cell, gmt = module_info
){
  object <- seurat_object
  string.to.evaluate <- paste0(cell_column,'==', '\'', cell, '\'')




  print(string.to.evaluate)

  subset_object <- subset(object, !!str2lang(string.to.evaluate))

  subset_pseudobulk <- AverageExpression(subset_object, return.seurat = TRUE,
                                         group.by = c(sample_column))

  exp <- subset_pseudobulk[['RNA']]$counts
  ### GSVA ###
  if (cell == 'other T'){
    cell_clean = 'other_T'
  }else{
    cell_clean = gsub( ' ', '', cell)
  }


  print(gmt)
  gmt <- dplyr::filter(gmt, gmt$cell_type == cell_clean)

  gmt <- dplyr::select(gmt, c(5,1))
  gmt <- split(gmt$gene, gmt$name, drop = TRUE)
  gsva_result <- GSVA::gsva(as.matrix(exp), gmt, kcdf = 'Gaussian',
                      verbose = F)

  return(gsva_result)
}




#' compare between two conditions
#'
#' @param seurat_object anotated seurat object
#' @param gsva GSVA results
#' @param sample_column Sample column
#' @param condition_column Condition column
#' @param condition1 First condition to compare
#' @param condition2 Second condition to compare
#' @return Differentially expressed modules
#' @export
compare_condition <- function(seurat_object,gsva, sample_column, condition_column,
                              condition1, condition2){

  meta <- seurat_object@meta.data
  meta <- meta[!duplicated(meta[,sample_column]), ]

  module_term <- rownames(gsva)

  cell_t <- tibble::rownames_to_column(as.data.frame(t(gsva)))

  meta_cell <-  merge(x = meta,y= cell_t,
                      by.x = sample_column, by.y = "rowname",all.x=TRUE)

  #############################################################################
  gsva_matrix <- gsva
  meta_cell <- meta_cell[meta_cell[,sample_column] %in% colnames(gsva_matrix),]

  sample <- meta_cell[,sample_column]
  #sample <- gsub( '_', '-', sample)


  print(sample)
  print(colnames(gsva_matrix))
  gsva_matrix <- gsva_matrix[, sample]

  Group <- factor(meta_cell[,condition_column])


  design_limma <- stats::model.matrix(~Group)
  colnames(design_limma) <- levels(Group)
  fit <- limma::lmFit(gsva_matrix, design_limma)
  x <- c(paste0(condition1, '-', condition2))
  print(x)
  contMatrix <- limma::makeContrasts(contrasts = x, levels = design_limma)

  fit <- limma::contrasts.fit(fit, contMatrix)
  fit <- limma::eBayes(fit)
  results <- limma::topTable(fit, coef = x[1], number = nrow(fit))

  results <- tibble::rownames_to_column(results)
  merge_info <- module_info[, c('name', 'merge_module')]
  merge_info <- merge_info[!duplicated(merge_info$name), ]
  results <- dplyr::left_join(results, merge_info, by = c('rowname' = 'name'))

  return(results)
}





#' compare between two conditions
#'
#' @param seurat_object anotated seurat object
#' @param sample_column Sample column
#' @param condition_column Condition column
#' @param condition1 First condition to compare
#' @param condition2 Second condition to compare
#' @return Differentially expressed modules
#' @export
do_dep <- function(seurat_object, sample_column, gmt = module_info, condition_column = 'Groups', condition1 = 'SLE',
                   condition2 = 'HC'){
  seurat_object@meta.data[,sample_column] <- gsub('_','-', seurat_object@meta.data[,sample_column])
  sample_data <- unique(seurat_object@meta.data[, c(sample_column, condition_column)])
  final_results <- data.frame()
  for (i in c('other T',"CD8 T", "B", "Mono", "CD4 T", "DC", "NK")){
    gsva_r <- gsva_cell_type(seurat_object, cell = i , sample_column = sample_column
                             )


    limma_result <- compare_condition(seurat_object, gsva_r, sample_column = sample_column,
                                      condition_column = condition_column, condition1 = condition1,
                                      condition2 = condition2)
    limma_result$cell_type <- i
    final_results <- rbind(final_results, limma_result)
    #   write.csv(limma_result, paste0('./', condition1, '/dep_', i, '.csv'))
    if (i == 'other T'){
      cell_clean = 'other_T'
    }else{
      cell_clean = gsub( ' ', '', i)
    }


  }
  final_results$disease <- condition1

  return(final_results)


}

#' compare between two conditions
#'
#' @param gene_list INPUT GENE LIST
#' @return gene_res Differentially expressed modules
#' @export
gene_enrichment_scimmuneco <- function(gene_list){
  gmt <- module_info
  gmt <- dplyr::select(gmt, c(5,1))
  colnames(gmt) <- c('term', 'gene')  # Rename columns to "term" and "gene"
  gmt <- split(gmt$gene, gmt$term, drop = TRUE)
  gmt_fixed <- utils::stack(gmt)  # Converts list to 2-column dataframe (ind, values)
  colnames(gmt_fixed) <- c("gene", "term")  # Rename columns to "term" and "gene"
  gmt_fixed <- gmt_fixed[, c("term", "gene")]  # Reorder columns

  gene_res <- clusterProfiler::enricher(gene = gene_list , TERM2GENE = gmt_fixed, pAdjustMethod = 'none', pvalueCutoff = 0.05)

  return(gene_res)

}



#' compare between two conditions
#'
#' @param results DEM results
#' @importFrom dplyr filter mutate summarise group_by n %>%
#' @return Differentially expressed modules
#' @export

draw_volcano <- function(comparison, p_cutoff = 0.05, logFC_cutoff = 0.25) {
  print(comparison)
  select_results <- comparison
  select_results <- filter(select_results, select_results$P.Value < p_cutoff)
  select_results$p_val_adj <- select_results$P.Value
  select_results$dir <- ifelse(select_results$logFC > 0, 'up', 'down')
  count_df <- select_results %>%
    group_by(cell_type, dir) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(count = ifelse(dir == "down", -count, count))  # Make "down" negative


  print(count_df)
  beautiful_colors <- c(
    'Mono' = "#E63946",  # Vibrant red
    'DC' = "#F4A261",  # Warm orange
    'CD8 T' = "#2A9D8F",  # Teal green
    'CD4 T' = "#264653",  # Dark blue-green
    'NK' = "#E9C46A",  # Golden yellow
    'other T' = "#6A4C93",  # Rich purple
    'B' = "#3A86FF"   # Bright blue
  )


  colnames(select_results)[2] <- 'avg_log2FC'
  colnames(select_results)[5] <- 'p_val'
  colnames(select_results)[9] <- 'cluster'
  colnames(select_results)[1] <- 'gene'



  count_df <- count_df %>%
    mutate(
      x_pos = ifelse(dir == "up",
                     max(select_results$avg_log2FC) * 0.9,
                     min(select_results$avg_log2FC) * 0.9),
      y_pos = max(-log10(select_results$p_val)) * 0.9
    )


  scRNAtoolVis::jjVolcano(diffData = select_results, topGeneN = 0, tile.col = beautiful_colors,
            back.col= 'white', log2FC.cutoff = logFC_cutoff) + ggplot2::xlab(''
            ) + ggplot2::ylab('') + ggplot2::theme(legend.position = "none")


}


#' Visualize Differential Expression by Module
#'
#' Creates a plot showing differentially expressed genes grouped by modules,
#' with modules ordered numerically (e.g., CD4T_1, CD4T_2, etc.).
#' Significant genes (P.Value < threshold) are colored by logFC while
#' non-significant genes appear in gray.
#'
#' @param deg_results A dataframe containing differential expression results.
#' Must include columns: `merge_module`, `logFC`, and `P.Value`.
#' @param pval_threshold P-value cutoff for significance (default = 0.05).
#' @param point_size Size of points in the plot (default = 3).
#' @param jitter_width Width of jitter for point positioning (default = 0.15).
#'
#' @return A ggplot object showing module-separated DEG visualization.
#'
#' @importFrom dplyr mutate group_by ungroup
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_color_gradient2
#' @importFrom ggplot2 theme_minimal theme labs position_jitter
#' @importFrom stats reorder
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data(example_deg_results)  # Your DEG dataframe
#' plot_module_deg(example_deg_results)
#' }
plot_module_deg <- function(deg_results,
                            pval_threshold = 0.05,
                            point_size = 3,
                            jitter_width = 0.15) {

  # Process data with explicit package references
  plot_data <- deg_results %>%
    dplyr::mutate(
      module_number = as.numeric(stringr::str_extract(.data$merge_module, "(?<=_)\\d+")),
      merge_module = stats::reorder(.data$merge_module, .data$module_number)
    ) %>%
    dplyr::group_by(.data$merge_module) %>%
    dplyr::mutate(module_position = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  # Create plot with explicit package references
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$merge_module, y = .data$logFC)) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = as.numeric(.data$merge_module) - 0.5,
        xmax = as.numeric(.data$merge_module) + 0.5,
        ymin = min(.data$logFC) - 0.5,
        ymax = max(.data$logFC) + 0.5
      ),
      fill = NA, color = "gray70", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = ifelse(.data$P.Value < pval_threshold, .data$logFC, NA)),
      size = point_size,
      position = ggplot2::position_jitter(width = jitter_width, height = 0)
    ) +
    ggplot2::scale_color_gradient2(
      low = "skyblue", mid = "white", high = "red",
      midpoint = 0,
      na.value = "gray80",
      name = paste0("log2FC\n(P < ", pval_threshold, ")")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Gene Modules", y = "log2 Fold Change")
  # Add y-axis limits based on logFC range
  p
  return(p)
}

