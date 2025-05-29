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
  colnames(select_results)[8] <- 'cluster'
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





