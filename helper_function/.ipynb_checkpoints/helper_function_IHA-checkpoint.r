#' Read Pseudobulk Expression Data
#'
#' This function reads a list of CSV files containing pseudobulk expression data
#' and returns a list of data frames. Each file is read in parallel using the specified 
#' number of cores.
#'
#' @param file_list A character vector containing the paths to the CSV files.
#' @param mc_cores An integer specifying the number of cores to use for parallel processing. Default is 60.
#'
#' @return A list of data frames, each corresponding to a CSV file in the `file_list`.
#'
#' @import parallel
#' @export
read_pseudobulk_expression <- function(file_list, mc_cores = 60) {
  
  # Measure the total time taken to read and process all files
  total_time <- system.time({
    df_list <- mclapply(file_list, function(x) {
      df <- read.csv(x, check.names = FALSE, row.names = 1)
      colnames(df) <- paste0(gsub("^.*/(.*)\\.csv$", "\\1", x), ":", colnames(df))
      return(df)
    }, mc.cores = mc_cores)
  })
  
  print(paste("Total reading time:", total_time["elapsed"], "seconds"))

  # Check if the length of the list is the same as the length of the input path
  if (length(df_list) == length(file_list)) {
    print("The length of the list matches the length of the input path.")
  } else {
    warning("The length of the list does not match the length of the input path.")
  }
    
  return(df_list)
}

#' Filter Genes and Cell Types of Interest
#'
#' This function filters genes and cell types of interest from a list of data frames.
#' It processes each data frame to select the specified genes and cell types, and
#' optionally converts the result to a long format.
#'
#' @param df_list A list of data frames containing expression data.
#' @param selected_genes A character vector of genes to select.
#' @param celltype_single A character string specifying the cell type to select.
#' @param longformat A logical value indicating whether to return the data in long format. Default is TRUE.
#'
#' @return If `longformat` is TRUE, a data frame in long format. Otherwise, a combined matrix of selected data.
#'
#' @import reshape2
#' @export

filter_genes_and_celltype <- function(df_list, selected_genes, celltype_single, longformat = TRUE) {
  # Process each dataframe in the list
  celltype_list <- lapply(df_list, function(df) {
    # Check if the cell type is in the column names
    if (any(grepl(celltype_single, names(df), fixed = TRUE))) {
      df <- df[selected_genes, grep(celltype_single, names(df), fixed = TRUE), drop = FALSE]
    } else {
      warning(paste("Cell type", celltype_single, "not found in the dataframe."))
      return(NULL)
    }
    return(df)
  })
  
  # Remove any NULL elements from the list
  celltype_list <- celltype_list[!sapply(celltype_list, is.null)]

  # Combine the processed dataframes into one matrix
  exp_matrix <- do.call(cbind, celltype_list)
  if (longformat) {
    # Convert the matrix to long format
    exp_matrix$RowNames <- rownames(exp_matrix)
    long_format <- melt(exp_matrix, id.vars = "RowNames")
    colnames(long_format) <- c('Gene', 'Samples', 'Mean_Expression')
    
    # Split the sample names to extract pbmc_sample_id and celltype
    split_strings <- strsplit(as.vector(long_format$Samples), ":")
    long_format$pbmc_sample_id <- sapply(split_strings, `[`, 1)
    long_format$celltype <- sapply(split_strings, `[`, 2)
    
    return(long_format)
  } else {
    return(exp_matrix)
  }
}


#' DESeq2 Differential Expression Analysis
#'
#' This function performs differential expression analysis using DESeq2.
#' It takes an expression matrix, metadata, a set of filtered genes, a design formula,
#' a cell type, and a list of comparisons to perform.
#'
#' @param exp_matrix A matrix of expression data with genes as rows and samples as columns.
#' @param meta_data A data frame of metadata with rows corresponding to samples in the expression matrix.
#' @param filtered_gene_set A vector containing the set of genes to be analyzed.
#' @param formula A formula specifying the design of the experiment.
#' @param celltype A character string specifying the cell type being analyzed.
#' @param comparisons A list of comparisons to perform. Each comparison is a list where the first element is the contrast name and the subsequent elements are the levels to be compared.
#'
#' @return A data frame containing the combined results of all comparisons, including gene names, log2 fold changes, adjusted p-values, direction of change, and the cell type.
# @exmaple     
# res=deseq2_analysis(exp_matrix,
#                      meta_data=meta_data_subset,
#                      filtered_gene_set=filtered_gene_set_filtered,
#                      formula= ~  cohort.cohortGuid+subject.biologicalSex+CMV,
#                      comparisons=list(c("subject.biologicalSex", "Male", "Female"),
#                                       c("cohort.cohortGuid", "BR2", "BR1"),
#                                       c("CMV", "Positive", "Negative")),
#                      celltype=celltype)
deseq2_analysis <- function(exp_matrix, meta_data, filtered_gene_set, formula, celltype, comparisons,...) {
  # Ensure colnames of exp_matrix match rownames of meta_data
  if (!all(colnames(exp_matrix) %in% rownames(meta_data))) {
    stop("Column names of the expression matrix must match row names of the metadata.")
  }
  
  # Create DESeqDataSet object from the expression matrix and metadata
  # check filtered_gene_set input 
  dds <- DESeqDataSetFromMatrix(
    exp_matrix[filtered_gene_set, ], 
    colData = meta_data[colnames(exp_matrix), ,drop=FALSE], 
    design = formula
  )
  
  # Perform differential expression analysis
  dds <- DESeq(dds, parallel = FALSE)
  
  # Initialize results list to store the results of each comparison
  results_list <- lapply(comparisons, function(comparison) {
    # Extract contrast name and levels
    contrast_name <- comparison[[1]]
    contrast_levels <- comparison[-1]
    
    # Perform differential expression analysis for the given contrast
    res <- data.frame(results(dds, contrast = c(contrast_name, 
                                                contrast_levels[1], 
                                                contrast_levels[2]),...))
    res$contrast <- contrast_name
    
    # Arrange results by adjusted p-value and add direction of change
    res <- res %>% 
      arrange(padj) %>% 
      mutate(Direction = case_when(
        log2FoldChange < 0 ~ paste0('HigherIn', contrast_levels[2]),
        log2FoldChange > 0 ~ paste0('HigherIn', contrast_levels[1])
      ))
    
    # Add gene names to the results
    res$gene <- rownames(res)
    rownames(res) <- NULL
    
    return(res)
  })
  
  # Combine results from all comparisons into a single data frame
  combined_results <- do.call(rbind, results_list)
  combined_results$celltype <- celltype
  return(combined_results)
}

#check if vector contain zeros
clr_transform <- function(x) {
  if (length(x) == 0) {
    return(NA)  # return NA for empty vectors
  }
  geom_mean <- exp(mean(log(x)))
  return(log(x / geom_mean))
}
