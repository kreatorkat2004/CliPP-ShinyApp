#' Calculate recommended smoothing factor based on total mutations
#' @param total_mutations The total number of mutations
#' @return Numeric value representing the recommended smoothing factor
recommended_smoothing_factor <- function(total_mutations) {
  min_smooth <- 0.2
  max_smooth <- 0.7
  min_mut <- 100
  max_mut <- 10000
  
  smoothing_factor <- max_smooth - (total_mutations - min_mut) * (max_smooth - min_smooth) / (max_mut - min_mut)
  smoothing_factor <- max(min_smooth, min(max_smooth, smoothing_factor))
  
  return(round(smoothing_factor, 2))
}

#' Read and validate input files for analysis
#' @param snvFile Path to SNV file
#' @param cnaFile Path to CNA file
#' @param purityFile Path to purity file
#' @return List containing validated data frames or NULL if validation fails
read_and_validate_files <- function(snvFile, cnaFile, purityFile) {
  if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) {
    return(NULL)
  }
  
  # Read files
  df_snv <- tryCatch(
    read.table(snvFile, header = TRUE),
    error = function(e) NULL
  )
  
  df_cna <- tryCatch(
    read.table(cnaFile, header = TRUE),
    error = function(e) NULL
  )
  
  purity <- tryCatch(
    as.numeric(readLines(purityFile, warn = FALSE)),
    error = function(e) NULL
  )
  
  # Validate SNV data columns
  required_cols <- c("chromosome_index", "position", "ref_count", "alt_count")
  if (!all(required_cols %in% colnames(df_snv))) {
    return(NULL)
  }
  
  list(
    df_snv = df_snv,
    df_cna = df_cna,
    purity = purity
  )
}

#' Process mutation data to calculate metrics
#' @param df_snv SNV data frame
#' @param df_cna CNA data frame
#' @param purity Purity value
#' @param df_assigned Assigned mutations data frame
#' @return List containing processed data and summary statistics
process_mutation_data <- function(df_snv, df_cna, purity, df_assigned) {
  # Process SNV and CNA data
  df_snvcna <- df_snv %>%
    mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
      inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
      if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
    }))
  
  # Calculate read depth metrics
  df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
  mean_rd <- mean(df_snvcna$sum_counts)
  min_rd <- min(df_snvcna$sum_counts)
  max_rd <- max(df_snvcna$sum_counts)
  
  # Merge with assigned mutations
  merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
  merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
  
  # Calculate mutation metrics
  merged_df$total_cn <- as.numeric(merged_df$total_cn)
  merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
  merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / 
    (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
  merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
  
  # Add clonality classification
  merged_df <- merged_df %>%
    mutate(
      Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
      ReadDepth = sum_counts
    )
  
  # Calculate subclonality counts
  subclonality_counts <- merged_df %>%
    count(Subclonality) %>%
    rename(n = n)
  
  list(
    merged_data = merged_df,
    subclonality_counts = subclonality_counts,
    read_depth = list(
      mean = mean_rd,
      min = min_rd,
      max = max_rd
    )
  )
}

#' Store uploaded files in a structured format
#' @param input Shiny input object
#' @param numSamples Number of samples to process
#' @return List containing organized file paths
storeUploadedFiles <- function(input, numSamples) {
  samples <- list()
  
  for (i in 1:numSamples) {
    sampleName <- input[[paste0("sampleName", i)]]
    
    if (input$uploadMode == "all_at_once") {
      files <- input[[paste0("fileInput", i)]]
      if (!is.null(files)) {
        uploadedFiles <- files$datapath
        uploadedNames <- files$name
        
        samples[[sampleName]] <- list(
          snvFile = uploadedFiles[grepl("snv.txt", uploadedNames)],
          cnaFile = uploadedFiles[grepl("cna.txt", uploadedNames)],
          purityFile = uploadedFiles[grepl("purity.txt", uploadedNames)]
        )
      }
    } else {
      snvFile <- input[[paste0("snvFile", i)]]
      cnaFile <- input[[paste0("cnaFile", i)]]
      purityFile <- input[[paste0("purityFile", i)]]
      
      if (!is.null(snvFile) && !is.null(cnaFile) && !is.null(purityFile)) {
        samples[[sampleName]] <- list(
          snvFile = snvFile$datapath,
          cnaFile = cnaFile$datapath,
          purityFile = purityFile$datapath
        )
      }
    }
  }
  
  return(samples)
}

#' Process driver mutation data
#' @param cancer_type Selected cancer type
#' @param driver_mutation_data Global driver mutation dataset
#' @param mode Either "gene" or "sample" for different analysis modes
#' @return Filtered and processed driver mutation data
process_driver_data <- function(cancer_type, driver_mutation_data, mode = "gene") {
  if (is.null(cancer_type) || cancer_type == "") {
    return(NULL)
  }
  
  filtered_data <- driver_mutation_data %>%
    filter(cancer == cancer_type)
  
  if (mode == "sample") {
    filtered_data <- filtered_data %>%
      group_by(gene) %>%
      summarise(
        protein = paste(unique(protein), collapse = ", "),
        consequence = paste(unique(consequence), collapse = ", "),
        Subclonality = paste(unique(Subclonality), collapse = ", "),
        VAF = paste(sprintf("%.3f", VAF), collapse = ", "),
        CP_unpenalized = paste(sprintf("%.3f", CP_unpenalized), collapse = ", "),
        .groups = 'drop'
      )
  }
  
  if (nrow(filtered_data) == 0) {
    return(NULL)
  }
  
  return(filtered_data)
}

#' Prepare data for plotting
#' @param merged_data Processed mutation data
#' @param color_by Color grouping variable ("clonality" or "depth")
#' @return List containing plot-ready data and aesthetics
prepare_plot_data <- function(merged_data, color_by = "clonality") {
  # Ensure no NA or NULL values in columns used for plotting
  plot_data <- merged_data %>%
    filter(!is.na(CP_unpenalized), !is.na(Subclonality), !is.na(ReadDepth))
  
  # Set color aesthetic based on user input
  color_aesthetic <- if (color_by == "clonality") "Subclonality" else "ReadDepth"
  
  list(
    plot_data = plot_data,
    color_aesthetic = color_aesthetic
  )
}