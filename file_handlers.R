library(shiny)

#' Create a download handler for result files
#' @param prefix Prefix for the file identifier (e.g., "TCGA", "PCAWG")
#' @param result_path Path to the result files
#' @return Shiny downloadHandler object
create_download_handler <- function(prefix, result_path) {
  downloadHandler(
    filename = function() {
      selected_file <- input[[paste0("selectedFile", prefix)]]
      if (grepl(".txt$", selected_file)) {
        return(selected_file)
      } else {
        return(paste0(selected_file, ".txt"))
      }
    },
    content = function(file) {
      filePath <- file.path(result_path, input[[paste0("selectedFile", prefix)]])
      if (!file.exists(filePath)) {
        stop("File not found.")
      }
      file.copy(from = filePath, to = file)
    },
    contentType = "text/plain"
  )
}

#' Prepare and run CliPP analysis
#' @param sampleDir Directory containing sample files
#' @param prefix Sample prefix
#' @param selectedSample Selected sample name
#' @return Boolean indicating success
prepareAndRunCliPP <- function(sampleDir, prefix, selectedSample) {
  # Check if mutation assignments already exist
  matching_file1 <- list.files(sampleDir, 
                               pattern = "mutation_assignments_lam.*\\.txt", 
                               full.names = TRUE)
  matching_file2 <- list.files(sampleDir, 
                               pattern = "subclonal_structure_lam.*\\.txt", 
                               full.names = TRUE)
  
  if (length(matching_file1) == 0 || length(matching_file2) == 0) {
    # Set up input files
    snvFile <- file.path(sampleDir, "sample.snv.txt")
    cnaFile <- file.path(sampleDir, "sample.cna.txt")
    purityFile <- file.path(sampleDir, "sample.purity.txt")
    
    # Construct and run CliPP command
    command <- sprintf("/usr/local/bin/python3 ./forCliPP/CliPP/run_clipp_main.py %s %s %s",
                       shQuote(snvFile), shQuote(cnaFile), shQuote(purityFile))
    
    result <- tryCatch({
      system(command, intern = TRUE)
      TRUE
    }, error = function(e) {
      showNotification(paste("Error running CliPP:", e$message), type = "error")
      FALSE
    })
    
    if (!result) return(FALSE)
    
    # Move results to sample directory
    resultFiles <- list.files("forCliPP/CliPP/sample_id/final_result/Best_lambda/",
                              pattern = "\\.txt$", 
                              full.names = TRUE)
    
    file.copy(resultFiles, sampleDir)
  }
  
  # Verify results
  matching_file1 <- list.files(sampleDir, 
                               pattern = "mutation_assignments_lam.*\\.txt", 
                               full.names = TRUE)
  matching_file2 <- list.files(sampleDir, 
                               pattern = "subclonal_structure_lam.*\\.txt", 
                               full.names = TRUE)
  
  if (length(matching_file1) == 0 || length(matching_file2) == 0) {
    showNotification("Failed to generate mutation assignments", type = "error")
    return(FALSE)
  }
  
  return(TRUE)
}

#' Read and validate input files
#' @param snvFile Path to SNV file
#' @param cnaFile Path to CNA file
#' @param purityFile Path to purity file
#' @return List containing validation results and data
validate_and_read_files <- function(snvFile, cnaFile, purityFile) {
  if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) {
    return(list(
      success = FALSE,
      message = "One or more required files are missing.",
      data = NULL
    ))
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
    missing_cols <- required_cols[!required_cols %in% colnames(df_snv)]
    return(list(
      success = FALSE,
      message = paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
      data = NULL
    ))
  }
  
  return(list(
    success = TRUE,
    message = "Files read successfully",
    data = list(
      df_snv = df_snv,
      df_cna = df_cna,
      purity = purity
    )
  ))
}

#' Store uploaded files in a structured format
#' @param input Shiny input object
#' @param numSamples Number of samples
#' @param uploadMode Upload mode ("all_at_once" or "individual")
#' @return List containing organized file paths
store_uploaded_files <- function(input, numSamples, uploadMode) {
  samples <- list()
  
  for (i in 1:numSamples) {
    sampleName <- input[[paste0("sampleName", i)]]
    
    if (uploadMode == "all_at_once") {
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

#' Create driver mutation download handler
#' @param cancer_type Selected cancer type
#' @param filtered_data Filtered driver mutation data
#' @return Shiny downloadHandler object
create_driver_download_handler <- function(cancer_type, filtered_data) {
  downloadHandler(
    filename = function() {
      paste0(
        "Driver_Mutations_",
        cancer_type, "_",
        format(Sys.Date(), "%Y%m%d"),
        ".txt"
      )
    },
    content = function(file) {
      if (!is.null(filtered_data) && nrow(filtered_data) > 0) {
        write.table(
          filtered_data,
          file,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
      } else {
        writeLines("No data available", file)
      }
    },
    contentType = "text/plain"
  )
}

#' Check if required directories exist and create them if needed
#' @param base_dir Base directory path
#' @param required_dirs Vector of required directory names
#' @return Boolean indicating success
setup_directories <- function(base_dir, required_dirs) {
  success <- TRUE
  
  for (dir in required_dirs) {
    dir_path <- file.path(base_dir, dir)
    if (!dir.exists(dir_path)) {
      dir_created <- dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      if (!dir_created) {
        warning(paste("Failed to create directory:", dir_path))
        success <- FALSE
      }
    }
  }
  
  return(success)
}

#' Process uploaded sample files
#' @param sampleFiles List containing file paths for SNV, CNA, and purity files
#' @param sampleName Name of the sample being processed
#' @return List containing processing success status, message, and data if successful
process_sample_files <- function(sampleFiles, sampleName) {
  # Validate input files
  if (is.null(sampleFiles$snvFile) || is.null(sampleFiles$cnaFile) || is.null(sampleFiles$purityFile)) {
    return(list(
      success = FALSE,
      message = "One or more required files are missing."
    ))
  }
  
  # Read and validate files
  file_data <- read_and_validate_files(sampleFiles$snvFile, sampleFiles$cnaFile, sampleFiles$purityFile)
  
  if (is.null(file_data)) {
    return(list(
      success = FALSE,
      message = "Error reading or validating input files."
    ))
  }
  
  df_snv <- file_data$df_snv
  df_cna <- file_data$df_cna
  purity <- file_data$purity
  
  # Create/ensure the directory exists
  directoryPath <- "forCliPP/CliPP/sample_id/"
  if (!dir.exists(directoryPath)) {
    dir.create(directoryPath, recursive = TRUE, showWarnings = TRUE)
  } else {
    system(paste("rm -rf", shQuote(directoryPath)), intern = TRUE)
    dir.create(directoryPath, recursive = TRUE, showWarnings = TRUE)
  }
  
  # Run CliPP analysis
  command <- sprintf("/usr/local/bin/python3 ./forCliPP/CliPP/run_clipp_main.py %s %s %s",
                     shQuote(sampleFiles$snvFile), shQuote(sampleFiles$cnaFile), shQuote(sampleFiles$purityFile))
  
  result <- tryCatch({
    system(command, intern = TRUE)
    TRUE
  }, error = function(e) {
    return(list(
      success = FALSE,
      message = paste("Error running CliPP:", e$message)
    ))
  })
  
  if (!is.logical(result) || !result) {
    return(list(
      success = FALSE,
      message = "Failed to run CliPP analysis."
    ))
  }
  
  # Path to result files
  resultPath <- "forCliPP/CliPP/sample_id/final_result/Best_lambda"
  
  # Check if result files were created
  matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
  matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
  
  if (length(matching_file1) == 0 || length(matching_file2) == 0) {
    return(list(
      success = FALSE,
      message = "CliPP analysis did not produce expected output files."
    ))
  }
  
  # Read mutation and subclonal structure assignments
  df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
  df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
  
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
  
  # Add clonality classification based on cluster_index
  # FIX: Use cluster_index, not Subclonality (which doesn't exist yet)
  merged_df <- merged_df %>%
    mutate(
      Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
      ReadDepth = sum_counts
    )
  
  # Make Subclonality a factor with defined levels to avoid color scale warnings
  merged_df$Subclonality <- factor(merged_df$Subclonality, levels = c("Clonal", "Subclonal"))
  
  # Calculate subclonality counts
  subclonality_counts <- merged_df %>%
    count(Subclonality) %>%
    rename(n = n)
  
  # Compute total mutations count
  total_SNVs <- sum(df_subc$num_SNV)
  
  # Calculate recommended smoothing factor
  recommended_smooth1 <- recommended_smoothing_factor(total_SNVs)
  recommended_smooth2 <- recommended_smoothing_factor(total_SNVs)
  
  # Update available result files for download
  result_files <- list.files(resultPath, pattern = "\\.txt$")
  if (exists("session") && is.function(updateSelectInput)) {
    updateSelectInput(session, "selectedFile", choices = result_files)
  }
  
  # Return success with processed data
  return(list(
    success = TRUE,
    message = "Sample processed successfully.",
    data = list(
      merged_data = merged_df,
      subclonality_counts = subclonality_counts,
      read_depth = list(
        mean = mean_rd,
        min = min_rd,
        max = max_rd
      ),
      purity = purity,
      sample_name = sampleName,
      subclonal_structure = df_subc,
      total_mutations = total_SNVs
    )
  ))
}

#' Helper function to get current analysis data
#' @param sampleFiles List containing file paths for the selected sample
#' @return Current analysis data or NULL if processing fails
get_current_analysis <- function(sampleFiles) {
  if (is.null(sampleFiles)) return(NULL)
  
  # Get sample name
  sample_name <- names(sampleFiles)[1]
  if (is.null(sample_name)) sample_name <- "Sample"
  
  # Process the sample files
  process_result <- process_sample_files(sampleFiles, sample_name)
  
  if (!process_result$success) {
    return(NULL)
  }
  
  # Return the processed data
  return(process_result$data)
}

#' Clean up temporary files
#' @param temp_dir Temporary directory path
#' @return Boolean indicating success
cleanup_temp_files <- function(temp_dir) {
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
  return(!dir.exists(temp_dir))
}