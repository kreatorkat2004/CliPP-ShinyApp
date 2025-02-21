library(reticulate)
library(ggplot2)
library(dplyr)
library(purrr)
library(plotly)
library(readr)
library(shinyjs)
library(DT)
library(tidyr)

#old_path <- Sys.getenv("PATH")
#Sys.setenv(PATH = paste("/usr/local/R/4.2.0/bin", old_path, sep = ":"))
#Sys.setenv(LD_LIBRARY_PATH="/usr/local/gcc/7.0.0/lib64")

# Data source directory
global_TCGA_PCAWG <- "dataSource"

# Sample data from PCAWG and TCGA datasets
sample_PCAWG <- read.csv(file.path(global_TCGA_PCAWG, "CliPP_PCAWG.tsv"), sep='\t')
sample_TCGA <- read.csv(file.path(global_TCGA_PCAWG, "CliPP_TCGA.tsv"), sep='\t')

# Dataset identifiers for PCAWG and TCGA
sample_PCAWG$dataset <- "PCAWG"
sample_TCGA$dataset <- "TCGA"

# Both datasets combined in a single dataframe
sample_data <- rbind(sample_PCAWG, sample_TCGA)

# Load driver mutation data
driver_mutation_data <- read.delim("./driver_mut_data.tsv")

driver_mutation_data$cancer <- gsub("RCCC", "KIRC", driver_mutation_data$cancer)
driver_mutation_data$cancer <- gsub("COREAD", "CRC", driver_mutation_data$cancer)
driver_mutation_data$cancer <- gsub("AML", "LAML", driver_mutation_data$cancer)
driver_mutation_data$cancer <- gsub("CM", "SKCM", driver_mutation_data$cancer)

# Recommended smoothing factor function
recommended_smoothing_factor <- function(total_mutations) {
  min_smooth <- 0.2
  max_smooth <- 0.7
  min_mut <- 100
  max_mut <- 10000
  smoothing_factor <- max_smooth - (total_mutations - min_mut) * (max_smooth - min_smooth) / (max_mut - min_mut)
  smoothing_factor <- max(min_smooth, min(max_smooth, smoothing_factor))
  return(round(smoothing_factor, 2))
}

# Plot annotations function
create_annotations <- function(subclonality_counts) {
  annotations <- list()
  
  if ("Clonal" %in% subclonality_counts$Subclonality) {
    n_clonal <- subclonality_counts$n[subclonality_counts$Subclonality == "Clonal"]
    if (length(n_clonal) > 0 && n_clonal > 0) {
      annotations <- c(
        annotations, 
        annotate("text", x = 0.85, y = -1.65, label = paste0("n clonal = ", n_clonal), 
                 size = 3, color = "black", vjust = -1.5)
      )
    }
  }
  
  if ("Subclonal" %in% subclonality_counts$Subclonality) {
    n_subclonal <- subclonality_counts$n[subclonality_counts$Subclonality == "Subclonal"]
    if (length(n_subclonal) > 0 && n_subclonal > 0) {
      annotations <- c(
        annotations, 
        annotate("text", x = 0.85, y = -2.65, label = paste0("n subclonal = ", n_subclonal), 
                 size = 3, color = "black", vjust = -1.5)
      )
    }
  }
  
  return(annotations)
}

# Function to allow for file handling 
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

# File download function
create_download_handler <- function(path) {
  downloadHandler(
    filename = function() {
      if (grepl(".txt$", input$selectedFile)) {
        return(input$selectedFile)
      } else {
        return(paste(input$selectedFile, ".txt", sep=""))
      }
    },
    content = function(file) {
      filePath <- file.path(path, input$selectedFile)
      if (!file.exists(filePath)) {
        stop("File not found.")
      }
      file.copy(from = filePath, to = file)
    },
    contentType = "text/plain"
  )
}

# Backup CliPP in case mutation files do not exist
prepareAndRunCliPP <- function(sampleDir, prefix, selectedSample) {
  # Check if mutation assignments already exist
  matching_file1 <- list.files(sampleDir, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
  matching_file2 <- list.files(sampleDir, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
  
  if (length(matching_file1) == 0 || length(matching_file2) == 0) {
    
    snvFile <- file.path(sampleDir, "sample.snv.txt")
    cnaFile <- file.path(sampleDir, "sample.cna.txt")
    purityFile <- file.path(sampleDir, "sample.purity.txt")
    
    # Construct and run CliPP command
    command <- sprintf("/usr/local/bin/python3 ./forCliPP/CliPP/run_clipp_main.py %s %s %s",
                       shQuote(snvFile), shQuote(cnaFile), shQuote(purityFile))
    
    result <- system(command, intern = TRUE)
    
    # Move results to sample directory
    resultFiles <- list.files("forCliPP/CliPP/sample_id/final_result/Best_lambda/", 
                              pattern = "\\.txt$", full.names = TRUE)
    file.copy(resultFiles, sampleDir)
  }
  
  # Check again for mutation assignment files
  matching_file1 <- list.files(sampleDir, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
  matching_file2 <- list.files(sampleDir, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
  
  if (length(matching_file1) == 0 || length(matching_file2) == 0) {
    showNotification("Failed to generate mutation assignments", type = "error")
    return(FALSE)
  }
  
  return(TRUE)
}

# Function to create download handlers
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

# Functions to write captions
get_vaf_plot_caption <- function() {
  HTML("<p> - <b>Variant Allele Frequency (VAF)</b> indicates the proportion of sequencing reads that support a variant allele.<br>
       - <b>Purity</b> reflects the proportion of tumor cells in the sample. For instance, a purity of 0.9 means 90% of the cells are cancerous.<br>
       - <b>Mean read depth</b> is the average number of sequencing reads per SNV, influencing data reliability.<br>
       - <b>Density</b> estimates the probability density function of VAFs based on kernel density estimation (KDE). <br>
  </p>")
}

get_cp_plot_caption <- function() {
  HTML("<p> - <b>Cellular Prevalence (CP)</b> the fraction of all cells (both tumor and admixed normal cells) from the sequenced tissue carrying a particular SNV.<br>
       - Annotations on the x-axis indicate whether mutations are <b>clonal</b> (present in all tumor cells) or <b>subclonal</b> (found in a subset).<br>
       - Counts of clonal and subclonal mutations are indicated as <b>n clonal</b> and <b>n subclonal</b> respectively.<br>
  </p>")
}

server <- function(input, output, session) {
  
  # Reactive values to track plot readiness for different sections
  plotsReadyCliPP <- reactiveVal(FALSE)
  plotsReadyTCGA <- reactiveVal(FALSE)
  plotsReadyPCAWG <- reactiveVal(FALSE)
  
  samplesUploaded <- reactiveVal(FALSE)
  uploadedSamples <- reactiveValues(samples = list())
  
  # Python path based on environment
  if (Sys.info()[['user']] == 'shiny') {
    Sys.setenv(RETICULATE_PYTHON = "./python3")
    #Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
  }
  
  # Output directories for analysis results
  outputDir <- getwd()
  resultPath <- file.path(outputDir, "forCliPP/CliPP/sample_id/final_result/Best_lambda")
  
  # Function to read and validate input files (SNV, CNA, and purity files)
  read_and_validate_files <- function(snvFile, cnaFile, purityFile) {
    if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) {
      showNotification("One or more required files are missing.", type = "error")
      return(NULL)
    }
    
    # Input data
    df_snv <- read.table(snvFile, header = TRUE)
    df_cna <- read.table(cnaFile, header = TRUE)
    purity <- as.numeric(readLines(purityFile, warn = FALSE))
    
    required_cols <- c("chromosome_index", "position", "ref_count", "alt_count")
    if (!all(required_cols %in% colnames(df_snv))) {
      missing_cols <- required_cols[!required_cols %in% colnames(df_snv)]
      showNotification(paste("The SNV file does not contain the required columns:", paste(missing_cols, collapse = ", ")), type = "error")
      return(NULL)
    }
    
    list(df_snv = df_snv, df_cna = df_cna, purity = purity)
  }
  
  # Dynamic file upload inputs based on selected number of samples
  observeEvent(input$setSamples, {
    numSamples <- input$numSamples
    output$dynamicFileInputs <- renderUI({
      lapply(1:numSamples, function(i) {
        tagList(
          h4(paste("Sample", i)),
          textInput(paste0("sampleName", i), "Sample Name:", value = paste0("Sample", i)),
          # File upload mode (all files at once or individual uploads)
          if (input$uploadMode == "all_at_once") {
            fileInput(paste0("fileInput", i), "Upload all three required files (snv.txt, cna.txt, purity.txt):",
                      multiple = TRUE, accept = c(".txt"))
          } else {
            tagList(
              fileInput(paste0("snvFile", i), "SNV file:", multiple = FALSE, accept = c(".txt")),
              fileInput(paste0("cnaFile", i), "CNA file:", multiple = FALSE, accept = c(".txt")),
              fileInput(paste0("purityFile", i), "Purity file:", multiple = FALSE, accept = c(".txt")),
              br()
            )
          }
        )
      })
    })
    shinyjs::show("uploadSamples")
    shinyjs::hide("setSamples")
    shinyjs::disable("numSamples")
  })
  
  # Hide plots when a new sample is selected for analysis
  observeEvent(input$selectedSample, {
    shinyjs::hide("plot1")
    shinyjs::hide("plot2")
    shinyjs::hide("caption1")
    shinyjs::hide("caption2")
  })
  
  # Toggle visibility of TCGA and PCAWG submit options
  observeEvent(input$submitTCGA, {
    shinyjs::hide("submitTCGA")
    plotsReadyTCGA(TRUE)
    observe({
      req(plotsReadyTCGA())
      shinyjs::show("colorByTCGA")
      shinyjs::show("clearTCGA")
      shinyjs::show(selector = ".controls-area")
    })
  })
  
  observeEvent(input$clearTCGA, {
    updateSelectInput(session, "cancer_type_TCGA", selected = character(0))
    updateSelectInput(session, "sample_name_TCGA", selected = character(0))
    
    elements_to_hide <- c(
      "plot1TCGA", "plot2TCGA",
      "caption1TCGA", "caption2TCGA",
      "smoothingFactor1_TCGA", "smoothingFactor2_TCGA",
      "smoothingExplanation1_TCGA", "smoothingExplanation2_TCGA",
      "colorByTCGA", "clearTCGA"
    )
    
    for (element in elements_to_hide) {
      shinyjs::hide(element)
    }
    
    output$plot1TCGA <- renderPlotly({ NULL })
    output$plot2TCGA <- renderPlotly({ NULL })
    output$caption1TCGA <- renderUI({ NULL })
    output$caption2TCGA <- renderUI({ NULL })
    
    shinyjs::show("submitTCGA")
    plotsReadyTCGA(FALSE)
    
    shinyjs::enable("cancer_type_TCGA")
    shinyjs::enable("sample_name_TCGA")
  })
  
  observe({
    req(input$cancer_type_TCGA)
    if (!is.null(input$cancer_type_TCGA) && input$cancer_type_TCGA != "") {
      shinyjs::enable("sample_name_TCGA")
      shinyjs::enable("submitTCGA")
    }
  })
  
  observeEvent(input$submitPCAWG, {
    shinyjs::hide("submitPCAWG")
    plotsReadyPCAWG(TRUE)
    observe({
      req(plotsReadyPCAWG())
      shinyjs::show("colorByPCAWG")
      shinyjs::show("clearPCAWG")
      shinyjs::show(selector = ".controls-area")
    })
  })
  
  observeEvent(input$clearPCAWG, {
    updateSelectInput(session, "cancer_type_PCAWG", selected = character(0))
    updateSelectInput(session, "sample_name_PCAWG", selected = character(0))
    
    elements_to_hide <- c(
      "plot1PCAWG", "plot2PCAWG",
      "caption1PCAWG", "caption2PCAWG",
      "smoothingFactor1_PCAWG", "smoothingFactor2_PCAWG",
      "smoothingExplanation1_PCAWG", "smoothingExplanation2_PCAWG",
      "colorByPCAWG", "clearPCAWG"
    )
    
    for (element in elements_to_hide) {
      shinyjs::hide(element)
    }
    
    output$plot1PCAWG <- renderPlotly({ NULL })
    output$plot2PCAWG <- renderPlotly({ NULL })
    output$caption1PCAWG <- renderUI({ NULL })
    output$caption2PCAWG <- renderUI({ NULL })
    
    shinyjs::show("submitPCAWG")
    plotsReadyPCAWG(FALSE)
    
    shinyjs::enable("cancer_type_PCAWG")
    shinyjs::enable("sample_name_PCAWG")
  })
  
  observe({
    req(input$cancer_type_PCAWG)
    if (!is.null(input$cancer_type_PCAWG) && input$cancer_type_PCAWG != "") {
      shinyjs::enable("sample_name_PCAWG")
      shinyjs::enable("submitPCAWG")
    }
  })
  
  # Function to store uploaded files
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
        
        # Stores the samples so that they can be uploaded at a later time
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
  
  # Store uploaded samples and show relevant UI elements for further processing
  observeEvent(input$uploadSamples, {
    numSamples <- input$numSamples
    samples <- storeUploadedFiles(input, numSamples)
    uploadedSamples$samples <- samples
    samplesUploaded(TRUE)
    updateSelectInput(session, "selectedSample", choices = names(samples))
    updateSelectInput(session, "circosSample", choices = names(samples)) 
    showNotification("Samples uploaded successfully.", type = "message")
    shinyjs::hide("uploadSamples")
    shinyjs::show("selectedSample")
    shinyjs::show("submit")
    shinyjs::show("colorBy")
  })
  
  # Clear current session data and reset UI for new uploads
  observeEvent(input$clear, {
    uploadedSamples$samples <- list()
    samplesUploaded(FALSE)
    updateNumericInput(session, "numSamples", value = 1)
    shinyjs::enable("numSamples")
    shinyjs::reset("dynamicFileInputs")
    shinyjs::reset("dynamicFileInputs")
    
    elements_to_hide <- c(
      "selectedSample",
      "submit",
      "colorBy",
      "clear",
      "plot1",
      "plot2",
      "caption1",
      "caption2",
      "smoothingFactor1",
      "smoothingFactor2",
      "smoothingExplanation1",
      "smoothingExplanation2"
    )
    
    for (element in elements_to_hide) {
      shinyjs::hide(element)
    }
    
    shinyjs::show("setSamples")
    shinyjs::show("uploadSamples")
    
    output$plot1 <- renderPlotly({ NULL })
    output$plot2 <- renderPlotly({ NULL })
    output$caption1 <- renderUI({ NULL })
    output$caption2 <- renderUI({ NULL })
    output$plotOutputArea <- renderUI({ NULL })
    
    plotsReadyCliPP(FALSE)
    
    updateSelectInput(session, "selectedFile", choices = NULL)
  })
  
  # Observer for main navigation tab changes
  observeEvent(input$`navbar-container`, {
    shinyjs::hide("smoothingFactor1")
    shinyjs::hide("smoothingFactor2")
    shinyjs::hide("smoothingFactor1_TCGA")
    shinyjs::hide("smoothingFactor2_TCGA")
    shinyjs::hide("smoothingFactor1_PCAWG")
    shinyjs::hide("smoothingFactor2_PCAWG")
    shinyjs::hide("smoothingFactor1_driver")
    shinyjs::hide("smoothingFactor2_driver")
    
    # Hide explanations too
    shinyjs::hide("smoothingExplanation1")
    shinyjs::hide("smoothingExplanation2")
    shinyjs::hide("smoothingExplanation1_TCGA")
    shinyjs::hide("smoothingExplanation2_TCGA")
    shinyjs::hide("smoothingExplanation1_PCAWG")
    shinyjs::hide("smoothingExplanation2_PCAWG")
    shinyjs::hide("smoothingExplanation1_driver")
    shinyjs::hide("smoothingExplanation2_driver")
  })
  
  # Observer for CliPP Data Resources tab panel changes
  observeEvent(input$`tabset`, {
    shinyjs::hide("smoothingFactor1_TCGA")
    shinyjs::hide("smoothingFactor2_TCGA")
    shinyjs::hide("smoothingFactor1_PCAWG")
    shinyjs::hide("smoothingFactor2_PCAWG")
    shinyjs::hide("smoothingExplanation1_TCGA")
    shinyjs::hide("smoothingExplanation2_TCGA")
    shinyjs::hide("smoothingExplanation1_PCAWG")
    shinyjs::hide("smoothingExplanation2_PCAWG")
  })
  
  # Generate plots and other outputs only after submitting all required files
  observeEvent(input$submit, {
    withProgress(message = 'Processing data...', value = 0, {
      plotsReadyCliPP(FALSE)
      
      selectedSample <- input$selectedSample
      sampleFiles <- uploadedSamples$samples[[selectedSample]]
      
      if (is.null(sampleFiles)) {
        showNotification("Please select a sample to analyze.", type = "error")
        return()
      }
      
      snvFile <- sampleFiles$snvFile
      cnaFile <- sampleFiles$cnaFile
      purityFile <- sampleFiles$purityFile
      
      if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) {
        showNotification("One or more required files are missing.", type = "error")
        return()
      }
      
      incProgress(0.1, detail = "Reading and validating files...")
      file_data <- read_and_validate_files(snvFile, cnaFile, purityFile)
      if (is.null(file_data)) return()
      
      df_snv <- file_data$df_snv
      df_cna <- file_data$df_cna
      purity <- file_data$purity
      
      incProgress(0.3, detail = "Preparing data...")
      
      directoryPath <- "forCliPP/CliPP/sample_id/"
      
      if (!dir.exists(directoryPath)) {
        dir.create(directoryPath, recursive = TRUE, showWarnings = TRUE)
      } else {
        system(paste("rm -rf", shQuote(directoryPath)), intern = TRUE)
      }
      
      incProgress(0.5, detail = "Executing analysis...")
      
      command <- sprintf("/usr/local/bin/python3 ./forCliPP/CliPP/run_clipp_main.py %s %s %s",
                         shQuote(snvFile), shQuote(cnaFile), shQuote(purityFile))
      
      result <- system(command, intern = TRUE)
      
      # Debugging
      cat("Command output:", result, "\n")
      
      # Increment progress
      incProgress(0.8, detail = "Finalizing...")
      
      # Update the file selection dropdown
      updateSelectInput(session, "selectedFile", choices = list.files(resultPath, pattern = "\\.txt$"))
      
      # Output message
      output$result <- renderText({
        "Processing complete. Select a file to download."
      })
      
      # Mark plots as ready for display
      plotsReadyCliPP(TRUE)
      # Clear button after submission
      shinyjs::show("clear")
      shinyjs::show("plotOutputArea")
      shinyjs::show("plot1")
      shinyjs::show("plot2")
      shinyjs::show("caption1")
      shinyjs::show("caption2")
      shinyjs::show("plotOutputArea")
      shinyjs::show("controls-area")
      
      # Increment progress to completion
      incProgress(1, detail = "Done")
    })
  })
  
  # Handlers to enable file downloads after plot generation
  output$downloadResult <- create_download_handler("", resultPath)
  
  output$warning <- renderText({
    if (is.null(input$submit)) return()
  })
  
  output$plotReadyCliPP <- reactive({ plotsReadyCliPP() })
  outputOptions(output, "plotReadyCliPP", suspendWhenHidden = FALSE)
  
  # Dummy text output to trigger the conditional panel
  output$plotReadyTextCliPP <- renderUI({
    if (plotsReadyCliPP()) {
      div(id = "plotReadyTextCliPP", "true")
    } else {
      div(id = "plotReadyTextCliPP", "false")
    }
  })
  
  output$samplesUploaded <- reactive({ samplesUploaded() })
  outputOptions(output, "samplesUploaded", suspendWhenHidden = FALSE)
  
  # Plot 1
  output$plot1 <- renderPlotly({
    req(plotsReadyCliPP())
    # Initialize the progress bar for Plot 1
    withProgress(message = 'Generating Plot 1...', value = 0, {
      # Ensure files are uploaded
      if (input$submit == 0 || !plotsReadyCliPP()) return()
      
      # Increment progress
      incProgress(0.1, detail = "Loading data...")
      
      selectedSample <- input$selectedSample
      sampleFiles <- uploadedSamples$samples[[selectedSample]]
      
      if (is.null(sampleFiles)) return()
      
      snvFile <- sampleFiles$snvFile
      cnaFile <- sampleFiles$cnaFile
      purityFile <- sampleFiles$purityFile
      
      if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) return()
      
      file_data <- read_and_validate_files(snvFile, cnaFile, purityFile)
      if (is.null(file_data)) return()
      
      # Increment progress
      incProgress(0.3, detail = "Processing data...")
      
      df_snv <- file_data$df_snv
      df_cna <- file_data$df_cna
      purity <- file_data$purity
      
      # Read mutation and subclonal structure assignments
      matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) > 0 && length(matching_file2) > 0) {
        df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
        df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
        
        df_snvcna <- df_snv %>%
          mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
            inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
            if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
          }))
        
        # Calculate the sum of ref_count and alt_count
        df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
        
        # Compute the mean of the sum_counts column
        mean_rd <- mean(df_snvcna$sum_counts)
        min_rd <- min(df_snvcna$sum_counts)
        max_rd <- max(df_snvcna$sum_counts)
        
        merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
        merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
        
        merged_df$total_cn <- as.numeric(merged_df$total_cn)
        merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
        merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
        merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
        
        total_SNVs <- sum(df_subc$num_SNV)
        
        recommended_smooth1 <- recommended_smoothing_factor(total_SNVs)
        recommended_smooth2 <- recommended_smoothing_factor(total_SNVs)
        
        output$smoothingExplanation1 <- renderUI({
          HTML(paste0("The smoothing factor should be set to ", recommended_smooth1, " based on ", total_SNVs, " mutations."))
        })
        
        output$smoothingExplanation2 <- renderUI({
          HTML(paste0("The smoothing factor should be set to ", recommended_smooth2, " based on ", total_SNVs, " mutations."))
        })
        
        # Increment progress
        incProgress(0.6, detail = "Calculating plot parameters...")
        
        # Smoothing slider
        smoothing_factor1 <- input$smoothingFactor1
        
        merged_df <- merged_df %>%
          mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
                 ReadDepth = sum_counts)
        
        # Rename the count column to 'n'
        subclonality_counts <- merged_df %>%
          count(Subclonality) %>%
          rename(n = n)
        
        # Color aesthetic based on the user input
        color_aesthetic <- if (input$colorBy == "clonality") "Subclonality" else "ReadDepth"
        
        # Initialize annotations
        annotations <- create_annotations(subclonality_counts)
        
        # Increment progress
        incProgress(0.8, detail = "Rendering plot...")
        
        # Create the plot
        gg1 <- ggplot(merged_df, aes(x = AF)) +
          geom_density(fill = "gray", alpha = .6, adjust = smoothing_factor1) +
          geom_point(
            shape = 142,
            size = 5,
            alpha = 0.3,
            aes(
              y = -0.75,
              color = !!sym(color_aesthetic)
            )
          ) +
          xlim(0, 1) +
          ggtitle(paste(
            selectedSample,
            " - Purity =", round(purity, 3),
            "      ", "Mean read depth =", round(mean_rd),
            "      ", "Range read depth =", min_rd, "-", max_rd
          )) +
          annotations +  
          ggpubr::theme_pubr() + 
          xlab("VAF") + ylab("Density") + 
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        # Increment progress to completion
        incProgress(1, detail = "Done")
        
        ggplotly(gg1)
      }
    })
  })
  
  # Captions for graph 1 
  output$caption1 <- renderUI({
    if (plotsReadyCliPP()) {
      get_vaf_plot_caption()
    }
  })
  
  # Plot 2
  output$plot2 <- renderPlotly({
    req(plotsReadyCliPP())
    # Initialize the progress bar for Plot 2
    withProgress(message = 'Generating Plot 2...', value = 0, {
      if (input$submit == 0 || !plotsReadyCliPP()) return()
      
      # Increment progress
      incProgress(0.1, detail = "Loading data...")
      
      selectedSample <- input$selectedSample
      sampleFiles <- uploadedSamples$samples[[selectedSample]]
      
      if (is.null(sampleFiles)) return()
      
      snvFile <- sampleFiles$snvFile
      cnaFile <- sampleFiles$cnaFile
      purityFile <- sampleFiles$purityFile
      
      if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) return()
      
      file_data <- read_and_validate_files(snvFile, cnaFile, purityFile)
      if (is.null(file_data)) return()
      
      # Increment progress
      incProgress(0.3, detail = "Processing data...")
      
      df_snv <- file_data$df_snv
      df_cna <- file_data$df_cna
      purity <- file_data$purity
      
      # Read mutation and subclonal structure assignments
      matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) > 0 && length(matching_file2) > 0) {
        df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
        df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
        
        df_snvcna <- df_snv %>%
          mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
            inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
            if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
          }))
        
        # Calculate the sum of ref_count and alt_count
        df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
        
        # Compute the mean of the sum_counts column
        mean_rd <- mean(df_snvcna$sum_counts)
        min_rd <- min(df_snvcna$sum_counts)
        max_rd <- max(df_snvcna$sum_counts)
        
        merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
        merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
        
        merged_df$total_cn <- as.numeric(merged_df$total_cn)
        merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
        merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
        merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
        
        total_SNVs <- sum(df_subc$num_SNV)
        
        recommended_smooth1 <- recommended_smoothing_factor(total_SNVs)
        recommended_smooth2 <- recommended_smoothing_factor(total_SNVs)
        
        output$smoothingExplanation1 <- renderUI({
          HTML(paste0("The smoothing factor should be set to ", recommended_smooth1, " based on ", total_SNVs, " mutations."))
        })
        
        output$smoothingExplanation2 <- renderUI({
          HTML(paste0("The smoothing factor should be set to ", recommended_smooth2, " based on ", total_SNVs, " mutations."))
        })
        
        # Increment progress
        incProgress(0.6, detail = "Calculating plot parameters...")
        
        # Smoothing slider
        smoothing_factor2 <- input$smoothingFactor2
        
        merged_df <- merged_df %>%
          mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
                 ReadDepth = sum_counts)
        
        # Rename the count column to 'n'
        subclonality_counts <- merged_df %>%
          count(Subclonality) %>%
          rename(n = n)
        
        # Ensures no NA or NULL values in columns used for plotting
        merged_df <- merged_df %>%
          filter(!is.na(CP_unpenalized), !is.na(Subclonality), !is.na(ReadDepth))
        
        # Determine the color aesthetic based on the user input
        color_aesthetic <- if (input$colorBy == "clonality") "Subclonality" else "ReadDepth"
        
        # Initialize annotations
        annotations <- create_annotations(subclonality_counts)
        
        # Increment progress
        incProgress(0.8, detail = "Rendering plot...")
        
        # Create plot
        gg2 <- ggplot(merged_df, aes(x = CP_unpenalized)) +
          geom_density(fill = "seagreen", alpha = .6, adjust = smoothing_factor2) +
          geom_point(
            # Lines
            shape = 142,
            # Size
            size = 5,
            # Transparency
            alpha = 0.3,
            aes(
              # Fixed y-position
              y = -0.75,
              # Dynamic color
              color = !!sym(color_aesthetic)
            )
          ) +
          xlim(0, 1.5) +
          ggtitle(paste(selectedSample, " - Estimated mutation-specific CP")) +
          annotations +  
          ggpubr::theme_pubr() + xlab("Estimated mutation-specific CP") + ylab("Density") +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          geom_segment(data = df_subc, aes(x = cellular_prevalence, xend = cellular_prevalence, y = 0, yend = 5.7), lty = "dashed")
        
        # Increment progress to completion
        incProgress(1, detail = "Done")
        
        ggplotly(gg2)
      }
    })
  })
  
  # Captions for graph 2 
  output$caption2 <- renderUI({
    if (plotsReadyCliPP()) {
      get_cp_plot_caption()
    }
  })
  
  # Output file does not appear until after graphs are generated
  output$plotOutputArea <- renderUI({
    if (plotsReadyCliPP()) {
      fluidPage(
        h4("Output:"),
        textOutput("result"),
        selectInput("selectedFile", "Choose a file to download:", choices = c()),
        downloadButton("downloadResult", "Download Result")
      )
    }
  })
  
  # Set up handlers for TCGA and PCAWG sections with helper functions for each
  createSectionHandlers <- function(prefix) {
    plotsReadyPrefix <- switch(prefix,
                               "TCGA" = plotsReadyTCGA,
                               "PCAWG" = plotsReadyPCAWG
    )
    
    observeEvent(plotsReadyPrefix(), {
      if (plotsReadyPrefix()) {
        selectedSample <- input[[paste0("sample_name_", prefix)]]
        sampleDir <- file.path(global_TCGA_PCAWG, prefix, selectedSample)
        resultPath <- sampleDir  
        
        # File selection dropdown for result files
        result_files <- list.files(resultPath, pattern = "\\.txt$", full.names = FALSE)
        
        output[[paste0("downloadUI_", prefix)]] <- renderUI({
          tagList(
            selectInput(paste0("selectedFile", prefix), "Choose a file to download:", choices = result_files),
            downloadButton(paste0("downloadResult", prefix), "Download Result")
          )
        })
      }
    })
    
    # Download handler for TCGA/PCAWG results
    output[[paste0("downloadResult", prefix)]] <- create_download_handler(
      prefix,
      file.path(global_TCGA_PCAWG, prefix, input[[paste0("sample_name_", prefix)]])
    )
    
    observeEvent(input[[paste0("submit", prefix)]], {
      plotsReadyPrefix(FALSE)
      
      selectedSample <- input[[paste0("sample_name_", prefix)]]
      sampleDir <- file.path(global_TCGA_PCAWG, prefix, selectedSample)
      
      if (is.null(selectedSample) || selectedSample == "") {
        showNotification("Please select a sample to analyze.", type = "error")
        return()
      }
      
      snvFile <- list.files(sampleDir, pattern = "*.snv.txt", full.names = TRUE)
      cnaFile <- list.files(sampleDir, pattern = "*.cna.txt", full.names = TRUE)
      purityFile <- list.files(sampleDir, pattern = "*.purity.txt", full.names = TRUE)
      
      if (length(snvFile) == 0 || length(cnaFile) == 0 || length(purityFile) == 0) {
        showNotification("One or more required files are missing.", type = "error")
        return()
      }
      
      file_data <- read_and_validate_files(snvFile[1], cnaFile[1], purityFile[1])
      
      if (is.null(file_data)) return()
      
      df_snv <- file_data$df_snv
      df_cna <- file_data$df_cna
      purity <- file_data$purity
      
      # Cut-off for PCAWG
      if (prefix == "PCAWG") {
        num_SNVs <- nrow(df_snv)
        if (num_SNVs > 20000) {
          showNotification(paste("The selected sample has", num_SNVs, "SNVs, which exceeds the limit of 20,000 SNVs for PCAWG samples."), type = "error")
          return()
        }
      }
      
      # File selection dropdown
      updateSelectInput(session, paste0("selectedFile", prefix), choices = list.files(resultPath, pattern = "\\.txt$"))
      
      # Output message
      output[[paste0("result", prefix)]] <- renderText({
        "Processing complete. Select a file to download."
      })
      plotsReadyPrefix(TRUE)
      # Clear button only after submission
      shinyjs::show(paste0("clear", prefix))
    })
    
    observe({
      req(plotsReadyPrefix())
      shinyjs::show(paste0("smoothingFactor1_", prefix))
      shinyjs::show(paste0("smoothingFactor2_", prefix))
      shinyjs::show(paste0("colorBy", prefix))
      shinyjs::show(selector = paste0(".controls-area-", prefix))
      
      # Plot 1
      output[[paste0("plot1", prefix)]] <- renderPlotly({
        selectedSample <- input[[paste0("sample_name_", prefix)]]
        
        sampleDir <- file.path(global_TCGA_PCAWG, prefix, selectedSample)
        
        # Prepare data and run CliPP if needed
        if (!prepareAndRunCliPP(sampleDir, prefix, selectedSample)) {
          return()
        }
        
        if (is.null(selectedSample) || selectedSample == "") return()
        
        snvFile <- list.files(sampleDir, pattern = "*.snv.txt", full.names = TRUE)
        cnaFile <- list.files(sampleDir, pattern = "*.cna.txt", full.names = TRUE)
        purityFile <- list.files(sampleDir, pattern = "*.purity.txt", full.names = TRUE)
        
        if (length(snvFile) == 0 || length(cnaFile) == 0 || length(purityFile) == 0) return()
        
        file_data <- read_and_validate_files(snvFile[1], cnaFile[1], purityFile[1])
        if (is.null(file_data)) return()
        
        df_snv <- file_data$df_snv
        df_cna <- file_data$df_cna
        purity <- file_data$purity
        
        resultPath <- sampleDir
        
        # Read mutation and subclonal structure assignments
        matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
        matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
        
        if (length(matching_file1) > 0 && length(matching_file2) > 0) {
          df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
          df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
          
          df_snvcna <- df_snv %>%
            mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
              inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
              if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
            }))
          
          # Calculate the sum of ref_count and alt_count
          df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
          
          # Compute the mean of the sum_counts column
          mean_rd <- mean(df_snvcna$sum_counts)
          min_rd <- min(df_snvcna$sum_counts)
          max_rd <- max(df_snvcna$sum_counts)
          
          merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
          merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
          
          merged_df$total_cn <- as.numeric(merged_df$total_cn)
          merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
          merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
          merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
          
          total_SNVs <- sum(df_subc$num_SNV)
          
          recommended_smooth1 <- recommended_smoothing_factor(total_SNVs)
          recommended_smooth2 <- recommended_smoothing_factor(total_SNVs)
          
          output[[paste0("smoothingExplanation1_", prefix)]] <- renderUI({
            HTML(paste0("The smoothing factor should be set to ", recommended_smooth1, " based on ", total_SNVs, " mutations."))
          })
          
          output[[paste0("smoothingExplanation2_", prefix)]] <- renderUI({
            HTML(paste0("The smoothing factor should be set to ", recommended_smooth2, " based on ", total_SNVs, " mutations."))
          })
          
          # Smoothing slider
          smoothing_factor1 <- input[[paste0("smoothingFactor1_", prefix)]]
          
          merged_df <- merged_df %>%
            mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
                   ReadDepth = sum_counts)
          
          # Rename the count column to 'n'
          subclonality_counts <- merged_df %>%
            count(Subclonality) %>%
            rename(n = n)
          
          # Determine the color aesthetic based on the user input
          color_aesthetic <- if (input[[paste0("colorBy", prefix)]] == "clonality") "Subclonality" else "ReadDepth"
          
          # Initialize annotations
          annotations <- create_annotations(subclonality_counts)
          
          gg1 <- ggplot(merged_df, aes(x = AF)) +
            geom_density(fill = "gray", alpha = .6, adjust = smoothing_factor1) +
            geom_point(
              shape = 142,
              size = 5,
              alpha = 0.3,
              aes(
                y = -0.75,
                color = !!sym(color_aesthetic)
              )
            ) +
            xlim(0, 1) +
            ggtitle(paste(selectedSample, " - Purity =", round(purity, 3),
                          "      ", "Mean read depth =", round(mean_rd),
                          "      ", "Range read depth =", min_rd, "-", max_rd)) +
            annotations + 
            ggpubr::theme_pubr() + xlab("VAF") + ylab("Density") +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
          
          ggplotly(gg1)
        }
      })
      
      # Captions for graph 1 
      output[[paste0("caption1", prefix)]] <- renderUI({
        if (plotsReadyPrefix()) {
          get_vaf_plot_caption()
        }
      })
      
      # Plot 2
      output[[paste0("plot2", prefix)]] <- renderPlotly({
        if (input[[paste0("submit", prefix)]] == 0 || !plotsReadyPrefix()) return()
        
        selectedSample <- input[[paste0("sample_name_", prefix)]]
        sampleDir <- file.path(global_TCGA_PCAWG, prefix, selectedSample)
        
        if (is.null(selectedSample) || selectedSample == "") return()
        
        snvFile <- list.files(sampleDir, pattern = "*.snv.txt", full.names = TRUE)
        cnaFile <- list.files(sampleDir, pattern = "*.cna.txt", full.names = TRUE)
        purityFile <- list.files(sampleDir, pattern = "*.purity.txt", full.names = TRUE)
        
        if (length(snvFile) == 0 || length(cnaFile) == 0 || length(purityFile) == 0) return()
        
        file_data <- read_and_validate_files(snvFile[1], cnaFile[1], purityFile[1])
        if (is.null(file_data)) return()
        
        df_snv <- file_data$df_snv
        df_cna <- file_data$df_cna
        purity <- file_data$purity
        
        resultPath <- sampleDir
        
        # Read mutation and subclonal structure assignments
        matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
        matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
        
        if (length(matching_file1) > 0 && length(matching_file2) > 0) {
          df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
          df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
          
          df_snvcna <- df_snv %>%
            mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
              inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
              if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
            }))
          
          # Calculate the sum of ref_count and alt_count
          df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
          mean_rd <- mean(df_snvcna$sum_counts)
          min_rd <- min(df_snvcna$sum_counts)
          max_rd <- max(df_snvcna$sum_counts)
          
          merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
          merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
          
          merged_df$total_cn <- as.numeric(merged_df$total_cn)
          merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
          merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
          merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
          
          total_SNVs <- sum(df_subc$num_SNV)
          
          recommended_smooth1 <- recommended_smoothing_factor(total_SNVs)
          recommended_smooth2 <- recommended_smoothing_factor(total_SNVs)
          
          output[[paste0("smoothingExplanation1_", prefix)]] <- renderUI({
            HTML(paste0("The smoothing factor should be set to ", recommended_smooth1, " based on ", total_SNVs, " mutations."))
          })
          
          output[[paste0("smoothingExplanation2_", prefix)]] <- renderUI({
            HTML(paste0("The smoothing factor should be set to ", recommended_smooth2, " based on ", total_SNVs, " mutations."))
          })
          
          # Smoothing slider
          smoothing_factor2 <- input[[paste0("smoothingFactor2_", prefix)]]
          
          merged_df <- merged_df %>%
            mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
                   ReadDepth = sum_counts)
          
          # Rename the count column to 'n'
          subclonality_counts <- merged_df %>%
            count(Subclonality) %>%
            rename(n = n)
          
          # Determine the color aesthetic based on the user input
          color_aesthetic <- if (input[[paste0("colorBy", prefix)]] == "clonality") "Subclonality" else "ReadDepth"
          
          # Initialize annotations
          annotations <- create_annotations(subclonality_counts)
          
          gg2 <- ggplot(merged_df, aes(x = CP_unpenalized)) +
            geom_density(fill = "seagreen", alpha = .6, adjust = smoothing_factor2) +
            suppressWarnings(geom_point(
              # Lines
              shape = 142,
              # Size
              size = 5,
              # Transparency
              alpha = 0.3,
              # Fixed y-position and dynamic color
              aes(y = -0.75, color = !!sym(color_aesthetic))
            )) +
            xlim(0, 1.5) +
            ggtitle("Estimated mutation-specific CP") +
            annotations +  
            ggpubr::theme_pubr() + xlab("Estimated mutation-specific CP") + ylab("Density") +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            geom_segment(data = df_subc, aes(x = cellular_prevalence, xend = cellular_prevalence, y = 0, yend = 5.7), lty = "dashed")
          
          ggplotly(gg2)
        }
      })
      
      # Captions for graph 2 
      output[[paste0("caption2", prefix)]] <- renderUI({
        if (plotsReadyPrefix()) {
          get_cp_plot_caption()
        }
      })
      
      # Output file does not appear until after the graphs are generated
      output[[paste0("plotOutputArea", prefix)]] <- renderUI({
        if (plotsReadyPrefix()) {
          fluidPage(
            h4("Output:"),
            textOutput(paste0("result", prefix)),
            selectInput(paste0("selectedFile", prefix), "Choose a file to download:", choices = c()),
            downloadButton(paste0("downloadResult", prefix), "Download Result")
          )
        } else {
          shinyjs::hide(paste0("plotOutputArea", prefix))
        }
      })
    })
    
    output[[paste0("plotReady", prefix)]] <- reactive({ plotsReadyPrefix() })
    outputOptions(output, paste0("plotReady", prefix), suspendWhenHidden = FALSE)
    
    output[[paste0("plotReadyText", prefix)]] <- renderUI({
      if (plotsReadyPrefix()) {
        div(id = paste0("plotReadyText", prefix), "true")
      } else {
        div(id = paste0("plotReadyText", prefix), "false")
      }
    })
  }
  
  createSectionHandlers("TCGA")
  createSectionHandlers("PCAWG")
  
  list_files <- function(path) {
    list.files(path = path, full.names = TRUE)
  }
  
  observe({
    files <- list.files(file.path(global_TCGA_PCAWG, "TCGA"), full.names = FALSE, recursive = FALSE)
    updateSelectizeInput(session, "selectedSampleTCGA", choices = files, server = TRUE)
  })
  
  observe({
    files <- list.files(file.path(global_TCGA_PCAWG, "PCAWG"), full.names = FALSE, recursive = FALSE)
    updateSelectizeInput(session, "selectedSamplePCAWG", choices = files, server = TRUE)
  })
  
  # Dynamic UI elements for selecting TCGA/PCAWG samples based on cancer type
  output$cancer_type_ui_TCGA <- renderUI({
    selectInput("cancer_type_TCGA", "Select Cancer Type:", choices = unique(sample_TCGA$cancer_type))
  })
  
  output$sample_name_ui_TCGA <- renderUI({
    selectInput("sample_name_TCGA", "Select Sample Name:", choices = NULL)
  })
  
  observeEvent(input$cancer_type_TCGA, {
    selected_cancer_type <- input$cancer_type_TCGA
    sample_names <- sample_TCGA$samplename[sample_TCGA$cancer_type == selected_cancer_type]
    updateSelectInput(session, "sample_name_TCGA", choices = sample_names)
  })
  
  output$cancer_type_ui_PCAWG <- renderUI({
    selectInput("cancer_type_PCAWG", "Select Cancer Type:", choices = unique(sample_PCAWG$cancer_type))
  })
  
  output$sample_name_ui_PCAWG <- renderUI({
    selectInput("sample_name_PCAWG", "Select Sample Name:", choices = NULL)
  })
  
  observeEvent(input$cancer_type_PCAWG, {
    selected_cancer_type <- input$cancer_type_PCAWG
    sample_names <- sample_PCAWG$samplename[sample_PCAWG$cancer_type == selected_cancer_type]
    updateSelectInput(session, "sample_name_PCAWG", choices = sample_names)
  })
  
  # Dynamic UI elements for driver mutation section based on selected cancer type and gene
  activeDriverTab <- reactiveVal("gene")
  
  observeEvent(input$driver_mutation_tabs, {
    activeDriverTab(input$driver_mutation_tabs)
    
    # Update input options based on active tab
    if (input$driver_mutation_tabs == "sample" && !is.null(input$cancer_type_driver)) {
      # When on sample tab, show samples for selected cancer type
      samples <- unique(driver_mutation_data$sample[driver_mutation_data$cancer == input$cancer_type_driver])
      updateSelectInput(session, "gene_name_driver", 
                        label = "Select Sample:",
                        choices = samples)
    } else {
      # When on gene tab, show genes for selected cancer type
      genes <- unique(driver_mutation_data$gene[driver_mutation_data$cancer == input$cancer_type_driver])
      updateSelectInput(session, "gene_name_driver", 
                        label = "Select Gene:",
                        choices = genes)
    }
  })
  
  # Update choices when cancer type changes
  observeEvent(input$cancer_type_driver, {
    if (activeDriverTab() == "sample") {
      samples <- unique(driver_mutation_data$sample[driver_mutation_data$cancer == input$cancer_type_driver])
      updateSelectInput(session, "gene_name_driver", 
                        label = "Select Sample:",
                        choices = samples)
    } else {
      genes <- unique(driver_mutation_data$gene[driver_mutation_data$cancer == input$cancer_type_driver])
      updateSelectInput(session, "gene_name_driver", 
                        label = "Select Gene:",
                        choices = genes)
    }
  })
  
  driver_data <- reactive({
    req(input$cancer_type_driver)
    
    selected_cancer_type <- input$cancer_type_driver
    
    if (is.null(selected_cancer_type) || selected_cancer_type == "") {
      return(NULL)
    }
    
    if (activeDriverTab() == "sample") {
      selected_sample <- input$gene_name_driver  
      selected_data <- driver_mutation_data %>% 
        filter(cancer == selected_cancer_type,
               sample == selected_sample)
    } else {
      selected_data <- driver_mutation_data %>% 
        filter(cancer == selected_cancer_type)
    }
    
    if (nrow(selected_data) == 0) {
      showNotification("No data available for the selected cancer type.", type = "error")
      return(NULL)
    }
    
    selected_data %>% mutate(cancer = selected_cancer_type)
    
    return(selected_data)
  })
  
  observeEvent(input$cancer_type_driver, {
    shinyjs::hide("smoothingControls")
    
    selected_cancer_type <- input$cancer_type_driver
    
    # Get sample for this cancer type (same as TCGA/PCAWG section)
    selected_sample <- sample_TCGA$samplename[sample_TCGA$cancer_type == selected_cancer_type][1]
    if (is.null(selected_sample) || selected_sample == "") {
      selected_sample <- sample_PCAWG$samplename[sample_PCAWG$cancer_type == selected_cancer_type][1]
    }
    
    if (is.null(selected_sample) || selected_sample == "") {
      showNotification("No sample found for selected cancer type.", type = "error")
      return()
    }
    
    # Determine if this is a TCGA or PCAWG sample
    prefix <- if (selected_sample %in% sample_TCGA$samplename) "TCGA" else "PCAWG"
    sampleDir <- file.path(global_TCGA_PCAWG, prefix, selected_sample)
    
    # Read files exactly like in TCGA/PCAWG section
    snvFile <- list.files(sampleDir, pattern = "*.snv.txt", full.names = TRUE)
    cnaFile <- list.files(sampleDir, pattern = "*.cna.txt", full.names = TRUE)
    purityFile <- list.files(sampleDir, pattern = "*.purity.txt", full.names = TRUE)
    
    if (length(snvFile) == 0 || length(cnaFile) == 0 || length(purityFile) == 0) {
      showNotification("Required files not found.", type = "error")
      return()
    }
    
    file_data <- read_and_validate_files(snvFile[1], cnaFile[1], purityFile[1])
    if (is.null(file_data)) return()
    
    df_snv <- file_data$df_snv
    df_cna <- file_data$df_cna
    purity <- file_data$purity
    
    resultPath <- sampleDir
    
    # Read mutation assignments (same as TCGA/PCAWG)
    matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
    matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
    
    if (length(matching_file1) > 0 && length(matching_file2) > 0) {
      df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
      df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
      
      # Process mutations as before
      df_snvcna <- df_snv %>%
        mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
          inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
          if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
        }))
      
      # Calculate metrics
      df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
      mean_rd <- mean(df_snvcna$sum_counts)
      min_rd <- min(df_snvcna$sum_counts)
      max_rd <- max(df_snvcna$sum_counts)
      
      merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
      merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
      
      # Process mutations
      merged_df$total_cn <- as.numeric(merged_df$total_cn)
      merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
      merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
      merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
      
      # Create processed mutations dataframe
      processed_mutations <- merged_df %>%
        mutate(
          Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
          Subclonality = factor(Subclonality, levels = c("Clonal", "Subclonal")),
          ReadDepth = sum_counts,
          IsDriver = FALSE  
        )
      
      # Get driver mutations and add flag
      driver_mutations <- driver_mutation_data %>%
        filter(cancer == selected_cancer_type) %>%
        mutate(IsDriver = TRUE)  
      
      all_mutations <- processed_mutations %>%
        bind_rows(driver_mutations)
      
      subclonality_counts <- processed_mutations %>%
        count(Subclonality) %>%
        rename(n = n)
      
      output$driverPlot1 <- renderPlotly({
        gg1 <- ggplot(all_mutations, aes(x = AF)) +
          geom_density(fill = "gray", alpha = .6, adjust = input$smoothingFactor1_driver) +
          geom_point(
            shape = 142,
            size = 5,
            alpha = 0.3,
            aes(
              y = -0.75,
              color = Subclonality
            )
          ) +
          geom_point(
            data = driver_mutations,
            aes(x = VAF, y = -0.75, fill = Subclonality),
            shape = 21,
            size = 8,
            color = "black",
            alpha = 0.7
          ) +
          geom_text(
            data = driver_mutations,
            aes(x = VAF, y = -1.25, label = gene),
            size = 3,
            angle = 45,
            hjust = 0
          ) +
          xlim(0, 1) +
          ggtitle(paste(
            selected_cancer_type,
            " - Mean read depth =", round(mean_rd),
            "      ", "Range read depth =", min_rd, "-", max_rd
          )) +
          ggpubr::theme_pubr() +
          xlab("VAF") + ylab("Density") +
          scale_color_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
          scale_fill_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
        ggplotly(gg1)
      })
      
      output$driverPlot2 <- renderPlotly({
        gg2 <- ggplot(all_mutations, aes(x = CP_unpenalized)) +
          geom_density(fill = "seagreen", alpha = .6, adjust = input$smoothingFactor2_driver) +
          geom_point(
            shape = 142,
            size = 5,
            alpha = 0.3,
            aes(
              y = -0.75,
              color = Subclonality
            )
          ) +
          geom_point(
            data = driver_mutations,
            aes(x = CP_unpenalized, y = -0.75, fill = Subclonality),
            shape = 21,
            size = 8,
            color = "black",
            alpha = 0.7
          ) +
          # Driver mutation labels
          geom_text(
            data = driver_mutations,
            aes(x = CP_unpenalized, y = -1.25, label = gene),
            size = 3,
            angle = 45,
            hjust = 0
          ) +
          xlim(0, 1.5) +
          ggtitle(paste("CP Distribution -", selected_cancer_type)) +
          ggpubr::theme_pubr() +
          xlab("CP") + ylab("Density") +
          scale_color_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
          scale_fill_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          geom_segment(
            data = df_subc,
            aes(x = cellular_prevalence, xend = cellular_prevalence, y = 0, yend = 5.7),
            lty = "dashed"
          )
        
        ggplotly(gg2)
      })
    }
  })
  
  output$driverMutationPlots <- renderUI({
    if (!is.null(input$cancer_type_driver)) {
      tagList(
        tabsetPanel(
          tabPanel("Plots",
                   uiOutput("smoothingExplanation1_driver"),
                   plotlyOutput("driverPlot1"),
                   uiOutput("driverCaption1"),
                   
                   uiOutput("smoothingExplanation2_driver"),
                   plotlyOutput("driverPlot2"),
                   uiOutput("driverCaption2")
          ),
          tabPanel("Data Table",
                   DT::DTOutput("driverMutationTable"),
                   shinyjs::hidden(downloadButton("downloadDriverMutation", "Download Data")),
                   tags$h4(tags$b(tags$span("Reference:"))),
                   p("Martínez-Jiménez, Francisco, et al. 'A compendium of mutational cancer driver genes.' Nature Reviews Cancer 20.10 (2020): 555-572.")
          )
        )
      )
    }
  })
  
  observe({
    req(input$cancer_type_driver)
    selected_data <- driver_data()
    
    if (is.null(selected_data)) return()
    
    output$driverMutationTable <- DT::renderDataTable({
      DT::datatable(
        selected_data %>% 
          select(gene, sample, protein, consequence, Subclonality, VAF, CP_unpenalized),
        selection = 'single',
        options = list(
          pageLength = 10,
          searchHighlight = TRUE,
          autoWidth = FALSE,
          columnDefs = list(
            list(width = '30px', targets = 0:6)
          )
        ),
        rownames = FALSE
      ) %>%
        DT::formatStyle(
          'Subclonality',
          backgroundColor = DT::styleEqual(
            c("Clonal", "Subclonal"), 
            c('#C9EDEF', '#FED9D7')
          ),
          color = 'black'
        )
    })
    
    observeEvent(input$driverMutationTable_rows_selected, {
      row_selected <- input$driverMutationTable_rows_selected
      if (is.null(row_selected)) {
        shinyjs::hide("smoothingFactor1_driver")
        shinyjs::hide("smoothingFactor2_driver")
        shinyjs::hide("smoothingExplanation1_driver")
        shinyjs::hide("smoothingExplanation2_driver")
        return()
      }
      
      shinyjs::show("smoothingFactor1_driver")
      shinyjs::show("smoothingFactor2_driver")
      shinyjs::show("smoothingExplanation1_driver")
      shinyjs::show("smoothingExplanation2_driver")
      
      selected_data <- driver_data()
      if (is.null(selected_data)) return()
      
      selected_row <- selected_data[row_selected, ]
      selected_cancer_type <- input$cancer_type_driver
      
      selected_sample <- selected_row$sample
      prefix <- if (selected_sample %in% sample_TCGA$samplename) "TCGA" else "PCAWG"
      sampleDir <- file.path(global_TCGA_PCAWG, prefix, selected_sample)
      
      if (!prepareAndRunCliPP(sampleDir, prefix, selected_sample)) {
        showNotification("Failed to prepare data for plotting", type = "error")
        return()
      }
      
      snvFile <- list.files(sampleDir, pattern = "*.snv.txt", full.names = TRUE)
      cnaFile <- list.files(sampleDir, pattern = "*.cna.txt", full.names = TRUE)
      purityFile <- list.files(sampleDir, pattern = "*.purity.txt", full.names = TRUE)
      
      if (length(snvFile) == 0 || length(cnaFile) == 0 || length(purityFile) == 0) return()
      
      file_data <- read_and_validate_files(snvFile[1], cnaFile[1], purityFile[1])
      if (is.null(file_data)) return()
      
      df_snv <- file_data$df_snv
      df_cna <- file_data$df_cna
      purity <- file_data$purity
      
      matching_file1 <- list.files(sampleDir, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(sampleDir, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) > 0 && length(matching_file2) > 0) {
        df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
        df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
        
        df_snvcna <- df_snv %>%
          mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
            inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
            if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
          }))
        
        df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
        mean_rd <- mean(df_snvcna$sum_counts)
        min_rd <- min(df_snvcna$sum_counts)
        max_rd <- max(df_snvcna$sum_counts)
        
        merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
        merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
        
        merged_df$total_cn <- as.numeric(merged_df$total_cn)
        merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
        merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
        merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
        
        merged_df <- merged_df %>%
          mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
                 ReadDepth = sum_counts)
        
        output$driverPlot1 <- renderPlotly({
          gg1 <- ggplot(merged_df, aes(x = AF)) +
            geom_density(fill = "gray", alpha = .6, adjust = input$smoothingFactor1_driver) +
            geom_vline(xintercept = selected_row$VAF, linetype = "dashed", color = "black") +
            geom_point(
              shape = 142,
              size = 5,
              alpha = 0.3,
              aes(
                y = -0.75,
                color = Subclonality
              )
            ) +
            xlim(0, 1) +
            ggtitle(paste(selected_sample, " - Purity =", round(purity, 3),
                          "      ", "Mean read depth =", round(mean_rd),
                          "      ", "Range read depth =", min_rd, "-", max_rd)) +
            ggpubr::theme_pubr() +
            xlab("VAF") + ylab("Density") +
            scale_color_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
          
          ggplotly(gg1)
        })
        
        output$driverPlot2 <- renderPlotly({
          gg2 <- ggplot(merged_df, aes(x = CP_unpenalized)) +
            geom_density(fill = "seagreen", alpha = .6, adjust = input$smoothingFactor2_driver) +
            geom_vline(xintercept = selected_row$CP_unpenalized, linetype = "dashed", color = "black") +
            geom_point(
              shape = 142,
              size = 5,
              alpha = 0.3,
              aes(
                y = -0.75,
                color = Subclonality
              )
            ) +
            xlim(0, 1.5) +
            ggtitle("CP Distribution - Selected Mutation") +
            ggpubr::theme_pubr() +
            xlab("CP") + ylab("Density") +
            scale_color_manual(values = c("Clonal" = "#87CEEB", "Subclonal" = "#FFA07A")) +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            geom_segment(data = df_subc, aes(x = cellular_prevalence, xend = cellular_prevalence, y = 0, yend = 5.7), lty = "dashed")
          
          ggplotly(gg2)
        })
      }
    })
    
    output$driverCaption1 <- renderUI({
      get_vaf_plot_caption()
    })
    
    output$driverCaption2 <- renderUI({
      get_cp_plot_caption()
    })
  })
  
  observe({
    selected_data <- driver_data()
    if (!is.null(selected_data)) {
      shinyjs::show("downloadDriverMutation")
    } else {
      shinyjs::hide("downloadDriverMutation")
    }
  })
  
  observeEvent(input$gene_name_driver, {
    selected_gene_name <- input$gene_name_driver
    sample_names <- unique(driver_mutation_data$sample[driver_mutation_data$gene == selected_gene_name])
    updateSelectInput(session, "sample_name_driver", choices = sample_names)
  })
  
  # Generate driver mutation data table and configure file download for driver mutation data
  observeEvent(input$submitDriver, {
    selected_sample <- input$sample_name_driver
    selected_cancer_type <- input$cancer_type_driver  
    
    if (is.null(selected_sample) || selected_sample == "") {
      showNotification("Please select a valid sample.", type = "error")
      return()
    }
    
    selected_data <- driver_mutation_data %>% filter(sample == selected_sample)
    
    if (nrow(selected_data) == 0) {
      showNotification("No data available for the selected sample.", type = "error")
      return()
    }
    
    selected_data <- selected_data %>%
      mutate(cancer = selected_cancer_type)
    
    output$driverMutationTable <- DT::renderDataTable({
      DT::datatable(selected_data %>% select(gene, protein, consequence, cancer, Subclonality, VAF, CP_unpenalized),
                    options = list(pageLength = 10)) %>%
        DT::formatStyle(
          'Subclonality',
          backgroundColor = DT::styleEqual(c("Clonal", "Subclonal"), c('#C9EDEF', '#FED9D7')),
          color = DT::styleEqual(c("Clonal", "Subclonal"), c('white', 'white'))
        )
    })
    shinyjs::show("downloadDriverMutation")
  })
  
  output$downloadDriverMutation <- downloadHandler(
    filename = function() {
      paste("Driver_Mutation_Data-", input$cancer_type_driver, "-", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      selected_cancer_type <- input$cancer_type_driver
      if (!is.null(selected_cancer_type)) {
        selected_data <- driver_data() %>%
          select(gene, protein, consequence, cancer, Subclonality, VAF, CP_unpenalized)
        
        if (!is.null(selected_data) && nrow(selected_data) > 0) {
          selected_data <- as.data.frame(selected_data)
          
          header <- list(colnames = colnames(selected_data))
          
          rows <- lapply(1:nrow(selected_data), function(i) {
            as.list(selected_data[i,])
          })
          
          output_data <- list(
            header = header,
            rows = rows
          )
          
          write.table(selected_data, file, sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          writeLines("No data available", file)
        }
      }
    },
    contentType = "text/plain"
  )
  
  # Observer to control button visibility
  observe({
    selected_data <- driver_data()
    if (!is.null(selected_data) && nrow(selected_data) > 0) {
      shinyjs::show("downloadDriverMutation")
    } else {
      shinyjs::hide("downloadDriverMutation")
    }
  })
  
  # Add observer to initially hide smoothing controls
  observe({
    if (is.null(input$driverMutationTable_rows_selected)) {
      shinyjs::hide("smoothingFactor1_driver")
      shinyjs::hide("smoothingFactor2_driver")
      shinyjs::hide("smoothingExplanation1_driver")
      shinyjs::hide("smoothingExplanation2_driver")
    }
  })
  
  # Create a reactive expression for the filtered data
  filtered_driver_data_by_sample <- reactive({
    req(input$cancer_type_driver_sample)
    req(input$sample_name_driver_sample)
    
    driver_mutation_data %>%
      filter(
        cancer == input$cancer_type_driver_sample,
        sample == input$sample_name_driver_sample
      ) %>%
      group_by(gene) %>%
      summarise(
        protein = paste(unique(protein), collapse = ", "),
        consequence = paste(unique(consequence), collapse = ", "),
        Subclonality = paste(unique(Subclonality), collapse = ", "),
        VAF = paste(sprintf("%.3f", VAF), collapse = ", "),
        CP_unpenalized = paste(sprintf("%.3f", CP_unpenalized), collapse = ", "),
        .groups = 'drop'
      )
  })
  
  # Update the cancer type UI with available options
  output$cancer_type_ui_driver_sample <- renderUI({
    cancer_types <- sort(unique(driver_mutation_data$cancer))
    selectInput("cancer_type_driver_sample",
                "Select Cancer Type:",
                choices = cancer_types)
  })
  
  # Update the sample name UI based on selected cancer type
  output$sample_name_ui_driver_sample <- renderUI({
    req(input$cancer_type_driver_sample)
    samples <- sort(unique(driver_mutation_data$sample[driver_mutation_data$cancer == input$cancer_type_driver_sample]))
    selectInput("sample_name_driver_sample",
                "Select Sample:",
                choices = samples)
  })
  
  # Render the data table using the reactive filtered data
  output$driverMutationTableBySample <- DT::renderDataTable({
    data <- filtered_driver_data_by_sample()
    
    DT::datatable(
      data,
      selection = 'single',
      options = list(
        pageLength = 10,
        searchHighlight = TRUE,
        autoWidth = FALSE,
        columnDefs = list(
          list(width = '30px', targets = 0:5)
        )
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        'Subclonality',
        backgroundColor = DT::styleEqual(
          c("Clonal", "Subclonal"),
          c('#C9EDEF', '#FED9D7')
        ),
        color = 'black'
      )
  })
  
  # Download handler for the by-sample data
  output$downloadDriverMutationBySample <- downloadHandler(
    filename = function() {
      paste0(
        "Driver_Mutations_",
        input$cancer_type_driver_sample, "_",
        input$sample_name_driver_sample, "_",
        Sys.Date(),
        ".txt"
      )
    },
    content = function(file) {
      data <- filtered_driver_data_by_sample()
      write.table(
        data,
        file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  )
  
  # Observer to handle table row selection and show/hide smoothing controls
  observeEvent(input$driverMutationTableBySample_rows_selected, {
    row_selected <- input$driverMutationTableBySample_rows_selected
    
    if (is.null(row_selected)) {
      shinyjs::hide("smoothingFactor1_driver_sample")
      shinyjs::hide("smoothingFactor2_driver_sample")
      shinyjs::hide("smoothingExplanation1_driver_sample")
      shinyjs::hide("smoothingExplanation2_driver_sample")
      return()
    }
    
    shinyjs::show("smoothingFactor1_driver_sample")
    shinyjs::show("smoothingFactor2_driver_sample")
    shinyjs::show("smoothingExplanation1_driver_sample")
    shinyjs::show("smoothingExplanation2_driver_sample")
  })
}
