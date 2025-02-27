server <- function(input, output, session) {
  # Initialize reactive values
  plots_ready_clipp <- reactiveVal(FALSE)
  plots_ready_TCGA <- reactiveVal(FALSE)
  plots_ready_PCAWG <- reactiveVal(FALSE)
  
  samples_uploaded <- reactiveVal(FALSE)
  uploaded_samples <- reactiveValues(samples = list())
  
  # Debug logging
  cat("Server function starting\n")
  
  # Null coalescing operator for safe default values
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # Reactive values for sample data
  tcga_data <- reactive({
    sample_TCGA
  })
  
  pcawg_data <- reactive({
    sample_PCAWG
  })
  
  # Call module servers to make them available
  tcga_output <- tcgaServer("tcga", tcga_data, global_TCGA_PCAWG)
  pcawg_output <- pcawgServer("pcawg", pcawg_data, global_TCGA_PCAWG)
  
  # Dynamic file upload inputs based on selected number of samples
  observeEvent(input$setSamples, {
    numSamples <- input$numSamples
    output$dynamicFileInputs <- renderUI({
      lapply(1:numSamples, function(i) {
        tagList(
          h4(paste("Sample", i)),
          textInput(paste0("sampleName", i), "Sample Name:", value = paste0("Sample", i)),
          if (input$uploadMode == "all_at_once") {
            fileInput(paste0("fileInput", i),
                      "Upload all three required files (snv.txt, cna.txt, purity.txt):",
                      multiple = TRUE, 
                      accept = c(".txt"))
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
  
  # Store uploaded files
  observeEvent(input$uploadSamples, {
    numSamples <- input$numSamples
    samples <- store_uploaded_files(input, numSamples, input$uploadMode)
    uploaded_samples$samples <- samples
    samples_uploaded(TRUE)
    
    updateSelectInput(session, "selectedSample",
                      choices = names(samples))
    
    showNotification("Samples uploaded successfully.", type = "message")
    shinyjs::hide("uploadSamples")
    shinyjs::show("selectedSample")
    shinyjs::show("submit")
    shinyjs::show("colorBy")
  })
  
  # Clear current session data
  observeEvent(input$clear, {
    uploaded_samples$samples <- list()
    samples_uploaded(FALSE)
    updateNumericInput(session, "numSamples", value = 1)
    shinyjs::enable("numSamples")
    shinyjs::reset("dynamicFileInputs")
    
    # Hide elements
    elements_to_hide <- c(
      "selectedSample", "submit", "colorBy", "clear",
      "plot1", "plot2", "caption1", "caption2",
      "smoothingFactor1", "smoothingFactor2",
      "smoothingExplanation1", "smoothingExplanation2"
    )
    
    for (element in elements_to_hide) {
      shinyjs::hide(element)
    }
    
    shinyjs::show("setSamples")
    shinyjs::show("uploadSamples")
    
    # Clear plots and outputs
    output$plot1 <- renderPlotly({ NULL })
    output$plot2 <- renderPlotly({ NULL })
    output$caption1 <- renderUI({ NULL })
    output$caption2 <- renderUI({ NULL })
    output$plotOutputArea <- renderUI({ NULL })
    
    plots_ready_clipp(FALSE)
    updateSelectInput(session, "selectedFile", choices = NULL)
  })
  
  # Process uploaded samples
  observeEvent(input$submit, {
    withProgress(message = 'Processing data...', value = 0, {
      plots_ready_clipp(FALSE)
      
      selectedSample <- input$selectedSample
      sampleFiles <- uploaded_samples$samples[[selectedSample]]
      
      if (is.null(sampleFiles)) {
        showNotification("Please select a sample to analyze.", type = "error")
        return()
      }
      
      # Process files
      process_result <- process_sample_files(sampleFiles, selectedSample)
      
      if (!process_result$success) {
        showNotification(process_result$message, type = "error")
        return()
      }
      
      # Update UI elements
      shinyjs::show("clear")
      shinyjs::show("plotOutputArea")
      shinyjs::show("plot1")
      shinyjs::show("plot2")
      shinyjs::show("caption1")
      shinyjs::show("caption2")
      shinyjs::show("controls-area")
      
      # Mark plots as ready
      plots_ready_clipp(TRUE)
    })
  })
  
  # Main plot outputs for uploaded samples
  output$plot1 <- renderPlotly({
    req(plots_ready_clipp())
    
    selected_sample <- input$selectedSample
    current_analysis <- get_current_analysis(uploaded_samples$samples[[selected_sample]])
    
    if (is.null(current_analysis)) return(NULL)
    
    create_vaf_plot(
      merged_df = current_analysis$merged_data,
      sample_name = selected_sample,
      purity = current_analysis$purity,
      read_depth = current_analysis$read_depth,
      smoothing_factor = input$smoothingFactor1,
      color_by = input$colorBy,
      subclonality_counts = current_analysis$subclonality_counts
    )
  })
  
  output$plot2 <- renderPlotly({
    req(plots_ready_clipp())
    
    selected_sample <- input$selectedSample
    current_analysis <- get_current_analysis(uploaded_samples$samples[[selected_sample]])
    
    if (is.null(current_analysis)) return(NULL)
    
    create_cp_plot(
      merged_df = current_analysis$merged_data,
      sample_name = selected_sample,
      smoothing_factor = input$smoothingFactor2,
      color_by = input$colorBy,
      subclonality_counts = current_analysis$subclonality_counts,
      df_subc = current_analysis$subclonal_structure
    )
  })
  
  # Plot captions
  output$caption1 <- renderUI({
    req(plots_ready_clipp())
    get_vaf_plot_caption()
  })
  
  output$caption2 <- renderUI({
    req(plots_ready_clipp())
    get_cp_plot_caption()
  })
  
  # Results download section
  output$plotOutputArea <- renderUI({
    req(plots_ready_clipp())
    
    fluidPage(
      h4("Output:"),
      textOutput("result"),
      selectInput("selectedFile", "Choose a file to download:", choices = c()),
      downloadButton("downloadResult", "Download Result")
    )
  })
  
  # UI for TCGA cancer type selector
  output$cancer_type_ui_TCGA <- renderUI({
    selectInput("cancer_type_TCGA", "Select Cancer Type:",
                choices = sort(unique(sample_TCGA$cancer_type)))
  })
  
  # UI for TCGA sample selector
  output$sample_name_ui_TCGA <- renderUI({
    req(input$cancer_type_TCGA)
    filtered_samples <- sample_TCGA$samplename[sample_TCGA$cancer_type == input$cancer_type_TCGA]
    selectInput("sample_name_TCGA", "Select Sample Name:",
                choices = sort(filtered_samples))
  })
  
  # UI for PCAWG cancer type selector
  output$cancer_type_ui_PCAWG <- renderUI({
    selectInput("cancer_type_PCAWG", "Select Cancer Type:",
                choices = sort(unique(sample_PCAWG$cancer_type)))
  })
  
  # UI for PCAWG sample selector
  output$sample_name_ui_PCAWG <- renderUI({
    req(input$cancer_type_PCAWG)
    filtered_samples <- sample_PCAWG$samplename[sample_PCAWG$cancer_type == input$cancer_type_PCAWG]
    selectInput("sample_name_PCAWG", "Select Sample Name:",
                choices = sort(filtered_samples))
  })
  
  # TCGA Submit button handler
  observeEvent(input$Submit, {
    cat("TCGA Submit button clicked\n")
    plots_ready_TCGA(FALSE)
    
    req(input$cancer_type_TCGA, input$sample_name_TCGA)
    cat("Processing TCGA sample:", input$sample_name_TCGA, "\n")
    
    # Get sample directory
    sample_dir <- file.path(global_TCGA_PCAWG, "TCGA", input$sample_name_TCGA)
    
    # Run CliPP analysis
    success <- prepareAndRunCliPP(sample_dir, "TCGA", input$sample_name_TCGA)
    
    if (!success) {
      showNotification("Failed to process sample", type = "error")
      return()
    }
    
    # Read and validate files
    file_data <- read_tcga_files(sample_dir)
    
    if (is.null(file_data)) {
      showNotification("Error reading sample files", type = "error")
      return()
    }
    
    # Process data
    analysis_result <- process_tcga_data_internal(file_data)
    
    if (is.null(analysis_result)) {
      showNotification("Error processing data", type = "error")
      return()
    }
    
    # Store the analysis result for later use
    tcga_analysis(analysis_result)
    
    # Show controls and enable download
    shinyjs::show("colorByTCGA")
    shinyjs::show("clearTCGA")
    shinyjs::show("plot1TCGA")
    shinyjs::show("plot2TCGA")
    shinyjs::show("caption1TCGA")
    shinyjs::show("caption2TCGA")
    shinyjs::show("smoothingFactor1_TCGA")
    shinyjs::show("smoothingFactor2_TCGA")
    
    # Mark plots as ready
    plots_ready_TCGA(TRUE)
    
    # Create download UI
    output$downloadUI_TCGA <- renderUI({
      result_files <- list.files(sample_dir, pattern = "\\.txt$")
      if (length(result_files) > 0) {
        tagList(
          selectInput("selectedFileTCGA", "Choose a file to download:", choices = result_files),
          downloadButton("downloadResultTCGA", "Download Result")
        )
      }
    })
    
    # Show download section
    shinyjs::show("downloadUI_TCGA")
  })
  
  # PCAWG Submit button handler
  observeEvent(input$Submit_PCAWG, {
    cat("PCAWG Submit button clicked\n")
    plots_ready_PCAWG(FALSE)
    
    req(input$cancer_type_PCAWG, input$sample_name_PCAWG)
    cat("Processing PCAWG sample:", input$sample_name_PCAWG, "\n")
    
    # Get sample directory
    sample_dir <- file.path(global_TCGA_PCAWG, "PCAWG", input$sample_name_PCAWG)
    
    # Check SNV count
    snv_file <- list.files(sample_dir, pattern = "*.snv.txt", full.names = TRUE)[1]
    if (!is.null(snv_file)) {
      snv_count <- count_snvs(snv_file)
      if (snv_count > 20000) {
        showNotification(
          paste("This sample has", snv_count, "SNVs, exceeding the 20,000 limit."),
          type = "error"
        )
        return()
      }
    }
    
    # Run CliPP analysis
    success <- prepareAndRunCliPP(sample_dir, "PCAWG", input$sample_name_PCAWG)
    
    if (!success) {
      showNotification("Failed to process sample", type = "error")
      return()
    }
    
    # Read and validate files
    file_data <- read_pcawg_files(sample_dir)
    
    if (is.null(file_data)) {
      showNotification("Error reading sample files", type = "error")
      return()
    }
    
    # Process data
    analysis_result <- process_pcawg_data_internal(file_data)
    
    if (is.null(analysis_result)) {
      showNotification("Error processing data", type = "error")
      return()
    }
    
    # Store the analysis result for later use
    pcawg_analysis(analysis_result)
    
    # Show controls and enable download
    shinyjs::show("colorByPCAWG")
    shinyjs::show("clearPCAWG")
    shinyjs::show("plot1PCAWG")
    shinyjs::show("plot2PCAWG")
    shinyjs::show("caption1PCAWG")
    shinyjs::show("caption2PCAWG")
    shinyjs::show("smoothingFactor1_PCAWG")
    shinyjs::show("smoothingFactor2_PCAWG")
    
    # Mark plots as ready
    plots_ready_PCAWG(TRUE)
    
    # Create download UI
    output$downloadUI_PCAWG <- renderUI({
      result_files <- list.files(sample_dir, pattern = "\\.txt$")
      if (length(result_files) > 0) {
        tagList(
          selectInput("selectedFilePCAWG", "Choose a file to download:", choices = result_files),
          downloadButton("downloadResultPCAWG", "Download Result")
        )
      }
    })
    
    # Show download section
    shinyjs::show("downloadUI_PCAWG")
  })
  
  # Clear TCGA data
  observeEvent(input$clearTCGA, {
    updateSelectInput(session, "cancer_type_TCGA", selected = character(0))
    updateSelectInput(session, "sample_name_TCGA", selected = character(0))
    
    elements_to_hide <- c(
      "plot1TCGA", "plot2TCGA",
      "caption1TCGA", "caption2TCGA",
      "smoothingFactor1_TCGA", "smoothingFactor2_TCGA",
      "smoothingExplanation1_TCGA", "smoothingExplanation2_TCGA",
      "colorByTCGA", "clearTCGA",
      "downloadUI_TCGA"
    )
    
    for (element in elements_to_hide) {
      shinyjs::hide(element)
    }
    
    output$plot1TCGA <- renderPlotly({ NULL })
    output$plot2TCGA <- renderPlotly({ NULL })
    output$caption1TCGA <- renderUI({ NULL })
    output$caption2TCGA <- renderUI({ NULL })
    
    plots_ready_TCGA(FALSE)
  })
  
  # Clear PCAWG data
  observeEvent(input$clearPCAWG, {
    updateSelectInput(session, "cancer_type_PCAWG", selected = character(0))
    updateSelectInput(session, "sample_name_PCAWG", selected = character(0))
    
    elements_to_hide <- c(
      "plot1PCAWG", "plot2PCAWG",
      "caption1PCAWG", "caption2PCAWG",
      "smoothingFactor1_PCAWG", "smoothingFactor2_PCAWG",
      "smoothingExplanation1_PCAWG", "smoothingExplanation2_PCAWG",
      "colorByPCAWG", "clearPCAWG",
      "downloadUI_PCAWG"
    )
    
    for (element in elements_to_hide) {
      shinyjs::hide(element)
    }
    
    output$plot1PCAWG <- renderPlotly({ NULL })
    output$plot2PCAWG <- renderPlotly({ NULL })
    output$caption1PCAWG <- renderUI({ NULL })
    output$caption2PCAWG <- renderUI({ NULL })
    
    plots_ready_PCAWG(FALSE)
  })
  
  # Storage for analysis results
  tcga_analysis <- reactiveVal(NULL)
  pcawg_analysis <- reactiveVal(NULL)
  
  # TCGA Plot renderers
  output$plot1TCGA <- renderPlotly({
    req(plots_ready_TCGA(), tcga_analysis())
    
    cat("Rendering TCGA VAF plot\n")
    
    analysis <- tcga_analysis()
    
    tryCatch({
      create_vaf_plot(
        merged_df = analysis$merged_data,
        sample_name = input$sample_name_TCGA,
        purity = analysis$purity,
        read_depth = list(
          mean = analysis$read_depth,
          min = min(analysis$merged_data$sum_counts),
          max = max(analysis$merged_data$sum_counts)
        ),
        smoothing_factor = input$smoothingFactor1_TCGA %||% 0.5,
        color_by = input$colorByTCGA %||% "clonality",
        subclonality_counts = analysis$subclonality_counts
      )
    }, error = function(e) {
      cat("Error in TCGA VAF plot:", e$message, "\n")
      plotly::plot_ly() %>% plotly::add_annotations(
        text = paste("Error:", e$message),
        showarrow = FALSE
      )
    })
  })
  
  output$plot2TCGA <- renderPlotly({
    req(plots_ready_TCGA(), tcga_analysis())
    
    cat("Rendering TCGA CP plot\n")
    
    analysis <- tcga_analysis()
    
    tryCatch({
      create_cp_plot(
        merged_df = analysis$merged_data,
        sample_name = input$sample_name_TCGA,
        smoothing_factor = input$smoothingFactor2_TCGA %||% 0.5,
        color_by = input$colorByTCGA %||% "clonality",
        subclonality_counts = analysis$subclonality_counts,
        df_subc = analysis$subclonal_structure
      )
    }, error = function(e) {
      cat("Error in TCGA CP plot:", e$message, "\n")
      plotly::plot_ly() %>% plotly::add_annotations(
        text = paste("Error:", e$message),
        showarrow = FALSE
      )
    })
  })
  
  # TCGA captions and explanations
  output$caption1TCGA <- renderUI({
    req(plots_ready_TCGA())
    get_vaf_plot_caption()
  })
  
  output$caption2TCGA <- renderUI({
    req(plots_ready_TCGA())
    get_cp_plot_caption()
  })
  
  output$smoothingExplanation1_TCGA <- renderUI({
    req(tcga_analysis())
    recommended <- recommended_smoothing_factor(tcga_analysis()$total_mutations)
    HTML(paste0("Recommended smoothing factor: ", recommended))
  })
  
  output$smoothingExplanation2_TCGA <- renderUI({
    req(tcga_analysis())
    recommended <- recommended_smoothing_factor(tcga_analysis()$total_mutations)
    HTML(paste0("Recommended smoothing factor: ", recommended))
  })
  
  # PCAWG Plot renderers
  output$plot1PCAWG <- renderPlotly({
    req(plots_ready_PCAWG(), pcawg_analysis())
    
    cat("Rendering PCAWG VAF plot\n")
    
    analysis <- pcawg_analysis()
    
    tryCatch({
      create_vaf_plot(
        merged_df = analysis$merged_data,
        sample_name = input$sample_name_PCAWG,
        purity = analysis$purity,
        read_depth = list(
          mean = analysis$read_depth,
          min = min(analysis$merged_data$sum_counts),
          max = max(analysis$merged_data$sum_counts)
        ),
        smoothing_factor = input$smoothingFactor1_PCAWG %||% 0.5,
        color_by = input$colorByPCAWG %||% "clonality",
        subclonality_counts = analysis$subclonality_counts
      )
    }, error = function(e) {
      cat("Error in PCAWG VAF plot:", e$message, "\n")
      plotly::plot_ly() %>% plotly::add_annotations(
        text = paste("Error:", e$message),
        showarrow = FALSE
      )
    })
  })
  
  output$plot2PCAWG <- renderPlotly({
    req(plots_ready_PCAWG(), pcawg_analysis())
    
    cat("Rendering PCAWG CP plot\n")
    
    analysis <- pcawg_analysis()
    
    tryCatch({
      create_cp_plot(
        merged_df = analysis$merged_data,
        sample_name = input$sample_name_PCAWG,
        smoothing_factor = input$smoothingFactor2_PCAWG %||% 0.5,
        color_by = input$colorByPCAWG %||% "clonality",
        subclonality_counts = analysis$subclonality_counts,
        df_subc = analysis$subclonal_structure
      )
    }, error = function(e) {
      cat("Error in PCAWG CP plot:", e$message, "\n")
      plotly::plot_ly() %>% plotly::add_annotations(
        text = paste("Error:", e$message),
        showarrow = FALSE
      )
    })
  })
  
  # PCAWG captions and explanations
  output$caption1PCAWG <- renderUI({
    req(plots_ready_PCAWG())
    get_vaf_plot_caption()
  })
  
  output$caption2PCAWG <- renderUI({
    req(plots_ready_PCAWG())
    get_cp_plot_caption()
  })
  
  output$smoothingExplanation1_PCAWG <- renderUI({
    req(pcawg_analysis())
    recommended <- recommended_smoothing_factor(pcawg_analysis()$total_mutations)
    HTML(paste0("Recommended smoothing factor: ", recommended))
  })
  
  output$smoothingExplanation2_PCAWG <- renderUI({
    req(pcawg_analysis())
    recommended <- recommended_smoothing_factor(pcawg_analysis()$total_mutations)
    HTML(paste0("Recommended smoothing factor: ", recommended))
  })
  
  # Download handlers
  output$downloadResultTCGA <- downloadHandler(
    filename = function() {
      selected_file <- input$selectedFileTCGA
      if (grepl(".txt$", selected_file)) {
        selected_file
      } else {
        paste0(selected_file, ".txt")
      }
    },
    content = function(file) {
      req(input$sample_name_TCGA)
      result_path <- file.path(global_TCGA_PCAWG, "TCGA", input$sample_name_TCGA)
      file_path <- file.path(result_path, input$selectedFileTCGA)
      
      if (!file.exists(file_path)) {
        stop("File not found")
      }
      
      file.copy(from = file_path, to = file)
    },
    contentType = "text/plain"
  )
  
  output$downloadResultPCAWG <- downloadHandler(
    filename = function() {
      selected_file <- input$selectedFilePCAWG
      if (grepl(".txt$", selected_file)) {
        selected_file
      } else {
        paste0(selected_file, ".txt")
      }
    },
    content = function(file) {
      req(input$sample_name_PCAWG)
      result_path <- file.path(global_TCGA_PCAWG, "PCAWG", input$sample_name_PCAWG)
      file_path <- file.path(result_path, input$selectedFilePCAWG)
      
      if (!file.exists(file_path)) {
        stop("File not found")
      }
      
      file.copy(from = file_path, to = file)
    },
    contentType = "text/plain"
  )
  
  output$downloadResult <- downloadHandler(
    filename = function() {
      selected_file <- input$selectedFile
      if (grepl(".txt$", selected_file)) {
        selected_file
      } else {
        paste0(selected_file, ".txt")
      }
    },
    content = function(file) {
      selected_sample <- input$selectedSample
      result_path <- file.path("forCliPP/CliPP/sample_id/final_result/Best_lambda")
      file_path <- file.path(result_path, input$selectedFile)
      
      if (!file.exists(file_path)) {
        stop("File not found")
      }
      
      file.copy(from = file_path, to = file)
    },
    contentType = "text/plain"
  )
  
  # Handle navigation changes
  observeEvent(input$`navbar-container`, {
    # Hide smoothing controls for all sections
    smoothing_elements <- c(
      "smoothingFactor1", "smoothingFactor2",
      "smoothingFactor1_TCGA", "smoothingFactor2_TCGA",
      "smoothingFactor1_PCAWG", "smoothingFactor2_PCAWG",
      "smoothingFactor1_driver", "smoothingFactor2_driver",
      "smoothingExplanation1", "smoothingExplanation2",
      "smoothingExplanation1_TCGA", "smoothingExplanation2_TCGA",
      "smoothingExplanation1_PCAWG", "smoothingExplanation2_PCAWG",
      "smoothingExplanation1_driver", "smoothingExplanation2_driver"
    )
    
    for (element in smoothing_elements) {
      shinyjs::hide(element)
    }
  })
  
  # HELPER FUNCTIONS
  
  # Function to count SNVs in a file
  count_snvs <- function(file_path) {
    if (!file.exists(file_path)) return(0)
    df <- read.table(file_path, header = TRUE)
    nrow(df)
  }
  
  # Helper function to read TCGA files
  read_tcga_files <- function(sample_dir) {
    snv_file <- list.files(sample_dir, pattern = "*.snv.txt", full.names = TRUE)
    cna_file <- list.files(sample_dir, pattern = "*.cna.txt", full.names = TRUE)
    purity_file <- list.files(sample_dir, pattern = "*.purity.txt", full.names = TRUE)
    
    if (length(snv_file) == 0 || length(cna_file) == 0 || length(purity_file) == 0) {
      return(NULL)
    }
    
    validate_and_read_files(snv_file[1], cna_file[1], purity_file[1])$data
  }
  
  # Helper function to read PCAWG files
  read_pcawg_files <- function(sample_dir) {
    snv_file <- list.files(sample_dir, pattern = "*.snv.txt", full.names = TRUE)
    cna_file <- list.files(sample_dir, pattern = "*.cna.txt", full.names = TRUE)
    purity_file <- list.files(sample_dir, pattern = "*.purity.txt", full.names = TRUE)
    
    if (length(snv_file) == 0 || length(cna_file) == 0 || length(purity_file) == 0) {
      return(NULL)
    }
    
    validate_and_read_files(snv_file[1], cna_file[1], purity_file[1])$data
  }
  
  # Process TCGA data for plotting
  process_tcga_data_internal <- function(file_data) {
    if (is.null(file_data)) return(NULL)
    
    df_snv <- file_data$df_snv
    df_cna <- file_data$df_cna
    purity <- file_data$purity
    
    tryCatch({
      # Look for existing result files
      result_path <- file.path(global_TCGA_PCAWG, "TCGA", input$sample_name_TCGA)
      matching_file1 <- list.files(result_path, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(result_path, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) == 0 || length(matching_file2) == 0) {
        cat("Required mutation assignment files not found\n")
        return(NULL)
      }
      
      # Read the necessary files
      df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
      df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
      
      # Process the data using the functions from your old server.R
      df_snvcna <- df_snv %>%
        mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
          inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
          if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
        }))
      
      # Calculate read depth
      df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
      mean_rd <- mean(df_snvcna$sum_counts)
      min_rd <- min(df_snvcna$sum_counts)
      max_rd <- max(df_snvcna$sum_counts)
      
      # Merge the data
      merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
      merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
      
      # Calculate parameters
      merged_df$total_cn <- as.numeric(merged_df$total_cn)
      merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
      merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
      merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
      
      # Calculate subclonality
      merged_df <- merged_df %>%
        mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
               ReadDepth = sum_counts)
      
      # Create subclonality counts
      subclonality_counts <- merged_df %>%
        count(Subclonality) %>%
        rename(n = n)
      
      # Total mutations count for smoothing recommendations
      total_SNVs <- sum(df_subc$num_SNV)
      
      # Return analysis result
      return(list(
        merged_data = merged_df,
        purity = purity,
        read_depth = mean_rd,
        subclonality_counts = subclonality_counts,
        subclonal_structure = df_subc,
        total_mutations = total_SNVs
      ))
    }, error = function(e) {
      cat("Error processing TCGA data:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Process PCAWG data for plotting
  process_pcawg_data_internal <- function(file_data) {
    if (is.null(file_data)) return(NULL)
    
    df_snv <- file_data$df_snv
    df_cna <- file_data$df_cna
    purity <- file_data$purity
    
    tryCatch({
      # Look for existing result files
      result_path <- file.path(global_TCGA_PCAWG, "PCAWG", input$sample_name_PCAWG)
      matching_file1 <- list.files(result_path, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(result_path, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) == 0 || length(matching_file2) == 0) {
        cat("Required mutation assignment files not found\n")
        return(NULL)
      }
      
      # Read the necessary files
      df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
      df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
      
      # Process the data using the functions from your old server.R
      df_snvcna <- df_snv %>%
        mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
          inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
          if (any(inds)) as.character(df_cna$total_cn[which.max(inds)]) else NA
        }))
      
      # Calculate read depth
      df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
      mean_rd <- mean(df_snvcna$sum_counts)
      min_rd <- min(df_snvcna$sum_counts)
      max_rd <- max(df_snvcna$sum_counts)
      
      # Merge the data
      merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
      merged_df['AF'] = merged_df['alt_count'] / (merged_df['sum_counts'])
      
      # Calculate parameters
      merged_df$total_cn <- as.numeric(merged_df$total_cn)
      merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
      merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
      merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
      
      # Calculate subclonality
      merged_df <- merged_df %>%
        mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"),
               ReadDepth = sum_counts)
      
      # Create subclonality counts
      subclonality_counts <- merged_df %>%
        count(Subclonality) %>%
        rename(n = n)
      
      # Total mutations count for smoothing recommendations
      total_SNVs <- sum(df_subc$num_SNV)
      
      # Return analysis result
      return(list(
        merged_data = merged_df,
        purity = purity,
        read_depth = mean_rd,
        subclonality_counts = subclonality_counts,
        subclonal_structure = df_subc,
        total_mutations = total_SNVs
      ))
    }, error = function(e) {
      cat("Error processing PCAWG data:", e$message, "\n")
      return(NULL)
    })
  }
}
