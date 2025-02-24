server <- function(input, output, session) {
  # Initialize reactive values for file upload section
  plots_ready_clipp <- reactiveVal(FALSE)
  samples_uploaded <- reactiveVal(FALSE)
  uploaded_samples <- reactiveValues(samples = list())
  
  # Reactive values for sample data
  tcga_data <- reactive({
    sample_TCGA
  })
  
  pcawg_data <- reactive({
    sample_PCAWG
  })
  
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
  
  # Call module servers
  tcgaServer("tcga", tcga_data, global_TCGA_PCAWG)
  pcawgServer("pcawg", pcawg_data, global_TCGA_PCAWG)
  driverMutationServer("driver", reactive(driver_mutation_data), global_TCGA_PCAWG)
  
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
  
  # Download handler
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
}
