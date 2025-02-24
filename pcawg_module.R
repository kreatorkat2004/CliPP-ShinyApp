library(shiny)
library(shinyjs)
library(plotly)
library(dplyr)

#' UI function for PCAWG module
#' @param id The module ID
#' @return A Shiny UI definition
pcawgUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    fluidRow(
      column(12,
             wellPanel(
               fluidRow(
                 column(6,
                        selectInput(
                          ns("cancer_type"),
                          "Select Cancer Type:",
                          choices = NULL
                        )
                 ),
                 column(6,
                        selectInput(
                          ns("sample_name"),
                          "Select Sample Name:",
                          choices = NULL
                        )
                 )
               ),
               fluidRow(
                 column(4,
                        actionButton(
                          ns("submit"),
                          "Analyze Sample",
                          class = "btn-primary"
                        )
                 ),
                 column(4,
                        hidden(
                          selectInput(
                            ns("colorBy"),
                            "Color points by:",
                            choices = c(
                              "Clonality" = "clonality",
                              "Read Depth" = "depth"
                            )
                          )
                        )
                 ),
                 column(4,
                        actionButton(
                          ns("clear"),
                          "Clear",
                          class = "btn-warning"
                        )
                 )
               )
             )
      )
    ),
    
    # Info box for SNV limit warning
    hidden(
      div(id = ns("snvWarning"),
          wellPanel(
            class = "alert alert-warning",
            "Note: PCAWG samples are limited to 20,000 SNVs for analysis."
          )
      )
    ),
    
    # Smoothing factor controls
    hidden(
      div(id = ns("smoothingControls"),
          fluidRow(
            column(6,
                   sliderInput(
                     ns("smoothingFactor1"),
                     "Adjust VAF plot smoothing:",
                     min = 0.1,
                     max = 1,
                     value = 0.5,
                     step = 0.1
                   ),
                   uiOutput(ns("smoothingExplanation1"))
            ),
            column(6,
                   sliderInput(
                     ns("smoothingFactor2"),
                     "Adjust CP plot smoothing:",
                     min = 0.1,
                     max = 1,
                     value = 0.5,
                     step = 0.1
                   ),
                   uiOutput(ns("smoothingExplanation2"))
            )
          )
      )
    ),
    
    # Plot outputs
    fluidRow(
      column(12,
             # VAF Plot
             plotlyOutput(ns("plot1")),
             uiOutput(ns("caption1")),
             
             # CP Plot
             plotlyOutput(ns("plot2")),
             uiOutput(ns("caption2"))
      )
    ),
    
    # Results download section
    hidden(
      div(id = ns("downloadSection"),
          wellPanel(
            h4("Analysis Results"),
            selectInput(
              ns("selectedFile"),
              "Choose file to download:",
              choices = NULL
            ),
            downloadButton(
              ns("downloadResult"),
              "Download"
            )
          )
      )
    )
  )
}

#' Server function for PCAWG module
#' @param id The module ID
#' @param sample_data Reactive containing PCAWG sample data
#' @param global_path Path to data directory
#' @return A module server function
pcawgServer <- function(id, sample_data, global_path) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values
    plots_ready <- reactiveVal(FALSE)
    current_analysis <- reactiveVal(NULL)
    
    # Initialize cancer type choices
    observe({
      req(sample_data())
      cancer_types <- unique(sample_data()$cancer_type)
      updateSelectInput(session, "cancer_type",
                        choices = sort(cancer_types))
    })
    
    # Update sample names when cancer type changes
    observeEvent(input$cancer_type, {
      req(input$cancer_type, sample_data())
      sample_names <- sample_data() %>%
        filter(cancer_type == input$cancer_type) %>%
        pull(samplename)
      updateSelectInput(session, "sample_name",
                        choices = sort(sample_names))
    })
    
    # Handle analysis submission
    observeEvent(input$submit, {
      req(input$cancer_type, input$sample_name)
      
      # Reset plots
      plots_ready(FALSE)
      
      withProgress(message = 'Processing PCAWG sample...', value = 0, {
        # Get sample directory
        sample_dir <- file.path(global_path, "PCAWG", input$sample_name)
        
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
        success <- prepareAndRunCliPP(sample_dir, "PCAWG", input$sample_name)
        
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
        analysis_result <- process_pcawg_data(file_data)
        
        if (is.null(analysis_result)) {
          showNotification("Error processing data", type = "error")
          return()
        }
        
        # Store current analysis
        current_analysis(analysis_result)
        
        # Show controls and enable download
        shinyjs::show("smoothingControls")
        shinyjs::show("colorBy")
        shinyjs::show("downloadSection")
        
        # Update available result files
        update_result_files(sample_dir)
        
        # Mark plots as ready
        plots_ready(TRUE)
      })
    })
    
    # Create VAF plot
    output$plot1 <- renderPlotly({
      req(plots_ready(), current_analysis())
      
      analysis <- current_analysis()
      
      create_vaf_plot(
        merged_df = analysis$merged_data,
        sample_name = input$sample_name,
        purity = analysis$purity,
        read_depth = analysis$read_depth,
        smoothing_factor = input$smoothingFactor1,
        color_by = input$colorBy,
        subclonality_counts = analysis$subclonality_counts
      )
    })
    
    # Create CP plot
    output$plot2 <- renderPlotly({
      req(plots_ready(), current_analysis())
      
      analysis <- current_analysis()
      
      create_cp_plot(
        merged_df = analysis$merged_data,
        sample_name = input$sample_name,
        smoothing_factor = input$smoothingFactor2,
        color_by = input$colorBy,
        subclonality_counts = analysis$subclonality_counts,
        df_subc = analysis$subclonal_structure
      )
    })
    
    # Plot captions
    output$caption1 <- renderUI({
      req(plots_ready())
      get_vaf_plot_caption()
    })
    
    output$caption2 <- renderUI({
      req(plots_ready())
      get_cp_plot_caption()
    })
    
    # Smoothing factor explanations
    output$smoothingExplanation1 <- renderUI({
      req(current_analysis())
      analysis <- current_analysis()
      recommended <- recommended_smoothing_factor(analysis$total_mutations)
      HTML(paste0("Recommended smoothing factor: ", recommended))
    })
    
    output$smoothingExplanation2 <- renderUI({
      req(current_analysis())
      analysis <- current_analysis()
      recommended <- recommended_smoothing_factor(analysis$total_mutations)
      HTML(paste0("Recommended smoothing factor: ", recommended))
    })
    
    # Handle download
    output$downloadResult <- downloadHandler(
      filename = function() {
        if (grepl(".txt$", input$selectedFile)) {
          input$selectedFile
        } else {
          paste0(input$selectedFile, ".txt")
        }
      },
      content = function(file) {
        req(input$sample_name)
        result_path <- file.path(global_path, "PCAWG", input$sample_name)
        file_path <- file.path(result_path, input$selectedFile)
        
        if (!file.exists(file_path)) {
          stop("File not found")
        }
        
        file.copy(from = file_path, to = file)
      },
      contentType = "text/plain"
    )
    
    # Clear analysis
    observeEvent(input$clear, {
      plots_ready(FALSE)
      current_analysis(NULL)
      
      # Reset inputs
      updateSelectInput(session, "cancer_type", selected = character(0))
      updateSelectInput(session, "sample_name", selected = character(0))
      
      # Hide elements
      shinyjs::hide("smoothingControls")
      shinyjs::hide("colorBy")
      shinyjs::hide("downloadSection")
      shinyjs::hide("snvWarning")
      
      # Clear plots
      output$plot1 <- renderPlotly({ NULL })
      output$plot2 <- renderPlotly({ NULL })
    })
    
    # Helper Functions
    
    #' Count SNVs in a file
    #' @param file_path Path to SNV file
    #' @return Number of SNVs
    count_snvs <- function(file_path) {
      if (!file.exists(file_path)) return(0)
      df <- read.table(file_path, header = TRUE)
      nrow(df)
    }
    
    #' Read PCAWG files
    #' @param sample_dir Sample directory path
    #' @return List of data frames or NULL
    read_pcawg_files <- function(sample_dir) {
      snv_file <- list.files(sample_dir, pattern = "*.snv.txt", full.names = TRUE)
      cna_file <- list.files(sample_dir, pattern = "*.cna.txt", full.names = TRUE)
      purity_file <- list.files(sample_dir, pattern = "*.purity.txt", full.names = TRUE)
      
      if (length(snv_file) == 0 || length(cna_file) == 0 || length(purity_file) == 0) {
        return(NULL)
      }
      
      validate_and_read_files(snv_file[1], cna_file[1], purity_file[1])$data
    }
    
    #' Process PCAWG data
    #' @param file_data List of input data frames
    #' @return Processed analysis results or NULL
    process_pcawg_data <- function(file_data) {
      if (is.null(file_data)) return(NULL)
      
      tryCatch({
        # Process data using shared data processing function
        analysis_result <- process_mutation_data(
          df_snv = file_data$df_snv,
          df_cna = file_data$df_cna,
          purity = file_data$purity
        )
        
        return(analysis_result)
      }, error = function(e) {
        showNotification(paste("Error processing data:", e$message), type = "error")
        return(NULL)
      })
    }
    
    #' Update result files list
    #' @param sample_dir Sample directory path
    update_result_files <- function(sample_dir) {
      result_files <- list.files(sample_dir, pattern = "\\.txt$")
      updateSelectInput(session, "selectedFile",
                        choices = result_files)
    }
  })
}