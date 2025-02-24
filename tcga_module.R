library(shiny)
library(shinyjs)
library(plotly)
library(dplyr)

#' UI function for TCGA module
#' @param id The module ID
#' @return A Shiny UI definition
tcgaUI <- function(id) {
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

#' Server function for TCGA module
#' @param id The module ID
#' @param sample_data Reactive containing TCGA sample data
#' @param global_path Path to data directory
#' @return A module server function
tcgaServer <- function(id, sample_data, global_path) {
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
      
      withProgress(message = 'Processing TCGA sample...', value = 0, {
        # Get sample directory
        sample_dir <- file.path(global_path, "TCGA", input$sample_name)
        
        # Run CliPP analysis
        success <- prepareAndRunCliPP(sample_dir, "TCGA", input$sample_name)
        
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
        analysis_result <- process_tcga_data(file_data)
        
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
        result_path <- file.path(global_path, "TCGA", input$sample_name)
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
      
      # Clear plots
      output$plot1 <- renderPlotly({ NULL })
      output$plot2 <- renderPlotly({ NULL })
    })
    
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
    
    # Helper function to update result files list
    update_result_files <- function(sample_dir) {
      result_files <- list.files(sample_dir, pattern = "\\.txt$")
      updateSelectInput(session, "selectedFile",
                        choices = result_files)
    }
  })
}