library(shiny)
library(shinyjs)
library(DT)
library(plotly)
library(dplyr)

#' UI function for driver mutation module
#' @param id The module ID
#' @return A Shiny UI definition
driverMutationUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    useShinyjs(),
    
    # Analysis mode tabs
    tabsetPanel(
      id = ns("driver_mutation_tabs"),
      
      # Gene Analysis Tab
      tabPanel(
        "Gene Analysis",
        value = "gene",
        wellPanel(
          fluidRow(
            column(6,
                   selectInput(
                     ns("cancer_type_gene"),
                     "Select Cancer Type:",
                     choices = NULL
                   )
            ),
            column(6,
                   selectInput(
                     ns("gene_name"),
                     "Select Gene:",
                     choices = NULL
                   )
            )
          )
        ),
        
        # Gene analysis results
        uiOutput(ns("gene_analysis_results"))
      ),
      
      # Sample Analysis Tab
      tabPanel(
        "Sample Analysis",
        value = "sample",
        wellPanel(
          fluidRow(
            column(6,
                   selectInput(
                     ns("cancer_type_sample"),
                     "Select Cancer Type:",
                     choices = NULL
                   )
            ),
            column(6,
                   selectInput(
                     ns("sample_name"),
                     "Select Sample:",
                     choices = NULL
                   )
            )
          )
        ),
        
        # Sample analysis results
        uiOutput(ns("sample_analysis_results"))
      )
    ),
    
    # Smoothing controls (hidden by default)
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
    
    # Plots and data table
    uiOutput(ns("analysis_output"))
  )
}

#' Server function for driver mutation module
#' @param id The module ID
#' @param driver_data Reactive containing driver mutation data
#' @param global_path Path to data directory
#' @return A module server function
driverMutationServer <- function(id, driver_data, global_path) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values
    active_tab <- reactiveVal("gene")
    plots_ready <- reactiveVal(FALSE)
    current_analysis <- reactiveVal(NULL)
    
    # Track active tab
    observeEvent(input$driver_mutation_tabs, {
      active_tab(input$driver_mutation_tabs)
    })
    
    # Initialize cancer type choices for both tabs
    observe({
      req(driver_data())
      cancer_types <- sort(unique(driver_data()$cancer))
      
      updateSelectInput(session, "cancer_type_gene",
                        choices = cancer_types)
      updateSelectInput(session, "cancer_type_sample",
                        choices = cancer_types)
    })
    
    # Update gene choices based on cancer type
    observeEvent(input$cancer_type_gene, {
      req(driver_data(), input$cancer_type_gene)
      filtered_data <- driver_data()
      genes <- filtered_data[filtered_data$cancer == input$cancer_type_gene, "gene"]
      genes <- sort(unique(genes))
      
      updateSelectInput(session, "gene_name",
                        choices = genes)
    })
    
    # Update sample choices based on cancer type
    observeEvent(input$cancer_type_sample, {
      req(driver_data(), input$cancer_type_sample)
      filtered_data <- driver_data()
      samples <- filtered_data[filtered_data$cancer == input$cancer_type_sample, "sample"]
      samples <- sort(unique(samples))
      
      updateSelectInput(session, "sample_name",
                        choices = samples)
    })
    
    # Process data based on active tab
    process_selection <- reactive({
      req(driver_data())
      
      cancer_type <- if(active_tab() == "gene") input$cancer_type_gene else input$cancer_type_sample
      selection <- if(active_tab() == "gene") input$gene_name else input$sample_name
      
      if (is.null(cancer_type) || is.null(selection)) return(NULL)
      
      # Filter data based on mode
      filtered_data <- if(active_tab() == "gene") {
        subset(driver_data(), 
               cancer == cancer_type & gene == selection)
      } else {
        subset(driver_data(),
               cancer == cancer_type & sample == selection)
      }
      
      if (nrow(filtered_data) == 0) return(NULL)
      
      # Process the data
      result <- process_driver_mutation_data(filtered_data, active_tab())
      
      if (is.null(result)) return(NULL)
      
      plots_ready(TRUE)
      current_analysis(result)
      
      result
    })
    
    # Create analysis output UI
    output$analysis_output <- renderUI({
      req(plots_ready(), current_analysis())
      
      tagList(
        # Plots
        fluidRow(
          column(12,
                 plotlyOutput(ns("plot1")),
                 uiOutput(ns("caption1")),
                 plotlyOutput(ns("plot2")),
                 uiOutput(ns("caption2"))
          )
        ),
        
        # Data table
        fluidRow(
          column(12,
                 wellPanel(
                   h4("Driver Mutations"),
                   DTOutput(ns("mutation_table")),
                   downloadButton(ns("download_data"), "Download Data"),
                   tags$br(),
                   tags$br(),
                   tags$h4(tags$b("Reference:")),
                   tags$p("Martínez-Jiménez, Francisco, et al. 'A compendium of mutational cancer driver genes.' Nature Reviews Cancer 20.10 (2020): 555-572.")
                 )
          )
        )
      )
    })
    
    # Create plots
    output$plot1 <- renderPlotly({
      req(plots_ready(), current_analysis())
      analysis <- current_analysis()
      
      create_driver_vaf_plot(
        merged_df = analysis$merged_data,
        driver_mutations = analysis$driver_mutations,
        sample_name = analysis$sample_name,
        purity = analysis$purity,
        read_depth = analysis$read_depth,
        smoothing_factor = input$smoothingFactor1
      )
    })
    
    output$plot2 <- renderPlotly({
      req(plots_ready(), current_analysis())
      analysis <- current_analysis()
      
      create_driver_cp_plot(
        merged_df = analysis$merged_data,
        driver_mutations = analysis$driver_mutations,
        sample_name = analysis$sample_name,
        smoothing_factor = input$smoothingFactor2,
        df_subc = analysis$subclonal_structure
      )
    })
    
    # Create data table
    output$mutation_table <- renderDT({
      req(current_analysis())
      analysis <- current_analysis()
      
      datatable(
        analysis$driver_mutations %>%
          select(gene, sample, protein, consequence, Subclonality, VAF, CP_unpenalized),
        selection = 'single',
        options = list(
          pageLength = 10,
          searchHighlight = TRUE,
          autoWidth = FALSE
        ),
        rownames = FALSE
      ) %>%
        formatStyle(
          'Subclonality',
          backgroundColor = styleEqual(
            c("Clonal", "Subclonal"),
            c('#C9EDEF', '#FED9D7')
          ),
          color = 'black'
        )
    })
    
    # Download handler
    output$download_data <- downloadHandler(
      filename = function() {
        cancer_type <- if(active_tab() == "gene") input$cancer_type_gene else input$cancer_type_sample
        selection <- if(active_tab() == "gene") input$gene_name else input$sample_name
        paste0("Driver_Mutations_", cancer_type, "_", selection, "_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(current_analysis())
        analysis <- current_analysis()
        
        write.table(
          analysis$driver_mutations %>%
            select(gene, sample, protein, consequence, Subclonality, VAF, CP_unpenalized),
          file,
          sep = "\t",
          row.names = FALSE,
          quote = FALSE
        )
      }
    )
    
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
    
    # Observer for table row selection
    observeEvent(input$mutation_table_rows_selected, {
      row_selected <- input$mutation_table_rows_selected
      if (is.null(row_selected)) {
        shinyjs::hide("smoothingControls")
      } else {
        shinyjs::show("smoothingControls")
      }
    })
  })
}

#' Process driver mutation data
#' @param filtered_data Filtered driver mutation data
#' @param mode Analysis mode ("gene" or "sample")
#' @return List containing processed data and analysis results
process_driver_mutation_data <- function(filtered_data, mode) {
  if (nrow(filtered_data) == 0) return(NULL)
  
  # Get sample information
  sample_name <- if(mode == "gene") filtered_data$sample[1] else filtered_data$sample
  
  # Determine if sample is from TCGA or PCAWG
  prefix <- if(grepl("TCGA", sample_name)) "TCGA" else "PCAWG"
  sample_dir <- file.path(global_path, prefix, sample_name)
  
  # Read and process sample files
  file_data <- read_sample_files(sample_dir)
  if (is.null(file_data)) return(NULL)
  
  # Process mutation data
  analysis_result <- process_mutation_data(
    file_data$df_snv,
    file_data$df_cna,
    file_data$purity
  )
  
  if (is.null(analysis_result)) return(NULL)
  
  # Add driver mutation information
  analysis_result$driver_mutations <- filtered_data
  analysis_result$sample_name <- sample_name
  
  return(analysis_result)
}

#' Create driver mutation VAF plot
#' @param merged_df Processed mutation data
#' @param driver_mutations Driver mutation data
#' @param sample_name Sample name
#' @param purity Purity value
#' @param read_depth Read depth statistics
#' @param smoothing_factor Smoothing factor for density plot
#' @return plotly object
create_driver_vaf_plot <- function(merged_df, driver_mutations, sample_name,
                                   purity, read_depth, smoothing_factor) {
  create_vaf_plot(
    merged_df = merged_df,
    sample_name = sample_name,
    purity = purity,
    read_depth = read_depth,
    smoothing_factor = smoothing_factor,
    color_by = "clonality",
    subclonality_counts = table(merged_df$Subclonality),
    driver_mutations = driver_mutations
  )
}

#' Create driver mutation CP plot
#' @param merged_df Processed mutation data
#' @param driver_mutations Driver mutation data
#' @param sample_name Sample name
#' @param smoothing_factor Smoothing factor for density plot
#' @param df_subc Subclonal structure data
#' @return plotly object
create_driver_cp_plot <- function(merged_df, driver_mutations, sample_name,
                                  smoothing_factor, df_subc) {
  create_cp_plot(
    merged_df = merged_df,
    sample_name = sample_name,
    smoothing_factor = smoothing_factor,
    color_by = "clonality",
    subclonality_counts = table(merged_df$Subclonality),
    df_subc = df_subc,
    driver_mutations = driver_mutations
  )
}