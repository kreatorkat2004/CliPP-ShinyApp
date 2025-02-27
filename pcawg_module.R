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
    
    # Simply maintain minimal functionality to not break dependencies
    cat("PCAWG module initialized\n")
  })
}
