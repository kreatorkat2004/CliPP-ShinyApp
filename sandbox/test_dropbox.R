library(shiny)

# Load the sample data
sample_PCAWG <- read.csv("CliPP_PCAWG.tsv", sep='\t')
sample_TCGA <- read.csv("CliPP_TCGA.tsv", sep='\t')

# Add dataset column to distinguish between PCAWG and TCGA
sample_PCAWG$dataset <- "PCAWG"
sample_TCGA$dataset <- "TCGA"

sample_data <- rbind(sample_PCAWG, sample_TCGA)

head(sample_data)

# UI
ui <- fluidPage(
  titlePanel("Test Cancer Type Selector"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select Dataset:", choices = c("PCAWG", "TCGA")),
      uiOutput("cancer_type_ui"),
      uiOutput("sample_name_ui")
    ),
    mainPanel(
      textOutput("selected_sample")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Generate the cancer type UI dynamically based on the selected dataset
  output$cancer_type_ui <- renderUI({
    if (is.null(input$dataset)) return(NULL)
    selectInput("cancer_type", "Select Cancer Type:", choices = NULL)
  })
  
  # Generate the sample name UI dynamically based on the selected cancer type
  output$sample_name_ui <- renderUI({
    if (is.null(input$cancer_type)) return(NULL)
    selectInput("sample_name", "Select Sample Name:", choices = NULL)
  })
  
  # Update cancer type choices based on the selected dataset
  observeEvent(input$dataset, {
    selected_dataset <- input$dataset
    cancer_types <- unique(sample_data$cancer_type[sample_data$dataset == selected_dataset])
    updateSelectInput(session, "cancer_type", choices = cancer_types, selected = NULL) # Reset cancer type selection
  })
  
  # Update sample name choices based on the selected cancer type
  observeEvent(input$cancer_type, {
    selected_dataset <- input$dataset
    selected_cancer_type <- input$cancer_type
    sample_names <- sample_data$samplename[sample_data$dataset == selected_dataset & sample_data$cancer_type == selected_cancer_type]
    updateSelectInput(session, "sample_name", choices = sample_names)
  })
  
  # Display the selected sample name
  output$selected_sample <- renderText({
    paste("Selected Sample Name:", input$sample_name)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
