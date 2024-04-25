library(reticulate)

server <- function(input, output, session) {
  if (Sys.info()[['user']] == 'shiny'){
    Sys.setenv(RETICULATE_PYTHON = "./python3")
  }
  
  outputDir <- getwd()  # Define the output directory globally
  resultPath <- file.path(outputDir, "forCliPP/CliPP/sample_id/final_result/Best_lambda")
  
  observeEvent(input$submit, {
    req(input$file1, input$file2, input$file3)  # Ensure all files are uploaded

    # Define the directory path
    directoryPath <- "forCliPP/CliPP/sample_id/"
    
    # Check if the directory exists, if not, create it
    if (!dir.exists(directoryPath)) {
      dir.create(directoryPath, recursive = TRUE, showWarnings = TRUE)
    } else {
      # Clear previous files in the directory
      system(paste("rm -rf", shQuote(directoryPath)), intern = TRUE)
    }
        
    # Save the uploaded files temporarily
    inFile1 <- input$file1$datapath
    inFile2 <- input$file2$datapath
    inFile3 <- input$file3$datapath
    
    # Construct the command
    command <- sprintf("./python3 ./forCliPP/CliPP/run_clipp_main.py %s %s %s",
                       shQuote(inFile1), shQuote(inFile2), shQuote(inFile3))
    
    # Execute the command
    system(command, intern = TRUE)
    
    # Update the file selection dropdown
    updateSelectInput(session, "selectedFile", choices = list.files(resultPath, pattern = "\\.txt$"))
    
    # Output message
    output$result <- renderText({
      "Processing complete. Select a file to download."
    })
  })
  
  # Download handler
  output$downloadResult <- downloadHandler(
    filename = function() {
      file_path <- file.path(resultPath, input$selectedFile)
      validate(need(file.exists(file_path), "File not found."))
      basename(file_path)
    },
    content = function(file) {
      file.copy(from = file.path(resultPath, input$selectedFile), to = file)
    }
  )
}
