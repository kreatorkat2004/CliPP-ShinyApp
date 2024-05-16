library(reticulate)
library(ggplot2)
library(dplyr)
library(purrr)
library(plotly)
library(readr)

# Note: must install pandas in order to run this code

server <- function(input, output, session) {
  plotsReady <- reactiveVal(FALSE)
  samplesUploaded <- reactiveVal(FALSE)
  # Define reactive values for the uploaded sample files
  uploadedSamples <- reactiveValues(samples = list())
  
  if (Sys.info()[['user']] == 'shiny') {
    # Sys.setenv(RETICULATE_PYTHON = "./python3")
    Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
  }
  
  outputDir <- getwd() # Define the output directory globally
  resultPath <- file.path(outputDir, "forCliPP/CliPP/sample_id/final_result/Best_lambda")
  
  # Generate input samples dynamically based on input amount
  observeEvent(input$setSamples, {
    numSamples <- input$numSamples
    output$dynamicFileInputs <- renderUI({
      lapply(1:numSamples, function(i) {
        tagList(
          h4(paste("Sample", i)),
          textInput(paste0("sampleName", i), "Sample Name:", value = paste0("Sample", i)),
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
  })
  
  # Helper function to store uploaded files
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
  
  # Store uploaded samples
  observeEvent(input$uploadSamples, {
    numSamples <- input$numSamples
    samples <- storeUploadedFiles(input, numSamples)
    uploadedSamples$samples <- samples
    samplesUploaded(TRUE)
    updateSelectInput(session, "selectedSample", choices = names(samples))
    showNotification("Samples uploaded successfully.", type = "message")
  })
  
  # This event allows for the plots and output to only be displayed after 
  # the entire code has been run
  observeEvent(input$submit, {
    plotsReady(FALSE)
    
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
    
    # Define the directory path
    directoryPath <- "forCliPP/CliPP/sample_id/"
    
    # Check if the directory exists, if not, create it
    if (!dir.exists(directoryPath)) {
      dir.create(directoryPath, recursive = TRUE, showWarnings = TRUE)
    } else {
      # Clear previous files in the directory
      system(paste("rm -rf", shQuote(directoryPath)), intern = TRUE)
    }
    
    # Construct the command
    command <- sprintf("/usr/local/bin/python3 ./forCliPP/CliPP/run_clipp_main.py %s %s %s",
                       shQuote(snvFile), shQuote(cnaFile), shQuote(purityFile))
    
    # Execute the command
    system(command, intern = TRUE)
    
    # Update the file selection dropdown
    updateSelectInput(session, "selectedFile", choices = list.files(resultPath, pattern = "\\.txt$"))
    
    # Output message
    output$result <- renderText({
      "Processing complete. Select a file to download."
    })
    plotsReady(TRUE)
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
  
  output$warning <- renderText({
    if (is.null(input$submit)) return()
  })
  
  output$plotReady <- reactive({ plotsReady() })
  outputOptions(output, "plotReady", suspendWhenHidden = FALSE)
  
  # Dummy text output to trigger the conditional panel
  output$plotReadyText <- renderUI({
    if (plotsReady()) {
      div(id = "plotReadyText", "true")
    } else {
      div(id = "plotReadyText", "false")
    }
  })
  
  # Observe plot readiness and render plots
  observe({
    req(plotsReady())
    output$plot1 <- renderPlotly({
      
      if (input$submit == 0 || !plotsReady()) return()  # Ensure files are uploaded
      
      selectedSample <- input$selectedSample
      sampleFiles <- uploadedSamples$samples[[selectedSample]]
      
      if (is.null(sampleFiles)) return()
      
      snvFile <- sampleFiles$snvFile
      cnaFile <- sampleFiles$cnaFile
      purityFile <- sampleFiles$purityFile
      
      if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) return()
      
      # Reading data
      df_snv <- read.table(snvFile, header = TRUE, sep = "\t")
      df_cna <- read.table(cnaFile, header = TRUE, sep = "\t")
      purity <- as.numeric(readLines(purityFile))
      
      # Read mutation and subclonal structure assignments
      matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) > 0 && length(matching_file2) > 0) {
        df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
        df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
        
        df_snvcna <- df_snv %>%
          mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
            inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
            if (any(inds)) df_cna$total_cn[which.max(inds)] else NA
          }))
        
        # Calculate the sum of ref_count and alt_count
        df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
        
        # Compute the mean of the sum_counts column
        mean_rd <- mean(df_snvcna$sum_counts)
        min_rd <- min(df_snvcna$sum_counts)
        max_rd <- max(df_snvcna$sum_counts)
        
        merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
        merged_df['AF'] = merged_df['alt_count']/(merged_df['sum_counts'])
        
        merged_df$total_cn <- as.numeric(merged_df$total_cn)
        merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
        merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
        merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
        
        total_SNVs <- sum(df_subc$num_SNV)
        # smoothing is adjustable with slider
        smoothing_factor1 <- input$smoothingFactor1
        
        merged_df <- merged_df %>%
          mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"))
        
        print(head(merged_df))
        
        subclonality_counts <- merged_df %>%
          count(Subclonality)
        
        # Rename the count column to 'n' as specified in your desired data frame
        subclonality_counts <- subclonality_counts %>%
          rename(n = n)
        
        # Create plot
        gg1 <- ggplot(merged_df, aes(x = AF)) + 
          geom_density(fill="gray",alpha=.6,adjust=smoothing_factor1) +
          geom_point(
            ## draw horizontal lines instead of points
            shape = 142, 
            size = 5,#size will need to be adjusted depending on how large the figure is 
            alpha = .3, aes(y=-.750,color=Subclonality,
                            text=paste("Chromosome:", chromosome_index, "<br>Position:", position))
          )+xlim(0,1)+
          ggtitle(paste(selectedSample, " - Purity =", round(purity, 3), "      ", "Mean read depth =", round(mean_rd,), "      ", "Range read depth =", min_rd, "-", max_rd)) +
          geom_text(data=data.frame(),aes(label=paste0("n clonal = ",subclonality_counts%>%filter(Subclonality=="Clonal")%>%select(n)),x=.85,y=-1.65))+
          geom_text(data=data.frame(),aes(label=paste0("n subclonal = ",subclonality_counts%>%filter(Subclonality=="Subclonal")%>%select(n)),x=.85,y=-2.65))+
          ggpubr::theme_pubr()+xlab("VAF")+ylab("Density")+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
        
        ggplotly(gg1)
      }
    })
    
    # here is where the captions for graph 1 are stored
    output$caption1 <- renderUI({
      if (plotsReady()) {
        HTML("<p> - <b>Variant Allele Frequency (VAF)</b> indicates the proportion of sequencing reads that support a variant allele. A VAF close to 0 indicates no variant alleles, while 1 indicates only variant alleles.<br>
         - <b>Purity</b> reflects the proportion of tumor cells in the sample. For instance, a purity of 0.9 means 90% of the cells are cancerous.<br>
         - <b>Mean read depth</b> is the average number of sequencing reads per SNV, influencing data reliability.<br>
         - <b>Density</b> estimates the probability density function of VAFs based on kernel density estimation (KDE). <br>
    </p>")
      }
    })
    
    output$plot2 <- renderPlotly({
      if (input$submit == 0 || !plotsReady()) return()  
      
      selectedSample <- input$selectedSample
      sampleFiles <- uploadedSamples$samples[[selectedSample]]
      
      if (is.null(sampleFiles)) return()
      
      snvFile <- sampleFiles$snvFile
      cnaFile <- sampleFiles$cnaFile
      purityFile <- sampleFiles$purityFile
      
      if (is.null(snvFile) || is.null(cnaFile) || is.null(purityFile)) return()
      
      # Reading data
      df_snv <- read.table(snvFile, header = TRUE, sep = "\t")
      df_cna <- read.table(cnaFile, header = TRUE, sep = "\t")
      purity <- as.numeric(readLines(purityFile))
      
      # Read mutation and subclonal structure assignments
      matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
      matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
      
      if (length(matching_file1) > 0 && length(matching_file2) > 0) {
        df_assigned <- read.table(matching_file1[1], header = TRUE, sep = "\t")
        df_subc <- read.table(matching_file2[1], header = TRUE, sep = "\t")
        
        df_snvcna <- df_snv %>%
          mutate(total_cn = map2_chr(position, chromosome_index, function(x, y) {
            inds = x >= df_cna$start_position & x <= df_cna$end_position & y == df_cna$chromosome_index
            if (any(inds)) df_cna$total_cn[which.max(inds)] else NA
          }))
        
        # Calculate the sum of ref_count and alt_count
        df_snvcna$sum_counts <- df_snvcna$ref_count + df_snvcna$alt_count
        
        # Compute the mean of the sum_counts column
        mean_rd <- mean(df_snvcna$sum_counts)
        min_rd <- min(df_snvcna$sum_counts)
        max_rd <- max(df_snvcna$sum_counts)
        
        merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
        merged_df['AF'] = merged_df['alt_count']/(merged_df['sum_counts'])
        
        merged_df$total_cn <- as.numeric(merged_df$total_cn)
        merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
        merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
        merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
        
        total_SNVs <- sum(df_subc$num_SNV)
        # smoothing is adjustable with slider
        smoothing_factor2 <- input$smoothingFactor2
        
        merged_df <- merged_df %>%
          mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"))
        
        print(head(merged_df))
        
        subclonality_counts <- merged_df %>%
          count(Subclonality)
        
        # Rename the count column to 'n' as specified in your desired data frame
        subclonality_counts <- subclonality_counts %>%
          rename(n = n)
        
        # Create plot
        gg2 <- ggplot(merged_df, aes(x = CP_unpenalized)) + 
          geom_density(fill="seagreen",alpha=.6,adjust=smoothing_factor2) +
          geom_point(
            ## draw horizontal lines instead of points
            shape = 142, 
            size = 5,#size will need to be adjusted depending on how large the figure is
            alpha = .3, aes(y=-.750,color=Subclonality,
                            text=paste("Chromosome:", chromosome_index, "<br>Position:", position))
          )+xlim(0,1.5)+
          ggtitle(paste(selectedSample, " - Estimated mutation-specific CP")) +
          geom_text(data=data.frame(),aes(label=paste0("n clonal = ",subclonality_counts%>%filter(Subclonality=="Clonal")%>%select(n)),x=1.3,y=-1.65))+
          geom_text(data=data.frame(),aes(label=paste0("n subclonal = ",subclonality_counts%>%filter(Subclonality=="Subclonal")%>%select(n)),x=1.3,y=-2.65))+
          ggpubr::theme_pubr() + xlab("Estimated mutation-specific CP") + ylab("Density") + 
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
          geom_segment(data=df_subc,aes(x = cellular_prevalence, xend=cellular_prevalence, y=0, yend=5.7), lty="dashed") 
        
        ggplotly(gg2)
      }
    })
    
    # here is where the captions for graph 2 are stored
    output$caption2 <- renderUI({
      if (plotsReady()) {
        HTML("<p> - <b>Cellular Prevalence estimated at each SNV without penalization (Estimated mutation-specific CP) </b> estimates gene copies in the tumor adjusted for purity without penalization, reflecting the observed vs. expected gene copy numbers.<br>
         - Annotations on the x-axis indicate whether mutations are <b>clonal</b> (present in all tumor cells) or <b>subclonal</b> (found in a subset).<br>
         - Counts of clonal and subclonal mutations are indicated as <b>n clonal</b> and <b>n subclonal</b> respectively.<br>
    </p>")
      }
    })
    
    # Makes it so that the output file does not appear until after the graphs
    # are generated
    output$plotOutputArea <- renderUI({
      if (plotsReady()) {
        fluidPage(
          h4("Output:"),
          textOutput("result"),
          selectInput("selectedFile", "Choose a file to download:", choices = c()),
          downloadButton("downloadResult", "Download Result")
        )
      }
    })
  })
}
