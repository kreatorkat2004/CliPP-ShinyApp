library(reticulate)
library(ggplot2)
library(dplyr)
library(purrr)
library(plotly)
library(readr)

server <- function(input, output, session) {
  if (Sys.info()[['user']] == 'shiny') {
    Sys.setenv(RETICULATE_PYTHON = "./python3")
  }
  
  outputDir <- getwd()  # Define the output directory globally
  resultPath <- file.path(outputDir, "forCliPP/CliPP/sample_id/final_result/Best_lambda")
  
  observeEvent(input$submit, {
    req(input$file1, input$file2, input$file3)  # Ensure all files are uploaded
    
    # Save the uploaded files temporarily
    inFile1 <- input$file1$datapath
    inFile2 <- input$file2$datapath
    inFile3 <- input$file3$datapath
    
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
  
  output$warning <- renderText({
    if (is.null(input$submit)) return()
  })
  
  output$plot1 <- renderPlotly({
    #browser()
    if (input$submit==0) return()
    
    req(input$file1, input$file2, input$file3)  # Ensure files are uploaded
    
    # Reading the uploaded files
    inFile1 <- input$file1$datapath
    inFile2 <- input$file2$datapath
    inFile3 <- input$file3$datapath
    
    # Checking for matching files
    matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
    matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
    
    if (input$submit && file.exists(inFile1) && file.exists(inFile2) && file.exists(inFile3) &&
        length(matching_file1) > 0 && length(matching_file2) > 0) {
      
      # Reading data
      df_snv <- read.table(inFile1, header = TRUE, sep = "\t")
      df_cna <- read.table(inFile2, header = TRUE, sep = "\t")
      purity <- as.numeric(readLines(inFile3))
      
      # Read mutation and subclonal structure assignments
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
      
      merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
      merged_df['AF'] = merged_df['alt_count']/(merged_df['sum_counts'])
      
      merged_df$total_cn <- as.numeric(merged_df$total_cn)
      merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
      merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
      merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
      
      total_SNVs <- sum(df_subc$num_SNV)
      smoothing_factor <- case_when(
        total_SNVs <= 100 ~ 0.2,
        total_SNVs <= 200 & total_SNVs > 100 ~ 0.3,
        total_SNVs <= 300 & total_SNVs > 100 ~ 0.4,
        total_SNVs <= 500 & total_SNVs > 300 ~ 0.5,
        TRUE ~ 0.6  
      )
      
      merged_df <- merged_df %>%
        mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"))
      
      print(head(merged_df))
      
      subclonality_counts <- merged_df %>%
        count(Subclonality)
      
      # Rename the count column to 'n' as specified in your desired data frame
      subclonality_counts <- subclonality_counts %>%
        rename(n = n)
      
      gg1 <- ggplot(merged_df#this is the dataframe with VAF values and CP values for each input SNV, along with their clonal assignments
             , aes(x = AF)) + 
        geom_density(fill="gray",alpha=.6,adjust=smoothing_factor)+
        geom_point(
          ## draw horizontal lines instead of points
          shape = 142,
          size = 5,#size will need to be adjusted depending on how large the figure is
          alpha = .3,aes(y=-.750,color=Subclonality)
        )+xlim(0,1)+
        ggtitle(paste0("Purity = ",round(purity,3),"      ","Mean read-depth = ",round(mean_rd,)))+
        geom_text(data=data.frame(),aes(label=paste0("n clonal = ",subclonality_counts%>%filter(Subclonality=="Clonal")%>%select(n)),x=.85,y=-1.65))+
        geom_text(data=data.frame(),aes(label=paste0("n subclonal = ",subclonality_counts%>%filter(Subclonality=="Subclonal")%>%select(n)),x=.85,y=-2.65))+
        ggpubr::theme_pubr()+xlab("VAF")+ylab("Density")+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
      
      ggplotly(gg1)
    }
  })
  
  output$plot2 <- renderPlotly({
    #browser()
    if (input$submit==0) return()
    
    req(input$file1, input$file2, input$file3)  # Ensure files are uploaded
    
    # Reading the uploaded files
    inFile1 <- input$file1$datapath
    inFile2 <- input$file2$datapath
    inFile3 <- input$file3$datapath
    
    # Checking for matching files
    matching_file1 <- list.files(resultPath, pattern = "mutation_assignments_lam.*\\.txt", full.names = TRUE)
    matching_file2 <- list.files(resultPath, pattern = "subclonal_structure_lam.*\\.txt", full.names = TRUE)
    
    if (input$submit && file.exists(inFile1) && file.exists(inFile2) && file.exists(inFile3) &&
        length(matching_file1) > 0 && length(matching_file2) > 0) {
      
      # Reading data
      df_snv <- read.table(inFile1, header = TRUE, sep = "\t")
      df_cna <- read.table(inFile2, header = TRUE, sep = "\t")
      purity <- as.numeric(readLines(inFile3))
      
      # Read mutation and subclonal structure assignments
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
      
      merged_df <- merge(df_snvcna, df_assigned, by = c("chromosome_index", "position"))
      merged_df['AF'] = merged_df['alt_count']/(merged_df['sum_counts'])
      
      merged_df$total_cn <- as.numeric(merged_df$total_cn)
      merged_df$b_i_V <- pmax(1, round(merged_df$AF * (1 / purity) * ((purity * merged_df$total_cn) + (2 * (1 - purity))), 0))
      merged_df$CP_unpenalized <- (merged_df$alt_count * ((1 - purity) * 2 + purity * merged_df$total_cn)) / (merged_df$b_i_V * (merged_df$alt_count + merged_df$ref_count))
      merged_df$CCF_unpenalized <- merged_df$CP_unpenalized / purity
      
      total_SNVs <- sum(df_subc$num_SNV)
      smoothing_factor <- case_when(
        total_SNVs <= 100 ~ 0.2,
        total_SNVs <= 200 & total_SNVs > 100 ~ 0.3,
        total_SNVs <= 300 & total_SNVs > 100 ~ 0.4,
        total_SNVs <= 500 & total_SNVs > 300 ~ 0.5,
        TRUE ~ 0.6  
      )
      
      merged_df <- merged_df %>%
        mutate(Subclonality = ifelse(cluster_index == 0, "Clonal", "Subclonal"))
      
      print(head(merged_df))
      
      subclonality_counts <- merged_df %>%
        count(Subclonality)
      
      # Rename the count column to 'n' as specified in your desired data frame
      subclonality_counts <- subclonality_counts %>%
        rename(n = n)
      
      gg2 <- ggplot(merged_df#this is the dataframe with VAF values and CP values for each input SNV, along with their clonal assignments
                    , aes(x = CP_unpenalized)) + 
        geom_density(fill="seagreen",alpha=.6,adjust=smoothing_factor)+
        geom_point(
          ## draw horizontal lines instead of points
          shape = 142,
          size = 5,#size will need to be adjusted depending on how large the figure is
          alpha = .3,aes(y=-.750,color=Subclonality)
        )+xlim(0,1.5)+
        #ggtitle(paste0("Purity = ",round(purity,3),"      ","Mean read depth = ",round(mean_rd,)))+
        geom_text(data=data.frame(),aes(label=paste0("n clonal = ",subclonality_counts%>%filter(Subclonality=="Clonal")%>%select(n)),x=1.3,y=-1.65))+
        geom_text(data=data.frame(),aes(label=paste0("n subclonal = ",subclonality_counts%>%filter(Subclonality=="Subclonal")%>%select(n)),x=1.3,y=-2.65))+
        ggpubr::theme_pubr()+xlab("Unpenalized CP")+ylab("Density")+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
        geom_segment(data=df_subc,aes(x = cellular_prevalence, xend=cellular_prevalence,y=0,yend=5.7),lty="dashed")
      
      ggplotly(gg2)
    }
  })
}
