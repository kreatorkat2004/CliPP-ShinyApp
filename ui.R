library(shiny)
library(plotly)
library(shinyjs)
library(shinyWidgets)

# Function to write UI plot elements
create_plot_container <- function(plot_output_id, smoothing_id, explanation_id, caption_id) {
  div(class = "plot-container",
      div(class = "plot-area",
          plotlyOutput(plot_output_id)
      ),
      hidden(
        div(class = "controls-area",
            sliderInput(smoothing_id, "Smoothing Factor:",
                        min = 0.2, max = 0.7, value = 0.45, width = "100%"
            ),
            div(class = "recommended-smoothing",
                uiOutput(explanation_id)
            )
        )
      ),
      div(class = "caption-area",
          uiOutput(caption_id)
      )
  )
}

# UI for the Shiny application
ui <- fluidPage(
  # shinyjs for additional JavaScript functionality
  tags$head(includeHTML("google-analytics.html")),
  shinyjs::useShinyjs(),
  # CSS for styling a progress bar overlay
  tags$head(
    tags$style(HTML("
      #progress {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        width: 50%;
        z-index: 9999;
      }
      .plot-container {
        position: relative;
        width: 100%;
        margin-bottom: 20px;
      }
      .plot-with-controls {
        display: flex;
        flex-direction: column;
        gap: 10px;
        margin-bottom: 10px;
      }
      .plot-area {
        width: 100%;
      }
      .controls-area {
        position: absolute;
        top: 28%;
        left: 86%;
        width: 20%;
        min-width: 100px;
        max-width: 200px;
        background: #f5f5f5;
        padding: 15px;
        border-radius: 8px;
        z-index: 100;
      } 
      .caption-area {
        width: 100%;
        margin-bottom: 30px;
      }
      .recommended-smoothing {
        font-size: 0.9em;
        color: #666;
        margin-top: 10px;
      }
    "))
  ),
  
  # Title and attribution for research purposes
  titlePanel("CliPP-on-Web: Clonal structure identification through penalizing pairwise differences"),
  div(p("Copyright ", HTML("&copy;"), " 2025 The University of Texas MD Anderson Cancer Center. All rights reserved.")),
  
  # Multi-page navigation structure
  navbarPage("",
             
             # Page 1: About section with application description and contact info
             tabPanel("About",
                      tags$div(
                        tags$h4(tags$b(tags$span("Description"))),
                        p("Intra-tumor heterogeneity is an important driver of tumor evolution and therapy response. Advances in precision cancer treatment will require understanding of 
                        mutation clonality and subclonal architecture. Currently the slow computational speed of subclonal reconstruction hinders large cohort studies. To overcome this 
                        bottleneck, we developed Clonal structure identification through Pairwise Penalization, or CliPP, which clusters subclonal mutations using a regularized likelihood 
                        model. CliPP reliably processed whole-genome and whole-exome sequencing data from over 12,000 tumor samples within 24 hours, thus enabling large-scale downstream 
                        association analyses between subclonal structures and clinical outcomes. Through a pan-cancer investigation of 7,827 tumors from 32 cancer types, we found that high 
                        subclonal mutational load (sML), a measure of latency time in tumor evolution, was significantly associated with better patient outcomes in 16 cancer types with low to 
                        moderate tumor mutation burden (TMB). In a cohort of prostate cancer patients participating in an immunotherapy clinical trial, high sML was indicative of favorable 
                        response to immune checkpoint blockade. This comprehensive study using CliPP underscores sML as a key feature of cancer. sML may be essential for linking mutation 
                        dynamics with immunotherapy response in the large population of non-high TMB cancers. Also see our paper:",
                          a(href = "https://www.biorxiv.org/content/10.1101/2024.07.03.601939v1",
                            "https://www.biorxiv.org/content/10.1101/2024.07.03.601939v1"),"."),
                        tags$br(),
                        p(tags$b(tags$span(style = "color: #D22B2B;", "For issues with the app, please contact:"))),
                        p("Aaron Wu - ", a(href = "mailto:aw80@rice.edu", "aw80@rice.edu")),
                        p("Quang Tran - ", a(href = "mailto:qmtran@mdanderson.org", "qmtran@mdanderson.org"))
                      )
             ),
             
             # Page 2: Tutorial section with an embedded video
             tabPanel("Tutorial",
                      tags$div(
                        tags$h4(tags$b("Video Tutorial")),
                        tags$p("Watch the video below to learn how to use the CliPP-on-Web application:"),
                        tags$iframe(src="https://www.youtube.com/embed/ZPX_hC-uiYA",
                                    style="width:80%; height:700px; border:none;")
                      )
             ),
             
             # Page 3: Main functionality for uploading and analyzing data with CliPP-on-Web
             tabPanel("CliPP-on-Web",
                      tags$div(
                        tags$h4(tags$b("Upload Files")),
                        sidebarLayout(
                          
                          # Sidebar for file upload and input controls
                          sidebarPanel(
                            tags$div(
                              p("The online tool requires three files as input: SNV File, CNA File and Cancer Purity File.
                Please check the three links below to download the example files:"),
                              p(a("SNV file", href="sample.snv.txt", download="sample.snv.txt"),
                                "TSV file with columns: chromosome_index, position, ref_count, alt_count", br(),
                                a("CNA file", href="sample.cna.txt", download="sample.cna.txt"),
                                "TSV file with columns: chromosome_index, start_position, end_position, major_cn, minor_cn, and total_cn", br(),
                                a("Purity file", href="sample.purity.txt", download="sample.purity.txt"),
                                "File storing a scalar purity value between 0 and 1"
                              )
                            ),
                            radioButtons("uploadMode", "Choose your upload method:",
                                         choices = list("Upload all files at once" = "all_at_once",
                                                        "Upload files individually" = "individually"),
                                         selected = "all_at_once"),
                            numericInput("numSamples", "Number of samples:", min = 1, max = 10, value = 1),
                            actionButton("setSamples", "Set Samples", class = "btn-primary"),
                            br(),
                            uiOutput("dynamicFileInputs"),
                            hidden(actionButton("uploadSamples", "Upload Samples", class = "btn-primary")),
                            br(),
                            hidden(selectInput("selectedSample", "Select a sample to analyze:", choices = NULL)),
                            hidden(actionButton("submit", "Submit", class = "btn-primary")),
                            hidden(selectInput("colorBy", "Color points by:", 
                                               choices = list("Clonality" = "clonality", "Read Depth" = "readDepth"), 
                                               selected = "clonality")),
                            hidden(actionButton("clear", "Clear", class = "btn-danger"))
                          ),
                          
                          # Main panel to display plots and analysis output
                          mainPanel(
                            # First plot container
                            create_plot_container("plot1", "smoothingFactor1", "smoothingExplanation1", "caption1"),
                            
                            # Second plot container
                            create_plot_container("plot2", "smoothingFactor2", "smoothingExplanation2", "caption2"),
                            
                            # Plot output area and progress bar
                            uiOutput("plotOutputArea"),
                            hidden(div(id = "plotReadyTextCliPP")),
                            shinyjs::hidden(
                              shinyWidgets::progressBar(
                                id = "progress",
                                value = 0,
                                total = 100,
                                display_pct = TRUE
                              )
                            )
                          )
                        )
                      )
             ),
             
             # Page 4: CliPP Data Resources with automatic data loading and visualization
             tabPanel("CliPP Data Resources",
                      tags$div(
                        tabsetPanel(
                          # Subsection for TCGA data selection and analysis
                          tabPanel("TCGA",
                                   tags$div(
                                     tags$h4(tags$b("TCGA Select")),
                                     sidebarLayout(
                                       sidebarPanel(
                                         uiOutput("cancer_type_ui_TCGA"),
                                         uiOutput("sample_name_ui_TCGA"),
                                         actionButton("submitTCGA", "Submit", class = "btn-primary"),
                                         hidden(actionButton("clearTCGA", "Clear", class = "btn-danger")),
                                         hidden(selectInput("colorByTCGA", "Color points by:", 
                                                            choices = list("Clonality" = "clonality", "Read Depth" = "readDepth"), 
                                                            selected = "clonality"))
                                       ),
                                       mainPanel(
                                         # First TCGA plot container
                                         create_plot_container(
                                           paste0("plot1TCGA"), 
                                           paste0("smoothingFactor1_TCGA"),
                                           paste0("smoothingExplanation1_TCGA"),
                                           paste0("caption1TCGA")
                                         ),
                                         
                                         # Second TCGA plot container
                                         create_plot_container(
                                           paste0("plot2TCGA"), 
                                           paste0("smoothingFactor2_TCGA"),
                                           paste0("smoothingExplanation2_TCGA"),
                                           paste0("caption2TCGA")
                                         ),
                                         
                                         hidden(uiOutput("plotOutputAreaTCGA")),
                                         hidden(div(id = "plotReadyTextTCGA")),
                                         uiOutput("downloadUI_TCGA")
                                       )
                                     )
                                   )
                          ),
                          
                          # Subsection for PCAWG data selection and analysis
                          tabPanel("PCAWG",
                                   tags$div(
                                     tags$h4(tags$b("PCAWG Select")),
                                     sidebarLayout(
                                       sidebarPanel(
                                         uiOutput("cancer_type_ui_PCAWG"),
                                         uiOutput("sample_name_ui_PCAWG"),
                                         hidden(actionButton("clearPCAWG", "Clear", class = "btn-danger")),
                                         actionButton("submitPCAWG", "Submit", class = "btn-primary"),
                                         hidden(selectInput("colorByPCAWG", "Color points by:", 
                                                            choices = list("Clonality" = "clonality", "Read Depth" = "readDepth"), 
                                                            selected = "clonality"))
                                       ),
                                       mainPanel(
                                         # First PCAWG plot container
                                         create_plot_container(
                                           paste0("plot1PCAWG"), 
                                           paste0("smoothingFactor1_PCAWG"),
                                           paste0("smoothingExplanation1_PCAWG"),
                                           paste0("caption1PCAWG")
                                         ),
                                         
                                         # Second PCAWG plot container
                                         create_plot_container(
                                           paste0("plot2PCAWG"), 
                                           paste0("smoothingFactor2_PCAWG"),
                                           paste0("smoothingExplanation2_PCAWG"),
                                           paste0("caption2PCAWG")
                                         ),
                                         
                                         hidden(uiOutput("plotOutputAreaPCAWG")),
                                         hidden(div(id = "plotReadyTextPCAWG")),
                                         uiOutput("downloadUI_PCAWG")
                                       )
                                     )
                                   )
                          )
                        )
                      )
             ),
             
             # Page 5: Driver Mutation analysis tool
             tabPanel("Driver Mutation",
                      sidebarLayout(
                        sidebarPanel(
                          shinyjs::hidden(
                            div(id = "active_driver_tab", "gene")
                          ),
                          
                          # Cancer type selection
                          uiOutput("cancer_type_ui_driver"),
                          
                          # Gene selection 
                          conditionalPanel(
                            condition = "input.active_driver_tab == 'gene'",
                            uiOutput("gene_name_ui_driver")
                          ),
                          
                          # Sample selection 
                          conditionalPanel(
                            condition = "input.active_driver_tab == 'sample'",
                            uiOutput("sample_name_ui_driver")
                          )
                        ),
                        mainPanel(
                          tabsetPanel(id = "driver_mutation_tabs",
                                      tabPanel("By Gene", value = "gene",
                                               DT::DTOutput("driverMutationTable"),
                                               shinyjs::hidden(
                                                 downloadButton("downloadDriverMutation", "Download Data")
                                               ),
                                               tags$h4(tags$b(tags$span("Reference:"))),
                                               p("Martínez-Jiménez, Francisco, et al. 'A compendium of mutational cancer driver genes.' Nature Reviews Cancer 20.10 (2020): 555-572.")
                                      ),
                                      tabPanel("By Sample", value = "sample",
                                               DT::DTOutput("driverMutationTableBySample"),
                                               br(),
                                               downloadButton("downloadDriverMutationBySample", "Download Sample Data")
                                      )
                          )
                        )
                      )
             ),
             
             # Page 6: Citations for related studies and resources
             tabPanel("Citations",
                      tags$div(
                        tags$h4(tags$b(tags$span("Manuscript under preparation:"))),
                        p("Yujie Jiang, Matthew D Montierth, Kaixian Yu, Shuangxi Ji,
           Quang Tran, Xiaoqian Liu, Jessica C Lal, Shuai Guo, Aaron Wu,
           Seung Jun Shin, Shaolong Cao, Ruonan Li, Yuxin Tang, Tom Lesluyes,
           Scott Kopetz, Pavlos Msaouel, Anil K. Sood, Christopher Jones, Jaffer Ajani,
           Sumit K Subudhi, Ana Aparicio, Padmanee Sharma, John Paul Shen, Marek Kimmel,
           Jennifer R. Wang, Maxime Tarabichi, Rebecca Fitzgerald, Peter Van Loo, Hongtu Zhu,
           Wenyi Wang. Subclonal mutational load predicts survival and response to immunotherapy
           in cancers with low to moderate TMB.",
                          a(href = "https://www.biorxiv.org/content/10.1101/2024.07.03.601939v1",
                            "https://www.biorxiv.org/content/10.1101/2024.07.03.601939v1"),".")
                      )
             )
  )
)
