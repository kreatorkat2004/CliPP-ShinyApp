library(shiny)
library(plotly)
library(shinyjs)
library(shinyWidgets)

# UI for the Shiny application
ui <- fluidPage(
  # shinyjs for additional JavaScript functionality
  tags$head(includeHTML("google-analytics.html")),
  shinyjs::useShinyjs(),
  # CSS for styling a progress bar overlay
  tags$head(
    tags$style(HTML("
      /* Progress bar styling */
      #progress {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        width: 50%;
        z-index: 9999;
      }
      
      /* Plot and controls layout */
      .plot-container {
        width: 100%;
        margin-bottom: 20px;
      }
      .plot-with-controls {
        display: flex;
        gap: 20px;
        margin-bottom: 10px;
      }
      .plot-area {
        flex: 3;
      }
      .controls-area {
        flex: 1;
        background: #f5f5f5;
        padding: 15px;
        border-radius: 8px;
        min-width: 200px;
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
  div(p("Copyright ", HTML("&copy;"), " 2024 Wang Lab at MD Anderson. All rights reserved.")),
  
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
                            #selectInput("penaltyType", "Select Penalty Type:", choices = c("Lasso", "MCP"), selected = "Lasso"),
                            hidden(actionButton("clear", "Clear", class = "btn-danger"))
                          ),
                          
                          # Main panel to display plots and analysis output
                          mainPanel(
                            # First plot container
                            div(class = "plot-container",
                                div(class = "plot-area",
                                    plotlyOutput("plot1")
                                ),
                                hidden(div(class = "controls-area",
                                           sliderInput("smoothingFactor1", "Smoothing Factor:",
                                                       min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                           div(class = "recommended-smoothing",
                                               uiOutput("smoothingExplanation1")
                                           )
                                )),
                                div(class = "caption-area",
                                    uiOutput("caption1")
                                )
                            ),
                            
                            # Second plot container
                            div(class = "plot-container",
                                div(class = "plot-area",
                                    plotlyOutput("plot2")
                                ),
                                hidden(div(class = "controls-area",
                                           sliderInput("smoothingFactor2", "Smoothing Factor:",
                                                       min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                           div(class = "recommended-smoothing",
                                               uiOutput("smoothingExplanation2")
                                           )
                                )),
                                div(class = "caption-area",
                                    uiOutput("caption2")
                                )
                            ),
                            
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
                                         hidden(selectInput("colorByTCGA", "Color points by:", 
                                                            choices = list("Clonality" = "clonality", "Read Depth" = "readDepth"), 
                                                            selected = "clonality"))
                                       ),
                                       #selectInput("penaltyType", "Select Penalty Type:", choices = c("Lasso", "MCP"), selected = "Lasso"),
                                       mainPanel(
                                         # First TCGA plot container
                                         div(class = "plot-container",
                                             div(class = "plot-area",
                                                 plotlyOutput("plot1TCGA") 
                                             ),
                                             hidden(div(class = "controls-area",
                                                        sliderInput("smoothingFactor1_TCGA", "Smoothing Factor:",
                                                                    min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                                        div(class = "recommended-smoothing",
                                                            uiOutput("smoothingExplanation1_TCGA")
                                                        )
                                             )),
                                             div(class = "caption-area",
                                                 uiOutput("caption1TCGA")
                                             )
                                         ),
                                         
                                         # Second TCGA plot container
                                         div(class = "plot-container",
                                             div(class = "plot-area",
                                                 plotlyOutput("plot2TCGA") 
                                             ),
                                             hidden(div(class = "controls-area",
                                                        sliderInput("smoothingFactor2_TCGA", "Smoothing Factor:",
                                                                    min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                                        div(class = "recommended-smoothing",
                                                            uiOutput("smoothingExplanation2_TCGA")
                                                        )
                                             )),
                                             div(class = "caption-area",
                                                 uiOutput("caption2TCGA")
                                             )
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
                                         actionButton("submitPCAWG", "Submit", class = "btn-primary"),
                                         hidden(selectInput("colorByPCAWG", "Color points by:", 
                                                            choices = list("Clonality" = "clonality", "Read Depth" = "readDepth"), 
                                                            selected = "clonality"))
                                       ),
                                       #selectInput("penaltyType", "Select Penalty Type:", choices = c("Lasso", "MCP"), selected = "Lasso"),
                                       mainPanel(
                                         # First PCAWG plot container
                                         div(class = "plot-container",
                                             div(class = "plot-area",
                                                 plotlyOutput("plot1PCAWG")
                                             ),
                                             hidden(div(class = "controls-area",
                                                        sliderInput("smoothingFactor1_PCAWG", "Smoothing Factor:",
                                                                    min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                                        div(class = "recommended-smoothing",
                                                            uiOutput("smoothingExplanation1_PCAWG")
                                                        )
                                             )),
                                             div(class = "caption-area",
                                                 uiOutput("caption1PCAWG")
                                             )
                                         ),
                                         
                                         # Second PCAWG plot container
                                         div(class = "plot-container",
                                             div(class = "plot-area",
                                                 plotlyOutput("plot2PCAWG") 
                                             ),
                                             hidden(div(class = "controls-area",
                                                        sliderInput("smoothingFactor2_PCAWG", "Smoothing Factor:",
                                                                    min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                                        div(class = "recommended-smoothing",
                                                            uiOutput("smoothingExplanation2_PCAWG")
                                                        )
                                             )),
                                             div(class = "caption-area",
                                                 uiOutput("caption2PCAWG")
                                             )
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
                          uiOutput("cancer_type_ui_driver")
                        ),
                        mainPanel(
                          DT::DTOutput("driverMutationTable"),
                          uiOutput("selectedMutationInfo"),
                          downloadButton("downloadDriverMutation", "Download Data"),
                          tags$h4(tags$b(tags$span("Reference:"))),
                          p("Martínez-Jiménez, Francisco, et al. 'A compendium of mutational cancer driver genes.' Nature Reviews Cancer 20.10 (2020): 555-572."),
                          
                          # First driver mutation plot container
                          div(class = "plot-container",
                              div(class = "plot-area",
                                  plotlyOutput("driverPlot1") 
                              ),
                              hidden(div(class = "controls-area",
                                         sliderInput("smoothingFactor1_driver", "Smoothing Factor:",
                                                     min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                         div(class = "recommended-smoothing",
                                             uiOutput("smoothingExplanation1_driver")
                                         )
                              )),
                              div(class = "caption-area",
                                  uiOutput("driverCaption1")
                              )
                          ),
                          
                          # Second driver mutation plot container
                          div(class = "plot-container",
                              div(class = "plot-area",
                                  plotlyOutput("driverPlot2") 
                              ),
                              hidden(div(class = "controls-area",
                                         sliderInput("smoothingFactor2_driver", "Smoothing Factor:",
                                                     min = 0.2, max = 0.7, value = 0.45, width = "100%"),
                                         div(class = "recommended-smoothing",
                                             uiOutput("smoothingExplanation2_driver")
                                         )
                              )),
                              div(class = "caption-area",
                                  uiOutput("driverCaption2")
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
