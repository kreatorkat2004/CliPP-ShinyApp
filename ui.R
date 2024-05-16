library(plotly)
library(shinyjs)
ui <- fluidPage(
  shinyjs::useShinyjs(), 
  titlePanel("CliPP-on-Web: Clonal structure identification through penalizing pairwise differences"),
  # Introduction
  tags$div(
    p("Subpopulations of tumor cells characterized by mutation profiles may confer differential fitness 
      and consequently influence prognosis of cancers. Understanding subclonal architecture has the potential 
      to provide biological insight in tumor evolution and advance precision cancer treatment. Recent methods 
      comprehensively integrate single nucleotide variants (SNVs) and copy number aberrations (CNAs) to 
      reconstruct subclonal architecture using whole-genome or whole-exome sequencing (WGS, WES) data 
      from bulk tumor samples. However, the commonly used Bayesian methods require a large amount of 
      computational resources, a prior knowledge of the number of subclones, and extensive post-processing. 
      Regularized likelihood modeling approach, never explored for subclonal reconstruction, can inherently 
      address these drawbacks. We therefore propose a model-based method, Clonal structure identification 
      through Pair-wise Penalization, or CliPP, for clustering subclonal mutations without prior knowledge 
      or post-processing. The CliPP model is applicable to genomic regions with or without CNAs. CliPP 
      demonstrates high accuracy in subclonal reconstruction through extensive simulation studies. 
      A penalized likelihood framework for subclonal reconstruction will help address intrinsic drawbacks 
      of existing methods and expand the scope of computational analysis for cancer evolution in large 
      cancer genomic studies. Also see our paper:",
      a(href = "https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1",
        "https://www.biorxiv.org/content/10.1101/2021.03.31.437383v1"),"."),
    p("Copyright ", HTML("&copy;"), " 2024 Wang Lab at MD Anderson. All rights reserved."),
    br(),
    h5("Manuscript under preparation:"),
    p("Yujie Jiang, Matthew D Montierth, Kaixian Yu, Shuangxi Ji, Shuai Guo, Quang Tran, Seung Jun Shin, Shaolong Cao, 
      Ruonan Li, Yuxin Tang, Tom Lesluyes, Scott Kopetz, Jaffer Ajani, Pavlos Msaouel, Sumit K Subudhi, Ana Aparicio, 
      Padmanee Sharma, John Paul Shen, Anil K. Sood, Maxime Tarabichi, Jennifer R. Wang, Marek Kimmel, Peter Van Loo, 
      Hongtu Zhu, Wenyi Wang, Aaron Wu. Pan-cancer subclonal mutation analysis of 7,827 tumors predicts clinical outcome."),
    br(),
    p("For issues with the app, please contact:"),
    p("Wenyi Wang - ", a(href = "wwang7@mdanderson.org", "wwang7@mdanderson.org")),
    p("Quang Tran - ", a(href = "qmtran@mdanderson.org", "qmtran@mdanderson.org")),
    p("Matthew Montierth - ", a(href = "montiert@bcm.edu", "montiert@bcm.edu")),
    p("Aaron Wu - ", a(href = "aw80@rice.edu", "aw80@rice.edu")),
    br(),
    p("The online tool requires three files as input: SNV File, CNA File and Cancer Purity File.
  Please check the three links below to download the example files:"),
    p(a("SNV file", href="sample.snv.txt", download="sample.snv.txt"), 
      "TSV file with columns: chromosome_index, position, ref_count, alt_count", br(),
      a("CNA file", href="sample.cna.txt", download="sample.cna.txt"), 
      "TSV file with columns: chromosome_index, start_position, end_position, major_cn, minor_cn, and total_cn", br(),
      a("Purity file", href="sample.purity.txt", download="sample.purity.txt"), 
      "File storing a scalar purity value between 0 and 1")
  ),
  sidebarLayout(
    sidebarPanel(
      # Two inputs that alternate upload mode
      radioButtons("uploadMode", "Choose your upload method:",
                   choices = list("Upload all files at once" = "all_at_once",
                                  "Upload files individually" = "individually"),
                   selected = "all_at_once"),
      
      # Number input to give amount of samples
      numericInput("numSamples", "Number of samples:", min = 1, max = 10, value = 1),
      actionButton("setSamples", "Set Samples", class = "btn-primary"),
      br(),
      uiOutput("dynamicFileInputs"),
      actionButton("uploadSamples", "Upload Samples", class = "btn-primary"),
      br(),
      selectInput("selectedSample", "Select a sample to analyze:", choices = NULL),
      actionButton("submit", "Submit", class = "btn-primary"),
      # Adds the smoothing slider
      conditionalPanel(
        condition = "output.plotReadyText == 'true'", 
        sliderInput("smoothingFactor1", "Smoothing Factor for Plot 1:",
                    min = 0.2, max = 0.7, value = 0.5, step = 0.1),
        sliderInput("smoothingFactor2", "Smoothing Factor for Plot 2:",
                    min = 0.2, max = 0.7, value = 0.5, step = 0.1)
      )
    ),
    # Produces the main plots
    mainPanel(
      plotlyOutput("plot1"),
      uiOutput("caption1"),
      plotlyOutput("plot2"),
      uiOutput("caption2"),
      uiOutput("plotOutputArea"),
      hidden(div(id = "plotReadyText")), 
    )
  )
)
