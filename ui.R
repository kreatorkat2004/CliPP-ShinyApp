ui <- fluidPage(
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
      Hongtu Zhu, Wenyi Wang. Pan-cancer subclonal mutation analysis of 7,827 tumors predicts clinical outcome. Nature Genetics 2024."),
    br(),
    p("For issues with the app, please contact:"),
    p("Wenyi Wang - ", a(href = "wwang7@mdanderson.org", "wwang7@mdanderson.org")),
    p("Quang Tran - ", a(href = "qmtran@mdanderson.org", "qmtran@mdanderson.org")),
    br(),
    p("Please choose files to upload using the form below:")
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "SNV file:"),
      fileInput("file2", "CNA file:"),
      fileInput("file3", "Purity file:"),
      actionButton("submit", "Submit", class = "btn-primary")
    ),
    mainPanel(
      h4("Output:"),
      textOutput("result"),
      selectInput("selectedFile", "Choose a file to download:", choices = c()),
      downloadButton("downloadResult", "Download Result")
    )
  )
)
