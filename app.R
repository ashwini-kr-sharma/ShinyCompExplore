#-------------------------------------------------------------------------------
# Required packages
#-------------------------------------------------------------------------------

options(repos = BiocManager::repositories()) 

library(shiny)
library(shinythemes)
library(shinyjs)
library(pheatmap)
library(plotly)
library(waiter)
library(DT)
library(immunedeconv)
library(data.table)
library(pdist)
library(parallel)
library(shinyWidgets)
library(fgsea)
library(cowplot)
library(matrixStats)
library(nFactors)

#library(deconica)
#library(httpuv)
#library(ica)
#library(fastICA)
#library(NMF)
#library(CellMix)

#-------------------------------------------------------------------------------
# Internal input for Unsupervised enrichment 
#-------------------------------------------------------------------------------

dbmarkers = readRDS("data/CellMatchCurated.RDS")
cancertypes = levels(as.factor(dbmarkers$cancerType))

#-------------------------------------------------------------------------------
# Functions for Unsupervised enrichment 
#-------------------------------------------------------------------------------

enrichplot = function(mydata, markerscelltypes, showCategory){
  enrich = function(x, markerscelltypes, showCategory){
    genes=sort(x,decreasing=T)
    #set.seed(1)
    fgseaRes = fgseaMultilevel(markerscelltypes, genes, minSize=2, maxSize=200)
    fgseaResDF = data.frame(fgseaRes)
    fgseaResDFSIGN = fgseaResDF[which(fgseaResDF[,"pval"]<0.05 & fgseaResDF[,"ES"]>0),c("pathway", "pval", "ES")]
    fgseaResDFSIGN = fgseaResDFSIGN[order(fgseaResDFSIGN[,"pval"], decreasing = TRUE),]
    if(nrow(fgseaResDFSIGN)==0){
      p = print("No significant cell-type enrichment")
    } else {
      fgseaResDFSIGN = fgseaResDFSIGN[1:min(showCategory, nrow(x)),]
      fgseaResDFSIGN = fgseaResDFSIGN[!is.na(fgseaResDFSIGN[,1]),]
      fgseaResDFSIGN[,"-log10(pval)"] = -log10(fgseaResDFSIGN[,"pval"])
      fgseaResDFSIGN$pathway <- factor(fgseaResDFSIGN$pathway, levels=fgseaResDFSIGN$pathway)
      p <- ggplot(fgseaResDFSIGN, aes(x=pathway, y=-log10(pval), fill =ES))+
        geom_bar(stat='identity') + scale_fill_continuous(low="white", high="red")+coord_flip()
    }
    return(p)
  }
  
  allplots <- apply(mydata, 2,  enrich, markerscelltypes, showCategory)
  #subplot(allplots, nrows = length(allplots), margin = 0.07)#, margin = c(0,0,0.07,0.07))
  cowplot::plot_grid(plotlist = allplots, ncol = 2, labels = colnames(mydata))
}

orient_funct <- function(S) {
  orient <-
    apply(S, 2, function(x) {
      if (min(x) < -3 & max(x) > 3) {
        ifelse (sum(x > 3)  < sum(x < -3), -1, 1)
      } else {
        ifelse (sum(x > 2)  < sum(x < -2), -1, 1)
      }
    })
  S <- as.matrix(S)  %*% diag(orient)
  return(S)
}

#-------------------------------------------------------------------------------
# Functions for identifying Unsupervised number of components 
#-------------------------------------------------------------------------------

findElbowPoint_PCAtools <- function(variance) {
  if (is.unsorted(-variance)) {
    stop("'variance' should be sorted in decreasing order")
  }
  # Finding distance from each point on the curve to the diagonal.
  dy <- -diff(range(variance))
  dx <- length(variance) - 1
  l2 <- sqrt(dx^2 + dy^2)
  dx <- dx/l2
  dy <- dy/l2
  dy0 <- variance - variance[1]
  dx0 <- seq_along(variance) - 1
  
  parallel.l2 <- sqrt((dx0 * dx)^2 + (dy0 * dy)^2)
  normal.x <- dx0 - dx * parallel.l2
  normal.y <- dy0 - dy * parallel.l2
  normal.l2 <- sqrt(normal.x^2 + normal.y^2)
  
  # Picking the maximum normal that lies below the line.
  # If the entire curve is above the line, we just pick the last point.
  below.line <- normal.x < 0 & normal.y < 0
  if (!any(below.line)) {
    length(variance)
  } else {
    which(below.line)[which.max(normal.l2[below.line])]
  }
}

plot_k = function(D){
  #X <- bigstatsr::as_FBM(t(D))
  #svd <- bigstatsr::big_SVD(X, bigstatsr::big_scale())
  svd1 = svd(t(D))$d[1:min(30, ncol(D))]
  nbf = findElbowPoint_PCAtools(svd1)
  ggplot(data.frame(val = svd1, idx = 1:min(30, ncol(D))), aes(x = idx, y = val)) +
    geom_line() + geom_vline(xintercept=nbf, linetype="dashed", color = "red", size=2) +
    geom_point() + scale_x_continuous(breaks = 1:min(30, ncol(D)))+
    labs(x = "PC index", y = "Eigenvalues") +
    theme_minimal(base_size = 16)
}

#-------------------------------------------------------------------------------
# UI.app
#-------------------------------------------------------------------------------

ui <- fluidPage(
  
  useShinyjs(),
  use_waiter(),
  #shinythemes::themeSelector(),
  
#-------------------------------------------------------------------------------
# Navbar
#-------------------------------------------------------------------------------
  
  navbarPage("DeconExplorer",
             theme = shinytheme("united"),
             
#-------------------------------------------------------------------------------
# Navbar : About
#-------------------------------------------------------------------------------

tabPanel(
  "About",
  
  h1("DeconExplorer"),
  
  tags$b("DeconpExplorer - Deconvolution explorer"),
  "is a visulatization tool to guide the user in the analysis and interpretation of the the results from",
  em("Supervised"),
  "and",
  em("Unsupervised"),
  "gene expresssion deconvolution algorithms.",
  
  tags$br(),
  tags$br(),
  
  img(src = "deconvolution.jpg", width = 800),
  
  tags$br(),
  tags$br(),
  
  h1("Terminology"),
  
  tags$ol(
    tags$li(
      tags$b("Gene expression matrix"),
      "- it is a",
      em("k x m"),
      "matrix with",
      em("k"),
      "rows of genes and",
      em("m"),
      "columns of samples.",
      "Each data point in this matrix represents the expression of a given gene in a given sample"
    ),
    tags$li(
      tags$b("Cell proportion matrix"),
      "- it is a",
      em("n x m"),
      "matrix with",
      em("n"),
      "rows of cell types and",
      em("m"),
      "columns of samples.",
      "Each data point in this matrix represents the proportion a given cell type in a given sample"
    ),
    tags$li(
      tags$b("Gene signature matrix"),
      "- it is a",
      em("k x n"),
      "matrix with",
      em("n"),
      "rows of genes and",
      em("m"),
      "columns of cell fraction",
      "Each data point in this matrix represents the contribution of a gene towards a cell type"
    ),
  ),
  
  h1("Input"),
  
  "The user simply has to upload the deconvolution algorithm outputs as plain comma separated file",
  em(".csv"),
  "file.",
  " For example, supported input files for -",
  
  tags$ol(
    tags$li(
      tags$b("Unsupervised algorithms"),
      "- Gene expression matrix (Determine number of components), Gene signature matrix (Cell type enrichment), Cell proportion matrix (For vizualisation)"
    ),
    tags$li(
      tags$b("Supervised algorithms"),
      "- Gene expression matrix (For deconvolution of immune cell types), Cell proportion matrix (For vizualisation)"
    ),
  ),
  
  h1("Output"),
  tags$ol(tags$li(
    tags$b("Unsupervised algorithms -"),
    tags$ol(
      tags$li(
        "Help for the determination of the number of cell types from the Gene expression matrix"
      ),
      tags$li(
        "Enrichment analysis of the Gene signature matrix to predict cell types"
      ),
      tags$li(
        "Visualization of the Cell proportion matrix as an interactive heatmap"
      ),
      tags$li("All numeric data to download as ", em(".csv"), "files")
    ),
  ),
  
  tags$li(
    tags$b("Supervised algorithms -"),
    tags$ol(
      tags$li(
        "Visualization of the Cell proportion matrix as an interactive heatmap"
      ),
      tags$li("All numeric data to download as a", em(".csv"), "files"),
    ),
  )),
  
  h1("CIBERSORT"),
  
  tags$a(href = "https://cibersort.stanford.edu/", "Download CIBERSORT.R and LM22.txt here"), tags$br(), "Do note registration is required before use !!",
  tags$h6("NOTE: Only for this workshop, a copy of these files has been made available along with the app"),
  tags$b("PLEASE DO NOT SHARE THESE FILES!!"),
  
  
  h1("Optimal k selection "),
  p(
    "This module provides guidance in the determination of the number of putative cell types from omics data (e.g. Gene expression matrix)."
  ),
  p(
    "Cattell’s rule is here suggested for choosing K [1]: it states that components corresponding to eigenvalues to the left of the straight line should be retained. When the actual number of different cell types is equal to K, we expect that there are (K-1) eigenvalues would correspond to the mixture of cell types and that other eigenvalues would correspond to noise (or other unaccounted for confounders). Indeed, one PCA axis is needed to separate two types, two axes for three types, etc. However, when not accounting for confounders, Cattell’s rule overestimates the number of PCs."
  ),
  
  h1("Enrichment analysis"),
  p(
    "This module helps in the biological interpretation of the components identified by unsupervised approaches. It proposes to perform Gene Set Enrichment analysis (fgsea R package [2]) on the cell-type gene marker database 'CellMatch' [3]. Do note that ICA-based methods are subjected to a reorientation of the components, as proposed in the deconica R package [4]."
  ),
  
  h1("References"),
  tags$a(href = "https://cancer-heterogeneity.github.io/cometh.html", "This app was developed within the COMETH program"),
  tags$br(),
  tags$br(),
  
  p(
    "[1] Cattell RB. The scree test for the number of factors. Multivariate Behav Res.2010;1:245–76.",
    tags$br(),
    "[2] Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi: 10.1101/060012.",
    tags$br(),
    "[3] Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, iScience, Volume 23, Issue 3, 27 March 2020. doi: 10.1016/j.isci.2020.100882" ,
    tags$br(),
    "[4] https://urszulaczerwinska.github.io/DeconICA/"
  ),
  br()
  
), #Tab - About

#-------------------------------------------------------------------------------
# Navbar : CSV-convertor
#-------------------------------------------------------------------------------

tabPanel(
  "CSV-convertor",

  sidebarPanel(
    fileInput(
      inputId = 'file05',
      label = h4('Your csv file'),
      accept = c('.csv')
    ),
    p(
      'Note that output from the cometh web-app are in the french-format (Separator = ";" Decimal = ",")'
    ),
    tags$hr(),
    radioButtons(
      inputId = 'sep05',
      label = 'Separator',
      choices = c(
        Semicolon = ';',
        Comma = ',',
        Tab = '\t',
        Space = ''
      ),
      selected = ';'
    ),
    br(),
    radioButtons(
      inputId = 'dec05',
      label = 'Decimal',
      choices = c(Comma = ',' , Dot = '.'),
      selected = ','
    ),
    br(),
    h4("Convert you csv file into:"),

    textInput(
      inputId = 'name05',
      label = 'Filename (without csv extension)',
      value = ""
    ),

    radioButtons(
      inputId = 'format05',
      label = 'Format',
      choices = c(English = 'English' , French =
                    'French'),
      selected = 'English'
    ),
    p('English -   Separator = "," Decimal = "."'),
    p('French -   Separator = ";" Decimal = ","'),

    br()
  ), #sidebarPanel

  mainPanel(
    p("Download your converted csv file"),
    downloadButton('downloadData1', 'Convert'),
    h5(
      'This module can be useful if your are, for instance, using an english version of excel: output of the cometh app
                           are csv files in the french format. Convert them into the english format will allow you to open them directly in excel.'
    )

  )
), # Tab - csv convertor

#-------------------------------------------------------------------------------
# Navbar : Optimal k - Deconvolution - Unsupervised
#-------------------------------------------------------------------------------

tabPanel(#"Deconvolution - Unsupervised",
  "Optimal k",
  
  sidebarPanel(
    fileInput(
      inputId = "file02",
      buttonLabel = "Upload",
      label = h4("Gene expression matrix"),
      accept = c(".csv")
    ),
    
    # actionButton('reset', 'Reset'),
    #
    # tags$hr(),
    #
    # sliderInput(
    #   "slider01",
    #   label = h4("Select appropriate", em("k")),
    #   min = 2,
    #   max = 15,
    #   step = 1,
    #   value = 5
    # ),
    #
    # tags$hr(),
    #
    # pickerInput(
    #   inputId = "select03",
    #   label = h4("Algorithm type"),
    #   choices = list(
    #     ICA = c("deconICA-1" = "MT1",
    #             "deconICA-2" = "MT14"),
    #
    #     NMF = c("ICA-nmf" = "MT2",
    #             "snmf" = "MT")
    #   ),
    # ),
    #
    # tags$hr(),
    #
    # actionButton("start", "Start deconvolution"),
  ), #Side bar panel
  
  mainPanel(
    fluidRow(
      plotOutput(outputId = "screeplot"),
      
      tags$br(),
      tags$br(),
      tags$br(),
      
      DT::DTOutput("cellpropmat", width = "60%"),
      
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      tags$br(),
      
      DT::DTOutput("genesigmat", width = "60%")
    )
  ),
), # Tab - Unsupervised deconvolution

#-------------------------------------------------------------------------------
# Navbar : Deconvolution - Supervised
#-------------------------------------------------------------------------------

tabPanel(
  #"Deconvolution - Supervised",
  "Deconvolution - Immune cells",
  
  
  sidebarPanel(
    
    # textInput("CBR", label = h4("CIBERSORT.R file"), value = "~/ShinyCompExplore/data/CIBERSORT.R"),
    # tags$hr(),
    # textInput(
    #   "LM22",
    #   label = h4("LM22 gene signature file"),
    #   value = "~/ShinyCompExplore/data/LM22.txt"
    # ),
    # tags$a(href = "https://cibersort.stanford.edu/", "Download CIBERSORT.R and LM22.txt here"), tags$br(), "Registration required.",
    # tags$h6(
    #   "NOTE: Only for this workshop, a copy has made available at - /ShinyCompExplore/data/CIBERSORT.R"
    # ),
    # tags$b("PLEASE DO NOT SHARE THESE FILES!!"),
    # 
    # tags$hr(),
    
    pickerInput(
      inputId = "select02",
      label = h4("Algorithm type"),
      choices = list(
        "CIBERSORT" = "cibersort",
        "CIBERSORT-ABS" = "cibersort_abs",
        "TIMER" = "timer",
        "MCPcounter" = "mcp_counter",
        "EPIC" = "epic",
        "xCell" = "xcell",
        "All" = "all"
      )
    ),
    
    tags$a(href = "https://cancer-heterogeneity.github.io/cometh_training.html", "Learn more about the methods"),
    
    tags$hr(),
    
    pickerInput(
      inputId = "select01",
      label = h4("Tumor indication"),
      
      choices = list(
        'None' = 'None',
        'Kidney Chromophobe (kich)' = 'kich',
        'Bladder Urothelial Carcinoma (blca)' = 'blca',
        'Breast invasive carcinoma (brca)' = 'brca',
        'Cervical squamous cell carcinoma and endocervical adenocarcinoma (cesc)' = 'cesc',
        'Glioblastoma multiforme (gbm)' = 'gbm',
        'Head and Neck squamous cell carcinoma (hnsc)' = 'hnsc',
        'Kidney renal papillary cell carcinoma (kirp)' = 'kirp',
        'Brain Lower Grade Glioma (lgg)' = 'lgg',
        'Liver hepatocellular carcinoma (lihc)' = 'lihc',
        'Lung adenocarcinoma (luad)' = 'luad',
        'Lung squamous cell carcinoma (lusc)' = 'lusc',
        'Prostate adenocarcinoma (prad)' = 'prad',
        'Sarcoma (sarc)' = 'sarc',
        'Pheochromocytoma and Paraganglioma (pcpg)' = 'pcpg',
        'Pancreatic adenocarcinoma (paad)' = 'paad',
        'Testicular Germ Cell Tumors (tgct)' = 'tgct',
        'Uterine Corpus Endometrial Carcinoma (ucec)' = 'ucec',
        'Ovarian serous cystadenocarcinoma (ov)'   = 'ov' ,
        'Skin Cutaneous Melanoma (skcm)' = 'skcm',
        'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (dlbc)' = 'dlbc',
        'Kidney renal clear cell carcinoma (kirc)' = 'kirc',
        'Adrenocortical carcinoma (acc)'  = 'acc',
        'Mesothelioma (meso)' = 'meso',
        'Thyroid carcinoma (thca)' = 'thca',
        'Uveal Melanoma (uvm)'  = 'uvm',
        'Uterine Carcinosarcoma (ucs)' = 'ucs',
        'Thymoma (thym)' = 'thym',
        'Esophageal carcinoma (esca)' = 'esca',
        'Stomach adenocarcinoma (stad)' = 'stad',
        'Rectum adenocarcinoma (read)' = 'read',
        'Colon adenocarcinoma (coad)' = 'coad',
        'Cholangiocarcinoma (chol)' = 'chol'
      )
    ),
    
    tags$b("Only needed if"),
    em("TIMER"),
    tags$b("or"),
    em("All"),
    tags$b("is selected as the algorithm type, else use the default"),
    em("None"),
    tags$b("option"),
    
    tags$br(),
    tags$br(),
    
    tags$a(href = "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations", "TCGA Tumor abbreviation source"),
    
    tags$hr(),
    
    fileInput(
      inputId = "file01",
      buttonLabel = "Upload",
      label = h4("Gene expression matrix"),
      accept = c(".csv")
    ),
    actionButton('reset', 'Reset')
  ),
  
  mainPanel(
    DT::DTOutput("supexpmat", width = "60%"),
    
    tags$br(), tags$br(), 
    
    h4("Download your deconvolution results as a .csv file"),
    downloadButton('downloadData2', 'Download')

    ),
  
), # Tab - Supervised deconvoliution

#-------------------------------------------------------------------------------
# Navbar : Visualization - Cell proportion matrix
#-------------------------------------------------------------------------------
             
tabPanel(
  "Cell proportions",
  
  sidebarPanel(
    fileInput(
      inputId = "file03",
      label = h4("Cell proportion matrix"),
      accept = c(".csv")
    ),
    
    tags$hr(),
    
    selectInput('samplist1', label = 'Select samples', choices = 'No choices here yet !!'),
    selectInput('celllist1', label = 'Select cell types', choices = 'No choices here yet !!'),
    
    tags$hr(),
    
    h4("NOTE - "),
    tags$br(),
    "(1)",
    tags$b("MCP-Counter, TIMER and xCell"),
    "can be used only for",
    tags$b("inter- sample comparisons"),
    "i.e. to compare one cell-type across multiple samples",
    
    tags$br(),
    tags$br(),
    
    "(2)",
    tags$b("CIBERSORT"),
    "can be used only for",
    tags$b("intra- sample comparisons"),
    "i.e. comparing different cell-types within each sample",
    
    tags$br(),
    tags$br(),
    
    "(3)",
    tags$b("CIBERSORT-ABS, EPIC and quanTIseq"),
    "can be used for both",
    tags$b("inter- and intra- sample comparisons"),
    "i.e. comparing one cell-type within one sample and across samples is possible"
  ),
  
  mainPanel(fluidRow(
    #h4("Heatmap - Celltypes vs Samples"),
    plotlyOutput(outputId = "viz1", height = "100%"),
    
    #h4("Cell types distribution in the selected sample"),
    plotlyOutput(outputId = "viz2", height = "100%"),
    
    #h4("Selected celltype proportion across samples"),
    plotlyOutput(outputId = "viz3", height = "100%"),
  )),
  
), # Visualization - Cell proportion matrix

#-------------------------------------------------------------------------------
# Navbar : Visualization - Enrichment
#-------------------------------------------------------------------------------

tabPanel(
  "Component enrichment",
  
  sidebarPanel(
    fileInput(
      inputId = "file04",
      label = h4("Gene signature matrix"),
      accept = c(".csv")
    ),
    
    radioButtons(
      inputId = "radio012",
      label = "Deconvolution method",
      choices = list("ICA-based" = 1,
                     "NMF-based" = 2)
    ),
    selectInput("select1", "Cancer type", c("ALL", cancertypes)),
    br(),
    h4("Download top markers"),
    p("Top100 gene markers for each component"),
    downloadButton('downloadData3', 'Download')
  ),
  
  mainPanel(fluidRow(plotOutput(outputId = "viz4"))),
) # Visualization - Enrichment

#-------------------------------------------------------------------------------
# End Navbar and fluidPage
#-------------------------------------------------------------------------------
  ) # NavBar
) # fluidPage




################################################################################
#----------------------------- UI.app end --------------------------------------
################################################################################





#-------------------------------------------------------------------------------
# SERVER.app
#-------------------------------------------------------------------------------

server <- function(input, output, session) {
  
options(shiny.maxRequestSize = 100*1024^2)
n_cores = detectCores()
if(n_cores > 7) {n_cores = 7}else{n_cores = c(n_cores - 1)}

#-------------------------------------------------------------------------------
# Render output csv convertor
#-------------------------------------------------------------------------------
   
inFile5 = reactive({
  req(input$file05)
  read.csv(
    input$file05$datapath,
    dec = input$dec05,
    sep = input$sep05,
    row.names = 1,
    header = TRUE
  )
})

output$downloadData1 <- downloadHandler(
  filename = function() {
    paste(input$name05, ".csv", sep = "")
  },
  content = function(file) {
   if (input$format05 == "English") {
        write.csv(inFile5(), file, row.names = TRUE)
      } else {
        write.csv2(inFile5(), file, row.names = TRUE)
      }
  }
)

#-------------------------------------------------------------------------------
# Unsupervised deconvolution : Optimal k : Read gene expression matrix
#-------------------------------------------------------------------------------

#----------------------------- REACTIVE values ---------------------------------

# Reactive for reset
rv2 <- reactiveValues(data = NULL)

# Input matrix
expmat_react2 <- reactiveVal()

observeEvent(input$file02, {
  #mytable <- data.frame(fread(input$file02$datapath, sep=";", dec="."), row.names = 1)
  mytable <- read.csv(input$file02$datapath, sep=";", dec=".", row.names = 1)
  rv2$data = mytable
  expmat_react2(mytable)
})

# # Number of cluster k choice
# kval_react <- reactiveVal()
# 
# observeEvent(input$slider01, {
#   kval_react(input$slider01)
# })
# 
# # Algorithm choice
# unsupalgo_react <- reactiveVal()
# 
# observeEvent(input$select03, {
#   unsupalgo_react(input$select03)
# })
# 
# # Supervised deconvolution output
# unsup_cellpropmat <- reactiveVal()
# unsup_genesigmat  <- reactiveVal()
# 
# Reset
observeEvent(input$reset, {
  rv2$data <- NULL
  reset('file02')
})

#---------------- RUN PCA for k val selection using Scree plot -----------------

# # Waiter waiting
# w_pca <- Waiter$new(html = tagList(
#   tags$strong(h1("Computing optimal k...")),
#   tags$br(),
#   html = spin_pulsar(),
#   tags$br(),
#   tags$strong(h1("Please wait for the results !!")),
#   tags$br()
# ))

observeEvent(input$file02, {
  output$screeplot <- renderPlot({
    mydata = expmat_react2()
    if (ncol(mydata) > 0)
    {
      mydata = mydata[order(matrixStats::rowMads(as.matrix(mydata)), decreasing = TRUE)[1:min(nrow(mydata), 10000)],]
      
      #w_pca$show()
      
      plot_k(mydata)
      
      #w_pca$hide()
    } else {
      text = c("Your input is empty. Check the separators in your csv file.")
      
    }
    
    # dat = as.matrix(t(expmat_react2()))
    # v = apply(dat, 2, var)
    # dat = dat[, v > quantile(v, 0.1)]
    # comps = 15
    # if(nrow(dat) < 15){comps = nrow(dat)}
    # p = prcomp(dat, center = T, scale = T)
    # stats::screeplot(p, npcs = comps, type = "lines", main = "")
    
    #shinyjs::reset('file02')
  })
})

#------------------------ Run UNSUPERVISED deconvolution -----------------------

# observeEvent(c(input$file02, input$slider01, input$select03, input$start), {
#
#   req(input$start)
#   req(input$file02)
#   req(input$slider01)
#   req(input$select03)
#
#   path = "/Users/ashwin/Documents/Projects/COMETH/src/compExplore/DeconCellTypesAPP/src/"
#   source(paste0(path,"Program_",input$select03(),"/program.R"))
#
# })

#------------------------ Render Cell proportion matrix ------------------------

# output$cellpropmat <- DT::renderDT({
#   datatable(
#     expmat_react2(),
#     extensions = 'Buttons',
#     class = "compact",
#     options = list(
#       dom = 'Bfrtip',
#       autoWidth = TRUE,
#       buttons = c('csv', 'excel')
#     )
#   )
# })

#------------------------- Render Gene signature matrix ------------------------

# output$genesigmat <- DT::renderDT({
#   datatable(
#     expmat_react2(),
#     extensions = 'Buttons',
#     class = "compact",
#     options = list(
#       dom = 'Bfrtip',
#       autoWidth = TRUE,
#       buttons = c('csv', 'excel')
#     )
#   )
# })

  
#-------------------------------------------------------------------------------
# Supervised deconvolution : Read gene expression matrix
#-------------------------------------------------------------------------------
  
#----------------------------- REACTIVE values ---------------------------------
  
# Reactive for filename reset
rv1 <- reactiveValues(data = NULL)

# Input matrix selection
expmat_react1 <- reactiveVal()

observeEvent(input$file01, {
  mytable <- data.frame(fread(input$file01$datapath, sep=";",  dec="."), row.names = 1)
  rv1$data = mytable
  expmat_react1(mytable)
})

# Algorithm choice
algo_react <- reactiveVal()

observeEvent(input$select02, {
  algo_react(input$select02)
})

# Tumor indication choice
canID_react <- reactiveVal()

observeEvent(input$select01, {
  canID_react(input$select01)
})

# Supervised deconvolution output
sup_decon_out <- reactiveVal(NULL)

# Filename reset
observeEvent(input$reset, {
  rv1$data <- NULL
  reset('file01')
})

# Waiter waiting
w <- Waiter$new(html = tagList(
  tags$strong(h1("Deconvolution in progress...")),
  tags$br(),
  html = spin_pulsar(),
  tags$br(),
  tags$strong(h1("Please wait for the results !!")),
  tags$br()
))

# Control flag to allow deconvolution processing
flag <- reactiveVal(NULL)

#----------------------------- SANITY checks -----------------------------------
  
  
observeEvent(c(algo_react(), canID_react()), {
  req(algo_react())
  req(canID_react())
  
  #----------------
  # Sanity checks 1
  #----------------
  
  if (is.element(algo_react(), c("timer", "all")) &
      canID_react() == "None")
  {
    sendSweetAlert(
      session = session,
      title = "Error: Tumor type for TIMER missing !!",
      text = "Please select the most matching tumor tissue type for your dataset",
      type = "error"
    )
    
    flag(NULL)
    
    #----------------
    # Sanity checks 2
    #----------------
    
  } else if (!is.element(algo_react(), c("timer", "all")) &
             !canID_react() == "None") {
    sendSweetAlert(
      session = session,
      title = "Error: Incorrect tumor indication !!",
      text = "Tumor type should be `None` when selecting methods other than TIMER or All !!",
      type = "error"
    )
    
    updatePickerInput(session = session,
                      inputId = "select01",
                      selected = "None")
    flag(NULL)
  } else {
    flag(TRUE)
  }
}, priority = 1000)

#------------------------ Run supervised deconvolution -------------------------
  
observeEvent(c(algo_react(), canID_react(), expmat_react1()), {
  # Control flag should be TRUE to start deconvolution !!
  req(flag())
  
  req(algo_react())
  req(canID_react())
  req(expmat_react1())
  
  #req(input$CBR)
  #req(input$LM22)
  
  w$show() # show the waiter
  
  # List of algorithms
  algo = c("quantiseq",
           "timer",
           "cibersort",
           "cibersort_abs",
           "mcp_counter",
           "xcell",
           "epic")
  
  set_cibersort_binary("/srv/shiny-server/data/CIBERSORT.R")
  set_cibersort_mat("/srv/shiny-server/data/LM22.txt")
  
  if (canID_react() == "None") {
    TIMER_type = NULL
  } else{
    TIMER_type = canID_react()
  }
  
  if (algo_react() == "all")
  {
    all_res = mclapply(algo, function(decon_method) {
      if (decon_method == "timer") {
        TIMER_type = rep(TIMER_type, ncol(expmat_react1()))
      }
      res = deconvolute(expmat_react1(), decon_method, indications = TIMER_type)
      res$cell_type = paste0(res$cell_type, "|", decon_method)
      return(res)
    }, mc.cores = n_cores)
    
    all_res = as.data.frame(do.call("rbind", all_res))
    
  } else{
    if (algo_react() == "timer") {
      TIMER_type = rep(TIMER_type, ncol(expmat_react1()))
    }
    all_res = deconvolute(expmat_react1(), algo_react(), indications = TIMER_type)
    all_res = as.data.frame(all_res)
  }
  
  all_res[, 2:ncol(all_res)] = signif(all_res[, 2:ncol(all_res)], 4)
  
  sup_decon_out(all_res)
  
  w$hide() # hide the waiter
  
  sendSweetAlert(
    session = session,
    title = "Deconvolution done !! " ,
    text = "Analysis complete",
    type = "info"
  )
  #shinyjs::reset('file01')
})

#------------------- DISPLAY supervised deconvolution table --------------------
  
output$supexpmat <- DT::renderDT({
  datatable(
    sup_decon_out(),
    extensions = 'Buttons',
    class = "compact",
    options = list(
      dom = 'Bfrtip',
      autoWidth = TRUE
      #buttons = c('csv', 'excel')
    )
  )
})

output$downloadData2 <- downloadHandler(
  filename = function() {
    paste("supervised_deconvolution_results_", algo_react() , ".csv", sep="")
  },
  content = function(file) {
      write.csv2(sup_decon_out(), file)
    }
) 

#-------------------------------------------------------------------------------
# Visualization : Cell proportions: Read sample names/ cell types
#-------------------------------------------------------------------------------
  
##-- Cell proportion matrix --##
inFile3 = reactiveVal(NULL)

# Update the input options based on the .csv file
observeEvent(input$file03, {
  mytable <- data.frame(fread(input$file03$datapath, sep=";",  dec=","))
  
  if(is.element("cell_type", colnames(mytable))){
    rownames(mytable) = mytable$cell_type
    mytable <- mytable[, 3:ncol(mytable)]
  }else{
    rownames(mytable) = mytable[,1]
    mytable = mytable[, -1]
  }
  
  updateSelectInput(session,
                    "samplist1",
                    label = "Select samples",
                    choices = colnames(mytable))
  
  updateSelectInput(session,
                    "celllist1",
                    label = "Select cell types",
                    choices = rownames(mytable))
  
  inFile3(mytable)
  
}, priority = 1000)

#-------------------------------------------------------------------------------
# Visualization : Cell proportions
#-------------------------------------------------------------------------------

# HEATMAP
observeEvent(c(inFile3(), input$samplist1, input$celllist1), {
  output$viz1 = renderPlotly({
    req(input$samplist1 != "No choices here yet !!")
    req(input$celllist1 != "No choices here yet !!")
    req(input$samplist1 %in% colnames(inFile3()))
    
    dat = as.matrix(inFile3())
    dat = log2(dat + 1)
    
    hm = pheatmap(dat, clustering_method = "ward.D2", silent = T)
    dat = dat[hm$tree_row$order, hm$tree_col$order]
    
    plot_ly(
      x = colnames(dat),
      y = rownames(dat),
      z = signif(dat, 4),
      type = "heatmap",
      height = 1000
    ) %>% layout(margin = c(0, 0, 100, 100), title = "Heatmap - Cell vs Samples")
    
  })
})

# Barplot: cells
observeEvent(c(inFile3(), input$samplist1, input$celllist1), {
  output$viz2 = renderPlotly({
    req(input$samplist1 != "No choices here yet !!")
    req(input$celllist1 != "No choices here yet !!")
    req(input$samplist1 %in% colnames(inFile3()))
    
    samp1 = inFile3()[, input$samplist1, drop = F]
    samp1 = log2(samp1 + 1)
    
    plot_ly(
      x = rownames(samp1),
      y = samp1[, 1],
      type = "bar",
      name = "Cells",
      height = 1000
    ) %>% layout(
      title = "Cell types distribution in the selected sample",
      yaxis = list(hoverformat = ".2f"),
      margin = c(0, 0, 100, 100)
    )
    
  })
})

# Barplot: samples
observeEvent(c(inFile3(), input$samplist1, input$celllist1), {
  output$viz3 = renderPlotly({
    req(input$samplist1 != "No choices here yet !!")
    req(input$celllist1 != "No choices here yet !!")
    req(input$samplist1 %in% colnames(inFile3()))
    
    cell1 = inFile3()[input$celllist1, , drop = F]
    cell1 = log2(cell1 + 1)
    
    plot_ly(
      x = colnames(cell1),
      y = as.numeric(cell1[1,]),
      type = "bar",
      name = "Samples",
      height = 1000
    ) %>% layout(
      title = "Selected celltype proportion across samples",
      yaxis = list(hoverformat = ".2f"),
      margin = c(0, 0, 100, 100)
    )
  })
})

#-------------------------------------------------------------------------------
# Unsupervised deconvolution - Visualization : Enrichment of cell types
#-------------------------------------------------------------------------------

# w_enr <- Waiter$new(html = tagList(
#   tags$strong(h1("Enrichment analysis in progress...")),
#   tags$br(),
#   html = spin_pulsar(),
#   tags$br(),
#   tags$strong(h1("Please wait for the results !!")),
#   tags$br()
# ))

inFile4 = reactiveVal()

# Read the gene signature matrix
observeEvent(input$file04, {
  gsig <- data.frame(fread(input$file04$datapath, sep=";",  dec=","), row.names = 1)
  inFile4(gsig)
})

observeEvent(c(inFile4(), input$radio012, input$select1), {
  output$viz4 = renderPlot({
    req(inFile4())
    req(input$radio012)
    req(input$select1)
    
    dat = as.matrix(inFile4())
    
    if (input$radio012 == "ICA-based") {
      datemp = orient_funct(dat)
      colnames(datemp) = colnames(dat)
      dat = datemp
    }
    
    if (input$select1 == "ALL") {
      markerscelltypes = tapply(dbmarkers$geneSymbol, dbmarkers$cellIDcancer, cbind)
    } else {
      dbmarkers = dbmarkers[dbmarkers$cancerType == input$select1, ]
      markerscelltypes = tapply(dbmarkers$geneSymbol, dbmarkers$cellID, cbind)
    }
    
    #w_enr$show()
    enrichplot(dat, markerscelltypes, showCategory = 10)
    #w_enr$hide()
  })
})

#
#
# inFile4 = reactive({
#   req(input$file04)
#   data.frame(fread(input$file04$datapath), row.names = 1)
# })
#
# output$viz4 = renderPlot({
#
#   if(! is.null(inFile4()) ){
#     dat = as.matrix(inFile4())
#
#     if(input$radio012 == "ICA-based"){
#       datemp = orient_funct(dat)
#       colnames(datemp) =  colnames(dat)
#       dat = datemp
#     }
#     if(input$select1 == "ALL"){
#       markerscelltypes = tapply(dbmarkers$geneSymbol,dbmarkers$cellIDcancer, cbind)
#     } else {
#       dbmarkers. = dbmarkers[dbmarkers$cancerType == input$select1,]
#       markerscelltypes = tapply(dbmarkers.$geneSymbol,dbmarkers.$cellID, cbind)
#     }
#
#     #w_enr$show()
#
#     enrichplot(dat, markerscelltypes, showCategory = 10)
#
#     #w_enr$hide()
#   }
#
#
# })

topgen = reactive({
  if (!is.null(inFile4())) {
    dat = as.matrix(inFile4())
    
    if (input$radio012 == "ICA-based") {
      datemp = orient_funct(dat)
      colnames(datemp) =  colnames(dat)
      dat = datemp
    }
    apply(dat, 2, function(x) {
      rownames(dat)[order(x, decreasing = TRUE)[1:100]]
    })
  }
  
})

output$downloadData3 <- downloadHandler(
  filename = function() {
    paste("results", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    write.csv(topgen(), file)
  }
)
} # Closing Server.app

################################################################################
#----------------------------- Server.app end ----------------------------------
################################################################################

#-------------------------------------------------------------------------------
# Create Shiny app
#-------------------------------------------------------------------------------
shinyApp(ui, server)
