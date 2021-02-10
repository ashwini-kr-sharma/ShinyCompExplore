#-------------------------------------------------------------------------------
# Required packages
#-------------------------------------------------------------------------------

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
library(deconica)
library(fgsea)
library(ica)
library(fastICA)
library(NMF)
library(CellMix)

library(cowplot)
library(matrixStats)
library(nFactors)

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

navbarPage("compExplorer",
theme = shinytheme("united"),
   
#-------------------------------------------------------------------------------
# Navbar : About
#-------------------------------------------------------------------------------
        
tabPanel(
  "About",
  
  h1("compExplorer"),
  
  tags$b("compExplore - Components explorer"),
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
  
  h1("References"),
  tags$a(href = "https://cancer-heterogeneity.github.io/cometh.html", "This is a contribution of the COMETH program"),
  br()
  
), #Tab - About

#-------------------------------------------------------------------------------
# Navbar : Deconvolution - Supervised
#-------------------------------------------------------------------------------

tabPanel(
  "Deconvolution - Supervised",
  sidebarPanel(
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
  
  mainPanel(DT::DTOutput("supexpmat", width = "60%")),
  
), # Tab - Supervised deconvoliution

#-------------------------------------------------------------------------------
# Navbar : Deconvolution - Unsupervised
#-------------------------------------------------------------------------------
    
tabPanel(
  "Deconvolution - Unsupervised",
  
  sidebarPanel(
    fileInput(
      inputId = "file02",
      buttonLabel = "Upload",
      label = h4("Gene expression matrix"),
      accept = c(".csv")
    ),
    actionButton('reset', 'Reset'),
    
    tags$hr(),
    
    sliderInput(
      "slider01",
      label = h4("Select appropriate", em("k")),
      min = 2,
      max = 15,
      step = 1,
      value = 5
    ),
    
    tags$hr(),
    
    pickerInput(
      inputId = "select03",
      label = h4("Algorithm type"),
      choices = list(
        ICA = c("deconICA-1" = "MT1",
                "deconICA-2" = "MT14"),
        
        NMF = c("ICA-nmf" = "MT2",
                "snmf" = "MT")
      ),
    ),
    
    tags$hr(),
    
    actionButton("start", "Start deconvolution"),
  ),
  
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
# Navbar : Visualization - Cell proportion matrix
#-------------------------------------------------------------------------------

tabPanel("Cell proportions",
    
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

tabPanel("Component enrichment",
  
  sidebarPanel(
    fileInput(
      inputId = "file04",
      label = h4("Gene signature matrix (csv file)"),
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
    downloadButton('downloadData', 'Download')
  ),
  
  mainPanel(fluidRow(
    plotlyOutput(outputId = "viz4", inline = T)
  )),
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
# Supervised deconvolution : Read gene expression matrix
#-------------------------------------------------------------------------------

#----------------------------- REACTIVE values ---------------------------------

# Reactive for filename reset
rv1 <- reactiveValues(data = NULL)

# Input matrix selection
expmat_react1 <- reactiveVal()

observeEvent(input$file01, {
  mytable <- read.csv(input$file01$datapath, row.names = 1)
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
  
  w$show() # show the waiter
  
  # List of algorithms
  algo = c("quantiseq",
           "timer",
           "cibersort",
           "cibersort_abs",
           "mcp_counter",
           "xcell",
           "epic")
  
  set_cibersort_binary(
    "/Users/ashwin/Documents/Projects/COMETH/src/compExplore/DeconCellTypesAPP/data/CIBERSORT.R"
  )
  set_cibersort_mat(
    "/Users/ashwin/Documents/Projects/COMETH/src/compExplore/DeconCellTypesAPP/data/LM22.txt"
  )
  
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
  
  shinyjs::reset('file01')

})

#------------------- DISPLAY supervised deconvolution table --------------------

output$supexpmat <- DT::renderDT({
  datatable(
    sup_decon_out(),
    extensions = 'Buttons',
    class = "compact",
    options = list(
      dom = 'Bfrtip',
      autoWidth = TRUE,
      buttons = c('csv', 'excel')
    )
  )
})

#-------------------------------------------------------------------------------
# Unsupervised deconvolution : Read gene expression matrix
#-------------------------------------------------------------------------------
      
#----------------------------- REACTIVE values ---------------------------------

# Reactive for reset
rv2 <- reactiveValues(data = NULL)

# Input matrix
expmat_react2 <- reactiveVal()

observeEvent(input$file02, {
  mytable <- read.csv(input$file02$datapath, row.names = 1, header = T)
  rv2$data = mytable
  expmat_react2(mytable)
})

# Number of cluster k choice
kval_react <- reactiveVal()

observeEvent(input$slider01, {
  kval_react(input$slider01)
})

# Algorithm choice
unsupalgo_react <- reactiveVal()

observeEvent(input$select03, {
  unsupalgo_react(input$select03)
})

# Supervised deconvolution output
unsup_cellpropmat <- reactiveVal()
unsup_genesigmat  <- reactiveVal()

# Reset
observeEvent(input$reset, {
  rv2$data <- NULL
  reset('file02')
})

#---------------- RUN PCA for k val selection using Scree plot -----------------

observeEvent(input$file02, {
  output$screeplot <- renderPlot({
    
    mydata = expmat_react2()
    if(ncol(mydata)>0)
    {
      mydata = mydata[order(matrixStats::rowMads(as.matrix(mydata)), decreasing = TRUE)[1:min(nrow(mydata), 10000)],]
      plot_k(mydata) } else {
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
}, priority = 1000)

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

output$cellpropmat <- DT::renderDT({
  datatable(
    expmat_react2(),
    extensions = 'Buttons',
    class = "compact",
    options = list(
      dom = 'Bfrtip',
      autoWidth = TRUE,
      buttons = c('csv', 'excel')
    )
  )
})

#------------------------- Render Gene signature matrix ------------------------

output$genesigmat <- DT::renderDT({
  datatable(
    expmat_react2(),
    extensions = 'Buttons',
    class = "compact",
    options = list(
      dom = 'Bfrtip',
      autoWidth = TRUE,
      buttons = c('csv', 'excel')
    )
  )
})

#-------------------------------------------------------------------------------
# Visualization : Cell proportions: Read sample names/ cell types
#-------------------------------------------------------------------------------
      
##-- Cell proportion matrix --##
inFile3 = reactiveVal(NULL)

# Update the input options based on the .csv file
observeEvent(input$file03, {
  mytable <- read.csv(input$file03$datapath, row.names = "cell_type")
  mytable <- mytable[, 2:ncol(mytable)]
  
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
      height = 1000) %>% layout(margin = c(0,0,100,100), title = "Heatmap - Cell vs Samples")
    
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
    ) %>% layout(title = "Cell types distribution in the selected sample",
                 yaxis = list(hoverformat = ".2f"), 
                 margin = c(0,0,100,100))
    
  })
})

# Barplot: samples
observeEvent(c(inFile3(), input$samplist1, input$celllist1), {
  output$viz3 = renderPlotly({
    req(input$samplist1 != "No choices here yet !!")
    req(input$celllist1 != "No choices here yet !!")
    req(input$samplist1 %in% colnames(inFile3()))
    
    cell1 = inFile3()[input$celllist1,, drop = F]
    cell1 = log2(cell1 + 1)
    
   plot_ly(
      x = colnames(cell1),
      y = as.numeric(cell1[1,]),
      type = "bar",
      name = "Samples",
      height = 1000
    ) %>% layout(title = "Selected celltype proportion across samples",
                 yaxis = list(hoverformat = ".2f"), 
                 margin = c(0,0,100,100))
  })
})

#-------------------------------------------------------------------------------
# Unsupervised deconvolution - Visualization : Enrichment of cell types
#-------------------------------------------------------------------------------

inFile4 = reactiveVal()

# Read the gene signature matrix
observeEvent(input$file04, {
  gsig <- read.csv(input$file04$datapath, row.names = 1)
  inFile4(gsig)
}, priority = 1000)

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
    enrichplot(dat, markerscelltypes, showCategory = 10)
  })
})

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

output$downloadData <- downloadHandler(
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
