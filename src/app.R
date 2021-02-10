library(shiny)
library(shinythemes)
library(pheatmap)
library(plotly)
library(cowplot)
library(ggplot2)
library(fgsea)
library(matrixStats)
library(nFactors)
#---------------------------
# Internal input for Unsupervised-enrichment 
#---------------------------

dbmarkers = readRDS("data/CellMatchCurated.RDS")
cancertypes = levels(as.factor(dbmarkers$cancerType))

#---------------------------
# Functions for Unspervised-enrichment 
#---------------------------

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


#---------------------------
# Functions for Unspervised-Number of components 
#---------------------------

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
  svd. = svd(t(D))$d[1:min(10, ncol(D))]
  nbf = findElbowPoint_PCAtools(svd.)
  ggplot(data.frame(val = svd., idx = 1:min(10, ncol(D))), aes(x = idx, y = val)) +
    geom_line() + geom_vline(xintercept=nbf, linetype="dashed", color = "red", size=2) +
    geom_point() + scale_x_continuous(breaks = 1:min(10, ncol(D)))+
    labs(x = "PC index", y = "Eigenvalues") +
    theme_minimal(base_size = 16)
}

#------------------------------------------------------#
#                   UI  
#-------------------------------------------------------#


# Define UI for data upload app ----
ui <- fluidPage(
  #shinythemes::themeSelector(),
  navbarPage("compExplorer",
             theme = shinytheme("united"),
             
             tabPanel("About",
                      
                      h1("compExplore"),
                      
                      tags$b("compExplore - Components explorer"),
                      "is a visualization tool to guide the user in the analysis and interpretation of the results from",
                      em("Supervised"),
                      "and",
                      em("Unsupervised"),
                      "gene expresssion deconvolution algorithms.",
                      br(),
                      
                      img(src = "deconvolution.jpg", width = 800),
                      
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
                        )
                      ),
                      
                      h1("Input"),
                      
                      "The user simply has to upload the deconvolution algorithm outputs as plain comma separated file",
                      em(".csv"),
                      "files",
                      " For example, supported input files for -",
                      
                      tags$ol(
                        tags$li(
                          tags$b("Unsupervised algorithms"),
                          "- Gene expression matrix (Determine number of components), Gene signature matrix (Cell type enrichment), Cell proportion matrix (For vizualisation) "
                        ),
                        tags$li(tags$b("Supervised algorithms"), "- Cell proportion matrix (For vizualisation)")
                      ),
                      
                      h1("Output"),
                      tags$ol(
                        tags$li(
                          tags$b("Unsupervised algorithms -"),
                          tags$ol(
                            tags$li(
                              "Help for the determination of the number of cell types from the Gene expression matrix."
                            ),
                            tags$li(
                              "Enrichment analysis of the Gene signature matrix to predict cell types."
                            ),
                            tags$li(
                              "Visualization of the Cell proportion matrix as an interactive heatmap."
                            ),
                            tags$li("All numeric data to download as ", em(".csv"), "files")
                          )
                        ),
                        tags$li(
                          tags$b("Supervised algorithms -"),
                          tags$ol(
                            tags$li(
                              "Visualization of the Cell proportion matrix as an interactive heatmap"
                            ),
                            tags$li("All numeric data to download as ", em(".csv"), "files")
                          )
                        )
                        
                      ),
                      
                      h1("References"),
                      tags$a(href="https://cancer-heterogeneity.github.io/cometh.html", "This is a contribution of the COMETH program."),
                      br()
                      
                      
             ), #Tab - About
             
             
             tabPanel("Unsupervised - number",
                      
                      sidebarPanel(
                        fileInput(inputId = 'file04', label = h4('Gene expression matrix (csv file)'),
                                  accept=c('.csv')),
                        tags$hr(),
                        radioButtons(inputId = 'sep', label = 'Separator', 
                                     choices = c(Comma=',' ,Semicolon=';'
                                                 ,Tab='\t', Space=''
                                     ), selected = ','),
                        br(),
                        h5("Only the 10,000 most variables features are kept for the analysis ")
                      ),#sidePanel
                      
                      mainPanel(
                        tabsetPanel(
                          
                          # tabPanel("Data table", value=3,tableOutput("table")),
                          tabPanel("Guide", value=2,
                                   br(),
                                   p("This module provides guidance in the determination of the number of putative cell types from omics data (e.g. Gene expression matrix)."),
                                   p("Cattell’s rule is here suggested for choosing K [1]: it states that components corresponding to eigenvalues to the left of the straight line should be retained. When the actual number of different cell types is equal to K, we expect that there are (K-1) 
                                     eigenvalues would correspond to the mixture of cell types and that other eigenvalues would correspond to noise (or other unaccounted for confounders). Indeed, one PCA axis is needed to separate two types, two axes for three types, etc. However, 
                                     when not accounting for confounders, Cattell’s rule overestimates the number of PCs."),
                                   br(),
                                   h4("References"),
                                   p("[1] Cattell RB. The scree test for the number of factors. Multivariate Behav Res.
                                     2010;1:245–76.")
                                   #br(),
                                   
                                   #p("This is part of the COMETH training course")
                          ),
                          tabPanel("Plot of eigenvalues",br(), span(textOutput("viz4warnings"), style="color:red"), br(), plotOutput("viz4")),
                          
                          id="tabselect"
                        )
                      )#main panel
                      
             ), # Tab - Unsupervised - number
             
            
             
             
             
              tabPanel("Unsupervised - enrich",
                      
                      sidebarPanel(
                        fileInput(inputId = "file03",
                                  label = h4("Gene signature matrix (csv file)"),
                                  accept = c(".csv")),
                        
                        radioButtons(
                          inputId = "radio012",
                          label = "Deconvolution method",
                          choices = list("ICA-based" = 1,
                                         "NMF-based" = 2)),
                        selectInput("select1", "Cancer type", c("ALL",cancertypes)),
                        br(),
                        h4("Download top markers"),
                        p("Top100 gene markers for each component"),
                        downloadButton('downloadData', 'Download')
                      ),
                      
                      #mainPanel(tags$b("Gene set enrichment analysis results using CellMatch DB"), plotOutput(outputId = "viz3"))
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Guide", value=2,
                                   br(),
                                   p("This module helps in the biological interpretation of the components identified by unsupervised approaches."),
                                   p("It proposes to perform Gene Set Enrichment analysis (fgsea R package [1]) on the cell-type gene marker database 'CellMatch' [2]", br(),
                                     "To be noted: ICA-based methods are subjected to a reorientation of the components, as proposed in the deconica R package [3]."),
                                   br(),
                                   h4("References"),
                                   p("[1] Korotkevich G, Sukhov V, Sergushichev A (2019). “Fast gene set enrichment analysis.” bioRxiv. doi: 10.1101/060012.", tags$br(),
                                     "[2] Shao et al., scCATCH:Automatic Annotation on Cell Types of Clusters from Single-Cell RNA Sequencing Data, iScience, Volume 23, Issue 3, 27 March 2020. doi: 10.1016/j.isci.2020.100882" , tags$br(),
                                     "[3] https://urszulaczerwinska.github.io/DeconICA/")
                          ),
                          tabPanel("Enrichment analysis",br(), tags$b("Gene set enrichment analysis results using CellMatch DB"), plotOutput(outputId = "viz3")),
                          
                          id="tabselect"
                        )
                      )#main panel
                      
                      
             ), # End Tab - Unsupervised - enrich
             
             tabPanel("Unsupervised - vizu",
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        fileInput(inputId = "file02",
                                  label = h4("Cell proportion matrix (csv file)"),
                                  accept = c(".csv")),
                        
                        tags$hr(),
                        
                        selectInput('samplist2', label = 'Select samples', choices = 'No choices here yet !!'),
                        selectInput('celllist2', label = 'Select cell types', choices = 'No choices here yet !!'),
                        
                        tags$hr()),
                      mainPanel(plotlyOutput(outputId = "viz2"))
                      
                      
                      
                      
             ), # Tab - Unsupervised - vizu
             
             tabPanel("Supervised - vizu",
                      
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        radioButtons(
                          inputId = "radio011",
                          label = h4("Immune cells - algorithm type"),
                          choices = list("CIBERSORT" = 1,
                                         "CIBERSORT-ABS" = 2,
                                         "TIMER" = 3,
                                         "MCPcounter" = 4,
                                         "EPIC" = 5,
                                         "xCell" = 6)
                        ),
                        tags$hr(),
                        
                        fileInput(inputId = "file01",
                                  label = h4("Cell proportion matrix (csv file)"),
                                  accept = c(".csv")
                        ),
                        
                        tags$hr(),
                        
                        selectInput('samplist1', label = 'Select samples', choices = 'No choices here yet !!'),
                        selectInput('celllist1', label = 'Select cell types', choices = 'No choices here yet !!')
                        
                      ),
                      
                      mainPanel(plotlyOutput(outputId = "viz1"))
                      
             ) # Tab Supervised
             
  ) # NavBar
) # fluidPage



#------------------------------------------------------#
#                   server  
#-------------------------------------------------------#

# Define server logic to read selected file ----
server <- function(input, output, session) {
  
  #-------------------------
  # Read input - Supervised
  #-------------------------
  
  ##-- Cell proportion matrix --##
  
  # Read the input .csv file
  inFile1 = reactive({
    req(input$file01)
    read.csv(input$file01$datapath, row.names = 1)
  })
  
  # Update the input options based on the .csv file
  observeEvent(input$file01, {
    mytable <- read.csv(input$file01$datapath, row.names = 1)
    updateSelectInput(session, "samplist1", label = "Select samples", choices = colnames(mytable))
    updateSelectInput(session, "celllist1", label = "Select cell types", choices = rownames(mytable))
  })
  
  # Catch user selection of sample
  samp1 = reactive({
    inFile1()[,input$samplist1, drop=F]
    #mytable[,input$samplist1, drop=F]
  })
  
  # Catch user selection of cell type
  cell1 =  reactive({
    inFile1()[input$celllist1,, drop=F]
    #mytable[input$celllist1,, drop=F]
  })
  
  #---------------------------
  # Read input - Unsupervised
  #---------------------------
  
  
  
  ##--  Cell proportion matrix --##
  inFile2 = reactive({
    req(input$file02)
    read.csv(input$file02$datapath, row.names = 1)
  })
  
  # Update the input options based on the .csv file
  observeEvent(input$file02, {
    mytable <- read.csv(input$file02$datapath, row.names = 1)
    updateSelectInput(session, "samplist2", label = "Select samples", choices = colnames(mytable))
    updateSelectInput(session, "celllist2", label = "Select cell types", choices = rownames(mytable))
  })
  
  # Catch user selection of sample
  samp2 = reactive({
    inFile2()[,input$samplist2, drop=F]
  })
  
  # Catch user selection of cell type
  cell2 =  reactive({
    inFile2()[input$celllist2,, drop=F]
  })
  
  #--  Gene signature matrix --##
  
  inFile3 = reactive({
    req(input$file03)
    read.csv(input$file03$datapath, row.names = 1)
  })
  
  #--  Gene expression matrix --##
  
  inFile4 = reactive({
    req(input$file04)
  })
  
  
  
  #---------------------------
  # Render output - Supervised
  #---------------------------
  
  output$viz1 = renderPlotly({
    
    dat = as.matrix(inFile1())
    
    hm = pheatmap(dat, clustering_method = "ward.D2", silent = T)
    dat = dat[hm$tree_row$order, hm$tree_col$order]
    
    p1 = plot_ly(x = colnames(dat), y = rownames(dat), z = round(dat, 2), type = "heatmap", height = 1000)
    p2 = plot_ly(x = rownames(samp1()), y = samp1()[,1], type="bar", name = "Cells", height = 1000) %>% layout(yaxis = list(hoverformat = ".2f"))
    p3 = plot_ly(x = colnames(cell1()), y = as.numeric(cell1()[1,]), type="bar", name = "Samples", height = 1000) %>% layout(yaxis = list(hoverformat = ".2f"))
    subplot(p1, p2, p3, nrows = 3, margin = c(0,0,0.07,0.07))
    
  })
  
  #-----------------------------
  # Render output - Unsupervised -vizualisation
  #-----------------------------
  
  output$viz2 = renderPlotly({
    
    if(! is.null(inFile2()) ){
      
      dat = as.matrix(inFile2())
      
      hm = pheatmap(dat, clustering_method = "ward.D2", silent = T)
      dat = dat[hm$tree_row$order, hm$tree_col$order]
      
      p1 = plot_ly(x = colnames(dat), y = rownames(dat), z = round(dat, 2), type = "heatmap", height = 1000)
      p2 = plot_ly(x = rownames(samp2()), y = samp2()[,1], type="bar", name = "Cells", height = 1000) %>% layout(yaxis = list(hoverformat = ".2f"))
      p3 = plot_ly(x = colnames(cell2()), y = as.numeric(cell2()[1,]), type="bar", name = "Samples", height = 1000) %>% layout(yaxis = list(hoverformat = ".2f"))
      subplot(p1, p2, p3, nrows = 3, margin = c(0,0,0.07,0.07))
    }
  })
  
  #------------------------------------------
  # Render output - Unsupervised - Enrichment
  #------------------------------------------
  #markerscelltypes = readRDS("./data/markerscelltypes.RDS")
  
  
  output$viz3 = renderPlot({
    if(! is.null(inFile3()) ){
      dat = as.matrix(inFile3())
      
      if(input$radio012 == "ICA-based"){
        datemp = orient_funct(dat)
        colnames(datemp) =  colnames(dat)
        dat = datemp
      }
      if(input$select1 == "ALL"){
        markerscelltypes = tapply(dbmarkers$geneSymbol,dbmarkers$cellIDcancer, cbind)
      } else {
        dbmarkers. = dbmarkers[dbmarkers$cancerType == input$select1,]
        markerscelltypes = tapply(dbmarkers.$geneSymbol,dbmarkers.$cellID, cbind)
      }
      
      enrichplot(dat, markerscelltypes, showCategory = 10)
    }
  })
  
  csvfile03 = reactive({
    if(! is.null(inFile3()) ){
      dat = as.matrix(inFile3())
      
      if(input$radio012 == "ICA-based"){
        datemp = orient_funct(dat)
        colnames(datemp) =  colnames(dat)
        dat = datemp
      }
      apply(dat, 2, function(x){rownames(dat)[order(x, decreasing = TRUE)[1:100]]})
    }
    
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("results",Sys.Date(),".csv", sep="")
    },
    content = function(file) {
      write.csv(csvfile03(), file)
    }
  )
  
  
  #  thedata <- reactive(iris)
  #  output$dto <- renderDataTable({thedata()})
  #  output$downloadData <- downloadHandler(
  #    filename = function(){"thename.csv"}, 
  #    content = function(fname){
  #      write.csv(thedata(), fname)
  #    }
  #  )
  
  #------------------------------------------
  # Render output - Unsupervised - Number
  #------------------------------------------
  #markerscelltypes = readRDS("./data/markerscelltypes.RDS")
  #--  Gene expression matrix --##
  output$viz4 = renderPlot({
    
    # input$file1 will be NULL initially. After the user selects and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can be found.
    
    inFile4 <- input$file04
    if (is.null(inFile4)) return(NULL)
    mydata <- read.csv(inFile4$datapath, row.names=1, header=T, sep=input$sep)
    if(ncol(mydata)>0)
    {
      mydata = mydata[order(matrixStats::rowMads(as.matrix(mydata)), decreasing = TRUE)[1:min(nrow(mydata), 10000)],]
      plot_k(mydata) } else {
        text = c("Your input is empty. Check the separators in your csv file.")
      }
    # package call
    # download results
    
  })
  
  output$viz4warnings <- renderText({
    inFile4 <- input$file04
    if (is.null(inFile4)) return(NULL)
    mydata <- read.csv(inFile4$datapath, row.names=1, header=T, sep=input$sep)
    if(ncol(mydata) == 0){
      print("Your input is empty. Check the separators in your csv file.")
    }
    
  })
  
}
# Create Shiny app ----
shinyApp(ui, server)
