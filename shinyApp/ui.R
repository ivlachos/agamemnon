source("common.R")
headerInfo = tags$h3("AGAMEMNON v0.1.0")
footerInfo = tags$small(tags$a("Visualization code", href = ""))

taxlist <- function(x)
{
  as.character(na.omit(unique(x)))
}

shinyUI(navbarPage("AGAMEMNON v0.1.0",
                   
  tabPanel("Microbial Abundances",
    tags$h6("Manhattan plots of microbial abundances where each bar represents the abundance in a sample."),
        sidebarLayout(position = "right",
          sidebarPanel(
            selectInput("level", "Taxonomic rank:", 
                        c("Phylum" = "Phylum", "Class" = "Class", "Order" = "Order", "Family" = "Family", 
                        "Genus" = "Genus", "Species" = "Species", "Scientific_Name" = "Scientific_Name"),
                        selected = c("Genus")),
            
              conditionalPanel(condition = "input.level == 'Scientific_Name'",
                selectInput('Scientific_Name', 'Microorganism:', taxlist(fData(MRobj)[, "Scientific_Name"]))),
            
              conditionalPanel(condition = "input.level == 'Species'",
                selectInput('Species', 'Microorganism:', taxlist(fData(MRobj)[, "Species"]))),
            
              conditionalPanel(condition = "input.level == 'Genus'",
                selectInput('Genus', 'Microorganism:', taxlist(fData(MRobj)[, "Genus"]))),
            
              conditionalPanel(condition = "input.level == 'Family'",
                selectInput('Family', 'Microorganism:', taxlist(fData(MRobj)[, "Family"]))),
            
              conditionalPanel(condition = "input.level == 'Order'",
                selectInput('Order', 'Microorganism:', taxlist(fData(MRobj)[, "Order"]))),
            
              conditionalPanel(condition = "input.level == 'Class'",
                selectInput('Class', 'Microorganism:', taxlist(fData(MRobj)[, "Class"]))),
            
              conditionalPanel(condition = "input.level == 'Phylum'",
                selectInput('Phylum', 'Microorganism:', taxlist(fData(MRobj)[, "Phylum"]))),
            
            selectInput("pd", "Phenotype column:", c(colnames(pData(MRobj)), "None"), selected = "Type"),
            
            checkboxInput("norm", "Normalize (CSS-normalization)", TRUE),
            
            checkboxInput("interactive", "Interactive plot", TRUE),
            
            conditionalPanel(condition = "input.interactive", checkboxInput("box", "Boxplots", FALSE))
          ),
          
          mainPanel(
            plotlyOutput("plotly1")
          )
      )
    )
  
  ,tabPanel("Heatmap",
      tags$h6("Heatmap of the top N microorganisms."),
      sidebarLayout(position = "left", sidebarPanel(
        selectInput('heatmap_level', 'Taxonomic rank', c('Phylum' = 'Phylum', 'Class' = 'Class', 'Order' = 'Order', 'Family' = 'Family',
                        'Genus' = 'Genus', 'Species' = 'Species', 'Scientific_Name' = 'Scientific_Name'), selected = 'Species'),
        
        selectInput("hpd", "Phenotype column:", colnames(pData(MRobj))),
        
        numericInput('heatNumber', 'Number of OTUs (rows) to display:', 100, min = 1, max = 1500),
        
        radioButtons("heat", "Choose microbes by:", c("Variability" = "sd", "Median Absolute Deviation" = "mad")),
        
        br(),
        
        tags$small("Warning: takes a few seconds")
        ),
        
      mainPanel(plotOutput("plotHeatmap", height = "800px")))
      )
  
  ,tabPanel("PCA / MDS",
    tags$h6("PCA or MDS to project samples, using the 200 most variable OTUs (Scientific_Names) onto 2 dimensions."),
        sidebarLayout(position = "right",
          sidebarPanel(
            radioButtons("pcaOrMds", "PCA or MDS:", c("PCA" = "TRUE", "MDS" = "FALSE"), selected = "TRUE"),
            
            radioButtons("useDist", "Make use of count distances:", c("False" = "FALSE", "True" = "TRUE"), selected = "TRUE"),
            
              conditionalPanel(condition = "input.useDist == 'TRUE'",
                selectInput("distance", "Distance:", 
                    choices = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", 
                      "horn", "raup", "binomial", "chao", "cao"), selected = "bray")
              ),
            
            selectInput("ppd", "Phenotype column:", colnames(pData(MRobj)), selected = "Type"),
            
            numericInput('dimensionx', 'X-axis dimension:', 1, min = 1, max = 4),
            
            numericInput('dimensiony', 'Y-axis dimension:', 2, min = 1, max = 4),
            
            br(),
            
            tags$small("Warning: takes a few seconds")
          ),
          
          mainPanel(
            plotlyOutput("pcaPlot", height = "650px")
          )
      )
    )
  
  ,tabPanel("Number of Taxa",
    tags$h6("Number of different taxa per sample."),
    sidebarLayout(position = "right", sidebarPanel(
      
      selectInput("spd", "Phenotype column:", colnames(pData(MRobj)), selected = "Type")),
      
      mainPanel(plotlyOutput("sparsity"), plotlyOutput("ntaxa",height="500px"), plotOutput("sparsityPCA"))
      )
    )
  
  
  ,tabPanel("Diversity",
      tags$h6("Boxplots of diversity indexes."),
      sidebarLayout(position = "right", sidebarPanel(
        selectInput("diversity", "Diversity index:", choices = c("shannon", "simpson", "invsimpson")),
        
        selectInput("dpd", "Phenotype column:", colnames(pData(MRobj)), selected = "Type")),
        
        mainPanel(plotlyOutput("diversity"), tableOutput("diversityTable")))
      )
  
  ,tabPanel("DE Analysis",
      tags$h6("DE Analysis between selected groups and Taxonomy tree level"),
      sidebarLayout(position = "left", 
        sidebarPanel(
         selectInput("phen", "Phenotype:", c(colnames(pData(MRobj)), "None"), selected = "None"),
         
         conditionalPanel(
           condition = "input.phen != 'None'",
              uiOutput("group1"),
              uiOutput("group2"),
              uiOutput("demethod"),
              actionButton("btn1", "Run DE Analysis")),
         conditionalPanel(
           condition = "input.btn1 > 0",
           br(),
           downloadButton("downloadData", "Export results to CSV"))),

        conditionalPanel(
           condition = "input.btn1 > 0",
              mainPanel(dataTableOutput("deanalysis")))
        )
      )
    
  
  ,tabPanel("Phenotype information", mainPanel(dataTableOutput("phenolist")))
  
  ,tabPanel("Reference Lineage", mainPanel(dataTableOutput("featurelist")))
  
))
