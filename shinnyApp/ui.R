source("common.R")
headerInfo = tags$h3("Microbiome")
footerInfo = tags$small(tags$a("Visualization code",href=""))

taxlist <- function(x){
  as.character(na.omit(unique(x)))
}

shinyUI(navbarPage("Microbiome",
  tabPanel("Bacterial Abundance",
    tags$h6("Manhattan plots of bacterial abundances where each bar represents the abundance in a sample."),
        sidebarLayout(position="right",
          sidebarPanel(
            selectInput("level","Level:",c(
                        "Phylum" = "Phylum",
                        "Class" = "Class","Order"="Order","Family"="Family",
                        "Genus" = "Genus","Species"="Species","OTU"="OTU"),
                        selected = c("Genus")),
              conditionalPanel(condition = "input.level == 'OTU'",
                # checkboxInput("otu","Index",TRUE),
                numericInput('feature', 'OTU index:', 1, min = 1, max = nrow(MRobj),value=1)
              ),
              conditionalPanel(condition = "input.level == 'Species'",
                selectInput('Species', 'Bacteria:',taxlist(fData(MRobj)[,"Species"]))
              ),
              conditionalPanel(condition = "input.level == 'Genus'",
                selectInput('Genus', 'Bacteria:',taxlist(fData(MRobj)[,"Genus"]))
              ),
              conditionalPanel(condition = "input.level == 'Family'",
                selectInput('Family', 'Bacteria:',taxlist(fData(MRobj)[,"Family"]))
              ),
              conditionalPanel(condition = "input.level == 'Order'",
                selectInput('Order', 'Bacteria:',taxlist(fData(MRobj)[,"Order"]))
              ),
              conditionalPanel(condition = "input.level == 'Class'",
                selectInput('Class', 'Bacteria:',taxlist(fData(MRobj)[,"Class"]))
              ),
              conditionalPanel(condition = "input.level == 'Phylum'",
                selectInput('Phylum', 'Bacteria:',taxlist(fData(MRobj)[,"Phylum"]))
              ),
            selectInput("pd","Phenotype color:",c(colnames(pData(MRobj)),"None"),selected="Type"),
            checkboxInput("norm", "Normalize (CSS-normalization)", TRUE),
            checkboxInput("interactive","Interactive plot",TRUE),
            conditionalPanel(condition = "input.interactive",
              checkboxInput("box", "Boxplots?", FALSE)
            )
          ),
          mainPanel(
            plotlyOutput("plotly1"),
            conditionalPanel(condition = "input.level == 'OTU'",
              dataTableOutput("table")
              # ,pre("OTU ID and sequence center for OTU above:",textOutput("clusterSequence"))
            )
            ,conditionalPanel(condition = "input.level != 'OTU'",
              pre("Prevalent (>20 samples) OTU ID and sequence centers for bacteria above:",textOutput("clusterSequences"))
            )
          )
      )
    ),
  tabPanel("Heatmap",
      tags$h6("Heatmap of the top N bacteria."),
      sidebarLayout(position="left",sidebarPanel(
        selectInput('heatmap_level','Tree level',c('Phylum' = 'Phylum',
                        'Class' = 'Class','Order'='Order','Family'='Family',
                        'Genus' = 'Genus','Species'='Species','OTU'='OTU')
                        ,selected='OTU'),
        selectInput("hpd","Phenotype columns:",colnames(pData(MRobj))),
        numericInput('heatNumber', 'Number of bacteria (rows) to display:', 150,
                 min = 5, max = 400),
        radioButtons("heat","Choose bacteria by:",c("Variability"="sd","Median Absolute Deviation"="mad")),
        br(),
        tags$small("Warning: takes a few seconds")
        ),
      mainPanel(plotOutput("plotHeatmap",height="800px")))
      )
  ,tabPanel("PCA / MDS",
    tags$h6("PCA or MDS to project samples, using the 200 most variable bacteria (OTUs) onto 2 dimensions."),
        sidebarLayout(position="right",
          sidebarPanel(
            radioButtons("pcaOrMds","PCA or MDS:",c("PCA"="TRUE","MDS"="FALSE"),selected="TRUE"),
            radioButtons("useDist","Make use of count distances:",c("False"="FALSE","True"="TRUE"),selected="TRUE"),
              conditionalPanel(condition = "input.useDist == 'TRUE'",
                selectInput("distance", "Distance:", 
                    choices=c("euclidean","manhattan","canberra","bray",
                      "kulczynski","jaccard","gower","altGower","morisita",
                      "horn","raup","binomial","chao","cao"),selected="bray")
              ),
            selectInput("ppd","Phenotype color:",colnames(pData(MRobj)),selected="Type"),
            numericInput('dimensionx', 'X-axis dimension:', 1,
                 min = 1, max = 4),
            numericInput('dimensiony', 'Y-axis dimension:', 2,
                 min = 1, max = 4),
            br(),
            tags$small("Warning: takes a few seconds")
          ),
          mainPanel(
            plotlyOutput("pcaPlot",height="650px")
          )
      )
    )
  ,tabPanel("Sequencing",
    tags$h6("Number of taxa per sample."),
    sidebarLayout(position="right",sidebarPanel(
      selectInput("spd","Phenotype:",colnames(pData(MRobj)),selected="Type")
      ),
      mainPanel(plotlyOutput("sparsity"),
                plotlyOutput("ntaxa",height="500px"),
                plotOutput("sparsityPCA"))
      )
    )
  ,tabPanel("Diversity",
      tags$h6("Boxplots of diversity indexes."),
      sidebarLayout(position="right",sidebarPanel(
        selectInput("diversity","Diversity index:",choices=c("shannon", "simpson", "invsimpson")),
        selectInput("dpd","Phenotype:",colnames(pData(MRobj)),selected="Type")
        ),
        mainPanel(plotlyOutput("diversity"),tableOutput("diversityTable")))
      )
  ,tabPanel("Rarefaction",
    tags$h6("The linear effect depth of coverage has on the number of bacteria detected."),
      sidebarLayout(position="right",
        sidebarPanel(
          selectInput("rpd","Phenotype color:",colnames(pData(MRobj)))
          ),
        mainPanel(
          plotlyOutput("plotRare",height="600px",width="600px")
          )
      )
    )
  ,tabPanel("Phenotype information",
        mainPanel(dataTableOutput("phenolist"))
      )
  ,tabPanel("OTU information",
        mainPanel(dataTableOutput("featurelist"))
      )
))
