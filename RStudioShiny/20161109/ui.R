library(shiny)
require("ggplot2")
require("reshape2")
require("graphics")

fluidPage(
    
  # Application title
  titlePanel("Permutation Test"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  sidebarLayout(
    sidebarPanel(
      radioButtons("GWAS", "GWAS_CVD_SNPs type:",
                   c("Pruned SNPs by 500Kb" = "GWAS_CVD_Pruned",
                     "Un-pruned SNPs" = "GWAS_CVD_UnPruned")),
      br(),
      
      radioButtons("eQTLs", "Transcribe eQTLs:",
                   c("Transcribe Baseline" = "Baseline",
                     "Transcribe Ischemia" = "Ischemia")),
      br(),
      
      
      sliderInput("n", 
                  "Bin Width:", 
                   value = .0003,
                   min = .0001, 
                   max = .0009)
    ),
    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
        tabPanel("Histogram Plot", plotOutput("plot")), 
        tabPanel("Summary", verbatimTextOutput("summary")), 
        tabPanel("Table", tableOutput("table"))
      )
    )
  )
)
