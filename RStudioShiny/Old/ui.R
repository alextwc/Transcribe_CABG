library(shiny)
require("ggplot2")
require("reshape2")
require("graphics")

shinyUI(pageWithSidebar(
  headerPanel('Permuting CVD associated GWAS SNPs to intersect with Transcribe cis-eQTLs'),
  sidebarPanel(
    selectInput('GWAS', 'GWAS_CVD_SNPs', c("GWAS_CVD_UnPruned","GWAS_CVD_Pruned")),
    selectInput('eQTLs', 'Transcribe_eQTLs', c("Baseline","Ischemia"),
                selected="Baseline")
   # numericInput('clusters', 'Pruning Distance', 500000,
   #              min = 500000, max = 1000000)
  ),
  mainPanel(
    plotOutput("plot1")
      # verbatimTextOutput("test1"),
      # verbatimTextOutput("test2")
      # verbatimTextOutput("test3")
      
  )
))
