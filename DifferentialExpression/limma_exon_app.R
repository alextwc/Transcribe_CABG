library(limma)
library(shiny)
library(DT)
library(ggplot2)
library(qqman)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(magrittr)
source('functions.R')

analyses <- make_analyses_list()
hgnc_table <- read_tsv('data/exon_to_gene_table.tsv')


#------------------------------------------------------------------------------

server <- function(input, output) {
    
    # reactives ---------------------------------------------------------------
    
    # general
    analysis_folder  <- reactive(get_analysis_folder(analyses, input$analysis))
    limma_object     <- reactive(read_rds_object(analysis_folder(), '/limma_object.rds'))
    limma_covariates <- reactive(get_limma_covariates(limma_object()))
    pheno_covariates <- reactive(get_pheno_covariates(analysis_folder()))
    limma_table      <- reactive(make_limma_exon_table(limma_object(), input$covariate))
    pheno_vector     <- reactive(make_pheno_vector(analysis_folder(), input$covariate2))
    pval_order_table <- reactive(select(limma_table(), ExonID, FDR))
    

    # tab 1
    filtered_table   <- reactive(filter_limma_exon_table(limma_table(), input$max_FDR))
    
    # tab 3
    norm_count_table  <- reactive(make_exon_count_table(analysis_folder(), '/norm_count_table.tsv', pval_order_table()))
    raw_count_table   <- reactive(make_exon_raw_count_table(analysis_folder(), '/raw_count_table.tsv', norm_count_table(), pval_order_table()))
    norm_count_vector <- reactive(make_count_vector(norm_count_table(), input$raw_count_table_row_last_clicked)) 
    raw_count_vector  <- reactive(make_count_vector(raw_count_table(),  input$raw_count_table_row_last_clicked))

    # tab 3
    MDS_object       <- reactive(read_rds_object(analysis_folder(), '/MDS_object.rds'))

    # output ------------------------------------------------------------------
    
    # sidebar
    output$covariates       <- renderUI(selectInput(
        "covariate", 
        label = h3("Choose analysis covariate"), 
        choices = limma_covariates()))
    
    output$covariates2       <- renderUI(selectInput(
        "covariate2", 
        label = h3("Choose grouping covariate"), 
        choices = pheno_covariates())) 
    
    # tab 1
    output$filtered_table   <- renderDataTable(filtered_table(), selection = 'single')
    output$exon_plot        <- renderPlot(create_exon_plot(limma_object(), input$covariate, input$filtered_table_row_last_clicked, filtered_table()))

    # tab 2
    output$qq_plot          <- renderPlot(qq(limma_table()$P.Value))
    output$histogram        <- renderPlot(create_histogram(limma_table()$P.Value))
    
    
    # tab 3
    output$raw_count_table  <- renderDataTable(raw_count_table(),  selection = 'single')
    output$norm_count_table <- renderDataTable(norm_count_table(), selection = 'single')
    output$raw_count        <- renderPlot(create_count_plot(pheno_vector(), raw_count_vector(), input$covariate2))
    output$norm_count       <- renderPlot(create_count_plot(pheno_vector(), norm_count_vector(), input$covariate2))

    # tab 4
    output$MDS_plot         <- renderPlot(create_mds_plot(MDS_object(), pheno_vector()))

    # tab 5
    output$MV_plot          <- renderImage(list(src = str_c(analysis_folder(), '/MV_plot.png')), deleteFile = FALSE)
    

}

#------------------------------------------------------------------------------
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput("analysis", 
                        label = h3("Choose analysis"),
                        choices = analyses$ui_name),
            conditionalPanel(condition = "input.conditionedPanels == 'Results table' | input.conditionedPanels == 'QQ/Histogram' | input.conditionedPanels == 'Count tables'" ,
                             uiOutput("covariates")),
            conditionalPanel(condition = "input.conditionedPanels == 'Count tables' | input.conditionedPanels == 'MDS plot'",
                             uiOutput("covariates2")),
            conditionalPanel(condition = "input.conditionedPanels == 'Results table' ",
                             numericInput("max_FDR", 
                                          label = h3("Choose max FDR value"), 
                                          value = .05,
                                          step  = .05)),
            conditionalPanel(condition = "input.conditionedPanels == 'Count tables'",
                             plotOutput("raw_count"),
                             plotOutput("norm_count"))),
        
        mainPanel(
            tabsetPanel(
                tabPanel("Results table",
                         dataTableOutput('filtered_table'),
                         plotOutput("exon_plot")),
                tabPanel("QQ/Histogram",
                         h4('QQ plot'),
                         plotOutput('qq_plot'),
                         h4('P value histogram'),
                         plotOutput('histogram')),
                tabPanel("Count tables",
                         h4("Raw counts"),
                         dataTableOutput('raw_count_table'),
                         h4("Normalized counts"),
                         dataTableOutput('norm_count_table')),
                tabPanel("MDS plot", plotOutput("MDS_plot")),
                tabPanel("MV plot", imageOutput("MV_plot")),
                id = "conditionedPanels"
            )
        
        
        
        )))
            


shinyApp(ui = ui, server = server)