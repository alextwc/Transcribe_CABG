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
hgnc_table <- read_tsv('data/ensemble_to_hgnc_gene_table.tsv')

#------------------------------------------------------------------------------

server <- function(input, output) {
    
    # reactives ---------------------------------------------------------------
    
    # general
    analysis_folder  <- reactive(get_analysis_folder(analyses, input$analysis))
    limma_object     <- reactive(read_rds_object(analysis_folder(), '/limma_object.rds'))
    limma_covariates <- reactive(get_limma_covariates(limma_object()))
    pheno_covariates <- reactive(get_pheno_covariates(analysis_folder()))
    limma_table      <- reactive(make_limma_table(limma_object(), input$covariate, hgnc_table))
    pheno_vector     <- reactive(make_pheno_vector(analysis_folder(), input$covariate2))
    pval_order_table <- reactive(select(limma_table(), ensemble, P.Value))

    # tab 1
    filtered_table   <- reactive(filter_limma_table(limma_table(), input$max_pvalue))
    
    # tab 3
    MDS_object       <- reactive(read_rds_object(analysis_folder(), '/MDS_object.rds'))
    
    # tab 4
    norm_count_table  <- reactive(make_count_table(analysis_folder(), '/norm_count_table.tsv', pval_order_table()))
    raw_count_table   <- reactive(make_raw_count_table(analysis_folder(), '/raw_count_table.tsv', norm_count_table(), pval_order_table()))
    norm_count_vector <- reactive(make_count_vector(norm_count_table(), input$raw_count_table_row_last_clicked)) 
    raw_count_vector  <- reactive(make_count_vector(raw_count_table(),  input$raw_count_table_row_last_clicked))
    
    # tab 6 
    go_table <- reactive(create_go_table(analysis_folder(), input$db, input$ontology, input$min_genes, input$max_genes, input$go_filter, input$min_go_level, input$max_go_level))
    
    # tab 7
    gage_table <- reactive(create_gage_table(analysis_folder(), input$gage_table, input$gage_min_genes, input$gage_max_genes))
    
    # tab 8
    roast_table <- reactive(create_roast_table(analysis_folder(), input$roast_db, input$roast_min_genes, input$roast_max_genes))
    
    # tab 9
    pheno_table <- reactive(read_tsv(str_c(analysis_folder(), '/pheno_file.tsv')))
    sva_object  <- reactive(read_rds_object(analysis_folder(), '/sva_object.rds'))
    
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
    
    output$surrogates  <- renderUI(radioButtons(
        "surrogate", 
        label = h3("Choose surrogate variable"), 
        choices = 1:sva_object()$n.sv))
    
    # tab 1
    output$filtered_table   <- renderDataTable(filtered_table())
    output$qq_plot          <- renderPlot(qq(limma_table()$P.Value))
    output$histogram        <- renderPlot(create_histogram(limma_table()$P.Value))
    
    # tab 2
    output$MA_plot          <- renderPlot(create_MA_plot(limma_table(), input$max_pvalue))
    
    # tab 3
    output$raw_count_table  <- renderDataTable(raw_count_table(),  selection = 'single')
    output$norm_count_table <- renderDataTable(norm_count_table(), selection = 'single')
    output$raw_count        <- renderPlot(create_count_plot(pheno_vector(), raw_count_vector(), input$covariate2))
    output$norm_count       <- renderPlot(create_count_plot(pheno_vector(), norm_count_vector(), input$covariate2))
    
    # tab 4
    output$MDS_plot         <- renderPlot(create_mds_plot(MDS_object(), pheno_vector()))
    
    # tab 5
    output$MV_plot          <- renderImage(list(src = str_c(analysis_folder(), '/MV_plot.png')), deleteFile = FALSE)
    
    # tab 6
    output$go_table        <- renderDataTable(go_table(), selection = 'single')
    
    # tab 7
    output$gage_table        <- renderDataTable(gage_table())
    
    # tab 8
    output$roast_table        <- renderDataTable(roast_table())
    
    # tab 9
    output$sva_plot        <- renderPlot(create_sva_plot(pheno_table(), sva_object(), input$covariate2, input$surrogate))
}

#------------------------------------------------------------------------------
ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput("analysis", 
                        label = h3("Choose analysis"),
                        choices = analyses$ui_name),
            conditionalPanel(condition = "input.conditionedPanels == 'Results table' | input.conditionedPanels == 'MA plot' | input.conditionedPanels == 'QQ/Histogram' | input.conditionedPanels == 'Count tables'" ,
                             uiOutput("covariates")),
            conditionalPanel(condition = "input.conditionedPanels == 'Count tables' | input.conditionedPanels == 'MDS plot' | input.conditionedPanels == 'SVA plot'",
                             uiOutput("covariates2")),
            conditionalPanel(condition = "input.conditionedPanels == 'Results table' | input.conditionedPanels == 'MA plot'",
                             numericInput("max_pvalue", 
                                          label = h3("Choose max padj value"), 
                                          value = .05,
                                          step  = .05)),
            conditionalPanel(condition = "input.conditionedPanels == 'Count tables'",
                             plotOutput("raw_count"),
                             plotOutput("norm_count")),
            conditionalPanel(condition = "input.conditionedPanels == 'GO table'",
                             selectInput(
                                 "ontology", 
                                 label = h3("Choose ontology"), 
                                 choices = c('all', 'MF', 'BP', 'CC')),
                             selectInput(
                                 "db", 
                                 label = h3("Choose DB"), 
                                 choices = c('all', 'KEGG', 'GO')),
                             numericInput(
                                 "min_genes", 
                                 label = h3("Min genes"), 
                                 value = 1), 
                             numericInput(
                                 "max_genes", 
                                 label = h3("Max genes"), 
                                 value = 20000),
                             h3("Filter by GO level"),
                             checkboxInput("go_filter", label = "filter", value = F),
                             numericInput(
                                 "min_go_level",
                                 label = h3("Min GO level"),
                                 value = 1),
                             numericInput(
                                 "max_go_level",
                                 label = h3("Max GO level"),
                                 value = 20)),
            conditionalPanel(condition = "input.conditionedPanels == 'Gage table'",
                             selectInput(
                                 "gage_table",
                                 label = h3("Choose table"),
                                 choices = c('kegg_up', 'kegg_down', 'go_up', 'go_down')),
                             numericInput(
                                 "gage_min_genes",
                                 label = h3("Min genes"),
                                 value = 1),
                             numericInput(
                                 "gage_max_genes",
                                 label = h3("Max genes"),
                                 value = 20000)),
            conditionalPanel(condition = "input.conditionedPanels == 'Roast table'",
                             selectInput(
                                 "roast_db", 
                                 label = h3("Choose DB"), 
                                 choices = c('all', 'KEGG', 'GO')),
                             numericInput(
                                 "roast_min_genes",
                                 label = h3("Min genes"),
                                 value = 1),
                             numericInput(
                                 "roast_max_genes",
                                 label = h3("Max genes"),
                                 value = 5000)),
            conditionalPanel(condition = "input.conditionedPanels == 'SVA plot'",
                             uiOutput("surrogates"))),
        mainPanel(
            tabsetPanel(
                tabPanel("Results table",dataTableOutput('filtered_table')),
                tabPanel("MA plot", plotOutput("MA_plot")),
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
                tabPanel("GO table", dataTableOutput('go_table')),
                tabPanel("Gage table", dataTableOutput('gage_table')),
                tabPanel("Roast table", dataTableOutput('roast_table')),
                tabPanel("SVA plot", plotOutput('sva_plot')),
                id = "conditionedPanels"
            ))))



shinyApp(ui = ui, server = server)