library(shiny)
require("ggplot2")
require("reshape2")
require("graphics")
load("./GWAS-CVD.Intersection.rda") 

# Define server logic for random distribution application
function(input, output) {
  
  # Reactive expression to generate the requested distribution.
  # This is called whenever the inputs change. The output
  # functions defined below then all use the value computed from
  # this expression
  GWAS_Data  <- reactive({
          get(input$GWAS)
      })
      
  eQTLs_Data <- reactive({
          get(input$eQTLs)
      })
      
  N          <- reactive({
  	      input$n
      })
      
  frac <- reactive({
      frac <- vector()
      print(length(frac))
      print(ncol(GWAS_Data()))
      print(length(GWAS_Data()[,1]))
      print(ncol(eQTLs_Data()))
         for (i in 1:ncol(GWAS_Data())){
           frac[i] <- length(levels(factor(GWAS_Data()[(GWAS_Data()[,i] %in% eQTLs_Data()$SNPsnap),i])))/length(levels(factor(GWAS_Data()[,i])))
           # print(eQTLs_Data()$SNPsnap )
         }
      pValue <- sum(frac >= frac[1])/length(frac)
      frac   <- as.data.frame(frac)
      frac$pValue <- pValue
      print(dim(frac))
      frac
  })
  # Generate a plot of the data. Also uses the inputs to build
  # the plot label. Note that the dependencies on both the inputs
  # and the data reactive expression are both tracked, and
  # all expressions are called in the sequence implied by the
  # dependency graph
  output$plot <- renderPlot({
   # hist(frac()$frac, breaks=1000, plot=TRUE)
   p <- ggplot(data=frac(), aes(frac())) + geom_histogram(binwidth=input$n, aes(x=frac()$frac, fill=..count..), col=I("blue"))
   p <- p + scale_fill_gradient("Counts", low = "cyan", high = "blue")
   p <- p + ggtitle("Percentage Distribution of CVD associated GWAS SNPs in Transcribe cis-eQTLs") + xlab("Percentage of cis-eQTLs in 1000 sets of matched GWAS_CVD SNPs") + ylab("Counts")
   p <- p + theme(axis.text=element_text(size=14), axis.title=element_text(size=12,face="bold"), axis.title.x=element_text(size=12,face="bold"), axis.title.y=element_text(size=12,face="bold"))
   p <- p + geom_vline(xintercept = frac()[1,1], color="red")
   p
  })
  
  # Generate a summary of the data
  output$summary <- renderPrint({
    summary(frac()$frac)
  })
  
  # Generate an HTML table view of the data
  output$table <- renderTable({
    data.frame(x=frac()$frac)
  })
}
