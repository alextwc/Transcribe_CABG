library(shiny)
require("ggplot2")
require("reshape2")
require("graphics")
load("./GWAS-CVD.Intersection.rda") 

# GWAS_CVD_UnPruned = data.frame(1:10, 1:10)
# myobj = "mytext"
# mytext = data.frame(1:10, 1:10)

shinyServer(function(input, output, session) {
      
      GWAS_Data = reactive({
          get(input$GWAS)
      })
      
      eQTLs_Data <- reactive({
          get(input$eQTLs)
      })
  
    # output$test1 = renderText(get(myobj))
    # output$test2 = renderText(ncol(GWAS_Data))
    # output$test3 = renderText(input$eQTLs)
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
  frac <- as.data.frame(frac)
  frac$pValue <- pValue
  print(dim(frac))
  frac
      })
      
  output$plot1 <- renderPlot({
   # hist(frac()$frac, breaks=1000, plot=TRUE)
   p <- ggplot(data=frac(), aes(frac())) + geom_histogram(binwidth=0.0003000, aes(x=frac()$frac, fill=..count..), col=I("blue"))
   p <- p + scale_fill_gradient("Counts", low = "cyan", high = "blue")
   p <- p + ggtitle("Percentage Distribution of CVD associated GWAS SNPs in Transcribe cis-eQTLs") + xlab("Percentage of cis-eQTLs in 1000 sets of matched GWAS_CVD SNPs") + ylab("Counts")
   p <- p + theme(axis.text=element_text(size=14), axis.title=element_text(size=12,face="bold"), axis.title.x=element_text(size=12,face="bold"), axis.title.y=element_text(size=12,face="bold"))
   p <- p + geom_vline(xintercept = frac()[1,1], color="red")
   p
  })
})
