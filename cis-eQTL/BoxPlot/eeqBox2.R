
eqBox2 <- function(gene, se, tf, snpgr, genome = "hg19", forceRs = TRUE, 
    ...)
{
    ans = gQTLstats:::prepEqData(gene, se, tf, snpgr, genome, forceRs = forceRs)
    ans$sex <- se@colData@listData$Sex
    bb  = boxplot(split(ans$ex, ans$gt), plot=FALSE)
    beeswarm::beeswarm (split(ans$ex, ans$gt), 
    pwcol = c('deeppink', 'blue')[as.numeric(factor(ans$sex))], pch = 16,
    #xlab = ans$coln, ylab = gene, 
    ...)
    bxp(bb, add=TRUE, boxfill="transparent")
}
