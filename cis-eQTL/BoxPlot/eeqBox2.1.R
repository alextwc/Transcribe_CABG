eqBox2 <- function(gene, se, tf, snpgr, genome = "hg19", forceRs = TRUE, 
    ...)
{
    ans = gQTLstats:::prepEqData(gene, se, tf, snpgr, genome, forceRs = forceRs)
    bb  = boxplot(split(ans$ex, ans$gt), plot=FALSE)
    beeswarm::beeswarm (split(ans$ex, ans$gt), pch = 21, col = 2:4, bg = "#00000050",
    #xlab = ans$coln, ylab = gene, 
    ...)
    bxp(bb, add=TRUE, boxfill="transparent")
}
