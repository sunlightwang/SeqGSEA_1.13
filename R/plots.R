# copyright: Xi Wang (xi.wang@newcastle.edu.au)
# plots for SeqGSEA

plotGeneScore <- function(score, perm.score=NULL, pdf=NULL, 
                          main=c("Overall", "Expression", "Splicing")) {
  main <- match.arg(main, c("Overall", "Expression", "Splicing"))
  main <- paste("Gene List", main, "Scoring Profile", sep=" ")
  ngene <- length(score)
  score.sorted <- sort(score, decreasing=TRUE)
  score.max <- max(score)
  score.avg <- mean(score)
  plot.mat <- matrix(rep(c(score.max, score.avg), each=2), 2)
  leg.txt <- c(paste("max. score:", signif(score.max, 5)), 
               paste("avg. score:", signif(score.avg, 5)))
  if(! is.null(perm.score) ) {
    score.perm.max <- max(perm.score)
    score.perm.avg <- mean(perm.score)
    plot.mat <- cbind(plot.mat, matrix(rep(c(score.perm.max, score.perm.avg), each=2), 2))
    leg.txt <- c(leg.txt,
                 paste("max. perm. score:", signif(score.perm.max, 5)),
                 paste("avg. perm. score:", signif(score.perm.avg, 5)))
  }
  loc <- 1:ngene
  if(! is.null(pdf)) 
    pdf(file = pdf, height = 6, width = 6)
  plot(loc, score.sorted, ylab = "Normlized Scores", xlab = "Gene List Location", 
       main = main, type = "l", lwd = 2, cex = 0.9, col = 1)
  for (i in seq(1, ngene, ceiling(ngene/500))) {
    lines(c(i, i), c(0, score.sorted[i]), lwd = 3, cex = 0.9, col = colors()[12])
  }
  points(loc, score.sorted, type = "l", lwd = 2, cex = 0.9, col = 1)
  lines(c(1, ngene), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1)
  #text(x=0.5 * ngene, y=0.95*score.max, adj = c(0, 1), labels=leg.txt, cex = 0.9)
  cols <- c(colors()[552], colors()[256], colors()[498], colors()[28])
  matlines(range(loc), plot.mat, type="l", lty=2, lwd=1, col=cols)
  legend(x=0.4 * ngene, y=0.95*score.max, bty="n", bg = "white", legend=leg.txt, 
         lty = 2, lwd = 1, col = cols, cex = 0.9)
  if(! is.null(pdf)) 
    dev.off()
}

plotES <- function(gene.set, pdf=NULL) {
  stopifnot( is( gene.set, "SeqGeneSet" ) )
  stopifnot( ncol(gene.set@GSEA.ES.perm) > 0 )
  stopifnot( gene.set@GSEA.normFlag )
  adjust.param <- 0.5
  ES.perm.density <- apply(gene.set@GSEA.ES.perm, 2, function(vec){
    d <- density(vec, adjust=adjust.param, n = 512, from=0, to=3.5)
    d$y / sum( d$y )
  } )
  ES.perm.density.mean <- rowMeans(ES.perm.density)
  d <- density(gene.set@GSEA.ES, adjust=adjust.param, n = 512, from=0, to=3.5)
  ES.density <- d$y / sum(d$y)
  x <- d$x
  y1 <- ES.perm.density.mean
  y2 <- ES.density
  x.plot.range <- range(x)
  y.plot.range <- c(-0.3 * max(c(y1,y2)), max(c(y1,y2)))
  if(! is.null(pdf)) 
    pdf(file = pdf, height = 6, width = 6)
  plot(x, y1, xlim = x.plot.range, ylim = 1.2*y.plot.range, type = "l", 
       lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", main = "Global Observed and Null Densities (Area Normalized)")
  temp <- apply(ES.perm.density, 2, function(vec) {
    points(x, vec, type = "l", lwd = 1, col = "#FFAA0011")
  })
  points(x, y1, type = "l", lwd = 2, col = colors()[555])
  points(x, y2, type = "l", lwd = 2, col = colors()[29])
  nset <- size(gene.set)
  for (i in 1:nset) {
    lines(c(gene.set@GSEA.ES[i], gene.set@GSEA.ES[i]), c(-0.2*max(c(y1, y2)), 0), lwd = 1, lty = 1, col = 1)
  }
  leg.txt <- c("Null Density", "Null Density Avg", "Observed Density", "Observed NES")
  c.vec <- c("#FFAA00", colors()[555], colors()[29], 1)
  lty.vec <- rep(1, 4)
  lwd.vec <- rep(2, 4)
  legend(x=2.0, y=1.2*y.plot.range[2], bty="n", bg = "white", legend=leg.txt, 
         lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 0.9)
  if(! is.null(pdf)) 
    dev.off()
}

plotSig <- function(gene.set, pdf=NULL) {
  stopifnot( is( gene.set, "SeqGeneSet" ) )
  stopifnot( length(gene.set@GSEA.FDR) > 0 )
  if(! is.null(pdf)) 
    pdf(file = pdf, height = 6, width = 6)
  plot(gene.set@GSEA.ES, gene.set@GSEA.FDR, ylim = c(0, 1), xlim = c(0, 1.5*max(gene.set@GSEA.ES)), 
       col = 1, bg = 1, type="p", pch = 22, cex = 0.75, xlab = "NES", 
       main = "p-values vs. NES", ylab ="p-val/q-val")
  points(gene.set@GSEA.ES, gene.set@GSEA.pval, type = "p", col = 2, bg = 2, pch = 22, cex = 0.75)
  points(gene.set@GSEA.ES, gene.set@GSEA.FWER, type = "p", col = colors()[577], bg = colors()[577], pch = 22, cex = 0.75)
  leg.txt <- c("p-value", "FWER", "FDR")
  cols <- c(2, colors()[577], 1)
  pchs <- c(22, 22, 22)
  legend(x=2, y=0.5, bty="n", bg = "white", legend=leg.txt, pch = pchs, col = cols, pt.bg = cols, cex = 0.9)
  #lines(range(gene.set@GSEA.ES), c(pval.threshold, pval.threshold), lwd = 1, lty = 2, col = 2)
  #lines(range(gene.set@GSEA.ES), c(fwer.threshold, fwer.threshold), lwd = 1, lty = 2, col = colors()[577])
  #lines(range(gene.set@GSEA.ES), c(fdr.threshold, fdr.threshold), lwd = 1, lty = 2, col = 1)
  if(! is.null(pdf)) 
    dev.off()
}

plotSigGeneSet <- function(gene.set, i, gene.score, pdf=NULL) {
  stopifnot( is( gene.set, "SeqGeneSet" ) )
  stopifnot( length(i) == 1 )
  gs <- gene.set[i, drop=FALSE]
  ngene <- length(gene.score)
  gene.score.indx <- order(gene.score, decreasing=TRUE)
  gene.score.sorted <- sort(gene.score, decreasing=TRUE) 
  RES <- gs@GSEA.score.cumsum  # running enrichment score
  min.RES <- min(RES)
  max.RES <- max(RES)
  if (max.RES < 0.3) max.RES <- 0.3
  if (min.RES > -0.3) min.RES <- -0.3
  delta <- (max.RES - min.RES)*0.50 
  min.plot <- min.RES - 2*delta
  max.plot <- max.RES
  
  max.genescore <- max(gene.score)
  min.genescore <- 0
  gene.score.line <- (gene.score.sorted - min.genescore)/(max.genescore - min.genescore) *1.25*delta + min.plot
  zero.genescore.line <- (- min.genescore/(max.genescore - min.genescore)) *1.25*delta + min.plot
  
  if(! is.null(pdf))
    pdf(file = pdf, height = 6, width = 10)
  def.par <- par(no.readonly = TRUE)
  nf <- layout(matrix(c(1,2), 1, 2, byrow=T), 1, c(1, 1), TRUE)
  
  # (1) Running enrichment plot
  main.string <- paste("Gene Set No.", i, ":", gs@GSNames)
  sub.string <- paste("Number of genes: ", ngene, " in list, ", gs@GSSize, " in gene set", sep = "", collapse="")
  plot(1:ngene, RES, main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", 
       xlim=c(1, ngene), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = 2)
  lines(c(1, ngene), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
  lines(c(gs@GSEA.ES.pos, gs@GSEA.ES.pos), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 2) # max enrichment vertical line
  temp <- sapply(match(gene.set[i], gene.score.indx), function(x) {
    lines(c(x, x), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
  })
  for (k in seq(1, ngene, ceiling(ngene/500))) {
    lines(c(k, k), c(zero.genescore.line, gene.score.line[k]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
  }
  lines(1:ngene, gene.score.line, type = "l", lwd = 1, cex = 1, col = 1)
  lines(c(1, ngene), c(zero.genescore.line, zero.genescore.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
  leg.txt <- paste("Peak at ", gs@GSEA.ES.pos, "\nnormalized ES ", signif(gs@GSEA.ES, digits = 3), sep="", collapse="")
  text(x=gs@GSEA.ES.pos, y= -0.1*delta, adj = c(0, 1), labels=leg.txt, cex = 1.0)
  
  # (2) p-val histogram
  adjust.param = 0.5
  sub.string <- paste("normES =", signif(gs@GSEA.ES, digits = 3), 
                      "p-val=", signif(gs@GSEA.pval, digits = 3),
                      "FDR=", signif(gs@GSEA.FDR, digits = 3), 
                      "FWER=", signif(gs@GSEA.FWER, digits = 3))
  temp <- density(gs@GSEA.ES.perm, adjust=adjust.param)
  x.plot.range <- range(temp$x)
  y.plot.range <- c(-0.125*max(temp$y), 1.5*max(temp$y))
  plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, ylim = y.plot.range, 
       lwd = 2, col = 2, main = "Gene Set Null Distribution", xlab = "ES", ylab="P(ES)")
  x.loc <- which.min(abs(temp$x - gs@GSEA.ES))
  lines(c(gs@GSEA.ES, gs@GSEA.ES), c(0, temp$y[x.loc]), lwd = 2, lty = 1, cex = 1, col = 1)
  lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)
  leg.txt <- c("Null Density", "Observed NES")
  legend(x=(0.5*x.plot.range[1]+0.4*x.plot.range[2]), y=y.plot.range[2], bty="n", bg = "white", legend=leg.txt, 
         lty = 1, lwd = 2, col = c(2,1), cex = 1.0)
  
  if(! is.null(pdf))
    dev.off()
  par(def.par)
}


writeSigGeneSet <- function(gene.set, i, gene.score, file="") {
  stopifnot( is( gene.set, "SeqGeneSet" ) )
  stopifnot( length(i) == 1 )
  gs <- gene.set[i, drop=FALSE]
  ngene <- length(gene.score)
  gene.score.indx <- order(gene.score, decreasing=TRUE)
  gene.score.sorted <- sort(gene.score, decreasing=TRUE) 
  
  # a table with genesetName, genesetSize, genesetDesc, NES, pos, p-value, FDR, FWER
  gs.res <- data.frame(rbind(genesetName=paste(gs@name, gs@GSNames, sep=":"), 
                       genesetSize=gs@GSSize, 
                       genesetDesc=gs@GSDescs, 
                       NES = gs@GSEA.ES,
                       Pos = gs@GSEA.ES.pos,
                       pvalue = gs@GSEA.pval,
                       FDR = gs@GSEA.FDR,
                       FWER = gs@GSEA.FWER ) )
  write.table(paste("GSEA result for gene set No. ", i, ":", sep=""), file=file, col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(gs.res, file=file, append=TRUE, quote=FALSE, sep="\t", col.names=FALSE)
  # leading gene set with scores
  leadingset <- match(gs@GS[[1]], gene.score.indx[1:gs@GSEA.ES.pos])
  leadingset <- leadingset[!is.na(leadingset)]
  write.table("", file=file, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table("Leading set: ", file=file, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(sort(gene.score.sorted[leadingset], decreasing=TRUE), file=file, quote=FALSE, sep="\t", col.names=FALSE, append=TRUE)
  # gene scores for this gene set
  write.table("", file=file, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table("Whole gene set: ", file=file, append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE)
  write.table(sort(gene.score[gs@GS[[1]]], decreasing=TRUE), file=file, quote=FALSE, sep="\t", col.names=FALSE, append=TRUE)
}

writeScores <- function(DEscore, DSscore, geneScore=NULL, geneScoreAttr=NULL, file="") {
  if(! is.null(DSscore))
    stopifnot(length(DEscore) == length(DSscore))
  if(! is.null(geneScore)) 
    stopifnot(length(DEscore) == length(geneScore))
  data <- data.frame(cbind(DEscore, DSscore, geneScore)) 
  if(! is.null(geneScore) && ! is.null(geneScoreAttr)) {
    colnames(data)[ncol(data)] <- paste("geneScore", "(", geneScoreAttr, ")", sep="")
  }
  write.table(data, file = file, quote = FALSE, sep = "\t")
}
