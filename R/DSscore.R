# copyright: Xi Wang (xi.wang@newcastle.edu.au)
# DS_NB.R: implement Weichen Wang's DS methods 

DSresultExonTable <- function(RCS)
{
  stopifnot( is( RCS, "ReadCountSet" ) )
  result <- data.frame(
    geneID=geneID(RCS),
    exonID=exonID(RCS),
    testable=fData(RCS)$testable,
    NBstat=fData(RCS)$NBstat,
    pvalue=fData(RCS)$pvalue,
    padjust=fData(RCS)$padjust,
    meanCounts=rowMeans(counts(RCS)))
  result
}

DSresultGeneTable <- function(RCS)
{
  stopifnot( is( RCS, "ReadCountSet" ) )
  featureData_gene <- RCS@featureData_gene
  result <- data.frame(
    geneID = rownames(featureData_gene),
    NBstat = featureData_gene$NBstat,
    pvalue = featureData_gene$pval,
    padjust = featureData_gene$padj 
    #meanCounts=rowMeans(counts(RCS))
    )
  result
}

topDSExons <- function(RCS, n = 20, sortBy = c("pvalue", "NBstat")) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  sortBy <- match.arg(sortBy, c("pvalue", "NBstat"))
  res <- DSresultExonTable(RCS)
  switch(sortBy, pvalue = {
    res <- res[order(res$pvalue)[1:n],] 
  }, NBstat = {
    res <- res[order(res$NBstat, decreasing = TRUE)[1:n],]
  })
  res 
}

topDSGenes <- function(RCS, n = 20, sortBy = c("pvalue", "NBstat")) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  sortBy <- match.arg(sortBy, c("pvalue", "NBstat"))
  res <- DSresultGeneTable(RCS)
  switch(sortBy, pvalue = {
    res <- res[order(res$pvalue)[1:n],] 
  }, NBstat = {
    res <- res[order(res$NBstat, decreasing = TRUE)[1:n],]
  })
  res 
}

exonTestability <- function(RCS, cutoff=5) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  fData(RCS)$testable <- rowSums(counts(RCS)) >= cutoff
  RCS
}

geneTestability <- function(RCS) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  if(any(is.na(fData(RCS)$testable))) stop("Please run exonTestability first.")
  tapply( 1:nrow(RCS), geneID(RCS), function(rows)
    any( fData(RCS)$testable[rows] ) )
}

estiExonProbVar <- function(dcounts, testable, label) {
  nExon <- nrow(dcounts)
  if(is.null(nExon)) { nExon <- 1 }
  prob_case <-rep( NA_real_, nExon )
  prob_ctrl <-rep( NA_real_, nExon )
  var_case <-rep( NA_real_, nExon )
  var_ctrl <-rep( NA_real_, nExon )

  if( sum(testable) <= 1) 
    return (data.frame(prob_case=prob_case, prob_ctrl=prob_ctrl, var_case=var_case, var_ctrl=var_ctrl))
  
  label <- as.factor(label)
  classes <- levels(label)
  cases <- (label == classes[2])  
  n_case <- sum(cases)
  ctrls <- (label == classes[1])
  n_ctrl <- sum(ctrls)
  dcnt <- dcounts[testable,] + .5  # add a dummy read to all counts in case of 0's
  dcnt_case <- dcnt[,cases] 
  dcnt_ctrl <- dcnt[,ctrls] 
  m <- colSums(dcnt)
  q_case <- estiPhi(dcnt_case)
  p_case <- q_case$prob                 ###
  v_case <- calVar(q_case, m[cases])    ###
  
  q_ctrl <- estiPhi(dcnt_ctrl)
  p_ctrl <- q_ctrl$prob                 ###
  v_ctrl <- calVar(q_ctrl, m[ctrls])    ###

  prob_case[testable] <- p_case
  var_case[testable] <- v_case
  prob_ctrl[testable] <- p_ctrl
  var_ctrl[testable] <- v_ctrl
  
  data.frame(prob_case=prob_case, prob_ctrl=prob_ctrl, var_case=var_case, var_ctrl=var_ctrl)
}

estiExonNBstat <- function(RCS) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  if(any(is.na(fData(RCS)$testable))) stop("Please run exonTestability first.")
  dcounts <- counts(RCS)
  geneIDs <- geneID(RCS)
  testable <- fData(RCS)$testable
  label <- pData(RCS)$label
  groupStat <- do.call(rbind, tapply(1:nrow(dcounts), geneIDs, 
                                     function(rows) { 
                                       estiExonProbVar(dcounts[rows,], testable[rows], label) } ) )
  #geneUniq <- levels(geneIDs)
  #groupStat <- foreach(i = geneUniq, .combine=rbind) %dopar% {
  #  rows <- which(geneIDs == i)
  #  estiExonProbVar(dcounts[rows,], testable[rows], label) 
  #} 
  NBStat <- (groupStat$prob_case - groupStat$prob_ctrl) * (groupStat$prob_case - groupStat$prob_ctrl) /
    (groupStat$var_case + groupStat$var_ctrl)
  fData(RCS)$prob_case <- groupStat$prob_case
  fData(RCS)$prob_ctrl <- groupStat$prob_ctrl
  fData(RCS)$var_case <- groupStat$var_case
  fData(RCS)$var_ctrl <- groupStat$var_ctrl
  fData(RCS)$NBstat <- NBStat
  RCS
}

estiGeneNBstat <- function(RCS) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  if(all(is.na(fData(RCS)$NBstat))) stop("Please run estiExonNBstat first.")
  geneIDs <- geneID(RCS) 
  n_exon <- length(fData(RCS)$exonIDs)
  geneNBstat <- tapply(1:n_exon, geneIDs, function(rows) {
    NBstats <- fData(RCS)$NBstat[rows]
    testables <- fData(RCS)$testable[rows]
    mean(NBstats[testables]) } )
  RCS@featureData_gene$NBstat <- as.numeric(geneNBstat)
  rownames(RCS@featureData_gene) <- rownames(geneNBstat)
  RCS
}

genpermuteMat <- function(RCS, times = 1000) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  label <- as.numeric(pData(RCS)$label) - 1
  n_sam <- length(label)
  permuteMat <- matrix(0, n_sam, times)
  for(i in 1:times) {
    permuteMat[,i] <- label[sample(n_sam,n_sam)]
  }
  permuteMat # rows for samples, and cols for every permutation 
}
  
DSpermutePval <- function(RCS, permuteMat) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  if(all(is.na(fData(RCS)$NBstat))) stop("Please run estiExonNBstat first.")
  if(all(is.na(RCS@featureData_gene$NBstat))) stop("Please run estiGeneNBstat first.")
  times <- ncol(permuteMat)
  dcounts <- counts(RCS)
  n_exon <- length(fData(RCS)$exonIDs)
  geneIDs <- geneID(RCS) 
  n_gene <- length(unique(geneIDs))
  testables <- fData(RCS)$testable
  permuteNBstatExon <- matrix(NA_real_, n_exon, times)
  permuteNBstatGene <- matrix(NA_real_, n_gene, times)
  for(i in 1:times) {
    groupStat <- do.call(rbind, tapply(1:n_exon, geneIDs, 
                                       function(rows) { 
                                         estiExonProbVar(dcounts[rows,], testables[rows], 
                                                         as.factor(permuteMat[,i])) } ) )
    permuteNBstatExon[,i] <- (groupStat$prob_case - groupStat$prob_ctrl) * (groupStat$prob_case - groupStat$prob_ctrl) /
      (groupStat$var_case + groupStat$var_ctrl)
    permuteNBstatGene[,i] <- tapply(1:n_exon, geneIDs, function(rows) {
      mean(permuteNBstatExon[rows,i][testables[rows]]) } )
  }
  permutePvalExon <- rowSums(fData(RCS)$NBstat <= permuteNBstatExon) / times
  permutePvalGene <- rowSums(RCS@featureData_gene$NBstat <= permuteNBstatGene) / times
  
  RCS@permute_NBstat_exon <- permuteNBstatExon
  RCS@permute_NBstat_gene <- permuteNBstatGene
  fData(RCS)$pvalue <- permutePvalExon
  RCS@featureData_gene$pval <- permutePvalGene
  fData(RCS)$padjust <- p.adjust(permutePvalExon, method="BH")
  RCS@featureData_gene$padj <- p.adjust(permutePvalGene, method="BH")
  RCS
}

DSpermute4GSEA <- function(RCS, permuteMat) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  if(any(is.na(fData(RCS)$testable))) stop("Please run exonTestability first.")
  times <- ncol(permuteMat)
  dcounts <- counts(RCS)
  n_exon <- length(fData(RCS)$exonIDs)
  geneIDs <- geneID(RCS) 
  n_gene <- length(unique(geneIDs))
  testables <- fData(RCS)$testable
  permuteNBstatGene <- matrix(NA_real_, n_gene, times)
  #for(i in 1:times) {
  permuteNBstatGene <- foreach(i = 1:times, .combine='cbind', .packages=c("SeqGSEA"))  %dopar% {
    groupStat <- do.call(rbind, tapply(1:n_exon, geneIDs, 
                                       function(rows) { 
                                         SeqGSEA:::estiExonProbVar(dcounts[rows,], testables[rows], 
                                                         as.factor(permuteMat[,i])) } ) )
    permuteNBstatExon <- (groupStat$prob_case - groupStat$prob_ctrl) * (groupStat$prob_case - groupStat$prob_ctrl) /
      (groupStat$var_case + groupStat$var_ctrl)
    #permuteNBstatGene[,i] <- tapply(1:n_exon, geneIDs, function(rows) {
    #  mean(permuteNBstatExon[rows][testables[rows]]) } )
    tapply(1:n_exon, geneIDs, function(rows) {
      mean(permuteNBstatExon[rows][testables[rows]]) } )
  }
  RCS@permute_NBstat_gene <- permuteNBstatGene
  RCS
}

getGeneCount <- function(RCS) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  do.call( rbind,
           tapply( 1:nrow(RCS), geneID(RCS), function(rows)
             colSums( counts(RCS)[rows,,drop=FALSE] ) ) )
}

## functions related to DS above ##
