# copyright: Xi Wang (xi.wang@newcastle.edu.au)
## DESeq pipeline w/ permutation for RNASeq_GSEA
#require(DESeq)
#require(locfit)
runDESeq <- function(geneCounts, label) {
  DEG <- newCountDataSet( geneCounts, label )
  DEG <- estimateSizeFactors(DEG)
  DEG <- estimateDispersions(DEG, method="per-condition", sharingMode="gene-est-only")
  DEG
}

# modified DESeq function to report variances for each group
DENBTest <- function (cds) 
{
  stopifnot(is(cds, "CountDataSet"))
  if (cds@multivariateConditions) 
    stop("For CountDataSets with multivariate conditions, only the GLM-based test can be used.")
  if (all(is.na(dispTable(cds)))) 
    stop("Call 'estimateDispersions' first.")
  
  colA <- conditions(cds) == levels(conditions(cds))[2]
  colB <- conditions(cds) == levels(conditions(cds))[1]
  bmv <- getBaseMeansAndVariances(counts(cds)[, colA | colB], 
                                  sizeFactors(cds)[colA | colB])
  rawScvA <- fData(cds)[, paste("disp", dispTable(cds)[ levels(conditions(cds))[2] ], 
                                sep = "_")]
  rawScvB <- fData(cds)[, paste("disp", dispTable(cds)[ levels(conditions(cds))[1] ], 
                                sep = "_")]
  pval <- nbinomTestForMatrices(counts(cds)[, colA], 
                                counts(cds)[, colB], 
                                sizeFactors(cds)[colA], 
                                sizeFactors(cds)[colB], 
                                rawScvA, rawScvB)
#   bmvA <- getBaseMeansAndVariances(counts(cds)[, colA], 
#                                    sizeFactors(cds)[colA])
#   bmvB <- getBaseMeansAndVariances(counts(cds)[, colB], 
#                                    sizeFactors(cds)[colB])

  sizeFactorsA <- sizeFactors(cds)[colA]
  sizeFactorsB <- sizeFactors(cds)[colB]
  dispsA <- pmax(fData(cds)[, paste("disp", dispTable(cds)[levels(conditions(cds))[2]], sep = "_")], 1e-8)
  dispsB <- pmax(fData(cds)[, paste("disp", dispTable(cds)[levels(conditions(cds))[1]], sep = "_")], 1e-8)
  musA <- rowMeans(t(t(counts(cds)[, colA])/sizeFactorsA))
  musB <- rowMeans(t(t(counts(cds)[, colB])/sizeFactorsB))
  VarsA <- musA * sum(1/sizeFactorsA) / sum(colA)^2 + dispsA * musA^2 / sum(colA)
  VarsB <- musB * sum(1/sizeFactorsB) / sum(colB)^2 + dispsB * musB^2 / sum(colB)
  
  data.frame(id = rownames(counts(cds)), baseMean = bmv$baseMean, 
             baseMeanA = musA, VarA = VarsA,
             baseMeanB = musB, VarB = VarsB,
             NBstat = (musA - musB) ^ 2 / (VarsA + VarsB), 
             foldChange = musB/musA, log2FoldChange = log2(musB/musA), 
             pval = pval, padj = p.adjust(pval, method = "BH"), 
             stringsAsFactors = FALSE)
}


DEpermutePval <- function(DEGres, permuteNBstat) {
  times <- ncol(permuteNBstat)
  permutePval <- rowSums(DEGres$NBstat <= permuteNBstat) / times
  data.frame(cbind(DEGres, perm.pval = permutePval, perm.padj = p.adjust(permutePval)), 
             stringsAsFactors = FALSE)  
}

## Added on Sep 12: new version for var 
DENBStat4GSEA <- function(cds) {
  stopifnot(is(cds, "CountDataSet"))
  stopifnot(length(levels(conditions(cds))) == 2)
  colA <- conditions(cds) == levels(conditions(cds))[2]
  colB <- conditions(cds) == levels(conditions(cds))[1]
  countsA <- counts(cds)[, colA]
  countsB <- counts(cds)[, colB]
  sizeFactorsA <- sizeFactors(cds)[colA]
  sizeFactorsB <- sizeFactors(cds)[colB]
  dispsA <- pmax(fData(cds)[, paste("disp", dispTable(cds)[levels(conditions(cds))[2]], sep = "_")], 1e-8)
  dispsB <- pmax(fData(cds)[, paste("disp", dispTable(cds)[levels(conditions(cds))[1]], sep = "_")], 1e-8)
  musA <- rowMeans(t(t(countsA)/sizeFactorsA))
  musB <- rowMeans(t(t(countsB)/sizeFactorsB))
  VarsA <- musA * sum(1/sizeFactorsA) / sum(colA)^2 + dispsA * musA^2 / sum(colA)
  VarsB <- musB * sum(1/sizeFactorsB) / sum(colB)^2 + dispsB * musB^2 / sum(colB)
  
  data.frame(id = rownames(counts(cds)), 
             baseMeanA = musA, VarA = VarsA,
             baseMeanB = musB, VarB = VarsB,
             NBstat = (musA - musB) ^ 2 / (VarsA + VarsB), 
             stringsAsFactors = FALSE)
}

## Added on Sep 12: new version for var
## permutate for GSEA
DENBStatPermut4GSEA <- function(DEG, permuteMat) {
  stopifnot( is( DEG, "CountDataSet" ) )
  times <- ncol(permuteMat)
  n_gene <- nrow(counts(DEG))
  #permuteNBstatGene <- matrix(NA_real_, n_gene, times)
  #for(i in 1:times) {
  foreach(i = 1:times, .combine='cbind', .packages = c("DESeq", "SeqGSEA"))  %dopar% {
    conditions(DEG) <- as.factor(permuteMat[,i])
    DEG <- estimateDispersions(DEG, method="per-condition", sharingMode="gene-est-only")
    DEGresPerm <- DENBStat4GSEA( DEG ) 
    #permuteNBstatGene[,i] <- DEGresPerm$NBstat
    DEGresPerm$NBstat
  }
  #permuteNBstatGene
}

topDEGenes <- function(DEGres, n = 20, sortBy = c("padj", "pval", "perm.pval", "perm.padj", "NBstat", "foldChange")) {
## DEGres is an output from DENBtest and DEpermutePval sequentially. 
  sortBy <- match.arg(sortBy, c("padj", "pval", "perm.pval", "perm.padj", "NBstat", "foldChange"))
  stopifnot(sortBy %in% colnames(DEGres))
  switch(sortBy, padj = {
    res <- DEGres[order(DEGres$padj)[1:n],] 
  }, pval = {
    res <- DEGres[order(DEGres$pval)[1:n],] 
  },  perm.padj = {
    res <- DEGres[order(DEGres$perm.padj)[1:n],] 
  }, perm.pval = {
    res <- DEGres[order(DEGres$perm.pval)[1:n],] 
  }, NBstat = {
    res <- DEGres[order(DEGres$NBstat, decreasing = TRUE)[1:n],]
  }, foldChange = {
    res <- DEGres[order(DEGres$foldChange, decreasing = TRUE)[1:n],]
  })
  res 
}

