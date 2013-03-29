\name{SeqGSEA-package}
\alias{SeqGSEA-package}
\alias{SeqGSEA}
\docType{package}
\title{
SeqGSEA: a Bioconductor package for gene set enrichment analysis of RNA-Seq data
}
\description{
SeqGSEA is an R package for gene set enrichment analysis of RNA-Seq data with the ability to integrate 
differential expression and differential splice in functional analysis. 
}
\details{
\tabular{ll}{
Package: \tab SeqGSEA\cr
Type: \tab Package\cr
License: \tab GPL (>= 3)\cr
}

A User's Guide is available as well as the usual help page documentation for each of the individual functions.


The most useful functions are listed below: 

* ReadCountSet class 
\itemize{
\item \code{\link{ReadCountSet-class}} 
\item \code{\link{ReadCountSet}}
\item \code{\link{exonID}}
\item \code{\link{geneID}}
\item \code{\link{counts-methods}}
\item \code{\link{label}}
\item \code{\link{subsetByGenes}}
}
 
* SeqGeneSet class 
\itemize{
\item \code{\link{SeqGeneSet-class}}
\item \code{\link{geneSetDescs}}
\item \code{\link{geneSetNames}}
\item \code{\link{geneSetSize}}
\item \code{\link{size}}
}

* Load data 
\itemize{
\item \code{\link{newReadCountSet}}
\item \code{\link{loadExonCountData}}
\item \code{\link{newGeneSets}}
\item \code{\link{loadGenesets}}
}

* DE analysis 
\itemize{
\item \code{\link{getGeneCount}}
\item \code{\link{runDESeq}}
\item \code{\link{DENBStat4GSEA}}
\item \code{\link{DENBStatPermut4GSEA}}
\item \code{\link{DENBTest}}
\item \code{\link{DEpermutePval}}
}

* DS analysis 
\itemize{
\item \code{\link{DSpermute4GSEA}}
\item \code{\link{DSpermutePval}}
\item \code{\link{exonTestability}}
\item \code{\link{geneTestability}}
\item \code{\link{estiExonNBstat}}
\item \code{\link{estiGeneNBstat}}
}

* GSEA main 
\itemize{
\item \code{\link{GSEnrichAnalyze}}
\item \code{\link{calES}}
\item \code{\link{calES.perm}}
\item \code{\link{genePermuteScore}}
\item \code{\link{geneScore}}
\item \code{\link{rankCombine}}
\item \code{\link{normES}}
\item \code{\link{normFactor}}
\item \code{\link{scoreNormalization}}
\item \code{\link{signifES}}
}

* Result tables
\itemize{
\item \code{\link{GSEAresultTable}}
\item \code{\link{DSresultExonTable}}
\item \code{\link{DSresultGeneTable}}
\item \code{\link{topDEGenes}}
\item \code{\link{topDSExons}}
\item \code{\link{topDSGenes}}
\item \code{\link{topGeneSets}}
}

* Result displays
\itemize{
\item \code{\link{plotES}}
\item \code{\link{plotGeneScore}}
\item \code{\link{plotSig}}
\item \code{\link{plotSigGeneSet}}
\item \code{\link{writeSigGeneSet}}
}

* Miscellaneous
\itemize{
\item \code{\link{genpermuteMat}}
\item \code{\link{convertEnsembl2Symbol}}
\item \code{\link{convertSymbol2Ensembl}}
}
}

\author{
Xi Wang and Murray J. Cairns

Maintainer: Xi Wang <xi.wang@newcastle.edu.au>
}

\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, in press
}
