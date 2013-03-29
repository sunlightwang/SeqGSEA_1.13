\name{GSEnrichAnalyze}
\alias{GSEnrichAnalyze}
\title{
Main function of gene set enrichment analysis 
}
\description{
The main function of gene set enrichment analysis 
}
\usage{
GSEnrichAnalyze(gene.set, gene.score, gene.score.perm, weighted.type = 1)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object.}
  \item{gene.score}{a vector of integrated gene scores in the same order as genes listed in the \code{geneList} slot of \code{gene.set}. }
  \item{gene.score.perm}{a matrix of integrated gene scores on permutation data sets; row: genes; col: permutation. }
  \item{weighted.type}{weight type for gene scores; default: 1. }
}
\value{
A SeqGeneSet object with many slots updated, such as \code{GSEA.ES} and \code{GSEA.pval}.
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, in press
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{normES}},
\code{\link{signifES}}
}
\examples{
data(DEscore, package="SeqGSEA")
data(DSscore, package="SeqGSEA")
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
data(DEscore.perm, package="SeqGSEA")
data(DSscore.perm, package="SeqGSEA")
gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, method="linear",  DEweight=0.3)
data(GS_example, package="SeqGSEA") 
GS_example <- GSEnrichAnalyze(GS_example, gene.score, gene.score.perm)
topGeneSets(GS_example, 5)
}
