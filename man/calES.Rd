\name{calES}
\alias{calES}
\title{
Calculate running enrichment scores of gene sets 
}
\description{
This is an internal function to calculate running enrichment scores of each gene set in the SeqGeneSet object specified
}
\usage{
calES(gene.set, gene.score, weighted.type = 1)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object.}
  \item{gene.score}{a vector of gene scores corresponding to the \code{geneList} slot of \code{gene.set}.}
  \item{weighted.type}{gene score weight type.}
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}}, 
\code{\link{calES.perm}},
}
\examples{
data(DEscore, package="SeqGSEA")
data(DSscore, package="SeqGSEA")
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
data(GS_example, package="SeqGSEA")
rES <- calES(GS_example, gene.score)
rES[1,]
}
