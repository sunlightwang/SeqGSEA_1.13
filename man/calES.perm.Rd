\name{calES.perm}
\alias{calES.perm}
\title{
Calculate enrichment scores for gene sets in the permutation data sets 
}
\description{
This is an internal function to calculate enrichment scores for gene sets in the permutation data sets.
}
\usage{
calES.perm(gene.set, gene.score.perm, weighted.type = 1)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object.}
  \item{gene.score.perm}{a matrix of gene scores on the permutation data sets.}
  \item{weighted.type}{gene score weight type.}
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}}, 
\code{\link{calES}},
}
\examples{
data(DEscore.perm, package="SeqGSEA")
data(DSscore.perm, package="SeqGSEA")
gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, method="linear",  DEweight=0.3)
data(GS_example, package="SeqGSEA")
ES.perm <- calES.perm(GS_example, gene.score.perm)
ES.perm[1:5,1:5]
}
