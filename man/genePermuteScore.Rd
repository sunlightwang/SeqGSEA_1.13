\name{genePermuteScore}
\alias{genePermuteScore}
\title{
Calculate gene scores on permutation data sets
}
\description{
Calculate gene scores on permutation data sets
}
\usage{
genePermuteScore(DEscoreMat, DSscoreMat = NULL, method = c("linear", "quadratic", "rank"), 
                 DEweight = 0.5)
}
\arguments{
  \item{DEscoreMat}{normalized DE scores on permutation data sets.}
  \item{DSscoreMat}{normalized DS scores on permutation data sets.}
  \item{method}{one of the integration methods: linear, quadratic, or rank; default: linear.}
  \item{DEweight}{any number between 0 and 1 (included), the weight of differential expression scores (the weight for differential splice is (1-DEweight)).}
}
\details{
The integration methods including "linear", "quadratic", and "rank" are detailed in Wang and Cairns (2013).
Here the rank method refers only to the method using data-set-specific ranks.

For DE-only analysis, just specify DEweight to be 1, and the DSscoreMat value can be NULL.
}
\value{
A gene score matrix. 
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, 14(Suppl 5):S16.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{geneScore}}
}
\examples{
data(DEscore.perm, package="SeqGSEA")
data(DSscore.perm, package="SeqGSEA")
# linear combination with weight for DE 0.3 
gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, method="linear",  DEweight=0.3)
# DE only analysis 
gene.score.perm <- genePermuteScore(DEscore.perm, DEweight=1)
}
