\name{geneScore}
\alias{geneScore}
\title{
Calculate gene scores by integrating DE and DS scores 
}
\description{
Calculate gene scores by integrating DE and DS scores 
}
\usage{
geneScore(DEscore, DSscore, method = c("linear", "quadratic", "rank"), DEweight = 0.5)
}
\arguments{
  \item{DEscore}{normalized DE scores. }
  \item{DSscore}{normalized DS scores. }
  \item{method}{one of the integration methods: linear, quadratic, or rank; default: linear.}
  \item{DEweight}{any number between 0 and 1 (included), the weight of differential expression scores (the weight for differential splice is (1-DEweight)).}
}
\details{
The integration methods are detailed in Wang and Cairns (2013). Here the rank method refers only to the 
method using data-set-specific ranks. 
}
\value{
A vector of gene scores. 
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, in press
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{genePermuteScore}}
}
\examples{
data(DEscore, package="SeqGSEA")
data(DSscore, package="SeqGSEA")
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
}
