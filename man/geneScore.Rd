\name{geneScore}
\alias{geneScore}
\title{
Calculate gene scores by integrating DE and DS scores 
}
\description{
Calculate gene scores by integrating DE and DS scores 
}
\usage{
geneScore(DEscore, DSscore = NULL, method = c("linear", "quadratic", "rank"), DEweight = 0.5)
}
\arguments{
  \item{DEscore}{normalized DE scores. }
  \item{DSscore}{normalized DS scores. }
  \item{method}{one of the integration methods: linear, quadratic, or rank; default: linear.}
  \item{DEweight}{any number between 0 and 1 (included), the weight of differential expression scores (the weight for differential splice is (1-DEweight)).}
}
\details{
The integration methods including "linear", "quadratic", and "rank" are detailed in Wang and Cairns (2013). 
Here the rank method refers only to the method using data-set-specific ranks. 

For DE-only analysis, just specify DEweight to be 1, and the DSscore value can be NULL. 
}
\value{
A vector of gene scores. 
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, 14(Suppl 5):S16.
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
# linear combination with weight for DE 0.3 
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
# DE only analysis 
gene.score <- geneScore(DEscore, DEweight = 1) 
}
