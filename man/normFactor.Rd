\name{normFactor}
\alias{normFactor}
\title{
Get normalization factors for normalization DE or DS scores
}
\description{
Get normalization factors from permutation scores for normalization DE or DS scores
}
\usage{
normFactor(permStat)
}
\arguments{
  \item{permStat}{a matrix of NB-statistics from permutation data sets, with row corresponding to genes and columns to permutations. }
}
\value{
A vector of normalization factors, each for one gene. 
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, in press
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{scoreNormalization}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
RCS_example <- exonTestability(RCS_example)
RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
RCS_example <- DSpermute4GSEA(RCS_example, permuteMat)
## (not run)
DSscore.normFac <- normFactor(RCS_example@permute_NBstat_gene)
DSscore <- scoreNormalization(RCS_example@featureData_gene$NBstat, DSscore.normFac)
DSscore.perm <- scoreNormalization(RCS_example@permute_NBstat_gene, DSscore.normFac)
## End (not run)
}
