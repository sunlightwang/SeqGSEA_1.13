\name{rankCombine}
\alias{rankCombine}
\title{
Integration of differential expression and differential splice scores with a rank-based strategy
}
\description{
Integration of differential expression and differential splice scores with a rank-based strategy, 
which simultaneously integrates observed scores and permutation scores using the same ranks. 
}
\usage{
rankCombine(DEscore, DSscore, DEscoreMat, DSscoreMat, DEweight = 0.5)
}
\arguments{
  \item{DEscore}{differential expression scores, normalized.}
  \item{DSscore}{differential splice scores, normalized.}
  \item{DEscoreMat}{differential expression scores in permuted data sets, normalized.}
  \item{DSscoreMat}{differential splice scores in permuted data sets, normalized.}
  \item{DEweight}{any number between 0 and 1 (included), the weight of differential expression scores (so the weight for differential splice is (1-DEweight)).}
}
\details{
This integration method is also known as integration with global ranks. See Wang and Cairns (2013) for details. 
}
\value{
A list with two elements \code{geneScore} and \code{genePermuteScore}.
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, 14(Suppl 5):S16.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{geneScore}}, 
\code{\link{genePermuteScore}}
}
\examples{
data(DEscore, package="SeqGSEA")
data(DSscore, package="SeqGSEA")
data(DEscore.perm, package="SeqGSEA")
data(DSscore.perm, package="SeqGSEA")
combine <- rankCombine(DEscore, DSscore, DEscore.perm, DSscore.perm, DEweight=0.3) 
gene.score <- combine$geneScore
gene.score.perm <- combine$genePermuteScore
}
