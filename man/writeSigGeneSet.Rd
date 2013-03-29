\name{writeSigGeneSet}
\alias{writeSigGeneSet}
\title{
Write gene set supporting information 
}
\description{
This function is to write the specified gene set (whose index is i) 
with significance information, including p-value and FDR, and gene 
scores for each gene in this set. 
}
\usage{
writeSigGeneSet(gene.set, i, gene.score, file = "")
}
\arguments{
  \item{gene.set}{an object of class SeqGeneSet with \code{\link{GSEnrichAnalyze}} done. }
  \item{i}{the i-th gene set in the SeqGeneSet object. \code{\link{topGeneSets}} is useful to find the most significantly overrepresented gene set. }
  \item{gene.score}{the vector of gene scores for running GSEA.}
  \item{file}{output file name, if not specified print to screen.}
}
\details{
See \code{\link{plotSigGeneSet}}, which shows graphic information of the gene set specified. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}},
\code{\link{topGeneSets}},
\code{\link{plotSigGeneSet}}
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
topGeneSets(GS_example, n=5)
writeSigGeneSet(GS_example, 9, gene.score) # 9th gene set is the most significant one.
}
