\name{plotSigGeneSet}
\alias{plotSigGeneSet}
\title{
Plot gene set details 
}
\description{
This function is to generate a two-panel plot showing detailed information of the gene set specified. 
One panel is showing the running enrichment scores and the position where the ES appear. 
The other panel shows the significance level of the ES, comparing with permutation ESs. 
}
\usage{
plotSigGeneSet(gene.set, i, gene.score, pdf = NULL)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object after running \code{\link{GSEnrichAnalyze}}.}
  \item{i}{the i-th gene set in the SeqGeneSet object. \code{\link{topGeneSets}} is useful to find the most significantly overrepresented gene set. }
  \item{gene.score}{the gene score vector containing gene scores for each gene. }
  \item{pdf}{whether to save the plot to PDF file; if yes, provide the name of the PDF file. }
}
\details{
See \code{\link{writeSigGeneSet}}, which writes the detailed gene set information to a file or to the screen. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}},
\code{\link{topGeneSets}},
\code{\link{plotSig}},
\code{\link{plotES}},
\code{\link{writeSigGeneSet}}
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
plotSigGeneSet(GS_example, 9, gene.score) # 9th gene set is the most significant one.
}
