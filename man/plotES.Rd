\name{plotES}
\alias{plotES}
\title{
Plot the distribution of enrichment scores
}
\description{
This function is to plot the distribution of enrichment scores, with comparison with permutation 
enrichment scores. 
}
\usage{
plotES(gene.set, pdf = NULL)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object after running \code{\link{GSEnrichAnalyze}}.}
  \item{pdf}{whether to save the plot to PDF file; if yes, provide the name of the PDF file.}
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}},
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
plotES(GS_example) 
}
