\name{GSEAresultTable}
\alias{GSEAresultTable}
\title{
Form a table for GSEA results
}
\description{
Form a table for GSEA results. 
}
\usage{
GSEAresultTable(gene.set, GSDesc = FALSE)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object after running \code{\link{GSEnrichAnalyze}}.}
  \item{GSDesc}{logical indicating whether to output gene set descriptions. default: FALSE}
}
\value{
A data frame containing columns of GSName, GSSize, ES, ES.pos, pval, FDR, and FWER.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}},
\code{\link{topGeneSets}}
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
head(GSEAresultTable(GS_example))
}
