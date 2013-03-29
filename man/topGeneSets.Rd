\name{topGeneSets}
\alias{topGeneSets}
\title{
Extract top significant gene sets
}
\description{
This function is to extract n top significant gene sets overrepresented in the samples 
studied, ranked by FDR, p-values, or FWER. 
}
\usage{
topGeneSets(gene.set, n = 20, sortBy = c("FDR", "pvalue", "FWER"), GSDesc = FALSE)
}
\arguments{
  \item{gene.set}{an object of class SeqGeneSet after GSEA runs.}
  \item{n}{the number of top gene sets.}
  \item{sortBy}{indicating which method to rank gene sets.}
  \item{GSDesc}{logical indicating whether or not to output gene set descriptions.}
}
\value{
A data frame for top n gene sets detected with respect to the ranking method specified. Information includes: GSName, GSSize, ES, ES.pos, pval, FDR, and FWER.

}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
 \code{\link{GSEnrichAnalyze}},
 \code{\link{GSEAresultTable}}
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
}
