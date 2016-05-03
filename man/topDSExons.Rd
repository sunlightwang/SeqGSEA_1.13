\name{topDSExons}
\alias{topDSExons}
\title{
Extract top differentially spliced exons 
}
\description{
This function is to extract top n differentially spliced exons, ranked by p-values or 
NB-stats. 
}
\usage{
topDSExons(RCS, n = 20, sortBy = c("pvalue", "NBstat"))
}
\arguments{
  \item{RCS}{a ReadCountSet object after running \code{\link{DSpermutePval}}.}
  \item{n}{the number of top genes.}
  \item{sortBy}{indicating whether p-value or NBstat to be used for ranking genes.}
}
\value{
A table for top n exons. Columns include: geneID, exonID, testable, NBstat, pvalue, padjust, and meanCounts.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
 \code{\link{topDSGenes}},
 \code{\link{DSpermutePval}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
RCS_example <- exonTestability(RCS_example)
RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
RCS_example <- DSpermutePval(RCS_example, permuteMat)
topDSExons(RCS_example, 10, "NB")
}
