\name{topDSGenes}
\alias{topDSGenes}
\title{
Extract top differentially spliced genes
}
\description{
This function to extract top n differentially spliced genes, ranked 
by p-values or NBstats. 
}
\usage{
topDSGenes(RCS, n = 20, sortBy = c("pvalue", "NBstat"))
}
\arguments{
  \item{RCS}{a ReadCountSet object after running \code{\link{DSpermutePval}}.}
  \item{n}{the number of top genes.}
  \item{sortBy}{indicating whether p-value or NBstat to be used for ranking genes.}
}
\value{
A table for top n genes. Columns include: geneID, NBstat, pvalue, and padjust.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
 \code{\link{topDSExons}},
 \code{\link{DSpermutePval}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
RCS_example <- exonTestability(RCS_example)
RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
RCS_example <- DSpermutePval(RCS_example, permuteMat)
topDSGenes(RCS_example, 10, "NB")
}
