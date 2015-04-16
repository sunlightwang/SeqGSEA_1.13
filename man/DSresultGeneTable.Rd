\name{DSresultGeneTable}
\alias{DSresultGeneTable}
\title{
Form a table for DS analysis results at the gene level
}
\description{
Form a table for differential splicing analysis results at the gene level.
}
\usage{
DSresultGeneTable(RCS)
}
\arguments{
  \item{RCS}{A ReadCountSet object with \code{\link{DSpermutePval}} done. }
}
\value{
A data frame containing each gene's NB-statistics, p-values and adjusted p-values for differential splicing analysis. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{DSresultExonTable}},
\code{\link{DSpermutePval}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
RCS_example <- exonTestability(RCS_example)
RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
RCS_example <- DSpermutePval(RCS_example, permuteMat)
head(DSresultGeneTable(RCS_example))
}
