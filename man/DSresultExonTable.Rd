\name{DSresultExonTable}
\alias{DSresultExonTable}
\title{
Form a table for DS analysis results at the Exon level
}
\description{
Form a table for differential splicing analysis results at the Exon level.
}
\usage{
DSresultExonTable(RCS)
}
\arguments{
  \item{RCS}{A ReadCountSet object with \code{\link{DSpermutePval}} done. }
}
\details{
A data frame containing each exon's NB-statistics, p-values and adjusted p-values for differential splicing analysis. 
}
\value{
A matrix containing exon DS analysis results, including testability, NBstats, p-values and adjusted p-values. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{DSresultGeneTable}},
\code{\link{DSpermutePval}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
RCS_example <- exonTestability(RCS_example)
RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
RCS_example <- DSpermutePval(RCS_example, permuteMat)
head(DSresultExonTable(RCS_example))
}
