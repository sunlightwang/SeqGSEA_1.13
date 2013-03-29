\name{getGeneCount}
\alias{getGeneCount}
\title{
Calculate read counts of genes from a ReadCountSet object
}
\description{
Calculate read counts of genes from a ReadCountSet object
}
\usage{
getGeneCount(RCS)
}
\arguments{
  \item{RCS}{a ReadCountSet object}
}
\details{
This function can be used to get gene read counts from exon read counts. 
}
\value{
a matrix of gene read counts for each gene (row) and each sample (col).
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{loadExonCountData}},
\code{\link{runDESeq}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneCounts <- getGeneCount(RCS_example)
}
