\name{estiExonNBstat}
\alias{estiExonNBstat}
\title{
Calculate NB-statistics quantifying differential splicing for individual exons
}
\description{
Calculate NB-statistics quantifying differential splicing for individual exons between two groups of samples compared.
}
\usage{
estiExonNBstat(RCS)
}
\arguments{
  \item{RCS}{a ReadCountSet object after running \code{exonTestability}.}
}
\note{
Please run \code{\link{exonTestability}} before you run this function. 
}
\value{
A ReadCountSet object with the slot \code{featureData} updated. 
}
\references{
Weichen Wang, Zhiyi Qin, Zhixing Feng, Xi Wang and Xuegong Zhang (2013). Identifying differentially spliced genes from two groups of RNA-seq samples. Gene, 518(1):164-170.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{exonTestability}},
\code{\link{estiGeneNBstat}}
}
\examples{
data(RCS_example, package="SeqGSEA")
RCS_example <- exonTestability(RCS_example, cutoff=5)
RCS_example <- estiExonNBstat(RCS_example)
head(fData(RCS_example))
}
