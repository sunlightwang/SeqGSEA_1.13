\name{estiGeneNBstat}
\alias{estiGeneNBstat}
\title{
Calculate NB-statistics quantifying differential splicing for each gene
}
\description{
Calculate NB-statistics quantifying differential splicing for each gene between two groups of samples compared. The results will be used for GSEA run (as DS-scores).
}
\usage{
estiGeneNBstat(RCS)
}
\arguments{
  \item{RCS}{a ReadCountSet object after running \code{estiExonNBstat}.}
}
\note{
Please run \code{\link{estiExonNBstat}} before run this function. 
}
\value{
A ReadCountSet object with slot \code{featureData_gene} updated. 
}
\references{
Weichen Wang, Zhiyi Qin, Zhixing Feng, Xi Wang and Xuegong Zhang (2013). Identifying differentially spliced genes from two groups of RNA-seq samples. Gene, 518(1):164-170.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{estiExonNBstat}}
}
\examples{
data(RCS_example, package="SeqGSEA")
RCS_example <- exonTestability(RCS_example, cutoff=5)
RCS_example <- estiExonNBstat(RCS_example)
RCS_example <- estiGeneNBstat(RCS_example)
head(RCS_example@featureData_gene)
}
