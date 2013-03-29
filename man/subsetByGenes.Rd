\name{subsetByGenes}
\alias{subsetByGenes}
\title{
Get a new ReadCountSet with specified gene IDs. 
}
\description{
Get a new ReadCountSet with specified gene IDs. 
}
\usage{
subsetByGenes(RCS, genes)
}
\arguments{
  \item{RCS}{a ReadCountSet object.}
  \item{genes}{a list of gene IDS.}
}
\value{
This function returns a new ReadCountSet object, with changes in slots \code{assayData}, \code{featureData}, \code{featureData_gene}, and \code{permute_NBstat_exon} and \code{permute_NBstat_gene} if they have been calculated.  
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{newReadCountSet}},
\code{\link{ReadCountSet}}
}
\examples{
data(RCS_example, package="SeqGSEA")
RCS_example
genes <- c("ENSG00000000938", "ENSG00000000005")
RCS_sub <- subsetByGenes(RCS_example, genes)
RCS_sub
}
