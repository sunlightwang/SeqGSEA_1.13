\name{geneSetSize}
\alias{geneSetSize}
\title{
Get the numbers of genes in each gene set in a SeqGeneSet object
}
\description{
Get the numbers of genes in each gene set in a SeqGeneSet object
}
\usage{
geneSetSize(GS)
}
\arguments{
  \item{GS}{a SeqGeneSet object.}
}
\details{
Gene sets with size less than \code{GSSizeMin} or more than \code{GSSizeMax} are not included.
}
\value{
A vector of integers indicating the number of genes in each gene set in this SeqGeneSet object. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{geneSetNames}},
\code{\link{geneSetDescs}},
\code{\link{SeqGeneSet-class}},
\code{\link{loadGenesets}}
}
\examples{
data(GS_example, package="SeqGSEA")
geneSetSize(GS_example)
}
