\name{geneSetDescs}
\alias{geneSetDescs}
\title{
Get the descriptions of gene sets in a SeqGeneSet object
}
\description{
Get the descriptions of gene sets in a SeqGeneSet object
}
\usage{
geneSetDescs(GS)
}
\arguments{
  \item{GS}{a SeqGeneSet object.}
}
\details{
Gene sets with size less than \code{GSSizeMin} or more than \code{GSSizeMax} are not included.
}
\value{
A vector of descriptions of each gene set in the SeqGeneSet object.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{geneSetNames}}, 
\code{\link{geneSetSize}}, 
\code{\link{SeqGeneSet-class}}, 
\code{\link{loadGenesets}}
}
\examples{
data(GS_example, package="SeqGSEA")
geneSetDescs(GS_example)
}
