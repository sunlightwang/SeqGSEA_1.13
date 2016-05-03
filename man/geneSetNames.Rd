\name{geneSetNames}
\alias{geneSetNames}
\title{
Get the names of gene set in a SeqGeneSet object
}
\description{
Get the names of gene set in a SeqGeneSet object
}
\usage{
geneSetNames(GS)
}
\arguments{
  \item{GS}{a SeqGeneSet object.}
}
\details{
Gene sets with size less than \code{GSSizeMin} or more than \code{GSSizeMax} are not included. 
}
\value{
A vector of gene set names in this SeqGeneSet object.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{geneSetDescs}},
\code{\link{geneSetSize}},
\code{\link{SeqGeneSet-class}},
\code{\link{loadGenesets}}
}
\examples{
data(GS_example, package="SeqGSEA")
geneSetNames(GS_example)
}
