\name{size}
\alias{size}
\title{
Number of gene sets in a SeqGeneSet object
}
\description{
This function to get the number of gene sets in a SeqGeneSet object.
}
\usage{
size(GS)
}
\arguments{
  \item{GS}{an object of class SeqGeneSet.}
}
\details{
Gene sets with size less than \code{GSSizeMin} or more than \code{GSSizeMax} are not included.
}
\value{
The number of gene sets in this SeqGeneSet object. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{SeqGeneSet-class}},
\code{\link{loadGenesets}}
}
\examples{
data(GS_example, package="SeqGSEA")
size(GS_example)
}
