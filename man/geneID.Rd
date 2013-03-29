\name{geneID}
\alias{geneID}
\alias{geneID<-}
\title{
Accessor to the geneID slot of ReadCountSet objects
}
\description{
Accessor to the geneID slot of ReadCountSet objects
}
\usage{
geneID(RCS)
geneID(RCS) <- value
}
\arguments{
  \item{RCS}{a ReadCountSet object}
  \item{value}{a vector of gene IDs}
}
\value{
A character vector of gene IDs, which can be duplicated; or a ReadCountSet object.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{newReadCountSet}},
\code{\link{exonID}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneID(RCS_example)
}
