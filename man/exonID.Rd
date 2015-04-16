\name{exonID}
\alias{exonID}
\alias{exonID<-}
\title{
Accessor to the exonID slot of ReadCountSet objects
}
\description{
Accessor to the exonID slot of ReadCountSet objects
}
\usage{
exonID(RCS)
exonID(RCS) <- value
}
\arguments{
  \item{RCS}{a ReadCountSet object}
  \item{value}{a vector of exon IDs}
}
\value{
A character vector of exon IDs; or a ReadCountSet object.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{newReadCountSet}},
\code{\link{geneID}}
}
\examples{
data(RCS_example, package="SeqGSEA")
exonID(RCS_example)
}
