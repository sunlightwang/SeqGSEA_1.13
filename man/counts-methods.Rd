\name{counts-methods}
\docType{methods}
\alias{counts-methods}
\alias{counts,ReadCountSet-method}
\alias{counts<-,ReadCountSet,matrix-method}
\alias{counts}
\title{ Accessors for the 'counts' slot of a ReadCountSet object. }
\description{
Accessors for the 'counts' slot of a ReadCountSet object.
}
\usage{
\S4method{counts}{ReadCountSet}(object)
\S4method{counts}{ReadCountSet,matrix}(object) <- value
}
\arguments{
  \item{object}{a ReadCountSet object}
  \item{value}{a matrix of read counts}
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\examples{
data(RCS_example, package="SeqGSEA")
readCounts <- counts(RCS_example)
head(readCounts)
}
