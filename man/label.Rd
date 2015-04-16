\name{label}
\alias{label}
\title{
Get the labels of samples in a ReadCountSet object
}
\description{
Get the labels of samples in a ReadCountSet object
}
\usage{
label(RCS)
}
\arguments{
  \item{RCS}{ a ReadCountSet object }
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{newReadCountSet}}
}
\examples{
data(RCS_example, package="SeqGSEA")
label(RCS_example)
}
