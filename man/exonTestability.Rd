\name{exonTestability}
\alias{exonTestability}
\title{
Check exon testability 
}
\description{
Check exon testability, filtering out exons with very few (default: 5) read counts 
}
\usage{
exonTestability(RCS, cutoff = 5)
}
\arguments{
  \item{RCS}{a ReadCountSet object.}
  \item{cutoff}{exons with read counts less than this cutoff are to be marked as untestable. }
}
\value{
a ReadCountSet object with slot \code{fData} updated. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{geneTestability}}
}
\examples{
data(RCS_example, package="SeqGSEA")
RCS_example <- exonTestability(RCS_example, cutoff=5)
head(fData(RCS_example))
}
