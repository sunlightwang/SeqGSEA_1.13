\name{geneTestability}
\alias{geneTestability}
\title{
Check gene testability
}
\description{
This function is to determine each gene's testability. A gene is testable if at least one of its exons are testable.
}
\usage{
geneTestability(RCS)
}
\arguments{
  \item{RCS}{a ReadCountSet object after exon testability checked, usually the output of \code{\link{exonTestability}}. }
}
\details{
This result can applied to filter out genes not expressed. 
}
\value{
A logical vector indicating which genes are testable, i.e., having at least one exon testable. 
}
\note{
Please run \code{\link{exonTestability}} before run this function. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{exonTestability}},
\code{\link{subsetByGenes}}
}
\examples{
data(RCS_example, package="SeqGSEA")
RCS_example <- exonTestability(RCS_example, cutoff=5)
geneTestable <- geneTestability(RCS_example)
head(geneTestable)
}
