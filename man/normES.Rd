\name{normES}
\alias{normES}
\title{
Normalize enrichment scores 
}
\description{
This is an internal function to normalize enrichment scores. For advanced users only.
}
\usage{
normES(gene.set)
}
\arguments{
  \item{gene.set}{a SeqGeneSet object after running \code{\link{calES}} and \code{\link{calES.perm}}. }
}
\value{
A SeqGeneSet object with ES scores normalized. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}},
\code{\link{signifES}}
}
