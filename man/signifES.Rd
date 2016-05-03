\name{signifES}
\alias{signifES}
\title{
Calculate significance of ESs
}
\description{
The is an internal function to calculate significance of ESs of each gene set. For advanced users only.
}
\usage{
signifES(gene.set)
}
\arguments{
  \item{gene.set}{a GeneSet object after running \code{\link{normES}}.} 
}
\value{
A SeqGeneSet object with gene set enrichment significance metrics calculated. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{GSEnrichAnalyze}},
\code{\link{normES}}
}
