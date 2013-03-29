\name{convertEnsembl2Symbol}
\alias{convertEnsembl2Symbol}
\title{
Convert ensembl gene IDs to gene symbols
}
\description{
Convert ensembl gene IDs to gene symbols
}
\usage{
convertEnsembl2Symbol(ensembl.genes)
}
\arguments{
  \item{ensembl.genes}{ensembl gene ID(s).}
}
\value{
A 2-column matrix showing the correspondence of ensembl gene IDs and gene symbols. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{convertSymbol2Ensembl}}
}
\examples{
convertEnsembl2Symbol("ENSG00000162946") #DISC1
}
