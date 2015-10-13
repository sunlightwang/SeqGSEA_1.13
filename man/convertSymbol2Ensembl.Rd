\name{convertSymbol2Ensembl}
\alias{convertSymbol2Ensembl}
\title{
Convert gene symbols to ensembl gene IDs
}
\description{
Convert gene symbols to ensembl gene IDs
}
\usage{
convertSymbol2Ensembl(symbols)
}
\arguments{
  \item{symbols}{gene symbol(s). }
}
\value{
A 2-column matrix showing the correspondence of gene symbols and ensembl gene IDs.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{convertEnsembl2Symbol}}
}
\examples{
convertSymbol2Ensembl("DISC1") #ENSG00000162946
}
