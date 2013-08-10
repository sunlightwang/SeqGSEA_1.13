\name{genpermuteMat}
\alias{genpermuteMat}
\title{
Generate permutation matrix 
}
\description{
Generate permutation matrix from ReadCountSet objects or from label vectors. 
}
\usage{
genpermuteMat(obj, times = 1000)
}
\arguments{
  \item{obj}{a ReadCountSet object or a label vector. This function needs the original sample label information to generate permutation matrix. }
  \item{times}{an integer indication the times of permutation. }
}
\value{
A sample label shuffled matrix, rows corresponding to samples and columns for each permutation. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{DSpermute4GSEA}},
\code{\link{DENBStatPermut4GSEA}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
RCS_example <- exonTestability(RCS_example)
RCS_example <- DSpermute4GSEA(RCS_example, permuteMat)
}
