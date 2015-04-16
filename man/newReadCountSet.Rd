\name{newReadCountSet}
\alias{newReadCountSet}
\title{
Generate a new ReadCountSet object
}
\description{
This is a internal function to generate a new ReadCountSet object. 
}
\usage{
newReadCountSet(readCounts, exonIDs, geneIDs)
}
\arguments{
  \item{readCounts}{a data frame, read counts for each exon of each samples. Must have colnames, which indicate the label of samples.}
  \item{exonIDs}{a character vector indicating exon IDs.}
  \item{geneIDs}{a character vector indicating gene IDs.}
}
\value{
A object of the ReadCountSet class.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{loadExonCountData}},
\code{\link{ReadCountSet-class}}
}
\examples{
rcounts <- cbind(t(sapply(1:10, function(x) {rnbinom(5, size=10, prob=runif(1))} ) ) , 
                 t(sapply(1:10, function(x) {rnbinom(5, size=10, prob=runif(1))} ) ) )
colnames(rcounts) <- c(paste("S", 1:5, sep=""), paste("C", 1:5, sep="")) 
geneIDs <- c(rep("G1", 4), rep("G2", 6))
exonIDs <- c(paste("E", 1:4, sep=""), paste("E", 1:6, sep=""))
## 
RCS <- newReadCountSet(rcounts, exonIDs, geneIDs)
RCS 
## End
}
