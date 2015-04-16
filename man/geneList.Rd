\name{geneList}
\alias{geneList}
\title{
Get the gene list in a SeqGeneSet object
}
\description{
Get the gene list in a SeqGeneSet object
}
\usage{
geneList(GS) 
}
\arguments{
  \item{GS}{A SeqGeneSet object.}
}
\details{
TBA
}
\value{
A vector of gene IDs. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{loadGenesets}},
\code{\link{SeqGeneSet-class}}
}
\examples{
##
gs <- newGeneSets(GS=list(1:10, 6:15, 11:20),
                  geneList=paste("Gene", 1:22, sep=""), 
                  GSNames=c("gs1","gs2","gs3"), 
                  GSDescs=c("test1","test2","test3"), 
                  name="gs examples")
geneList(gs)
## End 
}
