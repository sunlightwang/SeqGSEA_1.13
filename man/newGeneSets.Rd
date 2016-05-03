\name{newGeneSets}
\alias{newGeneSets}
\title{
Initialize a new SeqGeneSet object
}
\description{
This is an internal function to generate a new SeqGeneSet object.
}
\usage{
newGeneSets(GS, GSNames, GSDescs, geneList, name = NA_character_, 
            sourceFile = NA_character_, GSSizeMin = 5, GSSizeMax = 1000)
}
\arguments{
  \item{GS}{a list, each element is an integer vector, indicating the indexes of genes in each gene set. See \emph{Details} below.}
  \item{GSNames}{a character string vector, each is the name of each gene set.}
  \item{GSDescs}{a character string vector, each is the description of each gene set.}
  \item{geneList}{a character string vector of gene IDs. See \emph{Details} below. }
  \item{name}{the name of this category of gene sets.}
  \item{sourceFile}{the source file name of this category of gene sets. }
  \item{GSSizeMin}{the minimum number of genes in a gene set to be analyzed. Default: 5}
  \item{GSSizeMax}{the maximum number of genes in a gene set to be analyzed. Default: 1000}
}
\details{
TBA
}
\value{
A SeqGeneSet object. 
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
gs 
## End 
}
