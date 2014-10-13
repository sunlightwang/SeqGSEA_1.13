\name{loadGenesets}
\alias{loadGenesets}
\title{
Load gene sets from files 
}
\description{
This function is to load annotation of gene sets from files. The files are in the format of Molecular Signatures Database (MSigDB), and 
those files can be downloaded at \url{http://www.broadinstitute.org/gsea/msigdb/index.jsp}. 
}
\usage{
loadGenesets(geneset.file, geneIDs, geneID.type = c("gene.symbol", "ensembl"), 
             genesetsize.min = 5, genesetsize.max = 1000)
}
\arguments{
  \item{geneset.file}{the file containing the gene set annotation.}
  \item{geneIDs}{gene IDs that have expression values in the studied data set.}
  \item{geneID.type}{indicating the type of gene IDs, gene symbol or emsembl gene IDs.}
  \item{genesetsize.min}{the minimum number of genes in a gene set that will be treated in the analysis.}
  \item{genesetsize.max}{the maximum number of genes in a gene set that will be treated in the analysis.}
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
\code{\link{newGeneSets}},
\code{\link{SeqGeneSet-class}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneIDs <- geneID(RCS_example)
geneID.type <- "ensembl"
geneset.file <- system.file("extdata", "gs_symb.txt",  package="SeqGSEA", mustWork=TRUE)
##
GS <- loadGenesets(geneset.file, geneIDs, geneID.type = geneID.type)
GS
## end 
}
