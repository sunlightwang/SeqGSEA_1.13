\name{SeqGeneSet-class}
\Rdversion{1.1}
\docType{class}
\alias{SeqGeneSet-class}
\alias{SeqGeneSet}
\alias{[,SeqGeneSet,numeric-method}
\alias{show,SeqGeneSet-method}

\title{Class \code{"SeqGeneSet"}}
\description{
SeqGeneSet class
}
\section{Objects from the Class}{
Objects can be created by calls of the function \code{\link{newGeneSets}}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} the name of this gene set category }
    \item{\code{sourceFile}:}{Object of class \code{"character"} the source file of gene set category  }
    \item{\code{geneList}:}{Object of class \code{"character"} the gene ID list indicating genes involved in this GSEA }
    \item{\code{GS}:}{Object of class \code{"list"} a list of gene indexes corresponding to \code{geneList}, each element in the list 
                      indicating which genes are in each gene set of this SeqGeneSet object }
    \item{\code{GSNames}:}{Object of class \code{"character"}. Gene set names. }
    \item{\code{GSDescs}:}{Object of class \code{"character"}. Gene set descriptions. }
    \item{\code{GSSize}:}{Object of class \code{"numeric"}. Gene set sizes. }
    \item{\code{GSSizeMin}:}{Object of class \code{"numeric"}. The minimum gene set size to be analyzed. }
    \item{\code{GSSizeMax}:}{Object of class \code{"numeric"}. The maximum gene set size to be analyzed. }
    \item{\code{GS.Excluded}:}{Object of class \code{"list"}. Gene sets excluded to be analyzed. }
    \item{\code{GSNames.Excluded}:}{Object of class \code{"character"}. Gene set names excluded to be analyzed. }
    \item{\code{GSDescs.Excluded}:}{Object of class \code{"character"}. Gene set descriptions excluded to be analyzed. }
    \item{\code{GSEA.ES}:}{Object of class \code{"numeric"}. Enrichment scores. }
    \item{\code{GSEA.ES.pos}:}{Object of class \code{"numeric"}. The positions where enrichment scores appear. }
    \item{\code{GSEA.ES.perm}:}{Object of class \code{"matrix"}. The enrichment scores of the permutation data sets. }
    \item{\code{GSEA.score.cumsum}:}{Object of class \code{"matrix"}. Running enrichment scores. }
    \item{\code{GSEA.normFlag}:}{Object of class \code{"logical"}. Logical indicating whether \code{GSEA.ES} has been normalized. }
    \item{\code{GSEA.pval}:}{Object of class \code{"numeric"}. P-values of each gene set. }
    \item{\code{GSEA.FWER}:}{Object of class \code{"numeric"}. Family-wise error rate of each gene set. }
    \item{\code{GSEA.FDR}:}{Object of class \code{"numeric"}. False discovery rate of each gene set. }
    \item{\code{version}:}{Object of class \code{"Versions"}. Version information. }
  }
}
\section{Methods}{
  \describe{
    \item{[}{Get a sub-list of gene sets, and return a SeqGeneSet object. }
    \item{show}{Show basic information of the SeqGeneSet object. }
   }
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, 14(Suppl 5):S16.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{newGeneSets}},
\code{\link{size}},
\code{\link{geneSetNames}},
\code{\link{geneSetDescs}},
\code{\link{geneSetSize}}
}
\examples{
showClass("SeqGeneSet")
}
\keyword{classes}
