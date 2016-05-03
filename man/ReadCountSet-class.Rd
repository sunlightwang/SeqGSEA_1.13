\name{ReadCountSet-class}
\Rdversion{1.1}
\docType{class}
\alias{ReadCountSet-class}
\alias{ReadCountSet}

\title{Class \code{"ReadCountSet"}}
\description{
ReadCountSet class
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{newReadCountSet}}.
}
\section{Slots}{
  \describe{
    \item{\code{featureData_gene}:}{Object of class \code{"data.frame"}. Data for each genes. }
    \item{\code{permute_NBstat_exon}:}{Object of class \code{"matrix"}. NB statistics of exons on the permutation data sets. }
    \item{\code{permute_NBstat_gene}:}{Object of class \code{"matrix"}. NB statistics of genes on the permutation data sets. }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}. The read count data. }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}. Data for each samples. }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}. Data for each exons. }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}. Experiment data. }
    \item{\code{annotation}:}{Object of class \code{"character"}. Not used. }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}. Protocol information. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}. Version information. }
  }
}
\section{Methods}{
  \describe{
    \item{counts}{Get counts from a ReadCountSet object. See \code{\link{counts}}. }
    \item{counts<-}{Set counts to a ReadCountSet object. See \code{\link{counts}}. }
   }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, 14(Suppl 5):S16.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{newReadCountSet}},
\code{\link{loadExonCountData}},
\code{\link{exonID}},
\code{\link{geneID}},
\code{\link{counts-methods}},
\code{\link{label}},
\code{\link{subsetByGenes}}
}
\examples{
showClass("ReadCountSet")
}
\keyword{classes}
