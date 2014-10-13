\name{runDESeq}
\alias{runDESeq}
\title{
Run DESeq for differential expression analysis 
}
\description{
This function provides a wrapper to run DESeq for differential expression analysis. It includes two steps, 
\code{DESeq::estimateSizeFactors} and \code{DESeq::estimateDispersions}. 
}
\usage{
runDESeq(geneCounts, label)
}
\arguments{
  \item{geneCounts}{a matrix containing read counts for each gene, can be the output of \code{\link{getGeneCount}}.}
  \item{label}{the sample classification labels.}
}
\value{
A CountDataSet object with size factors and dispersion parameters been estimated. 
}
\references{
Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data, Genome Biol, 11, R106.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{getGeneCount}},
\code{\link{DENBTest}},
\code{\link{DENBStat4GSEA}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneCounts <- getGeneCount(RCS_example)
label <- label(RCS_example)
DEG <- runDESeq(geneCounts, label)
}
