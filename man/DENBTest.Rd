\name{DENBTest}
\alias{DENBTest}
\title{
Perform negative binomial exact test for differential expression 
}
\description{
Perform negative binomial exact test for differential expression - a modified version of nbinomTest in DESeq package. 
}
\usage{
DENBTest(cds)
}
\arguments{
  \item{cds}{A CountDataSet object with size factors and dispersion parameters estimated. Recommended to take the output of \code{\link{runDESeq}}.}
}
\value{
A data frame of the test results. Information contains mean expression values, NB-statistics, (log) fold-changes, p-values, and adjusted p-values. 
}
\references{
Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data, Genome Biol, 11, R106.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{runDESeq}},
\code{\link{DENBStat4GSEA}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneCounts <- getGeneCount(RCS_example)
label <- label(RCS_example)
DEG <- runDESeq(geneCounts, label)
DEGres <- DENBTest(DEG)
head(DEGres)
}
