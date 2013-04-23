\name{DENBStat4GSEA}
\alias{DENBStat4GSEA}
\title{
Calculate NB-statistics quantifying differential expression for each gene
}
\description{
Calculate NB-statistics quantifying differential expression between two groups of samples compared. 
The results will be used for GSEA run. Comparing with \code{\link{DENBTest}}, this function will not 
calculate NB test p-values. 

This function only works with two-group comparison. 
}
\usage{
DENBStat4GSEA(cds)
}
\arguments{
  \item{cds}{A CountDataSet object with size factors and dispersion parameters estimated. Recommended to take the output of \code{\link{runDESeq}}. }
}
\note{
The results with the output of \code{\link{DENBStatPermut4GSEA}} can also be used to run \code{\link{DEpermutePval}}.
}
\value{
A data frame containing each gene's expression means and variances in each group, and each gene's DE NB-statistics. 
}
\references{
Xi Wang and Murray J. Cairns (2013). Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing. BMC Bioinformatics, 14(Suppl 5):S16.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{DENBTest}},
\code{\link{runDESeq}},
\code{\link{DENBStatPermut4GSEA}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneCounts <- getGeneCount(RCS_example)
label <- label(RCS_example)
DEG <- runDESeq(geneCounts, label)
DEGres <- DENBStat4GSEA(DEG)
head(DEGres)
}
