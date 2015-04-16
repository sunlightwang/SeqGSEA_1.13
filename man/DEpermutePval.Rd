\name{DEpermutePval}
\alias{DEpermutePval}
\title{
Permutation for p-values in differential expression analysis 
}
\description{
Calculate permutation p-values in differential expression analysis for each genes. 
}
\usage{
DEpermutePval(DEGres, permuteNBstat)
}
\arguments{
  \item{DEGres}{the output of \code{\link{DENBStat4GSEA}}.}
  \item{permuteNBstat}{the output of \code{\link{DENBStatPermut4GSEA}}.}
}
\value{
A data frame containing the expression means and variances for each gene in each group compared, and NB-stats, permutation p-values and adjusted p-values for each gene. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{runDESeq}},
\code{\link{DENBStat4GSEA}},
\code{\link{DENBStatPermut4GSEA}},
\code{\link{DENBTest}}
}
\examples{
data(RCS_example, package="SeqGSEA")
permuteMat <- genpermuteMat(RCS_example, times=10)
geneCounts <- getGeneCount(RCS_example)
label <- label(RCS_example)
DEG <- runDESeq(geneCounts, label)
DEGres <- DENBStat4GSEA(DEG)
DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat)
DEGres <- DEpermutePval(DEGres, DEpermNBstat) 
head(DEGres)
}
