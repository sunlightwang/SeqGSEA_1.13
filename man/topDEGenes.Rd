\name{topDEGenes}
\alias{topDEGenes}
\title{
Extract top differentially expressed genes. 
}
\description{
This function is to extract top n differentially expressed genes, ranked by either DESeq p-values, 
DESeq adjusted p-values, permutation p-values, permutation adjusted p-values, or NB-statistics. 
}
\usage{
topDEGenes(DEGres, n = 20, 
           sortBy = c("padj", "pval", "perm.pval", "perm.padj", "NBstat", "foldChange"))}
\arguments{
  \item{DEGres}{DE analysis results. }
  \item{n}{the number of top DE genes.}
  \item{sortBy}{indicating which method to rank genes.}
}
\details{
If the sortBy method is not among the column names, the function will result in an error. 
}
\value{
A table for top n DE genes with significance metrics.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{topDSGenes}},
\code{\link{topDSExons}}
}
\examples{
data(RCS_example, package="SeqGSEA")
geneCounts <- getGeneCount(RCS_example)
label <- label(RCS_example)
DEG <- runDESeq(geneCounts, label)
permuteMat <- genpermuteMat(RCS_example, times=10)
DEGres <- DENBTest(DEG)
DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat)
DEGres <- DEpermutePval(DEGres, DEpermNBstat) 
topDEGenes(DEGres, n = 10, sortBy = "NBstat")
}
