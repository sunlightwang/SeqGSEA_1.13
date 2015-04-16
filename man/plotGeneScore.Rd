\name{plotGeneScore}
\alias{plotGeneScore}
\title{
Plot gene (DE/DS) scores 
}
\description{
This function is to plot gene scores, as well as DE scores and DS scores 
}
\usage{
plotGeneScore(score, perm.score = NULL, pdf = NULL, 
              main = c("Overall", "Expression", "Splicing"))
}
\arguments{
  \item{score}{the gene/DE/DS score vector.}
  \item{perm.score}{a matrix of the corresponding gene/DE/DS scores on the permutation data sets.}
  \item{pdf}{if a PDF file name provided, plot will be save to that file. }
  \item{main}{the key words representing the type of scores that will be shown in the plot main title.}
}
\details{
The plot shows the ranked scores from the largest to the smallest. Lines also show the maximum and average scores, values shown on the top left. 
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\examples{
data(DEscore, package="SeqGSEA")
plotGeneScore(DEscore, main="Expression")
data(DSscore, package="SeqGSEA")
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
plotGeneScore(gene.score)
}
