\name{writeScores}
\alias{writeScores}
\title{
Write DE/DS scores and gene scores
}
\description{
This function is to write DE and DS scores, and optionally gene scores. 
}
\usage{
writeScores(DEscore, DSscore, geneScore=NULL, geneScoreAttr=NULL, file="")
}
\arguments{
  \item{DEscore}{normalized DE scores. } 
  \item{DSscore}{normalized DS scores. }
  \item{geneScore}{gene scores integrated from DE and DS scores. }
  \item{geneScoreAttr}{the parameters for integrating DE and DS scores. }
  \item{file}{output file name, if not specified print to screen.}
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
\code{\link{DEscore}},
\code{\link{geneScore}}
}
\examples{
data(DEscore, package="SeqGSEA")
data(DSscore, package="SeqGSEA")
gene.score <- geneScore(DEscore, DSscore, method="linear", DEweight = 0.3)
writeScores(DEscore, DSscore) # without gene scores
writeScores(DEscore, DSscore, geneScore = gene.score, 
            geneScoreAttr = "linear,0.3") # gene scores with attr.
}
