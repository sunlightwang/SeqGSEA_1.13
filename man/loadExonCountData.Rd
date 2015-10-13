\name{loadExonCountData}
\alias{loadExonCountData}
\title{
Load Exon Count Data
}
\description{
  This function is used to load (sub-)exon count data. Exon count data can 
be got by the Python script \code{count_in_exons.py}. 
}
\usage{
loadExonCountData(case.files, control.files)
}
\arguments{
  \item{case.files}{a character vector containing the exon count file names for case samples}
  \item{control.files}{a character vector containing the exon count file names for control samples}
}
\details{
You may need the Python script count_in_exons.py (released with this package) 
to generate your exon count files from read mapping results (say BAM files). 
The detailed usage can be obtained by simply typing 
\code{ python \\path\\to\\count_in_exons.py}. Users can also use other scripts or 
software for exon read counting.

The format of the exon count file is: 

\preformatted{
GeneName1:001[tab]Count11
GeneName1:002[tab]Count12
...
GeneName1:00N[tab]Count1N
GeneName2:001[tab]Count21
...
}
}
\value{
This function returns a ReadCountSet object.
}
\author{
Xi Wang, xi.wang@newcastle.edu.au
}
\seealso{
 \code{\link{newReadCountSet}},
 \code{\link{ReadCountSet-class}}
}
\examples{
library(SeqGSEA)
dat.dir = system.file("extdata", package="SeqGSEA", mustWork=TRUE)
case.pattern <- "^SC"
ctrl.pattern <- "^SN"
case.files <- dir(dat.dir, pattern=case.pattern, full.names = TRUE)
control.files <- dir(dat.dir, pattern=ctrl.pattern, full.names = TRUE)

## Not run: 

RCS <- loadExonCountData(case.files, control.files)
RCS 

## End (Not run)
}
