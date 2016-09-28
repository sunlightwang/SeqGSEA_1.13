# SeqGSEA

The package generally provides methods for gene set enrichment analysis of high-throughput RNA-Seq data by integrating differential expression and splicing. It uses negative binomial distribution to model read count data, which accounts for sequencing biases and biological variation. Based on permutation tests, statistical significance can also be achieved regarding each gene's differential expression and splicing, respectively.


Author: Xi Wang (Xi.Wang at newcastle.edu.au)

Maintainer: Xi Wang (Xi.Wang at dkfz-heidelberg.de)


##Citation
Wang X and Cairns MJ (2013). “Gene Set Enrichment Analysis of RNA-Seq Data: Integrating Differential Expression and Splicing.” BMC Bioinformatics, 14(Suppl 5), pp. S16.

Wang X and Cairns MJ (2014). “SeqGSEA: a Bioconductor package for gene set enrichment analysis of RNA-Seq data integrating differential expression and splicing.” Bioinformatics, 30(12), pp. 1777-9.

* from within R, enter citation("SeqGSEA")

##Installation
To install this package, start R and enter:

source("https://bioconductor.org/biocLite.R")

biocLite("SeqGSEA")


##Documentation
To view documentation for the version of this package installed in your system, start R and enter:
browseVignettes("SeqGSEA")
