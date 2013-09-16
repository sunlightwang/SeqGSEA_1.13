# copyright: Xi Wang (xi.wang@newcastle.edu.au)
# integrating DE and DS for SeqGSEA
normFactor <- function(permStat) {
  rowMeans(permStat)
}

scoreNormalization <- function(scores, norm.factor) {
  stopifnot( nrow(as.matrix(scores)) == length(norm.factor) )
  scores <- scores / norm.factor
  scores[ is.na(scores) | is.infinite(scores) ] <- 0
  scores
}

geneScore <- function (DEscore, DSscore = NULL,
                       method = c("linear", "quadratic", "rank"),
                       DEweight = 0.5) {
# DEscore: differential expression score 
# DSscore: differential splicing score
# Method: linear: weigth * a + (1 - weight) * b
#         quadratic:  sqrt(weigth * a ^ 2 + (1 - weight) * b ^ 2)
#         rank: (rank_a * a + rank_b * b) / (rank_a + rank_b), where ranks are in ascending order
# DEweight: the weight for DE, (1-DEweight) for DS
  DEscore[is.na(DEscore) | is.infinite(DEscore)] <- 0
  if(DEweight == 1) return(DEscore)
  stopifnot(!is.null(DSscore))
  stopifnot(length(DEscore) == length(DSscore))
  DSscore[is.na(DSscore) | is.infinite(DSscore)] <- 0
  if(DEweight == 0) return(DSscore)

  method <- match.arg(method, c("linear", "quadratic", "rank"))
  switch(method, linear = {
    DEscore * DEweight + DSscore * (1 - DEweight)
  }, quadratic = {
    sqrt(DEscore^2 * DEweight + DSscore^2 * (1 - DEweight))
  }, rank = {
    DErank <- rank(DEscore, ties.method = "min")
    DSrank <- rank(DSscore, ties.method = "min")
    (DEscore * DErank * DEweight + DSscore * DSrank *
       (1 - DEweight))/(DErank * DEweight + DSrank * (1 - DEweight))
  })
}

genePermuteScore <- function (DEscoreMat, DSscoreMat = NULL,
                              method = c("linear", "quadratic", "rank"),
                              DEweight = 0.5) {
# parameters as function `geneScore`
  DEscoreMat[is.na(DEscoreMat)] <- 0
  if(DEweight == 1) return(DEscoreMat)
  stopifnot(!is.null(DSscoreMat))
  stopifnot(all(dim(DEscoreMat) == dim(DSscoreMat)))
  DSscoreMat[is.na(DSscoreMat)] <- 0
  if(DEweight == 0) return(DSscoreMat)

  method <- match.arg(method, c("linear", "quadratic", "rank"))
  switch(method, linear = {
    DEscoreMat * DEweight + DSscoreMat * (1 - DEweight)
  }, quadratic = {
    sqrt(DEscoreMat^2 * DEweight + DSscoreMat^2 * (1 - DEweight))
  }, rank = {
    DErankMat <- apply(DEscoreMat, 2, rank, ties.method = "min")
    DSrankMat <- apply(DSscoreMat, 2, rank, ties.method = "min")
    (DEscoreMat * DErankMat * DEweight + DSscoreMat * DSrankMat *
       (1 - DEweight))/(DErankMat * DEweight + DSrankMat *
                          (1 - DEweight))
  })
}

rankCombine <- function(DEscore, DSscore, DEscoreMat, DSscoreMat, DEweight=0.5) {
# combining DE and DS with the same weight according to rank in the observe data 
  stopifnot( length(DEscore) == length(DSscore))
  DEscore[ is.na(DEscore) | is.infinite(DEscore) ] <- 0
  DSscore[ is.na(DSscore) | is.infinite(DSscore) ] <- 0
  stopifnot(all(dim(DEscoreMat) == dim(DSscoreMat)))
  DEscoreMat[is.na(DEscoreMat)] <- 0
  DSscoreMat[is.na(DSscoreMat)] <- 0  
  DErank <- rank(DEscore, ties.method="min")
  DSrank <- rank(DSscore, ties.method="min")
  geneScore <- ( DEscore * DErank * DEweight + DSscore * DSrank * (1 - DEweight) ) / (DErank * DEweight + DSrank * (1 - DEweight) ) 
  geneScoreMat <- ( DEscoreMat * DErank * DEweight + DSscoreMat * DSrank * (1 - DEweight) ) / (DErank * DEweight + DSrank * (1 - DEweight) ) 
  list(geneScore=geneScore, genePermuteScore=geneScoreMat)
}

### dealing with geneset file ### 
convertEnsembl2Symbol <- function(ensembl.genes) {
  #require(biomaRt)
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  getBM(values = ensembl.genes, attributes = c('ensembl_gene_id','hgnc_symbol'), 
        filters = 'ensembl_gene_id', mart = ensembl, bmHeader=FALSE )
}

convertSymbol2Ensembl <- function(symbols) {
  #require(biomaRt)
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  getBM(values = symbols, attributes = c('hgnc_symbol', 'ensembl_gene_id'), 
        filters = 'hgnc_symbol', mart = ensembl, bmHeader=FALSE )
}

loadGenesets <- function(geneset.file, geneIDs, geneID.type=c("gene.symbol","ensembl"), 
                         genesetsize.min = 5, genesetsize.max = 1000) {
  # geneIDs can contain more than one genes, splited by '+' (because of HTSeq counting)
  # geneIDs can be in either "gene symbols" or "ensembl gene names"
  geneset.name <- basename(geneset.file)
  geneIDs <- unique(as.character(geneIDs))
  nGeneID <- length(geneIDs)
  geneID.type <- match.arg(geneID.type, c("gene.symbol","ensembl"))
  splitGeneIDs <- strsplit(as.character(geneIDs), "+", fixed=TRUE)
  uniGenes <- unlist(splitGeneIDs, use.names=FALSE)
  idxGenes <- rep(seq_along(splitGeneIDs), sapply(splitGeneIDs, length))
  stopifnot( length(uniGenes) == length(idxGenes))
  if (geneID.type == "ensembl") {
    temp = convertEnsembl2Symbol(uniGenes)
    idx <- data.frame(idx=idxGenes, ensembl = uniGenes, 
                      symbol = rep(NA_character_, length(uniGenes)),
                      row.names=uniGenes)
    idx$symbol <- temp$hgnc_symbol [ match( uniGenes, temp$ensembl_gene_id ) ]
    
    # solving duplicated mapping (say, one ensembl ID maps to multiply gene symbols)
    multi.map <- temp[duplicated(temp$ensembl_gene_id),]     # using 'duplicated' for identifying duplicated matches
    addi.idx <- idx[match( multi.map$ensembl_gene_id, uniGenes ), ]
    addi.idx$symbol <- multi.map$hgnc_symbol
    idx <- rbind(idx, addi.idx)
  } else {
    idx <- data.frame(idx=idxGenes, symbol = uniGenes, row.names=uniGenes)
  }
  gs.lines <- readLines(geneset.file)
  nGS <- length(gs.lines)
  gs <- vector("list", length=nGS)
  gs.name <- vector("character", length=nGS)
  gs.descs <- vector("character", length=nGS)
  for (i in 1:nGS) {
    strs <- noquote(unlist(strsplit(gs.lines[[i]], "\t")))
    gs.name[i] <- strs[1]
    gs.descs[i] <- strs[2] 
    temp.genes <- do.call(c, as.list(sapply(strs[3:length(strs)], function(x) {
      x <- gsub(" ", "", x)
      unlist(strsplit(as.character(x),"\\///")) # multiple genes are separated by "///" in one entry (str[j])
    }, USE.NAMES = FALSE)))
    gs[[i]] <- unique( idx$idx[ idx$symbol %in% temp.genes ] )
  } 
  gene.set <- newGeneSets(name = geneset.name, sourceFile = geneset.file, 
                  geneList=geneIDs, GS=gs, GSNames=gs.name, GSDescs=gs.descs, 
                  GSSizeMin=genesetsize.min, GSSizeMax=genesetsize.max) 
}

calES <- function(gene.set, gene.score, weighted.type=1) {
  ngene <- length(gene.score)
  nset <- size(gene.set)
  gene.set.size <- geneSetSize(gene.set)
  sort.idx <- order(gene.score, decreasing=TRUE)
  gene.score.sorted <- gene.score[sort.idx]
  cumsum.score <- matrix(0, nrow = ngene, ncol = nset)
  for (i in 1:nset) {
    ngene.hit <- gene.set.size[i]
    ngene.miss <- ngene - ngene.hit
    sort.idx.hit <- match(gene.set[i], sort.idx)
    cumsum.score[,i] <- - 1.0 / ngene.miss
    if( weighted.type == 0 ) {
      cumsum.score[sort.idx.hit, i] <- 1.0 / ngene.hit
    } else if ( weighted.type == 1 ) {
      gene.score.sorted.hit <- gene.score.sorted[sort.idx.hit]
      cumsum.score[sort.idx.hit, i] <- gene.score.sorted.hit / sum(gene.score.sorted.hit)
    } else {
      gene.score.sorted.hit <- gene.score.sorted[sort.idx.hit] ** weighted.type
      cumsum.score[sort.idx.hit, i] <- gene.score.sorted.hit / sum(gene.score.sorted.hit)      
    }
  }
  cumsum.score <- apply(cumsum.score, 2, cumsum)
  t(cumsum.score) # nset rows * ngene cols
}

calES.perm <- function(gene.set, gene.score.perm, weighted.type=1) {
  # only return ES (nset * times)
  ngene <- nrow(gene.score.perm)
  times <- ncol(gene.score.perm)
  nset <- size(gene.set)
  gene.set.size <- geneSetSize(gene.set)
  #ES.perm <- matrix(0, nrow = nset, ncol = times) # nset rows, times cols
  #ES.perm.vec <- matrix(0, nrow = nset, ncol = 1) # nset rows, 1 col
  #for (i in 1:times) {
  foreach (i = 1:times, .combine='cbind') %dopar% {
    sort.idx <- order(gene.score.perm[,i], decreasing=TRUE)
    gene.score.sorted <- gene.score.perm[sort.idx,i]
    #for (j in 1:nset) {
    foreach (j = 1:nset, .combine='c') %do% {
      ngene.hit <- gene.set.size[j]
      ngene.miss <- ngene - ngene.hit
      sort.idx.hit <- match(gene.set[j], sort.idx)
      score <- as.vector( rep( - 1.0 / ngene.miss, ngene ) )
      if( weighted.type == 0 ) {
        score[sort.idx.hit] <- 1.0 / ngene.hit
      } else if ( weighted.type == 1 ) {
        gene.score.sorted.hit <- gene.score.sorted[sort.idx.hit]
        score[sort.idx.hit] <- gene.score.sorted.hit / sum(gene.score.sorted.hit)
      } else {
        gene.score.sorted.hit <- gene.score.sorted[sort.idx.hit] ** weighted.type
        score[sort.idx.hit] <- gene.score.sorted.hit / sum(gene.score.sorted.hit)      
      }
      #ES.perm[j,i] <- max(cumsum(score))
      #ES.perm.vec[j,1] <- max(cumsum(score))
      max(cumsum(score))
    }
    #ES.perm.vec
  }
  #ES.perm
}

normES <- function(gene.set) {
  stopifnot(is(gene.set, "SeqGeneSet"))
  if(gene.set@GSEA.normFlag) 
    return(gene.set)
  ES.mean <- rowMeans(gene.set@GSEA.ES.perm) 
  gene.set@GSEA.ES <- gene.set@GSEA.ES / ES.mean
  gene.set@GSEA.ES.perm <- gene.set@GSEA.ES.perm / ES.mean
  gene.set@GSEA.normFlag <- TRUE
  gene.set
}

signifES <- function(gene.set) {
  stopifnot(is(gene.set, "SeqGeneSet"))
  if(! gene.set@GSEA.normFlag)
    gene.set <- normES(gene.set)
  times <- ncol(gene.set@GSEA.ES.perm)
  # pval
  pval <- rowSums(gene.set@GSEA.ES <= gene.set@GSEA.ES.perm) / times
  # FWER 
  each.perm.max.ES <- apply(gene.set@GSEA.ES.perm, 2, max)
  FWER <- sapply(gene.set@GSEA.ES, function(x) {
    sum(x <= each.perm.max.ES) / times
  })
  FWER <- ifelse(FWER > 1, 1, FWER)
  # FDR mean
  #FWER <- sapply(gene.set@GSEA.ES, function(x) {
  # mean( apply(x <= gene.set@GSEA.ES.perm, 2, sum) / 
  #   sum( x <= gene.set@GSEA.ES) )
  #})
  # FDR median
  FDR <- sapply(gene.set@GSEA.ES, function(x) {
    median( apply(x <= gene.set@GSEA.ES.perm, 2, sum) / 
      sum( x <= gene.set@GSEA.ES) )
  }) 
  FDR <- ifelse(FDR > 1, 1, FDR)
  
  gene.set@GSEA.pval <- signif(pval, 5)
  gene.set@GSEA.FWER <- signif(FWER, 5)
  gene.set@GSEA.FDR <- signif(FDR, 5)
  gene.set
}

GSEnrichAnalyze <- function(gene.set, gene.score, gene.score.perm, weighted.type=1) {
  stopifnot(is(gene.set, "SeqGeneSet"))
  stopifnot(all(names(gene.score) == gene.set@geneList))
  GSEA.score.cumsum <- calES(gene.set, gene.score, weighted.type=weighted.type)
  GSEA.ES <- apply(GSEA.score.cumsum, 1, max)
  GSEA.ES.pos <- apply(GSEA.score.cumsum, 1, which.max)
  GSEA.ES.perm <- calES.perm(gene.set, gene.score.perm, weighted.type=weighted.type)
  stopifnot(length(GSEA.ES) == nrow(GSEA.ES.perm))
  stopifnot(length(GSEA.ES) == length(GSEA.ES.pos))
  gene.set@GSEA.score.cumsum <- GSEA.score.cumsum
  gene.set@GSEA.ES <- GSEA.ES
  gene.set@GSEA.ES.pos <- GSEA.ES.pos
  gene.set@GSEA.ES.perm <- GSEA.ES.perm
  gene.set <- normES(gene.set)
  gene.set <- signifES(gene.set)
  gene.set
}

GSEAresultTable <- function(gene.set, GSDesc = FALSE)
{
  stopifnot( is( gene.set, "SeqGeneSet" ) )
  if(length(gene.set@GSEA.ES) == 0) stop("Please run GSEnrichAnalyze first.")
  result <- data.frame(
    GSName = gene.set@GSNames,
    GSSize = gene.set@GSSize,
    ES = gene.set@GSEA.ES,
    ES.pos = gene.set@GSEA.ES.pos,
    pval = gene.set@GSEA.pval, 
    FDR = gene.set@GSEA.FDR, 
    FWER = gene.set@GSEA.FWER
    )
  if(! GSDesc)
    return(result)
  data.frame(cbind(result, 
                   GSDesc = gene.set@GSDescs
                   ))
}

topGeneSets <- function(gene.set, n = 20, sortBy = c("FDR", "pvalue", "FWER"), GSDesc = FALSE) {
  stopifnot( is( gene.set, "SeqGeneSet" ) )
  sortBy <- match.arg(sortBy, c("FDR", "pvalue", "FWER"))
  res <- GSEAresultTable(gene.set, GSDesc)
  switch(sortBy, FDR = {
    res <- res[order(res$FDR)[1:n],]
  }, pvalue = {
    res <- res[order(res$pvalue)[1:n],] 
  }, FWER = {
    res <- res[order(res$FWER)[1:n],]
  })
  res
}

runSeqGSEA <- function(data.dir, case.pattern, ctrl.pattern, geneset.file, 
                       output.prefix, topGS=10, 
                       geneID.type=c("gene.symbol", "ensembl"), 
                       nCores=1, perm.times=1000, seed=NULL, minExonReadCount=5, 
                       integrationMethod=c("linear", "quadratic", "rank"),
                       DEweight=c(0.5), DEonly=FALSE, 
                       minGSsize=5, maxGSsize=1000, GSEA.WeightedType=1) 
## Assuming starting with exon reads counts, even for DEonly analysis
{
  # 0 # prepairation
  # input count data files
  case.files <- dir(data.dir, pattern=case.pattern, full.names = TRUE)
  control.files <- dir(data.dir, pattern=ctrl.pattern, full.names = TRUE)
  stopifnot (length(case.files)> 0)
  stopifnot (length(control.files)> 0)
  
  # setup parallel backend
  if (nCores > 1) {
    cl <- makeCluster(nCores)
    registerDoParallel(cl) # parallel backend registration
  }
  
  # 1 # DS analysis
  # load exon read count data
  RCS <- loadExonCountData(case.files, control.files)
  # remove genes with low exprssion
  RCS <- exonTestability(RCS, cutoff=minExonReadCount)
  geneTestable <- geneTestability(RCS)
  RCS <- subsetByGenes(RCS, unique(geneID(RCS))[ geneTestable ])
  # get gene IDs, which will be used in initialization of gene set
  geneIDs <- unique(geneID(RCS))
  permuteMat <- genpermuteMat(RCS, times=perm.times, seed=seed)
  
  if(! DEonly) {
    # calculate DS NB statistics
    RCS <- estiExonNBstat(RCS)
    RCS <- estiGeneNBstat(RCS)
    # calculate DS NB statistics on the permutation data sets
    RCS <- DSpermute4GSEA(RCS, permuteMat)    
  }

  # 2 # DE analysis
  # get gene read counts
  geneCounts <- getGeneCount(RCS)
  # calculate DE NB statistics
  label <- label(RCS)
  DEG <-runDESeq(geneCounts, label)
  DEGres <- DENBStat4GSEA(DEG)
  # calculate DE NB statistics on the permutation data sets
  DEpermNBstat <- DENBStatPermut4GSEA(DEG, permuteMat) # permutation
  
  # 3 #  score normalization
  # DE score normalization
  DEscore.normFac <- normFactor(DEpermNBstat)
  DEscore <- scoreNormalization(DEGres$NBstat, DEscore.normFac)
  DEscore.perm <- scoreNormalization(DEpermNBstat, DEscore.normFac)
  if(DEonly) {
    DSscore <- NULL
    DSscore.perm <- NULL
  } else {
    # DS score normalization
    DSscore.normFac <- normFactor(RCS@permute_NBstat_gene)
    DSscore <- scoreNormalization(RCS@featureData_gene$NBstat, DSscore.normFac)
    DSscore.perm <- scoreNormalization(RCS@permute_NBstat_gene, DSscore.normFac)
  }
  
  # visilization of DE & DS scores
  if(! DEonly) {
    plotGeneScore(DEscore, DEscore.perm, pdf=paste(output.prefix,".DEScore.pdf",sep=""), main="Expression")
    plotGeneScore(DSscore, DSscore.perm, pdf=paste(output.prefix,".DSScore.pdf",sep=""), main="Splicing")
  }
  # output DE and DS scores
  writeScores(DEscore, DSscore, file = paste(output.prefix,".DEDS_Score.txt",sep=""))
  
  if(DEonly) 
    DEweight <- 1
  GSEAres.list <- vector("list", length(DEweight))
  names(GSEAres.list) <- paste0("weight",sub(".", "_", DEweight, fixed =TRUE))
  for(i in 1:length(DEweight)) {
    if (DEweight[i] == 1) {
      output.prefix0 <- paste0(output.prefix, ".DEonly")
    } else if (DEweight[i] == 0) {
      output.prefix0 <- paste0(output.prefix, ".DSonly")
    } else 
      output.prefix0 <- paste0(output.prefix, ".weight", sub(".", "_", DEweight[i], fixed =TRUE))
    
    # 4 # score integration
    gene.score <- geneScore(DEscore, DSscore, method=integrationMethod, DEweight=DEweight[i])
    gene.score.perm <- genePermuteScore(DEscore.perm, DSscore.perm, method=integrationMethod, DEweight=DEweight[i])
    # visilization of gene scores
    plotGeneScore(gene.score, gene.score.perm, pdf=paste(output.prefix0,".GeneScore.pdf",sep=""))
    # output gene score
    writeScores(DEscore, DSscore, geneScore = gene.score, 
                geneScoreAttr = paste(integrationMethod, DEweight[i], sep=","),
                file = paste(output.prefix0,".GeneScore.txt",sep=""))
    
    # 5 #  main GSEA
    # load gene set data
    gene.set <- loadGenesets(geneset.file, geneIDs, geneID.type=geneID.type,
                             genesetsize.min = minGSsize, genesetsize.max = maxGSsize)
    # enrichment analysis
    gene.set <- GSEnrichAnalyze(gene.set, gene.score, gene.score.perm, weighted.type=GSEA.WeightedType)
    # format enrichment analysis results
    GSEAres <- GSEAresultTable(gene.set, TRUE)
    GSEAres.list[[i]] <- GSEAres
    
    # 6 #output results
    plotES(gene.set, pdf=paste(output.prefix0,".SeqGSEA.EnrichScore.pdf",sep=""))
    plotSig(gene.set, pdf=paste(output.prefix0,".SeqGSEA.FDR.pdf",sep="")) 
    write.table(GSEAres, paste(output.prefix0,".SeqGSEA.result.txt",sep=""),
                quote=FALSE, sep="\t", row.names=FALSE)
    topList <- order(GSEAres$FDR, GSEAres$pval)
    for(j in 1:min(topGS,length(topList))) {
      output.prefix00 <- paste0(output.prefix0, ".top_", j, "_GS_detail")
      plotSigGeneSet(gene.set, topList[j], gene.score, pdf=paste0(output.prefix00, ".pdf"))
      writeSigGeneSet(gene.set, topList[j], gene.score, file=paste0(output.prefix00, ".txt")) 
    }
  }
  GSEAres.list # return a list of GSEA results for meta analysis
}

