#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input sce file list")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-s", "--id", type="character", default="AID", help="a id")
parser$add_argument("-n", "--ncore", type="integer",default=12L, help="[default %(default)s]")
parser$add_argument("-m", "--mode", type="character",default="svm", help="one of ['svm','rf'] [default %(default)s]")
args <- parser$parse_args()
print(args)

args.infile <- args$inFile
args.outprefix <- args$outPrefix
args.ncore <- args$ncore
args.mode <- args$mode
args.id <- args$id

##library("scReClassify")
library("tictoc")
#library("ggpubr")
#library("ggplot2")
#library("ComplexHeatmap")
#library("RColorBrewer")
#library("circlize")
library("data.table")
library("plyr")
library("SingleCellExperiment")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

#out.prefix <- "OUT.dev.startrac2/test"
dir.create(dirname(args.outprefix),F,T)

mcrf <- function(ntree, ncore, ...)
{
    ntrees <- rep(ntree%/%ncore, ncore) + c(rep(0, ncore-ntree%%ncore), rep(1, ntree%%ncore))
    rfs <- mclapply(ntrees, function(nn) {
			a.mod <- randomForest::randomForest(ntree = nn, ...)
			return(a.mod)
    }, mc.cores = ncore)
    do.call(randomForest::combine, rfs)
}

my.multiAdaSampling <- function(data, label, seed=1, classifier="svm", percent=1, L=10,
				prob=FALSE, balance=FALSE, iter=3,
				rf.ntree=1000,ncore=12)
{
  models <- list()
  for(l in 1:L) {
    set.seed(seed+l)
    X <- data
    Y <- label

    model <- c()
    prob.mat <- c()
    for (i in 1:iter) {

      if (classifier == "rf") {
        ###model <- randomForest::randomForest(t(X), factor(Y), ntree = 100)
        model <- mcrf(ntree=rf.ntree,ncore=ncore, t(X), factor(Y))
        prob.mat <- predict(model, newdata=t(data), type="prob")
      }
      if (classifier == "svm") {
        tmp <- t(X)
        rownames(tmp) <- NULL
        model <- e1071::svm(tmp, factor(Y), probability = TRUE)
        prob.mat <- attr(predict(model, t(data), decision.values = FALSE, probability = TRUE), "probabilities")
      }

      X <- c()
      Y <- c()
      for(j in 1:ncol(prob.mat)) {
        voteClass <- prob.mat[label==colnames(prob.mat)[j],,drop=F]
        idx <- c()
        if (balance == FALSE) {
	  ####
          idx <- sample(1:nrow(voteClass), size=nrow(voteClass)*percent, replace = TRUE, prob=voteClass[,j])
        } else {
          sampleSize <- round(median(table(label)))
          if (nrow(voteClass) > sampleSize) {
            idx <- sample(1:nrow(voteClass), size=sampleSize*percent, replace = TRUE, prob=voteClass[,j])
          } else {
            idx <- sample(1:nrow(voteClass), size=nrow(voteClass)*percent, replace = TRUE, prob=voteClass[,j])
          }
        }

	### update X and Y
        X <- cbind(X, data[, rownames(voteClass)[idx]])
        Y <- c(Y, label[rownames(voteClass)[idx]])
      }
    }
    models[[l]] <- model
  }

  predictMat <- matrix(0, nrow=ncol(data), ncol=length(table(label)))
  final <- c()
  for (l in 1:L) {
    if (classifier == "svm") {
      tmp <- attr(predict(models[[l]], newdata=t(data), probability = TRUE), "prob")[,names(table(label))]
      predictMat <- predictMat + tmp
    } else {
      tmp <- predict(models[[l]], newdata=t(data), type="prob")[,names(table(label))]
      predictMat <- predictMat + tmp
    }
  }

  if(prob==TRUE) {
    final <- apply(predictMat, 1, max)
    names(final) <- names(table(label))[apply(predictMat, 1, which.max)]
  } else {
    final <- names(table(label))[apply(predictMat, 1, which.max)]
  }

  # return (final)
  return(
    list(
      final = final,
      models = models,
      prob = predictMat
    )
  )
}

{

    sce <- readRDS(args.infile)
    dat.pc <- t(reducedDim(sce,"seurat.pca"))[1:30,]
    rownames(dat.pc) <- gsub("_","",rownames(dat.pc))
    dim(dat.pc)
    dat.pc[1:4,1:3]
    cellLabel.ori <- as.character(sce$meta.cluster)
    names(cellLabel.ori) <- colnames(sce)

    tic("my.multiAdaSampling")
    cellLable.reclassify <- my.multiAdaSampling(data=dat.pc, label=cellLabel.ori, seed = 1,
						classifier = args.mode, percent = 1, L = 10,
						iter=3,
						ncore=args.ncore)
    toc()

    print(cellLable.reclassify$prob[1:2,])
    rowSums(cellLable.reclassify$prob) %>% str %>% print
    prob.mat <- sweep(cellLable.reclassify$prob,1,rowSums(cellLable.reclassify$prob),"/")
    prob.mat[1:2,1:4] %>% print
    saveRDS(cellLable.reclassify,sprintf("%s.scReclassify.%s.%s.rds",args.outprefix,args.mode,args.id))
    saveRDS(prob.mat,sprintf("%s.prob.mat.%s.%s.rds",args.outprefix,args.mode,args.id))

}

