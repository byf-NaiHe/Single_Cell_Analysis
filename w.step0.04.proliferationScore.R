#!/usr/bin/env Rscript

library("sscClust")
library("Seurat")
library("tictoc")
library("plyr")
library("dplyr")
library("tibble")
library("doParallel")
library("sscClust")
library("Matrix")
library("data.table")
library("R.utils")
library("gplots")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("cowplot")
library("limma")
library("reticulate")
library("ggrastr")
options(stringsAsFactors = FALSE)


out.prefix <- "../data/metaInfo/panC.proliferation"

dir.create(dirname(out.prefix),F,T)
list.file <- "list/obj.T.list"
list.tb <- fread(list.file,head=F)
colnames(list.tb) <- c("dataset","measurement","platform","seu.file","sce.file","opt.res")
gene.prol.tb <- fread("../data/external/gene.ProliferationScore.list")

colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 1)

######################
if(!file.exists(sprintf("%s.txt.gz",out.prefix))){

    out.tb <- as.data.table(ldply(seq_len(nrow(list.tb)),function(i){
        seu.file <- list.tb$seu.file[i]
        sce.file <- list.tb$sce.file[i]
        measurement <- list.tb$measurement[i]
        dataset <- list.tb$dataset[i]

        if(grepl("\\.rds$",seu.file)){
            seu <- readRDS(seu.file)
            sce <- readRDS(sce.file)
        }else{
            env.a <- loadToEnv(seu.file)
            env.b <- loadToEnv(sce.file)
            obj.name.a <- names(env.a)[1]
            obj.name.b <- names(env.b)[1]
            seu <- env.a[[obj.name.a]]
            sce <- env.b[[obj.name.b]]
            rm(env.a)
            rm(env.b)
        }

        if("percent.mito" %in% colnames(seu[[]])){
            if(max(seu$percent.mito) <1){
                seu  <- subset(seu, subset = percent.mito<0.1)
            }else{
                seu  <- subset(seu, subset = percent.mito<10)
            }
        }

        f.cell <- intersect(colnames(seu),colnames(sce))
        seu <- seu[,f.cell]
        sce <- sce[,colnames(seu)]

        if(measurement=="TPM"){
            #### TPM
            assay.name <- "log2TPM"
        }else{
            #### counts, cpm
            assay.name <- "norm_exprs"
        }

        tic("calProliferationScore")

        dat.out <- calProliferationScore(sce,assay.name,gene.prol=gene.prol.tb$geneSymbol,
                         out.prefix=sprintf("%s.%s",out.prefix,dataset),
                         method="mean")
        o.tb <- dat.out$out.tb
        o.tb$dataset <- dataset
        toc()

        return(o.tb)
    },.parallel=T))

    colnames(out.tb) <- c("cellID","proliferationScore.bin","proliferationScore","classification","dataset")

    conn <- gzfile(sprintf("%s.txt.gz",out.prefix),"w")
    write.table(out.tb,file=conn,row.names = F,quote = F,sep = "\t")
    close(conn)

}else{
    cat(sprintf("file %s.txt.gz exists!",out.prefix))
}

