#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input file")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-c", "--curDir", type="character", default=".", help="absolute to current directory")
parser$add_argument("-d", "--dataID", type="character",help="dataset id")
args <- parser$parse_args()
print(args)
out.prefix <- args$outPrefix
inFile <- args$inFile
data.id <- args$dataID
curDir <- args$curDir

library("sscVis")
library("Startrac")
library("tictoc")
library("ggpubr")
library("ggplot2")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("data.table")
library("plyr")
library("ggpubr")
#source("./func.R")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

dir.create(dirname(out.prefix),F,T)

{
    gset.prol.file <- sprintf("%s/../data/external/gene.ProliferationScore.list",curDir)
    gset.file <- sprintf("%s/../data/external/c2.cp.v7.0.symbols.gmt",curDir)
    gset.tb <- as.data.table(clusterProfiler::read.gmt(gset.file))
    gset.prol.tb <- fread(gset.prol.file)

    fname.TCR.signaling <- c("REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
			     "REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE",
			     "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
			     "REACTOME_DOWNSTREAM_TCR_SIGNALING",
			     "REACTOME_TCR_SIGNALING")

    gset.list <- llply(fname.TCR.signaling,function(x){
				gset.tb[term==x,][["gene"]]
			     })
    names(gset.list) <- fname.TCR.signaling
    gset.list[["proliferation"]] <- gset.prol.tb$geneSymbol
}

##### GSEA (by clusterProfiler)
{
    library("clusterProfiler")
    sce.x <- readRDS(sprintf("%s",inFile))
    gset.df <- ldply(names(gset.list),function(x){ data.table(gs_name=x,geneID=gset.list[[x]])  })
    
    dat.n <- assay(sce.x,"norm_exprs")
    rownames(dat.n) <- rowData(sce.x)$display.name
    #method.trans <- "zscore"
    #method.trans <- "rank"
    #method.trans <- "minusAvg"
    method.trans <- "zscore.mini"
    if(method.trans=="rank"){
        dat.z <- apply(dat.n,2,function(x){ rank(x,ties.method="random") })
    }else if(method.trans=="minusAvg"){
        dat.z <- t(scale(t(dat.n),scale=F))
    }else if(method.trans=="zscore"){
        dat.z <- t(scale(t(dat.n)))
    }else{
        #dat.z <- t(scale(t(dat.n)))
        sce.x.z <- ssc.scale(sce.x,gene.id=rownames(sce.x),assay.name="norm_exprs",do.scale=T)
        sce.x.mini <- ssc.average.cell(sce.x.z,assay.name="norm_exprs.scale",column="miniCluster",ret.type="sce")
        dat.z <- assay(sce.x.mini,"norm_exprs.scale")
        rownames(dat.z) <- rowData(sce.x.mini)$display.name
    }
    print(dat.z[1:4,1:5])
    do.gsea.parallel <- T

    tic("GSEA per cell ...")
    gsea.out.list <- llply(seq_len(ncol(dat.z)),function(i){
	gene.vec <- dat.z[,i]
	gene.vec <- sort(gene.vec,decreasing=T)

	gsea.out <- GSEA(gene.vec, TERM2GENE = gset.df,pvalueCutoff=1,
			 #by="DOSE",nPerm=1000,
			 verbose=F)
	
	return(gsea.out@result)
		       },.parallel=do.gsea.parallel)
    names(gsea.out.list) <- colnames(dat.z)
    toc()

    saveRDS(gsea.out.list,file=sprintf("%s.gsea.out.list.%s.%s.rds",out.prefix,method.trans,data.id))

    #### print some result
    gsea.act.NES.tb <- do.call(cbind,llply(seq_len(length(gsea.out.list)),function(i){
			     out.tb <- gsea.out.list[[i]][,"NES"]
			     names(out.tb) <- rownames(gsea.out.list[[i]])
			     return(out.tb)
		       }))
    colnames(gsea.act.NES.tb) <- names(gsea.out.list)
    print(gsea.act.NES.tb[1:4,1:5])

    gsea.act.pvalue.tb <- do.call(cbind,llply(seq_len(length(gsea.out.list)),function(i){
			     out.tb <- gsea.out.list[[i]][,"pvalue"]
			     names(out.tb) <- rownames(gsea.out.list[[i]])
			     return(out.tb)
		       }))
    colnames(gsea.act.pvalue.tb) <- names(gsea.out.list)
    print(gsea.act.pvalue.tb[1:4,1:5])

}

