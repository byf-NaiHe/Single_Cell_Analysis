#!/usr/bin/env Rscript

library("ggpubr")
library("ggplot2")
library("data.table")
library("sscVis")
library("plyr")
library("ComplexHeatmap")
library("Startrac")
source("./func.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)

out.prefix <- "OUT_Fig3/Fig3"
dir.create(dirname(out.prefix),F,T)

### pTrans (fig. 25B, 26C)
{

    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    out <- readRDS("../data/tcr/startrac/CD4/CD4.out.startrac.nperm1000.rds")

    pIndex.sig.tran <- as.data.table(out@pIndex.sig.tran)

    dat.plot <- pIndex.sig.tran
    dat.plot[,index:=as.character(index)]
    dat.plot[,index.name:=fetchMetaClusterID2CusterFullName("cluster.name")[index] ]
    g.colSet$cluster.name <- g.colSet$meta.cluster
    names(g.colSet$cluster.name) <- fetchMetaClusterID2CusterFullName("cluster.name")[names(g.colSet$cluster.name)]
    f.col <- !is.na(names(g.colSet$cluster.name))
    g.colSet$cluster.name <- g.colSet$cluster.name[f.col]

    mcls.moi <- c("CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9")
    #mcls.moi <- c("CD8.c12.Tex.CXCL13")

    l_ply(mcls.moi,function(x){
        ###dat.tmp <- dat.plot[majorCluster==x & aid=="panC" & !is.na(value),]
        dat.tmp <- dat.plot[majorCluster==x & aid=="panC" &
                    !is.na(value),]
        dat.med.tmp <- dat.tmp[order(value),]
        dat.tmp[,index:=factor(index,levels=dat.med.tmp$index)]
        dat.tmp[,index.name:=factor(index.name,levels=dat.med.tmp$index.name)]
        p <- ggplot(dat.tmp, aes(index.name,value)) +
            geom_col(fill="steelblue",col=NA,width=0.8) +
            ##geom_hline(yintercept=0.10,linetype="dashed",color="black",alpha=0.8) +
            xlab("") + ylab(sprintf("pTrans Of %s",fetchMetaClusterID2CusterFullName("cluster.name")[x])) +
            theme_pubr() +
            theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))
        ##ggsave(sprintf("%s.pTrans.Fig.barplot.byPatientF.%s.pdf",out.prefix,x),width=5.5,height=4.5)
        ggsave(sprintf("%s.pTrans.Fig.barplot.byPatientF.%s.pdf",out.prefix,x),width=5.2,height=3.5)
        
	},.parallel=T)

}

