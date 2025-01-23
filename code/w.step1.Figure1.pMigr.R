#!/usr/bin/env Rscript

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
library("cowplot")
source("./func.R")
##source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/TCR.ana.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

#### CD8+ T cells
{
    startrac.file <- "../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds"
    startrac.byCancerType.file <- "../data/tcr/startrac/CD8/CD8.out.startrac.byCancerType.rds"
    #tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
    out.prefix <- "./OUT_Fig1/pMigr.CD8"
    dir.create(dirname(out.prefix),F,T)

    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

    ### load data
    out <- readRDS(startrac.file)
    res.byCancerType <- readRDS(startrac.byCancerType.file)
    cancerType.vec <- names(res.byCancerType)

    ######## pMigr, panC, (fig. S11A)
    {

        ###
        dat.plot <- as.data.table(out@pIndex.sig.migr)[aid!="panC" & index %in% c("N-P","P-T"),]
        dat.plot <- dat.plot[order(majorCluster),]
        dat.plot[,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name")[as.character(majorCluster)]]
        dat.plot[,cluster.name:=factor(cluster.name,levels=unique(cluster.name))]
        #dat.plot[is.na(value),value:=0]
        dat.plot <- dat.plot[!is.na(value),]
        dat.plot[index=="N-P",index:="P-N"]
        dat.plot[,index:=factor(index,levels=c("P-N","P-T"))]
        dat.sig.tb <- as.data.table(ldply(unique(dat.plot$index),function(x){
                          aa <- dat.plot[index==x & majorCluster=="CD8.c07.Temra.CX3CR1",][["value"]]
                          dat.plot[majorCluster!="CD8.c07.Temra.CX3CR1",
                                {
                                bb <- .SD$value
                                p.value <- as.double(NA)
                                if(length(aa)>3 && length(bb)>3){
                                    res <- wilcox.test(aa,bb)
                                    p.value <- as.double(res$p.value)
                                }
                                p.char <- if(is.na(p.value)) "na" else if(p.value<0.001) "***" else if(p.value<0.01) "**" else if(p.value < 0.05) "*" else "ns"
                                .(N=.N,N.aa=length(aa),N.bb=length(bb),p.value=p.value,
                                  p.char=p.char)
                                }, by=c("majorCluster","cluster.name","index")]
                               }))
        write.table(dat.sig.tb,file=sprintf("%s.Fig01.S05.01.boxplot.test.txt",out.prefix),row.names=F,sep="\t",quote=F)
        
        #dat.plot[,majorCluster:=factor(majorCluster,levels=sort(unique(majorCluster)))]
        p <- ggboxplot(dat.plot, x="cluster.name",y="value",
               xlab="",ylab="pMigr", add = "jitter",outlier.shape=NA,ylim=c(0,1),
               group="index",color="index",fill="index",alpha=0.5) +
            facet_wrap(~index,ncol=1) +
            scale_color_manual(values=c("P-N"="#4DBBD5FF","P-T"="#E64B35FF"))+
            scale_fill_manual(values=c("P-N"="#4DBBD5FF","P-T"="#E64B35FF"))+
            geom_text(data=dat.sig.tb,aes(label=p.char,y=0.9)) +
            geom_hline(yintercept=0.1,linetype="dashed",color="gray") +
            #stat_compare_means(ref.group="Temra/Teff",aes(label = ..p.signif..)) +
            theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
              legend.position="none",
              strip.text=element_text(size=12),
              strip.background=element_blank())
        ggsave(sprintf("%s.Fig01.S05.01.boxplot.pdf",out.prefix),width=6,height=6,useDingbats=F)
        saveRDS(p,file=sprintf("%s.Fig01.S05.01.boxplot.rds",out.prefix))


    }

    ##### pMigr, individual cancer types (fig. S11B)
    {

        z.max <- 3

        ######## pIndex (migr)
        dat.plot.tb <- as.data.table(ldply(cancerType.vec,function(x){
                           if(nrow(as.data.table(res.byCancerType[[x]]@pIndex.sig.migr))>0){
                               as.data.table(res.byCancerType[[x]]@pIndex.sig.migr)[aid==x,]
                           }}))
        dat.plot.tb[,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name")[majorCluster]]

        #### P-N
        dat.plot.dcast.tb <- dcast(data=dat.plot.tb[index=="N-P",],aid~cluster.name,value.var="value",fill=0)
        dat.plot.mtx <- as.matrix(dat.plot.dcast.tb[,-1])
        rownames(dat.plot.mtx) <- dat.plot.dcast.tb[[1]]
        dat.plot.mtx[is.na(dat.plot.mtx)] <- 0
        ###sscVis::plotMatrix.simple(dat.plot.mtx,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.NP.ht.scaleF",out.prefix),
        ###					      clust.row=T,clust.column=T,show.dendrogram=T,
        ###					      pdf.width = 8, pdf.height = 5,exp.name="pIndex(migr.NP)")
        dat.plot.mtx.tmp <- t(scale(t(dat.plot.mtx)))
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp > z.max ] <- z.max
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp < -z.max ] <- -z.max
        dat.plot.mtx.tmp[ abs(dat.plot.mtx) < 0.1 ] <- 0
        sscVis::plotMatrix.simple(dat.plot.mtx.tmp,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.NP.ht.scaleT",out.prefix),
                                  clust.row=T,clust.column=T,
                                  show.dendrogram=T,
                                  z.lo=-z.max,z.hi=z.max,z.len=50,
                                  palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                                  par.legend=list(at=seq(-z.max,z.max,1)),
                                  par.heatmap=list(row_names_gp=gpar(fontsize = 12),
                                           column_names_gp=gpar(fontsize = 12)),
                                  pdf.width = 8, pdf.height = 5,exp.name="pIndex(migr.NP)")

        ### P-T
        dat.plot.dcast.tb <- dcast(data=dat.plot.tb[index=="P-T",],aid~cluster.name,value.var="value",fill=0)
        dat.plot.mtx <- as.matrix(dat.plot.dcast.tb[,-1])
        rownames(dat.plot.mtx) <- dat.plot.dcast.tb[[1]]
        dat.plot.mtx[is.na(dat.plot.mtx)] <- 0
        ###sscVis::plotMatrix.simple(dat.plot.mtx,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.PT.ht.scaleF",out.prefix),
        ###					      clust.row=T,clust.column=T,show.dendrogram=T,
        ###					      pdf.width = 8, pdf.height = 5.5,exp.name="pIndex(migr.PT)")
        dat.plot.mtx.tmp <- t(scale(t(dat.plot.mtx)))
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp > z.max ] <- z.max
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp < -z.max ] <- -z.max
        dat.plot.mtx.tmp[ abs(dat.plot.mtx) < 0.1 ] <- 0
        sscVis::plotMatrix.simple(dat.plot.mtx.tmp,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.PT.ht.scaleT",out.prefix),
                                  clust.row=T,clust.column=T,show.dendrogram=T,
                                  z.lo=-z.max,z.hi=z.max,z.len=50,
                                  palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                                  par.legend=list(at=seq(-z.max,z.max,1)),
                                  par.heatmap=list(row_names_gp=gpar(fontsize = 12),
                                           column_names_gp=gpar(fontsize = 12)),
                                  pdf.width = 8, pdf.height = 5.5,exp.name="pIndex(migr.PT)")

    }

}

#### CD4+ T cells
{
    startrac.file <- "../data/tcr/startrac/CD4/CD4.out.startrac.nperm1000.rds"
    startrac.byCancerType.file <- "../data/tcr/startrac/CD4/CD4.out.startrac.byCancerType.rds"
    #tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
    out.prefix <- "./OUT_Fig1/pMigr.CD4"
    dir.create(dirname(out.prefix),F,T)

    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

    ### load data
    out <- readRDS(startrac.file)
    res.byCancerType <- readRDS(startrac.byCancerType.file)
    cancerType.vec <- names(res.byCancerType)

    ######## pMigr, panC, (fig. S11A)
    {

        ###
        dat.plot <- as.data.table(out@pIndex.sig.migr)[aid!="panC" & index %in% c("N-P","P-T"),]
        dat.plot <- dat.plot[order(majorCluster),]
        dat.plot[,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name")[as.character(majorCluster)]]
        dat.plot[,cluster.name:=factor(cluster.name,levels=unique(cluster.name))]
        #dat.plot[is.na(value),value:=0]
        dat.plot <- dat.plot[!is.na(value),]
        dat.plot[index=="N-P",index:="P-N"]
        dat.plot[,index:=factor(index,levels=c("P-N","P-T"))]
        dat.sig.tb <- as.data.table(ldply(unique(dat.plot$index),function(x){
                              aa <- dat.plot[index==x & majorCluster=="CD4.c13.Temra.CX3CR1",][["value"]]
                              dat.plot[majorCluster!="CD4.c13.Temra.CX3CR1",
                                {
                                    bb <- .SD$value
                                    p.value <- as.double(NA)
                                    if(length(aa)>3 && length(bb)>3){
                                    res <- wilcox.test(aa,bb)
                                    p.value <- as.double(res$p.value)
                                    }
                                    p.char <- if(is.na(p.value)) "na" else if(p.value<0.001) "***" else if(p.value<0.01) "**" else if(p.value < 0.05) "*" else "ns"
                                    .(N=.N,N.aa=length(aa),N.bb=length(bb),p.value=p.value,
                                      p.char=p.char)
                                }, by=c("majorCluster","cluster.name","index")]
                               }))
        write.table(dat.sig.tb,file=sprintf("%s.Fig01.S05.01.boxplot.test.txt",out.prefix),row.names=F,sep="\t",quote=F)

        #dat.plot[,majorCluster:=factor(majorCluster,levels=sort(unique(majorCluster)))]
        p <- ggboxplot(dat.plot, x="cluster.name",y="value",
                   xlab="",ylab="pMigr", add = "jitter",outlier.shape=NA,ylim=c(0,1),
                   group="index",color="index",fill="index",alpha=0.5) +
            facet_wrap(~index,ncol=1) +
            scale_color_manual(values=c("P-N"="#4DBBD5FF","P-T"="#E64B35FF"))+
            scale_fill_manual(values=c("P-N"="#4DBBD5FF","P-T"="#E64B35FF"))+
            geom_text(data=dat.sig.tb,aes(label=p.char,y=0.9)) +
            geom_hline(yintercept=0.1,linetype="dashed",color="gray") +
            theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
                  legend.position="none",
                  strip.text=element_text(size=12),
                  strip.background=element_blank())
        ggsave(sprintf("%s.Fig01.S05.01.boxplot.pdf",out.prefix),width=6,height=6,useDingbats=F)
        saveRDS(p,file=sprintf("%s.Fig01.S05.01.boxplot.rds",out.prefix))

    }

    ##### pMigr, individual cancer types (fig. S11B)
    {

        z.max <- 3

        ######## pIndex (migr)
        dat.plot.tb <- as.data.table(ldply(cancerType.vec,function(x){
                               if(nrow(as.data.table(res.byCancerType[[x]]@pIndex.sig.migr))>0){
                               as.data.table(res.byCancerType[[x]]@pIndex.sig.migr)[aid==x,]
                               }}))
        dat.plot.tb[,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name")[majorCluster]]

        #### P-N
        dat.plot.dcast.tb <- dcast(data=dat.plot.tb[index=="N-P",],aid~cluster.name,value.var="value",fill=0)
        dat.plot.mtx <- as.matrix(dat.plot.dcast.tb[,-1])
        rownames(dat.plot.mtx) <- dat.plot.dcast.tb[[1]]
        dat.plot.mtx[is.na(dat.plot.mtx)] <- 0
        #f.aid <- apply(dat.plot.mtx,1,function(x){ all(x==0) })
        #dat.plot.mtx <- dat.plot.mtx[!f.aid,]
        ###sscVis::plotMatrix.simple(dat.plot.mtx,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.NP.ht.scaleF",out.prefix),
        ###                          clust.row=T,clust.column=T,show.dendrogram=T,
        ###                          pdf.width = 8, pdf.height = 5,exp.name="pIndex(migr.NP)")
        dat.plot.mtx.tmp <- t(scale(t(dat.plot.mtx)))
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp > z.max ] <- z.max
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp < -z.max ] <- -z.max
        dat.plot.mtx.tmp[ abs(dat.plot.mtx) < 0.1 ] <- 0
        sscVis::plotMatrix.simple(dat.plot.mtx.tmp,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.NP.ht.scaleT",out.prefix),
                                  clust.row=T,clust.column=T,show.dendrogram=T,
                                  z.lo=-z.max,z.hi=z.max,z.len=50,
                                  palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                                  par.legend=list(at=seq(-z.max,z.max,1)),
                                  par.heatmap=list(row_names_gp=gpar(fontsize = 12),
                                           column_names_gp=gpar(fontsize=12)),
                                  pdf.width = 8, pdf.height = 5,exp.name="pIndex(migr.NP)")

        #### P-T
        dat.plot.dcast.tb <- dcast(data=dat.plot.tb[index=="P-T",],aid~cluster.name,value.var="value",fill=0)
        dat.plot.mtx <- as.matrix(dat.plot.dcast.tb[,-1])
        rownames(dat.plot.mtx) <- dat.plot.dcast.tb[[1]]
        dat.plot.mtx[is.na(dat.plot.mtx)] <- 0
        #f.aid <- apply(dat.plot.mtx,1,function(x){ all(x==0) })
        #dat.plot.mtx <- dat.plot.mtx[!f.aid,]
        ###sscVis::plotMatrix.simple(dat.plot.mtx,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.PT.ht.scaleF",out.prefix),
        ###                          clust.row=T,clust.column=T,show.dendrogram=T,
        ###                          pdf.width = 8, pdf.height = 5.5,exp.name="pIndex(migr.PT)")
        dat.plot.mtx.tmp <- t(scale(t(dat.plot.mtx)))
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp > z.max ] <- z.max
        dat.plot.mtx.tmp[ dat.plot.mtx.tmp < -z.max ] <- -z.max
        dat.plot.mtx.tmp[ abs(dat.plot.mtx) < 0.1 ] <- 0
        sscVis::plotMatrix.simple(dat.plot.mtx.tmp,out.prefix=sprintf("%s.CmpCancerType.pIndex.migr.PT.ht.scaleT",out.prefix),
                                  clust.row=T,clust.column=T,show.dendrogram=T,
                                  z.lo=-z.max,z.hi=z.max,z.len=50,
                                  palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                                  par.legend=list(at=seq(-z.max,z.max,1)),
                                  par.heatmap=list(row_names_gp=gpar(fontsize = 12),
                                           column_names_gp=gpar(fontsize=12)),
                                  pdf.width = 8, pdf.height = 5.5,exp.name="pIndex(migr.PT)")

    }

}

#### merged (fig. S11A)
{
    out.prefix <- "./OUT_Fig1/pMigr.merged"
    a.file <- "./OUT_Fig1/pMigr.CD8.Fig01.S05.01.boxplot.rds"
    b.file <- "./OUT_Fig1/pMigr.CD4.Fig01.S05.01.boxplot.rds"
    obj.a <- readRDS(a.file)
    obj.b <- readRDS(b.file)
    pp <- plot_grid(obj.a,obj.b,align="hv",nrow=1,rel_widths=c(0.8,1),
		    label_y=1,label_x=-0.15,
		    labels=c("CD8+ compartment","CD4+ compartment"))
    ggsave(sprintf("%s.startrac.pMigr.boxplot.01.pdf",out.prefix),width=12,height=6,useDingbats=F)
}


###############################


