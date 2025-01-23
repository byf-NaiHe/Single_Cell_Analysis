#!/usr/bin/env Rscript

library("ggpubr")
library("ggplot2")
library("data.table")
library("sscVis")
library("plyr")
library("Startrac")
library("RColorBrewer")
source("./func.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)


out.prefix <- "OUT_Fig2/Fig2"
dir.create(dirname(out.prefix),F,T)


# load data
{
    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    in.startrac.file.list <- list("CD8"="../data/tcr/startrac/CD8/CD8.out.startrac.byCancerType.rds",
                                  "CD8.onlyT"="../data/tcr/startrac/CD8/CD8.out.startrac.byCancerType.onlyT.rds",
                                  "CD4"="../data/tcr/startrac/CD4/CD4.out.startrac.byCancerType.rds",
                                  "CD4.onlyT"="../data/tcr/startrac/CD4/CD4.out.startrac.byCancerType.onlyT.rds")

}

#######  startrac pTrans across cancerTypes ########
## version 1 (used in paper) (Fig. 2D, Fig. 3C ??)
{
    z.max <- 3
    doit <- function(in.startrac.file,out.prefix,
                     stype="CD8",
                     mcls.x="CD8.c12.Tex.CXCL13",pdf.width=6.5,pdf.height=6,
                     th.pTrans=0.1,th.mcls.sd=0.00,do.clust.col=F)
    {
        dat.startrac <- readRDS(in.startrac.file)

        ####### pIndex (tran)
        {

            #mcls.x <- "CD8.c12.Tex.CXCL13"
            cancerType.vec <- names(dat.startrac)
            dat.plot.tb <- as.data.table(ldply(cancerType.vec,function(x){
                                   as.data.table(dat.startrac[[x]]@pIndex.sig.tran)[aid==x & majorCluster==mcls.x,]
                                                }))
            dat.plot.tb <- dat.plot.tb[!is.na(value) & index!="CD8.c03.Tm.RPS12",]
            dat.plot.dcast.tb <- dcast(data=dat.plot.tb,aid~index,value.var="value",fill=0)
    #        if(stype=="CD8"){
    #            dat.plot.dcast.tb <- dat.plot.dcast.tb[order(-(CD8.c11.Tex.PDCD1),-(CD8.c10.Trm.ZNF683)),]
    #        }

            dat.plot.mtx <- as.matrix(dat.plot.dcast.tb[,-1])
            rownames(dat.plot.mtx) <- dat.plot.dcast.tb[[1]]
            colnames(dat.plot.mtx)[colMaxs(dat.plot.mtx) < 0.1]
            colnames(dat.plot.mtx) <- fetchMetaClusterID2CusterFullName("cluster.name")[colnames(dat.plot.mtx)]

            mcls.info.tb <- fread("../data/metaInfo/name.conversion.txt")
            mcls.info.tb <- mcls.info.tb[grepl(sprintf("^%s",stype),meta.cluster,perl=T),]
            setkey(mcls.info.tb,"cluster.name")
            mcls.info.tb <- mcls.info.tb[colnames(dat.plot.mtx),]
            mcls.info.tb[,mcls.group:=factor(mcls.group,levels=c("Start","PathO","PathISG","PathTrm","PathTem","PathTex",
                                                                 "PathTfh","PathTreg"))]

            dat.plot.mtx.tmp <- t(scale(t(dat.plot.mtx)))
            dat.plot.mtx.tmp[ dat.plot.mtx.tmp > z.max ] <- z.max
            dat.plot.mtx.tmp[ dat.plot.mtx.tmp < -z.max ] <- -z.max

            #pTrans.hclust.row <- run.cutree(dat.plot.mtx,k=2,method.distance="cosine",method.hclust="ward.D2")
            f.mcls <- colSds(dat.plot.mtx) > th.mcls.sd
            colnames(dat.plot.mtx)[f.mcls]
            dat.plot.mtx.tmp[ abs(dat.plot.mtx) < th.pTrans] <- 0
            f.cancerType <- rowSds(dat.plot.mtx.tmp)==0
            dat.plot.mtx <- dat.plot.mtx[!f.cancerType,,drop=F]
            dat.plot.mtx.tmp <- dat.plot.mtx.tmp[!f.cancerType,,drop=F]
            pTrans.hclust.row <- run.cutree(dat.plot.mtx.tmp[,f.mcls],
                            k=4,method.distance="cosine",method.hclust="ward.D2")
            if(do.clust.col){
                hclust.col <- run.cutree(t(dat.plot.mtx.tmp),k=2,method.distance="cosine",method.hclust="ward.D2")
                hclust.col <- hclust.col$branch
            }else{
                hclust.col=F
            }
            #pTrans.hclust.row <- run.cutree(dat.plot.mtx,k=2,method.distance="",method.hclust="complete")
            #pTrans.hclust.row <- run.cutree(dat.plot.mtx[,mcls.info.tb[mcls.group %in% c("PathTem","PathTrm"),][["meta.cluster"]]]>0.1,
            #				    k=2,method.distance="spearman",method.hclust="ward.D2")

            sscVis::plotMatrix.simple(dat.plot.mtx.tmp,out.prefix=sprintf("%s.startrac.%s.ht.scaleT.c",
                                              out.prefix,mcls.x),
                          clust.row=pTrans.hclust.row$branch,
                          clust.column=hclust.col,
                          #clust.row=T,
                          #clust.column=F,
                          show.dendrogram=T,
                          z.lo=-z.max,z.hi=z.max,z.len=50,
                          #palatte=((scales::viridis_pal(option = "cividis"))(9)),
                          palatte=rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")),
                          par.legend=list(at=seq(-z.max,z.max,1)),
                          column.split=if(do.clust.col) NULL else mcls.info.tb$mcls.group,
                          par.heatmap=list(cex.row=1.2,
                                   border=F,
                                   row_names_gp = gpar(fontsize = 12),
                                   column_names_gp = gpar(fontsize = 12),
                                   column_title_gp=gpar(fontsize=0),
                                   column_gap=unit(0,"mm")),
                          pdf.width = pdf.width, pdf.height = pdf.height,exp.name="row Z Score(pTran)")

            if(F){
                sscVis::plotMatrix.simple(dat.plot.mtx,out.prefix=sprintf("%s.startrac.%s.ht.scaleF.c",
                                                  out.prefix,mcls.x),
                              clust.row=pTrans.hclust.row$branch,
                              clust.column=hclust.col,
                              #clust.row=T,
                              #clust.column=F,
                              show.dendrogram=T,
                              z.lo=0,z.hi=0.2,z.len=50,
                              ####palatte=((scales::viridis_pal(option = "cividis"))(9)),
                              palatte=rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")),
                              par.legend=list(at=seq(0,0.2,0.05)),
                              column.split=if(do.clust.col) NULL else mcls.info.tb$mcls.group,
                              par.heatmap=list(cex.row=1.2,
                                       border=F,
                                       row_names_gp = gpar(fontsize = 12),
                                       column_names_gp = gpar(fontsize = 12),
                                       column_title_gp=gpar(fontsize=0),
                                       column_gap=unit(0,"mm")),
                              pdf.width = pdf.width, pdf.height = pdf.height,exp.name="pTran")
            }

        }
    }

    {

        doit(in.startrac.file.list[["CD8"]],out.prefix=sprintf("%s.clustColF",out.prefix),
             stype="CD8",mcls.x="CD8.c12.Tex.CXCL13")
        doit(in.startrac.file.list[["CD8"]],out.prefix,stype="CD8",mcls.x="CD8.c12.Tex.CXCL13",
             do.clust.col=T)
        doit(in.startrac.file.list[["CD4"]],out.prefix,stype="CD4",mcls.x="CD4.c17.TfhTh1.CXCL13",
             pdf.width=7.5,do.clust.col=T)
        doit(in.startrac.file.list[["CD4"]],out.prefix,stype="CD4",mcls.x="CD4.c20.Treg.TNFRSF9",
             pdf.width=7.5,th.pTrans=0.01,do.clust.col=T)

    #    doit(in.startrac.file.list[["CD8.onlyT"]],out.prefix=sprintf("%s.onlyT.clustColF",out.prefix),
    #         stype="CD8",mcls.x="CD8.c12.Tex.CXCL13")
    #    doit(in.startrac.file.list[["CD8.onlyT"]],out.prefix=sprintf("%s.onlyT",out.prefix),
    #         stype="CD8",mcls.x="CD8.c12.Tex.CXCL13",do.clust.col=T)
    #    doit(in.startrac.file.list[["CD4.onlyT"]],out.prefix=sprintf("%s.onlyT",out.prefix),
    #         stype="CD4",mcls.x="CD4.c17.TfhTh1.CXCL13",pdf.width=7.5,do.clust.col=T)
    #    doit(in.startrac.file.list[["CD4.onlyT"]],out.prefix=sprintf("%s.onlyT",out.prefix),
    #         stype="CD4",mcls.x="CD4.c20.Treg.TNFRSF9",pdf.width=7.5,th.pTrans=0.01,do.clust.col=T)

    }
}

### pTrans of terminal Tex (Fig. 2C etc.)
{

    out <- readRDS("../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds")

    pIndex.sig.tran <- as.data.table(out@pIndex.sig.tran)

    dat.plot <- pIndex.sig.tran
    dat.plot[,index:=as.character(index)]
    dat.plot[,index.name:=fetchMetaClusterID2CusterFullName("cluster.name")[index] ]
    g.colSet$cluster.name <- g.colSet$meta.cluster
    names(g.colSet$cluster.name) <- fetchMetaClusterID2CusterFullName("cluster.name")[names(g.colSet$cluster.name)]
    f.col <- !is.na(names(g.colSet$cluster.name))
    g.colSet$cluster.name <- g.colSet$cluster.name[f.col]

    mcls.moi <- c("CD8.c12.Tex.CXCL13","CD8.c15.ISG.IFIT1","CD8.c14.Tex.TCF7","CD8.c17.Tm.NME1")
    #mcls.moi <- c("CD8.c12.Tex.CXCL13")

    l_ply(mcls.moi,function(x){
        ###dat.tmp <- dat.plot[majorCluster==x & aid=="panC" & !is.na(value),]
        dat.tmp <- dat.plot[majorCluster==x & aid=="panC" &
                    ##### uncharacterized cluster
                    !(index %in% c("CD8.c03.Tm.RPS12")) &
                    !is.na(value),]
        dat.med.tmp <- dat.tmp[order(value),]
        dat.tmp[,index:=factor(index,levels=dat.med.tmp$index)]
        dat.tmp[,index.name:=factor(index.name,levels=dat.med.tmp$index.name)]
        p <- ggplot(dat.tmp, aes(index.name,value)) +
            geom_col(fill="steelblue",col=NA,width=0.8) +
            geom_hline(yintercept=0.10,linetype="dashed",color="black",alpha=0.8) +
            xlab("") + ylab(sprintf("pTrans Of %s",fetchMetaClusterID2CusterFullName("cluster.name")[x])) +
            theme_pubr() +
            theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1))
        ggsave(sprintf("%s.pTrans.Fig.barplot.byPatientF.%s.pdf",out.prefix,x),width=5.5,height=4.5)
        if(mcls.moi=="CD8.c12.Tex.CXCL13"){
            #ggsave(sprintf("%s.pTrans.Fig.barplot.byPatientF.%s.Fig2C.pdf",out.prefix,x),width=3.5,height=3.5)
            ggsave(sprintf("%s.pTrans.Fig.barplot.byPatientF.%s.Fig2C.v2.pdf",out.prefix,x),width=3.8,height=3.8)
        }
	},.parallel=T)

}

### pTrans of terminal Tex, per tumor (fig. S15 A)
{

    out <- readRDS("../data/tcr/startrac/CD8/CD8.out.onlyT.startrac.nperm1000.rds")
    
    pIndex.sig.tran <- as.data.table(out@pIndex.sig.tran)
    dat.plot <- pIndex.sig.tran[aid!="panC" &
		    ##!(index %in% c("CD8.c12.Tex.CXCL13",mcls.exclude)) &
            ##### uncharacterized cluster
		    !(index %in% c("CD8.c12.Tex.CXCL13","CD8.c03.Tm.RPS12")) &
		    majorCluster=="CD8.c12.Tex.CXCL13",]
    dat.plot[,index:=as.character(index)]
    dat.plot[,index.clusterName:=fetchMetaClusterID2CusterFullName()[as.character(index)]]

    pTrans.mcls.med <- dat.plot[,.(pTrans.med=median(value)),
				by=c("majorCluster","index","index.clusterName")][order(-pTrans.med),]
    dat.plot[,index.fac:=factor(as.character(index),
				levels=pTrans.mcls.med$index)]
    dat.plot[,index.clusterName.fac:=factor(as.character(index.clusterName),
					   levels=pTrans.mcls.med$index.clusterName)]
    head(dat.plot,n=3)

    p <- ggboxplot(dat.plot,x="index.clusterName.fac",y="value",
		   color = "index", legend="none",title="",
		   fill = "index",
		   alpha=0.8,
		   xlab="",ylab="pTrans Of CD8.c12(terminal Tex)",
		   add = "jitter",outlier.shape=NA,
		   add.params=list(size=0.5)) +
	    scale_color_manual(values=g.colSet$meta.cluster) +
	    scale_fill_manual(values=g.colSet$meta.cluster) +
	    geom_hline(yintercept=c(0.1),linetype="dashed",
		       color="black",alpha=0.8)+
	    #stat_compare_means(label="p.format",comparisons=list(c("P","N"),c("N","T"),c("P","T"))) +
	    coord_cartesian(clip="off") +
	    theme(plot.title = element_text(hjust = 0.5,size=14),
		  axis.text.x=element_text(angle=60,hjust=1))
    ggsave(sprintf("%s.pTrans.Tex.boxplot.byTumor.v2.pdf", out.prefix), width=4.3,height=5,useDingbats=F)

}


### barplots showing clonotypes spanning more than two clusters and containing terminal Tex, cells from tumors only (fig. S15 B ?)
{

    ### clonotypes
    {
        tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
        in.dat <- readRDS(tcr.file)
        in.dat$majorCluster <- as.character(in.dat$meta.cluster)
        in.dat$clone.id <- in.dat$cloneID
        #### only T cells in tumors
        in.dat <- in.dat[loc %in% c("T"),]
        #### filter out patient.cluster with number of cell < 10
        ncell.patient.cluster <- sort(unclass(in.dat[,table(sprintf("%s.%s",patient,majorCluster))]))
        ##ncell.patient.cluster <- sort(unclass(in.dat[stype=="CD8", table(sprintf("%s.%s",patient,majorCluster))]))
        in.dat <- in.dat[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
        ####
    }

    ### probability matrix
    {
        prob.svm.mat.list <- llply(unique(in.dat$dataset.old),function(x){
        prob.svm.mat <- readRDS(sprintf("../data/expression/CD8/integration/reClassify.v1/T.CD8.%s/T.CD8.prob.mat.svm.%s.rds",x,x))
                        })
        names(prob.svm.mat.list) <- unique(in.dat$dataset.old)
    }

    #### sce data
    {
        sce.list <- list()
        for(data.id in unique(in.dat$dataset))
        {
            if(data.id=="NSCLC.XinyiGuo2018") { t.id <- "LC.XinyiGuo2018" 
            }else if(data.id=="STAD.BoxiKang2020"){ t.id <- "STAD.BoxiKang2019"
            }else{ t.id <- data.id }
            sce.list[[data.id]] <- readRDS(sprintf("../data/expression/CD8/byDataset/%s.sce.rds",t.id))
        }
    }

    ### clonotype level pIndex.tran
    {

        obj.star <- new("Startrac",in.dat[stype=="CD8",],aid="panC")
        
        object <- obj.star

        ##### main 2 paths (used in paper fig. S15B ??)
        {

            #### 
            l.colSet <- g.colSet
            l.colSet$meta.cluster <- l.colSet$meta.cluster[c("CD8.c02.Tm.IL7R","CD8.c10.Trm.ZNF683","CD8.c12.Tex.CXCL13")]
            out.Path.Trm <- ana.clonotypeAcrossMcls.moreThanTwo(obj.star,
                    in.dat,
                    out.prefix=sprintf("%s/PathTrm.newSize/%s",dirname(out.prefix),basename(out.prefix)),
                    aid="PathTrm",
                    par.barplot=list(pdf.width.byPatientF=7,pdf.height.byPatientF=4.5,
                             pdf.width.byPatientT=5.0,pdf.height.byPatientT=5),
                    show.N=T,
                    l.colSet=l.colSet,
                    mcls.moi=c("CD8.c02.Tm.IL7R","CD8.c10.Trm.ZNF683","CD8.c12.Tex.CXCL13"))
            out.Path.Trm$clone.info.flt.tb[,-c("G.obs")]
            
            #### 
            l.colSet <- g.colSet
            l.colSet$meta.cluster <- l.colSet$meta.cluster[c("CD8.c06.Tem.GZMK","CD8.c11.Tex.PDCD1", "CD8.c12.Tex.CXCL13")]
            out.Path.Tem <- ana.clonotypeAcrossMcls.moreThanTwo(obj.star,
                               in.dat,
                               out.prefix=sprintf("%s/PathTem.newSize/%s",dirname(out.prefix),basename(out.prefix)),
                               aid="PathTem",
                               par.barplot=list(pdf.width.byPatientF=7,pdf.height.byPatientF=4.5,
                                    pdf.width.byPatientT=5,pdf.height.byPatientT=5),
                               show.N=T,
                               l.colSet=l.colSet,
                               mcls.moi=c("CD8.c06.Tem.GZMK","CD8.c11.Tex.PDCD1", "CD8.c12.Tex.CXCL13"))
            out.Path.Tem$clone.info.flt.tb[,-c("G.obs")]
            #### 

        }

        ### connection with ISG (fig. S16 C ?? )
        {
           
            l.colSet <- g.colSet
            l.colSet$meta.cluster <- l.colSet$meta.cluster[c("CD8.c10.Trm.ZNF683","CD8.c15.ISG.IFIT1","CD8.c12.Tex.CXCL13")]
            out.Path.ISG.Trm <- ana.clonotypeAcrossMcls.moreThanTwo(obj.star,
                    in.dat,
                    out.prefix=sprintf("%s/PathISG.Trm/%s",dirname(out.prefix),basename(out.prefix)),
                    aid="PathISG",
                    par.barplot=list(pdf.width.byPatientF=7,pdf.height.byPatientF=6.5,
                             pdf.width.byPatientT=7,pdf.height.byPatientT=7),
                    l.colSet=l.colSet,show.N=T,
                    mcls.moi=c("CD8.c10.Trm.ZNF683","CD8.c15.ISG.IFIT1","CD8.c12.Tex.CXCL13"))
            out.Path.ISG.Trm$clone.info.flt.tb[,-c("G.obs")]

            l.colSet <- g.colSet
            l.colSet$meta.cluster <- l.colSet$meta.cluster[c("CD8.c11.Tex.PDCD1","CD8.c15.ISG.IFIT1","CD8.c12.Tex.CXCL13")]
            out.Path.ISG.GZMKTex <- ana.clonotypeAcrossMcls.moreThanTwo(obj.star,
                    in.dat,
                    out.prefix=sprintf("%s/PathISG.GZMKTex/%s",dirname(out.prefix),basename(out.prefix)),
                    aid="PathISG",
                    par.barplot=list(pdf.width.byPatientF=7,pdf.height.byPatientF=6.5,
                             pdf.width.byPatientT=7,pdf.height.byPatientT=7),
                    l.colSet=l.colSet,show.N=T,
                    mcls.moi=c("CD8.c11.Tex.PDCD1","CD8.c15.ISG.IFIT1","CD8.c12.Tex.CXCL13"))
            out.Path.ISG.GZMKTex$clone.info.flt.tb[,-c("G.obs")]

        }
        
    }

}


### the two paths in individual tumors (fig. S15C ?)
{

    {
        tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
        in.dat <- readRDS(tcr.file)

        in.dat$majorCluster <- as.character(in.dat$meta.cluster)
        in.dat$clone.id <- in.dat$cloneID

        #### filter out patient.cluster with number of cell < 10
        ncell.patient.cluster <- sort(unclass(in.dat[,table(sprintf("%s.%s",patient,majorCluster))]))
        #ncell.patient.cluster <- sort(unclass(in.dat[stype=="CD8",table(sprintf("%s.%s",patient,majorCluster))]))
        in.dat <- in.dat[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
        ####
        dim(in.dat)
        ### [1] 164271     33

        out <- readRDS(sprintf("../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds"))
        out.loc.list <- llply(c("P","N","T"),function(a.loc) {
            readRDS(sprintf("../data/tcr/startrac/CD8/CD8.out.only%s.startrac.nperm1000.rds",a.loc))
        })
        names(out.loc.list) <- c("P","N","T")
    }

    loc.i <- "T"
    cex.point <- 0.5
    in.dat.flt <- in.dat[stype=="CD8" & loc==loc.i,]
    dim(in.dat.flt)
    ## [1] 52268    33
    ncell.patient.cluster <- sort(unclass(in.dat.flt[,table(sprintf("%s.%s",patient,majorCluster))]))
    in.dat.flt <- in.dat.flt[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
    dim(in.dat.flt)
    ## [1] 51447    33

    patient2cancerType.tb <- in.dat.flt[,.N,by=c("patient","cancerType")][order(cancerType,patient),]
    TexSize.tb <- in.dat.flt[majorCluster=="CD8.c12.Tex.CXCL13",.(N.Tex=.N),
                 by=c("patient","cancerType")][order(cancerType,patient),]
    patient2cancerType.tb <- merge(patient2cancerType.tb,TexSize.tb,by=c("patient","cancerType"))

    do.it <- function(out.prefix,obj.startrac.out,mcls.pairs,vtype="cancerType")
    {
        require("ComplexHeatmap")
        require("cowplot")
        plot.list <- llply(seq_along(mcls.pairs),function(i){
            mcls.1 <- mcls.pairs[[i]][1]
            mcls.2 <- mcls.pairs[[i]][2]
            #mcls.1 <- "CD8.c10.Trm.ZNF683"
            #mcls.2 <- "CD8.c11.Tex.PDCD1"
            get_density <- function(x, y, ...) {
              dens <- MASS::kde2d(x, y, ...)
              ix <- findInterval(x, dens$x)
              iy <- findInterval(y, dens$y)
              ii <- cbind(ix, iy)
              return(dens$z[ii])
            }

            dat.Tex.pTran.2path.a.tb <- as.data.table(obj.startrac.out)[aid!="panC",
                                  ][majorCluster=="CD8.c12.Tex.CXCL13" &
                                  index %in% c(mcls.1,mcls.2),]
            dat.Tex.pTran.2path.b.tb <- dcast(data=dat.Tex.pTran.2path.a.tb,
                              formula=aid+majorCluster~index,value.var="value")
            colnames(dat.Tex.pTran.2path.b.tb)[c(1,3,4)] <- c("patient","mcls.1","mcls.2")
            dat.Tex.pTran.2path.b.tb[is.na(mcls.1),mcls.1:=0]
            dat.Tex.pTran.2path.b.tb[is.na(mcls.2),mcls.2:=0]
            dat.Tex.pTran.2path.b.tb <- merge(dat.Tex.pTran.2path.b.tb,
                              patient2cancerType.tb,by="patient")
            dat.Tex.pTran.2path.b.tb[,size:=(log10(N.Tex))*cex.point]
            dat.Tex.pTran.2path.b.tb[,cor.test(mcls.1,mcls.2)]
            dat.Tex.pTran.2path.b.tb[,density2D:=get_density(mcls.1,mcls.2,n=100)]
            setkey(dat.Tex.pTran.2path.b.tb,"patient")

            #### heatmap
            ##> head(dat.Tex.pTran.2path.b.tb)
            ##        patient       majorCluster     mcls.1     mcls.2 cancerType   N N.Tex      size density2D
            ##1: BC.P20190123 CD8.c12.Tex.CXCL13 0.14599497 0.04651163       BRCA 893    26 0.7074867  8.456591
            ##2:    CRC.P0701 CD8.c12.Tex.CXCL13 0.04445360 0.00000000        CRC  85    58 0.8817140  7.559410
            dat.plot.mat <- as.matrix(dat.Tex.pTran.2path.b.tb[,c("mcls.1","mcls.2"),with=F])
            rownames(dat.plot.mat) <- dat.Tex.pTran.2path.b.tb$patient
            ## pattern sort
            dat.plot.mat.pattern <- as.data.table(dat.plot.mat > 0.1)
            dat.plot.mat.pattern$rid <- rownames(dat.plot.mat)
            dat.plot.mat.pattern$V.mcls.1 <- dat.plot.mat[,"mcls.1"]
            dat.plot.mat.pattern$V.mcls.2 <- dat.plot.mat[,"mcls.2"]
            dat.plot.mat.pattern <- dat.plot.mat.pattern[order(-mcls.1,-mcls.2,-V.mcls.1,-V.mcls.2),]
            dat.plot.mat.pattern[,Group:=sprintf("%s%s",as.integer(mcls.1),as.integer(mcls.2))]
            dat.plot.mat.pattern[,Group:=factor(Group,levels=c("11","10","01","00"))]
            dat.plot.mat <- dat.plot.mat[dat.plot.mat.pattern$rid,,drop=F]
            colnames(dat.plot.mat) <- c(mcls.1,mcls.2)

            write.table(dat.plot.mat,file=sprintf("%s.mat.dump.txt",out.prefix),row.names=T,sep="\t",quote=F)

            dat.fisher <- matrix(dat.plot.mat.pattern[,.N,by=c("Group")][["N"]],nrow=2)
            res.fisher <- fisher.test(dat.fisher)

            ann.cancerType <- rowAnnotation(cate=dat.plot.mat.pattern$Group,
                            cancerType=dat.Tex.pTran.2path.b.tb[dat.plot.mat.pattern$rid,][["cancerType"]],
                                  col=c(g.colSet["cancerType"],
                        list("cate"=structure(RColorBrewer::brewer.pal(4,"Set2"),names=c("11","10","01","00")))
                        ),
                                  #annotation_legend_param=list(comb.ES=list(at = seq(ES.range[1],ES.range[2],ES.step),
                                  #                                     grid_width = unit(0.4,"cm"),
                                  #                                     grid_height = unit(0.4, "cm"),
                                  #                                     legend_height = unit(6, "cm"))),
                                  border=F)

            sscVis:::plotMatrix.simple(dat.plot.mat,out.prefix=sprintf("%s.mat.heatmap",out.prefix),
                           z.lo = 0, z.hi = 0.2,
                           row.split=dat.plot.mat.pattern$Group,
                           palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                           par.legend=list(at=seq(0,0.2,0.05)),
                           par.heatmap=list(
                               cex.row=0.5,
                                                   row_title_rot=0,border=T,
                               #row_names_gp=gpar(fontsize = 10),
                               column_names_gp=gpar(fontsize=12),
                                                   #raster_device="png",
                                                   #raster_quality = 5,
                                                   #right_annotation =  ha,
                                                   left_annotation = ann.cancerType
                                                               ),
                           exp.name="pTrans",pdf.width=4,pdf.height=8)


            #### scatter plot
            p <- ggscatter(dat.Tex.pTran.2path.b.tb,x="mcls.1",y="mcls.2",
                   size="size",
                   #color="cancerType",
                   #color="density2D",
                   color=vtype,
                   alpha=0.75) +
                labs(x=sprintf("pTrans(%s)",mcls.1),
                     y=sprintf("pTrans(%s)",mcls.2),
                     title="") +
                #geom_density_2d_filled(alpha = 0.5,contour_var="density") +
                #geom_density_2d() +
                geom_hline(yintercept=c(0.2),linetype="dashed",
                       color="gray",alpha=0.8) +
                geom_vline(xintercept=c(0.2),linetype="dashed",
                       color="gray",alpha=0.8) +
                geom_smooth(method='lm') +
                scale_size(breaks=log10(c(10,30,50,100))*cex.point,
                       range=c(0.2,6),labels=c(10,30,50,100),
                       limits=c(0,10)*cex.point) +
                stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
                stat_regline_equation(aes(label =  paste(..eq.label..)),vjust=2.5) +
                coord_cartesian(clip="off") +
                theme(legend.position = "right",
                  plot.title = element_text(hjust = 0.5))

            if(vtype=="cancerType"){
                p <- p + guides(size = guide_legend(order=2,direction="horizontal",
                                    title.position = "top",
                                    label.position = "bottom"),
                       color=guide_legend(order=1,override.aes=list(size=5),ncol=2),
                       legend.direction="horizontal") +
                     scale_color_manual(values=g.colSet[["cancerType"]])
            }else{
                p <- p + guides(size = guide_legend(order=2,direction="horizontal",
                                       title.position = "top",
                                       label.position = "bottom"),
                        legend.direction="horizontal") +
                       viridis::scale_color_viridis()
            }
            #pe <- ggExtra::ggMarginal(p)
            pp <- plot_grid(ggExtra::ggMarginal(p+theme(legend.position = "none"),
                            colour = '#FF0000',
                            fill = '#FFA500',
                            size=3),
                    get_legend(p),
                    ncol=1,
                    rel_heights=c(3.75,ifelse(vtype=="cancerType",3.75,5)),
                    align="v")
#            ggsave(sprintf("%s.pTrans.onlyT.%s.%s.%s.pdf",out.prefix,mcls.1,mcls.2,vtype),
#                   width=4.0,height=ifelse(vtype=="cancerType",7.5,8.75),useDingbats=F)

            return(p)
        })
        names(plot.list) <- names(mcls.pairs)

#        saveRDS(plot.list,sprintf("%s.pTrans.onlyT.%s.%s.rds",out.prefix,vtype,"plot.list"))

    }

    do.it(out.prefix=sprintf("%s",out.prefix),
      obj.startrac.out=out.loc.list$T@pIndex.sig.tran,
      mcls.pairs=list("pair1"=c("CD8.c10.Trm.ZNF683","CD8.c11.Tex.PDCD1")))

}

##########################################

