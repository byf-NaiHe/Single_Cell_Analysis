#!/usr/bin/env Rscript

# import libraries
install.packages("magrittr")
install.packages("sscVis")
install.packages("data.table")
install.packages("R.utils")
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("plyr")
install.packages("grid")
install.packages("cowplot")
install.packages("ggrepel")


library("magrittr")
library("sscVis")
library("data.table")
library("R.utils")
library("ggpubr")
library("ggplot2")
library("plyr")
library("grid")
library("cowplot")
library("ggrepel")
source("./func.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)

out.prefix <- "OUT_Fig1/Fig1"
dir.create(dirname(out.prefix),F,T)

# load data
{
    dir.metaInfo <- "../data/metaInfo"
    meta.tb <- readRDS(sprintf("%s/panC.freq.all.meta.tb.rds",dir.metaInfo))
    freq.all.ht.tb <- readRDS("../data/metaInfo/panC.freq.all.ht.tb.rds")
    freq.CD8.ht.tb <- readRDS("../data/metaInfo/panC.freq.CD8.ht.tb.rds")
    freq.CD4.ht.tb <- readRDS("../data/metaInfo/panC.freq.CD4.ht.tb.rds")

    dir.int.list <- list("CD8"="../data/expression/CD8/integration",
                         "CD4"="../data/expression/CD4/integration")
    #sce.merged.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds"),
    #                        "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.sce.merged.rds"))
    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    colSet.CD8 <- g.colSet["meta.cluster"]
    colSet.CD8$meta.cluster <- colSet.CD8$meta.cluster[ grepl("^CD8\\.c",names(colSet.CD8$meta.cluster))]
    colSet.CD4 <- g.colSet["meta.cluster"]
    colSet.CD4$meta.cluster <- colSet.CD4$meta.cluster[ grepl("^CD4\\.c",names(colSet.CD4$meta.cluster))]
    colSet.list <- list("CD8"=colSet.CD8,"CD4"=colSet.CD4)

    #gene.desc.tb.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds"),
    #                          "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.gene.tb.rds"))
}


############ cluster composition (fig. S2F, fig. S3F ??)
{
    ### cancer types (CD8)
    dat.plot <- meta.tb[stype=="CD8",.(N=.N),by=c("meta.cluster","cancerType")]
    dat.plot <- dat.plot[,.(cancerType=.SD$cancerType,N=.SD$N,freq=.SD$N/sum(.SD$N)),by="meta.cluster"]
    dat.plot[,meta.cluster:=factor(meta.cluster,levels=sort(unique(meta.cluster),decreasing=T))]
    dat.plot[,cancerType:=factor(cancerType,levels=sort(unique(cancerType),decreasing=F))]
    p <- ggbarplot(dat.plot,x="meta.cluster",y="freq",fill="cancerType",color=NA,
			       position=position_stack(reverse=T),
			       ylab="Frequency",xlab="") +
		    scale_fill_manual(values=g.colSet$cancerType) +
		    theme(legend.position="right") +
		    coord_flip()
    ggsave(sprintf("%s.clusterComp.cancerType.CD8.pdf",out.prefix),width=8,height=5)

    ### cancer types (CD4)
    dat.plot <- meta.tb[stype=="CD4",.(N=.N),by=c("meta.cluster","cancerType")]
    dat.plot <- dat.plot[,.(cancerType=.SD$cancerType,N=.SD$N,freq=.SD$N/sum(.SD$N)),by="meta.cluster"]
    dat.plot[,meta.cluster:=factor(meta.cluster,levels=sort(unique(meta.cluster),decreasing=T))]
    dat.plot[,cancerType:=factor(cancerType,levels=sort(unique(cancerType),decreasing=F))]
    p <- ggbarplot(dat.plot,x="meta.cluster",y="freq",fill="cancerType",color=NA,
			       position=position_stack(reverse=T),
			       ylab="Frequency",xlab="") +
		    scale_fill_manual(values=g.colSet$cancerType) +
		    theme(legend.position="right") +
		    coord_flip()
    ggsave(sprintf("%s.clusterComp.cancerType.CD4.pdf",out.prefix),width=8,height=5)

}

### data statistics (fig. S01)
{
    
    meta.tb[,table(cancerType)]
    meta.tb[,.N,by=c("cancerType","dataset")][order(cancerType),]
    meta.tb[,table(dataset,dataSource)]
    unique(meta.tb[,c("cancerType","dataset","patient"),with=F])[,.N,by=c("cancerType","dataset")][order(cancerType),]

    #### cell number
    meta.tb[1:2,]
    meta.tb[dataset=="NSCLC.XinyiGuo2018" & patient=="LUNG.P0512",pub:="thisStudy"]
    meta.tb[,table(pub)]
    meta.tb[,table(pub)/.N]
    ### pub
    ### published thisStudy
    ### 0.5354239 0.4645761

    makePlotDataStat <- function(m.tb,out.prefix,plot.mode="numOfCell",...)
    {

        ##dat.plot.stat <- meta.tb[,.(numOfCell=.N),by=c("cancerType","pub","patient.uid")]
        dat.plot.stat <- m.tb[,.(numOfCell=.N),by=c("cancerType","pub","patient.uid")]
        dat.plot.stat[,pub:=factor(pub,levels=c("thisStudy","published"))]
        dat.plot.stat.nDonor <- dat.plot.stat[,.(numOfDonor=.N),by=c("cancerType","pub")]
        dat.plot.stat.nCell <- dat.plot.stat[,.(numOfCell=sum(.SD$numOfCell)),by=c("cancerType","pub")]
        dat.plot.stat.nDonor.sum <- dat.plot.stat.nDonor[,.(numOfDonor=sum(.SD$numOfDonor)),by=c("cancerType")
                                 ][order(-numOfDonor),]
        dat.plot.stat.nCell.sum <- dat.plot.stat.nCell[,.(numOfCell=sum(.SD$numOfCell)),by=c("cancerType")
                                 ][order(-numOfCell),]
        dat.plot.stat.nDonor[,cancerType:=factor(cancerType,levels=dat.plot.stat.nDonor.sum$cancerType)]
        dat.plot.stat.nCell[,cancerType:=factor(cancerType,levels=dat.plot.stat.nCell.sum$cancerType)]

        
        
        if(plot.mode=="numOfCell"){
            p <- plotNightingaleRose(dat.plot.stat.nCell,...)
                         ##y.pretty=NULL,y.lim.min=NULL,
                         ##empty_bar=2
            ggsave(sprintf("%s.nCell.NRose.pdf",out.prefix),width=8,height=7)
        }else if(plot.mode=="numOfDonor"){
            p <- plotNightingaleRose(dat.plot.stat.nDonor,...)
                         ##empty_bar=2,y.colum="numOfDonor",
                         ##y.pretty=NULL,y.lim.min=NULL
            ggsave(sprintf("%s.nDonor.NRose.pdf",out.prefix),width=8,height=7)
        }

    }
   
    ### all 
    makePlotDataStat(m.tb=meta.tb,out.prefix=sprintf("%s.dataStat",out.prefix),
		     plot.mode="numOfCell",my.title="Number Of Cells (all)")
    makePlotDataStat(m.tb=meta.tb,out.prefix=sprintf("%s.dataStat",out.prefix),
		     plot.mode="numOfDonor",y.colum="numOfDonor",my.title="Number Of Donors (all)",
		     exp.tick=0.3,y.pretty=NULL,y.lim.min=NULL)

    #### Droplet-based
    makePlotDataStat(m.tb=meta.tb[tech.cate=="Droplet",],out.prefix=sprintf("%s.dataStat.Droplet",out.prefix),
		     plot.mode="numOfCell",my.title="Number Of Cells (Droplet)",y.pretty=NULL,y.lim.min=NULL)
    makePlotDataStat(m.tb=meta.tb[tech.cate=="Droplet",],out.prefix=sprintf("%s.dataStat.Droplet",out.prefix),
		     plot.mode="numOfDonor",y.colum="numOfDonor",my.title="Number Of Donors (Droplet)",
		     exp.tick=0.3,y.pretty=NULL,y.lim.min=NULL)

    #### SmartSeq2-based
    makePlotDataStat(m.tb=meta.tb[tech.cate=="SmartSeq2",],out.prefix=sprintf("%s.dataStat.SmartSeq2",out.prefix),
		     plot.mode="numOfCell",my.title="Number Of Cells (SmartSeq2)",
		     exp.tick=0.20,y.pretty=NULL,y.lim.min=NULL)
    makePlotDataStat(m.tb=meta.tb[tech.cate=="SmartSeq2",],out.prefix=sprintf("%s.dataStat.SmartSeq2",out.prefix),
		     plot.mode="numOfDonor",my.title="Number Of Donors (SmartSeq2)",
		     y.colum="numOfDonor",exp.tick=0.08,y.pretty=NULL,y.lim.min=NULL)

}

#### tissue distribution by OR (Fig. 1F ??)
{

    doPlotORHT <- function(OR.all.mtx,out.prefix,note.str="",
			   method.distance="",k=3,
			   do.hclust=T,pdf.width = 5.5, pdf.height = 10)
    {

        OR.all.mtx.tmp <- OR.all.mtx
        OR.all.mtx.tmp[OR.all.mtx.tmp > 3] <- 3

        OR.hclust.row <- NULL

        if(do.hclust){
            OR.hclust.row <- run.cutree(OR.all.mtx.tmp,k=k,method.distance=method.distance,method.hclust="ward.D2")
            ##OR.hclust.row <- run.cutreeDynamic(OR.all.mtx,deepSplit=1, minClusterSize=2,method.distance="")
            OR.hclust.row$branch <- dendextend::set(OR.hclust.row$branch,"branches_lwd", 2)
        }else{
            th.OR <- 1.5
            OR.pattern <- OR.all.mtx > th.OR
            colnames(OR.pattern) <- sprintf("bin.%s",colnames(OR.pattern))
            ##OR.ext.tb <- cbind(data.table(cluster.name=rownames(OR.all.mtx)),
            OR.ext.tb <- cbind(data.table(cluster.name=rownames(OR.all.mtx.tmp)),
                       OR.pattern,OR.all.mtx)

            OR.tb.list <- list()
            cname.used <- c()
            OR.tb.list[["enrichT"]] <- OR.ext.tb[bin.T==TRUE,][order(-T),]
            cname.used <- c(cname.used,OR.tb.list[["enrichT"]][["cluster.name"]])
            OR.tb.list[["enrichP"]] <- OR.ext.tb[!(cluster.name %in% cname.used) & bin.P==TRUE,][order(-P),]
            cname.used <- c(cname.used,OR.tb.list[["enrichP"]][["cluster.name"]])
            OR.tb.list[["enrichN"]] <- OR.ext.tb[!(cluster.name %in% cname.used) & bin.N==TRUE,][order(-N),]
            cname.used <- c(cname.used,OR.tb.list[["enrichN"]][["cluster.name"]])
            OR.tb.list[["enrichO"]] <- OR.ext.tb[!(cluster.name %in% cname.used),]

            OR.tb.order.list <- llply(names(OR.tb.list),function(x){
            #x <- "enrichT"
            x.hclust <- run.cutree(as.matrix(OR.tb.list[[x]][,c("P","N","T")]),
                           k=1,method.distance=method.distance,method.hclust="ward.D2")
            ##OR.tb.list[[x]][rev(x.hclust$hclust$order),]
            OR.tb.list[[x]][(x.hclust$hclust$order),]
            #setorderv()
                       })
            names(OR.tb.order.list) <- names(OR.tb.list)

            OR.order.mtx <- do.call(rbind,llply(names(OR.tb.order.list),function(x){
                          a.mtx <- as.matrix(OR.tb.order.list[[x]][,c("P","N","T")])
                          rownames(a.mtx) <- OR.tb.order.list[[x]][["cluster.name"]]
                          return(a.mtx)
                       }))
            OR.all.mtx <- OR.all.mtx[rownames(OR.order.mtx),]
        }

        mapping.col <- g.colSet$meta.cluster
        names(mapping.col) <- mcls2Name[names(mapping.col)]
        #col.row.man <- mapping.col[rownames(OR.all.mtx.tmp)[OR.hclust.row$hclust$order]]
        col.row.man <- mapping.col[rownames(OR.all.mtx)]

        sscVis:::plotMatrix.simple(OR.all.mtx,
                       col.ht=circlize::colorRamp2(c(0, 1, 3), viridis::viridis(3)),
                       out.prefix=sprintf("%s.OR.dist.rClust.withDend.a%s",out.prefix,note.str),
                       show.number=F, show.dendrogram=do.hclust,
                       clust.row=if(do.hclust) OR.hclust.row$branch else FALSE,
                       row_dend_width = unit(1.5, "cm"),
                       #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
                       #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                       exp.name=expression(italic(OR)),
                       z.hi=3,
                       #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                       #palatte=viridis::viridis(7),
                       par.heatmap=list(cex.row=1.5,
                                row_names_gp=gpar(col=col.row.man,fontsize=10)),
                       pdf.width = pdf.width, pdf.height = pdf.height)

        return(OR.hclust.row)

    }

    #### tissue Dist (baseline only)
    {

        cluster.name.tb <- fread("../data/metaInfo/name.conversion.txt",head=T)
        mcls2Name <- structure(factor(cluster.name.tb$cluster.name.full,
                      levels=cluster.name.tb$cluster.name.full),
                   names=cluster.name.tb$meta.cluster)

        OR.CD8.list <- do.tissueDist(cellInfo.tb=meta.tb[stype=="CD8" &  treatment=="baseline",],
                        out.prefix=sprintf("%s.STARTRAC.dist.T.baseline.CD8",out.prefix),
                        pdf.width=4,pdf.height=6,verbose=1)
        OR.CD8.mtx <- OR.CD8.list$OR.dist.mtx
        OR.CD4.list <- do.tissueDist(cellInfo.tb=meta.tb[stype=="CD4" &  treatment=="baseline",],
                        out.prefix=sprintf("%s.STARTRAC.dist.T.baseline.CD4",out.prefix),
                        pdf.width=4,pdf.height=6,verbose=1)
        OR.CD4.mtx <- OR.CD4.list$OR.dist.mtx
        OR.all.mtx <- rbind(OR.CD8.mtx,OR.CD4.mtx)
        rownames(OR.all.mtx) <- mcls2Name[rownames(OR.all.mtx)]
        count.dist.melt.ext.tb <- rbind(OR.CD8.list$count.dist.melt.ext.tb,
                        OR.CD4.list$count.dist.melt.ext.tb)
        count.dist.melt.ext.tb <- count.dist.melt.ext.tb[order(cid,adj.p.value,-OR),]
        write.table(count.dist.melt.ext.tb,file=sprintf("%s.OR.count.dist.ext.baseline.txt",
                                out.prefix),row.names=F,sep="\t",quote=F)
        count.dist.melt.ext.tb[OR> 1.2,summary(adj.p.value)]
        count.dist.melt.ext.tb[OR < 0.8,summary(adj.p.value)]
        count.dist.melt.ext.tb[OR> 1.5,summary(adj.p.value)]
        count.dist.melt.ext.tb[OR < 0.5,summary(adj.p.value)]

        write.table(OR.all.mtx,file=sprintf("%s.OR.all.baseline.csv",out.prefix),row.names=T,sep=",",quote=F)
        
        #tmp <- doPlotORHT(OR.all.mtx,out.prefix,note.str=".all.baseline")
        ### use this for paper (Fig. 1F ??)
        OR.hclust.cosine <- doPlotORHT(OR.all.mtx,out.prefix,note.str=".all.baseline.cosine",
                           k=4,method.distance="cosine")
        saveRDS(OR.hclust.cosine,file=sprintf("%s.OR.all.baseline.hclust.cosine.rds",out.prefix))
        #doPlotORHT(OR.all.mtx,out.prefix,note.str=".all.baseline.hclust.F",do.hclust=F,pdf.width=4.8)

    }
    
}

##### P/N/T comparison of panCancer (Fig. 1 D and E)
{

    dat.plot <- freq.all.ht.tb
    ###dat.plot[,loc:=factor(loc,levels=c("P","N","T"),labels=c("Blood","Normal","Tumor"))]
    dat.plot[,loc:=factor(loc,levels=c("P","N","T"))]
    loc.full <- c("P"="Blood","N"="Normal","T"="Tumor")
    dat.plot[,loc.full:=factor(loc.full[as.character(loc)],levels=loc.full)]
    dat.plot$mcls.cate <- "noChange"
    dat.plot[group.var %in% c("CD4.c22.ISG.IFIT1","CD8.c11.Tex.PDCD1","CD8.c13.Tex.myl12a",
						      "CD4.c20.Treg.TNFRSF9","CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13",
						      "CD8.c15.ISG.IFIT1","CD4.c21.Treg.OAS1"),mcls.cate:="up"]
    dat.plot[group.var %in% c("CD4.c13.Temra.CX3CR1","CD8.c07.Temra.CX3CR1",
						      "CD4.c01.Tn.TCF7","CD4.c03.Tn.ADSL"),mcls.cate:="down"]
    mcls.cate.tb <- unique(dat.plot[,c("group.var","mcls.cate"),with=F])[order(group.var),]
    col.axis.text.x <- c("noChange"="black","up"="#E41A1C","down"="#377EB8")
    out.prefix.loc.cmp <- sprintf("%s/loc.cmp/%s",dirname(out.prefix),basename(out.prefix))
    dir.create(dirname(out.prefix.loc.cmp),F,T)
    dat.plot[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(group.var)]]

    #### fig. S7 C and D ??
    p.list.perMcls <- llply(unique(sort(dat.plot$cluster.name)),function(mcls){
	    p <- ggboxplot(dat.plot[cluster.name==mcls,],x="loc",y="freq",
					       color = "loc", legend="none",title=mcls,
					       #fill = "loc",alpha=0.8,
					       #add = "none",outlier.shape=NA) +
					       xlab="",ylab="frequency",
					       add = "jitter",outlier.shape=NA) +
			    scale_color_manual(values=g.colSet$loc) +
			    #scale_color_manual(values=c("P"="#E41A1C","N"="#377EB8","T"="#4DAF4A")) +
			    #scale_fill_manual(values=c("P"="#E41A1C","N"="#377EB8","T"="#4DAF4A")) +
			    stat_compare_means(label="p.format",comparisons=list(c("P","N"),c("N","T"),c("P","T"))) +
			    coord_cartesian(clip="off") +
			    theme(plot.title = element_text(hjust = 0.5,size=14))
	    #ggsave(sprintf("%s.compostion.PNT.%s.png",out.prefix.loc.cmp,mcls),width=2,height=4.5)
	    ggsave(sprintf("%s.compostion.PNT.%s.pdf",out.prefix.loc.cmp,make.names(mcls)),width=2.5,height=4,useDingbats=F)
	    return(p)
						      },.parallel=T)
    names(p.list.perMcls) <- unique(sort(dat.plot$cluster.name))
    saveRDS(p.list.perMcls,file=sprintf("%s.compostion.PNT.rds",out.prefix.loc.cmp))

    dat.plot.avg.tb <- dat.plot[,{
        N <- sum(.SD$N)
        NTotal <- sum(.SD$NTotal)
        .(N=N,
          NTotal=NTotal,
          freq=N/NTotal,
          freq.avg=mean(.SD$freq),
          freq.sd=sd(.SD$freq),
          freq.sem=sd(.SD$freq)/sqrt(length(.SD$freq))
          )
    }, by=c("stype","group.var","loc","mcls.cate")]

    dat.plot.avg.tb[loc=="P" &stype=="CD4",cor(freq,freq.avg)]
    dat.plot.avg.tb[loc=="N" &stype=="CD4",cor(freq,freq.avg)]
    dat.plot.avg.tb[loc=="T" &stype=="CD4",cor(freq,freq.avg)]
    dat.plot.avg.tb[loc=="P" &stype=="CD8",cor(freq,freq.avg)]
    dat.plot.avg.tb[loc=="N" &stype=="CD8",cor(freq,freq.avg)]
    dat.plot.avg.tb[loc=="T" &stype=="CD8",cor(freq,freq.avg)]
    dat.plot.avg.tb[loc=="T" &stype=="CD8",]

    diversity.norm <- function(x)
    {
	    -sum(ifelse(x>0,x*log2(x),0))/log2(length(x))
    }

    dat.diversity.norm.perSample.tb <- dat.plot[,.(NMCls=.N,
						   NTotal=NTotal[1],
						   diversity=diversity.norm(freq)),
	  by=c("stype","cmp.var","donor.var","loc")]
    saveRDS(dat.diversity.norm.perSample.tb,file=sprintf("%s.diversity.norm.perSample.rds",out.prefix))
    write.table(dat.diversity.norm.perSample.tb[,.(diversity.avg=mean(diversity)), by=c("stype","loc")],
		    file=sprintf("%s.diversity.norm.avg.txt",out.prefix),row.names=F,sep="\t",quote=F)

    ##### diversity comparison (fig. S7, A and B ??)
    {

        p <- ggboxplot(dat.diversity.norm.perSample.tb,x="loc",y="diversity",xlab="",
                   add = "jitter",outlier.shape=NA,legend="none",
                           color="loc") +
                scale_color_manual(values=g.colSet[["loc"]]) +
                stat_compare_means(comparisons=list(c("P","N"),c("P","T"),c("N","T")),method="wilcox.test") +
                facet_wrap(~stype,nrow=1) +
                theme(strip.background=element_blank(),strip.text=element_text(size=12))
        ggsave(file=sprintf("%s.diversity.norm.pdf",out.prefix),width=4,height=4)

        p <- ggboxplot(dat.diversity.norm.perSample.tb[stype=="CD8",],x="loc",y="diversity",xlab="",
                   add = "jitter",outlier.shape=NA,legend="none",
                           color="loc") +
                scale_color_manual(values=g.colSet[["loc"]]) +
                stat_compare_means(comparisons=list(c("P","N"),c("P","T"),c("N","T")),method="wilcox.test") +
                facet_wrap(~stype,nrow=1) +
                theme(strip.background=element_blank(),strip.text=element_text(size=12))
        ggsave(file=sprintf("%s.diversity.norm.CD8.pdf",out.prefix),width=2.5,height=4)

        p <- ggboxplot(dat.diversity.norm.perSample.tb[stype=="CD4",],x="loc",y="diversity",xlab="",
                   add = "jitter",outlier.shape=NA,legend="none",
                           color="loc") +
                scale_color_manual(values=g.colSet[["loc"]]) +
                stat_compare_means(comparisons=list(c("P","N"),c("P","T"),c("N","T")),method="wilcox.test") +
                facet_wrap(~stype,nrow=1) +
                theme(strip.background=element_blank(),strip.text=element_text(size=12))
        ggsave(file=sprintf("%s.diversity.norm.CD4.pdf",out.prefix),width=2.5,height=4)

    }

    dat.diversity.norm <- data.table(index="Shannon.norm",
				     CD8.P=dat.plot.avg.tb[loc=="P" & stype=="CD8", diversity.norm(freq) ],
				     CD8.N=dat.plot.avg.tb[loc=="N" & stype=="CD8", diversity.norm(freq) ],
				     CD8.T=dat.plot.avg.tb[loc=="T" & stype=="CD8", diversity.norm(freq) ],
				     CD4.P=dat.plot.avg.tb[loc=="P" & stype=="CD4", diversity.norm(freq) ],
				     CD4.N=dat.plot.avg.tb[loc=="N" & stype=="CD4", diversity.norm(freq) ],
				     CD4.T=dat.plot.avg.tb[loc=="T" & stype=="CD4", diversity.norm(freq) ]
				     )
    write.table(dat.diversity.norm,file=sprintf("%s.diversity.norm.txt",out.prefix),row.names=F,sep="\t",quote=F)

    makePNTPlotFig1 <- function(out.prefix.loc.cmp,mcls.name.use="cluster.name")
    {

        #mcls.name.use <- "cluster.name"
        ### CD8
        dat.plot.tmp <- dat.plot[stype=="CD8",]
        dat.plot.avg.sort.tb <- dat.plot.avg.tb[loc=="T" & stype=="CD8",][order(freq.avg,decreasing=T),]
        dat.plot.avg.sort.tb[,cluster.fullName:=fetchMetaClusterID2CusterFullName(mcls.name.use)[as.character(group.var)]]
        dat.plot.CD8.tb <- dat.plot.avg.tb[stype=="CD8",]
        dat.plot.CD8.tb[,group.var:=factor(as.character(group.var),
                           levels=as.character(dat.plot.avg.sort.tb$group.var))]
        dat.plot.CD8.tb[,cluster.fullName:=factor(fetchMetaClusterID2CusterFullName(mcls.name.use)[as.character(group.var)],
                                                 levels=as.character(dat.plot.avg.sort.tb$cluster.fullName))]
        dat.plot.tmp[,group.var:=factor(as.character(group.var),
                        levels=as.character(dat.plot.avg.sort.tb$group.var))]
        dat.plot.tmp[,cluster.fullName:=factor(fetchMetaClusterID2CusterFullName(mcls.name.use)[as.character(group.var)],
                               levels=as.character(dat.plot.avg.sort.tb$cluster.fullName))]

        p1 <- ggbarplot(dat.plot.CD8.tb,x="cluster.fullName",y="freq.avg",
        #p1 <- ggbarplot(dat.plot.tmp,x="cluster.fullName",y="freq",
                       legend="none",xlab="",ylab="Frequency",
                       #width=0.8,
                       #add = "mean_sd",
                       #add.params = list(width = 1.6),
                       #error.plot = "upper_errorbar",
                       #position = position_dodge2(0.8),
                       color="loc",fill="loc") +
                geom_errorbar(aes(ymin = freq.avg, ymax = freq.avg+freq.sd,color=loc), width = 0.2) +
                scale_fill_manual(values=g.colSet$loc) +
                scale_color_manual(values=g.colSet$loc) +
                facet_wrap(~loc,ncol=1,scales="free_y") +
                theme(axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),
                      strip.text=element_text(size=12),
                      strip.background=element_blank())
        #ggsave(sprintf("%s.CD8.barplot.01.pdf",out.prefix.loc.cmp),width=4,height=6)
        #ggsave(sprintf("%s.CD8.barplot.02.pdf",out.prefix.loc.cmp),width=5,height=6)

        ### CD4
        dat.plot.tmp <- dat.plot[stype=="CD4",]
        dat.plot.avg.sort.tb <- dat.plot.avg.tb[loc=="T" & stype=="CD4",][order(freq.avg,decreasing=T),]
        dat.plot.avg.sort.tb[,cluster.fullName:=fetchMetaClusterID2CusterFullName(mcls.name.use)[as.character(group.var)]]
        dat.plot.CD4.tb <- dat.plot.avg.tb[stype=="CD4",]
        dat.plot.CD4.tb[,group.var:=factor(as.character(group.var),
                           levels=as.character(dat.plot.avg.sort.tb$group.var))]
        dat.plot.CD4.tb[,cluster.fullName:=factor(fetchMetaClusterID2CusterFullName(mcls.name.use)[as.character(group.var)],
                                                 levels=as.character(dat.plot.avg.sort.tb$cluster.fullName))]
        dat.plot.tmp[,group.var:=factor(as.character(group.var),
                        levels=as.character(dat.plot.avg.sort.tb$group.var))]
        dat.plot.tmp[,cluster.fullName:=factor(fetchMetaClusterID2CusterFullName(mcls.name.use)[as.character(group.var)],
                               levels=as.character(dat.plot.avg.sort.tb$cluster.fullName))]

        p2 <- ggbarplot(dat.plot.CD4.tb,x="cluster.fullName",y="freq.avg",
        #p2 <- ggbarplot(dat.plot.tmp,x="cluster.fullName",y="freq",
                       legend="none",xlab="",ylab="Frequency",
                       #width=0.8,
                       #add = "mean_sd",
                       #add.params = list(color = ""),
                       #error.plot = "upper_errorbar",
                       #position = position_dodge2(0.8),
                       color="loc",fill="loc") +
                geom_errorbar(aes(ymin = freq.avg, ymax = freq.avg+freq.sd,color=loc), width = 0.2) +
                scale_fill_manual(values=g.colSet$loc) +
                scale_color_manual(values=g.colSet$loc) +
                facet_wrap(~loc,ncol=1,scales="free_y") +
                theme(axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),
                      strip.text=element_text(size=12),
                      strip.background=element_blank())
        #ggsave(sprintf("%s.CD4.barplot.01.pdf",out.prefix.loc.cmp),width=5.5,height=6)

        pp <- plot_grid(p1,p2,align="hv",nrow=1,rel_widths=c(1,1.3))
        ########ggsave(sprintf("%s.merged.barplot.01.pdf",out.prefix.loc.cmp),width=9,height=6)
        ##ggsave(sprintf("%s.merged.barplot.02.pdf",out.prefix.loc.cmp),width=10,height=6.5)
        ggsave(sprintf("%s.merged.barplot.02.pdf",out.prefix.loc.cmp),width=9,height=5.8)

    }

    makePNTPlotFig1(out.prefix.loc.cmp,mcls.name.use="cluster.name")

}

#### expansion index and proliferation index (Fig. 1G)
{
    prol.tb <- fread(sprintf("%s/panC.proliferation.txt.gz",dir.metaInfo))
    prol.vec <- structure(prol.tb$proliferationScore.bin, names=sprintf("%s.%s",prol.tb$dataset,prol.tb$cellID))

    cellInfo.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")
    cellInfo.tb[,prol.bin:=prol.vec[sprintf("%s.%s",dataset.old,cellID)]]
    cellInfo.tb[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]

    tcr.dat <- readRDS("../data/tcr/byCell/tcr.zhangLab.comb.flt.rds")
    ncell.patient.cluster <- sort(unclass(tcr.dat[,table(sprintf("%s.%s",patient,meta.cluster))]))
    #ncell.patient.cluster <- sort(unclass(in.dat[stype=="CD8",table(sprintf("%s.%s",patient,majorCluster))]))
    tcr.dat <- tcr.dat[ncell.patient.cluster[sprintf("%s.%s",patient,meta.cluster)]>=10,]
    tcr.dat[,cellID.uniq:=sprintf("%s.%s",dataset.old,Cell_Name)]

    tcr.dat.list <- list()
    tcr.dat.list[["all"]] <- tcr.dat
    for(a.loc in c("P","N","T"))
    {
        in.dat.flt <- tcr.dat[loc==a.loc,]
        ncell.patient.cluster <- sort(unclass(in.dat.flt[,table(sprintf("%s.%s",patient,meta.cluster))]))
        in.dat.flt <- in.dat.flt[ncell.patient.cluster[sprintf("%s.%s",patient,meta.cluster)]>=10,]
        tcr.dat.list[[a.loc]] <- in.dat.flt
    }

    startrac.CD4 <- readRDS("../data/tcr/startrac/CD4/CD4.out.startrac.nperm1000.rds")
    startrac.CD8 <- readRDS("../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds")

    startrac.dat.list <- list()
    startrac.dat.list[["all"]] <- as.data.table(rbind(startrac.CD8@cluster.sig.data,
                                                      startrac.CD4@cluster.sig.data))[index=="expa" & aid=="panC",]
    for(a.loc in c("P","N","T"))
    {
        ss.CD8 <- readRDS(sprintf("../data/tcr/startrac/CD8/CD8.out.only%s.startrac.nperm1000.rds",a.loc))
        ss.CD4 <- readRDS(sprintf("../data/tcr/startrac/CD4/CD4.out.only%s.startrac.nperm1000.rds",a.loc))
        startrac.dat.list[[a.loc]] <- as.data.table(rbind(ss.CD8@cluster.sig.data,
                                                          ss.CD4@cluster.sig.data))[index=="expa" & aid=="panC",]
    }


    ##### overall
    cluster.name.use <- "cluster.id"

    dat.ret.tb <- as.data.table(ldply(c("all","P","N","T"),function(a.loc){

        startrac.tb <- startrac.dat.list[[a.loc]]

        my.seed <- 123456
        cex.point <- 1.0
        dat.plot.tb <- cellInfo.tb[cellID.uniq %in% tcr.dat.list[[a.loc]]$cellID.uniq,
                       .(N=.N,prol.n=sum(.SD$prol.bin),freq=sum(.SD$prol.bin)/.N)
                       , by=c("meta.cluster")][N>=30,]
        dat.plot.tb <- as.data.table(ldply(seq_len(nrow(dat.plot.tb)),
                         function(i){
                             N.i <- dat.plot.tb$N[i]
                             prol.n.i <- dat.plot.tb$prol.n[i]
                             NTotal <- sum(dat.plot.tb$N)
                             ProlNTotal <- sum(dat.plot.tb$prol.n)
                             mat.t <- matrix(c(prol.n.i,N.i-prol.n.i,
                                       ProlNTotal-prol.n.i,NTotal-N.i-ProlNTotal+prol.n.i),ncol=2)
                             res.fisher <- fisher.test(mat.t)
                             return(cbind(dat.plot.tb[i,],
                                  data.table(fisher.OR=res.fisher$estimate,
                                         fisher.p.value=res.fisher$p.value)))
                         }))
        dat.plot.tb <- merge(dat.plot.tb,startrac.tb,
                             by.x=c("meta.cluster"),by.y=c("majorCluster"))
        dat.plot.tb[,stype:=sapply(strsplit(meta.cluster,"\\.",perl=T),"[",1)]
        dat.plot.tb[,size:=(log10(N))*cex.point]
        dat.plot.tb <- dat.plot.tb[order(-value),]
        dat.plot.tb[,cluster.name:=fetchMetaClusterID2CusterFullName(cluster.name.use)[as.character(meta.cluster)]]
        dat.plot.tb[,loc:=a.loc]

        l.colSet <- g.colSet
        l.colSet$cluster.name <- l.colSet$meta.cluster
        names(l.colSet$cluster.name) <- fetchMetaClusterID2CusterFullName(cluster.name.use)[names(l.colSet$cluster.name)]
        l.colSet$cluster.name <- l.colSet$cluster.name[!is.na(names(l.colSet$cluster.name))]
        l.colSet$cluster.name <- l.colSet$cluster.name[dat.plot.tb$cluster.name]

        ###for(a.x in c("freq","fisher.OR"))
        for(a.x in c("freq"))
        {

            if(a.x=="freq"){
                dat.plot.hl.tb <- dat.plot.tb[p.value < 0.01 & value>0.02 & freq > 0.05]
            }else{
                dat.plot.hl.tb <- dat.plot.tb[p.value < 0.01 & value>0.02 & fisher.p.value < 0.01 & fisher.OR > 1]
            }
            dat.plot.hl.tb[,cate:="cur.prol"]
            dat.plot.hl.sup.tb <- head(dat.plot.tb[p.value < 0.01 & value>0.10 &
                           !(meta.cluster %in% dat.plot.hl.tb$meta.cluster),],
                       n=100)
            dat.plot.hl.sup.tb[,cate:="mem"]
            dat.plot.hl.tb <- rbind(dat.plot.hl.tb,dat.plot.hl.sup.tb)

            p <- ggplot(dat.plot.tb,aes_string(a.x,"value")) +
                geom_hline(yintercept=0.02,linetype=2,size=1.0,alpha=0.3) +
                geom_vline(xintercept=if(a.x=="freq") 0.05 else 1,
                       linetype=2,size=1.0,alpha=0.3) +
                geom_point(aes(
                           fill=cluster.name,
                           shape=stype,
                           size=size
                           ),
                       ##color="transparent",
                       color="NA",
                       alpha=0.8) +
                scale_fill_manual(values=l.colSet[["cluster.name"]]) +
                scale_shape_manual(values=c("CD4"=22,"CD8"=21),guide="none") +
                #scale_size(breaks=c(3,3.5,4,4.5)*cex.point,range=1.2*c(0.2,6),labels=c(3,3.5,4,4.5)) +
                scale_size(breaks=c(2.5,3,3.5,4)*cex.point,range=1.2*c(0.2,6),labels=c(2.5,3,3.5,4)) +
                xlab(if(a.x=="freq") "Freq. Of Proliferative Cells" else "Odd Ratio Of Proliferative Cells") +
                ylab("Expansion Index")
                if(a.x=="freq"){
                    p <- p + xlim(c(0,0.5)) + ylim(c(0,0.3))
                }
                p <- p + geom_text_repel(data=dat.plot.hl.tb,aes(label=cluster.name,color=cate),parse=T,
                                    segment.alpha=0.3,seed=my.seed, show.legend=F,
                                    size=4,force=10) +
                scale_color_manual(values=c("cur.prol"="#F8766D","mem"="#00BFC4")) +
                #facet_wrap(~stype,ncol=1,scales="free_y") +
                guides(
                       size = guide_legend(title="log10(size)",direction="horizontal",
                               override.aes=list(shape=c(19),color="black",fill="black")),
                       fill=guide_legend(override.aes=list(
                                       shape=c(rep(22,sum(dat.plot.tb$stype=="CD4")), rep(21,sum(dat.plot.tb$stype=="CD8"))),
                                       size=5
                                       ),
                             label.theme=element_text(size=10),ncol=3)
                       ) +
                theme_pubr() +
                theme(legend.position="right")
            ggsave(sprintf("%s.FreqProl.ExpaIndex.overall.%s.thFreq0.5.%s.pdf",out.prefix,a.loc,a.x),
                   width=7.5,height=3.75,useDingbats=F)

            #p <- p + facet_wrap(~stype,ncol=1,scales="free_y") +
            #ggsave(sprintf("%s.FreqProl.ExpaIndex.overall.%s.thFreq0.5.%s.splitStype.pdf",out.prefix,a.loc,a.x),
            #       width=7.5,height=6,useDingbats=F)


        }

        return(dat.plot.tb)
    }))
    write.table(dat.ret.tb,file=sprintf("%s.FreqProl.ExpaIndex.overall.txt",out.prefix),
            sep="\t",quote=F,row.names=F)

}

