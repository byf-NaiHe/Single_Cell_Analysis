#!/usr/bin/env Rscript

library("sscVis")
library("ggplot2")
library("plyr")
library("cowplot")
library("ggpubr")
source("./func.R")

meta.tb.file <- "../data/metaInfo/panC.freq.all.meta.tb.rds"
out.prefix <- "./OUT_Fig4/freq/panC.freq.all"
colSet.file <- "../data/metaInfo/panC.colSet.list.rds"
dir.create(dirname(out.prefix),F,T)

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

meta.tb <- readRDS(meta.tb.file)
colSet <- readRDS(colSet.file)
colSet$cluster.name.full <- colSet$meta.cluster
names(colSet$cluster.name.full) <- fetchMetaClusterID2CusterFullName()[names(colSet$cluster.name.full)]
f.na <- is.na(names(colSet$cluster.name.full))
colSet$cluster.name.full <- colSet$cluster.name.full[!f.na]

###############################


#### meta.cluster comparison among cancer types
makePlot.cmp.cancerType <- function(meta.tb,cmp.loc="T",method.use=NULL)
{
    p.list.meta.cluster <- list()
    p.list.meta.cluster.coarse <- list()
    loc.full <- c("P"="Blood","N"="Normal","T"="Tumor")
    for(cmp.stype in c("CD8","CD4"))
    {
        out.prefix.plot <- sprintf("%s/all.%s.%s%s/%s",dirname(out.prefix),cmp.loc,cmp.stype,
                       if(is.null(method.use)) "" else sprintf(".%s",method.use),
                       basename(out.prefix))
        dir.create(dirname(out.prefix.plot),F,T)
        mcls2fullname <- fetchMetaClusterID2CusterFullName()
        ii2mcls <- structure(names(mcls2fullname), names=mcls2fullname)

        dat.plot <- meta.tb[usedForFreq=="Y" & treatment=="baseline" & stype==cmp.stype & loc==cmp.loc,]
        dat.plot[,cluster.name:=mcls2fullname[as.character(meta.cluster)]]
        #dat.plot.debug <<- dat.plot
        ##p.names <- unique(dat.plot[["meta.cluster"]])
        p.names <- unique(dat.plot[["cluster.name"]])
        p.list <- llply(p.names,function(ii){
            p <- plotDistFromCellInfoTable(dat.plot, out.prefix=NULL,
                               plot.type="boxplot2",legend="none",
                               facet.ncol=3,
                               fill="cmp.var",alpha=0.5,
                               group.filter=ii,
                               sort.freq=T,
                               par.stat=list(
                                     method=method.use,
                                     label.x.npc=0.05,label.y.npc=0.9),min.NTotal=30,
                               ###cmp.var="cancerType",group.var="meta.cluster",donor.var="patient") +
                               cmp.var="cancerType",group.var="cluster.name",donor.var="patient") +
                geom_hline(yintercept=0.10,linetype=2) +
                labs(x="",y=sprintf("Frequency (%s)",loc.full[cmp.loc])) +
                theme(strip.background=element_blank(),
                      strip.text=element_text(size=14))+
                scale_color_manual(values=colSet$cancerType) +
                scale_fill_manual(values=colSet$cancerType)
            ggsave(sprintf("%s.cancerType.meta.cluster.dist.all.%s.%s.pdf",out.prefix.plot,cmp.loc,ii2mcls[ii]),
                   width=6,height=3,useDingbats=F)
            return(p)
                                   },.parallel=T)
        names(p.list) <- ii2mcls[p.names]
        p.list.meta.cluster <- c(p.list.meta.cluster,p.list)

#####        p.names <- unique(dat.plot[["meta.cluster.coarse"]])
#####        p.list <- llply(p.names,function(ii){
#####            p <- plotDistFromCellInfoTable(dat.plot, out.prefix=NULL,
#####                               plot.type="boxplot2",legend="none",
#####                               facet.ncol=3,
#####                               group.filter=ii,sort.freq=T,
#####                               par.stat=list(
#####                                     method = method.use,
#####                                     label.x.npc=0.05,label.y.npc=0.9),min.NTotal=30,
#####                               cmp.var="cancerType",group.var="meta.cluster.coarse",
#####                               donor.var="patient") +
#####                    geom_hline(yintercept=0.10,linetype=2) +
#####                    labs(x="",y=sprintf("Frequency (%s)",loc.full[cmp.loc])) +
#####                    theme(strip.background=element_blank())+
#####                    scale_color_manual(values=colSet$cancerType)
#####            ggsave(sprintf("%s.cancerType.meta.cluster.coarse.dist.all.%s.%s.pdf",
#####                       out.prefix.plot,cmp.loc,ii), width=6,height=3,useDingbats=F)
#####                                   },.parallel=T)
#####        names(p.list) <- p.names
#####        p.list.meta.cluster.coarse <- c(p.list.meta.cluster.coarse,p.list)

    }
    return(list("p.list.meta.cluster"=p.list.meta.cluster,
                "p.list.meta.cluster.coarse"=p.list.meta.cluster.coarse))
}

p.cmp.loc.T <- makePlot.cmp.cancerType(meta.tb,cmp.loc="T",method.use="anova")
p.cmp.loc.N <- makePlot.cmp.cancerType(meta.tb,cmp.loc="N",method.use="anova")
p.cmp.loc.P <- makePlot.cmp.cancerType(meta.tb,cmp.loc="P",method.use="anova")
save(p.cmp.loc.P,p.cmp.loc.N,p.cmp.loc.T,file=sprintf("%s.cancerType.cmp.plot.list.anova.RData",out.prefix))

lname <- load(file=sprintf("%s.cancerType.cmp.plot.list.anova.RData",out.prefix))

### compare freq among cancer types (Fig. 4A, fig. S29A, fig. S30)
{
    ### Fig.SXX (fig. S30 ??)
    pp <- plot_grid(plotlist=p.cmp.loc.T$p.list.meta.cluster[c("CD4.c08.Tm.CREM", "CD4.c15.Th17.IL23R", 
                                   "CD8.c09.Tk.KIR2DL4", "CD4.c10.Tm.CAPG",
                                   #"CD8.c03.Tm.RPS12",
                                   "CD8.c07.Temra.CX3CR1", "CD4.c01.Tn.TCF7",
                                   "CD8.c02.Tm.IL7R", "CD8.c11.Tex.PDCD1",
                                   "CD8.c06.Tem.GZMK", "CD4.c11.Tm.GZMA",
                                   "CD8.c15.ISG.IFIT1","CD8.c04.Tm.CD52")],
                    align="hv",ncol=3)
    ggsave(sprintf("%s.cancerType.cmp.meta.cluster.merge.Fig.SXX.T.var.anova.v3.pdf",out.prefix),
           width=15,height=2.625*4,useDingbats=F)

    ### Fig. 4A ???
    pp <- plot_grid(plotlist=p.cmp.loc.T$p.list.meta.cluster[c("CD8.c12.Tex.CXCL13", "CD8.c10.Trm.ZNF683",
                                   "CD4.c20.Treg.TNFRSF9", "CD8.c14.Tex.TCF7")],
                    align="hv",ncol=2)
    ggsave(sprintf("%s.cancerType.cmp.meta.cluster.merge.Fig.SXX.T.var.main.anova.v3.pdf",out.prefix),
           width=5*2,height=2.625*2,useDingbats=F)

    ### fig. S29A ???
    pp <- plot_grid(plotlist=p.cmp.loc.N$p.list.meta.cluster[c("CD8.c16.MAIT.SLC4A10","CD4.c16.Tfh.CXCR5", 
                                   "CD8.c09.Tk.KIR2DL4", "CD8.c07.Temra.CX3CR1")],
                    align="hv",ncol=2)
    ggsave(sprintf("%s.cancerType.cmp.meta.cluster.merge.Fig.SXX.N.var.anova.v2.pdf",out.prefix),
           width=10,height=5.15,useDingbats=F)

}

##################### 
getFreqTable <- function(tb,stype,type.return="mtx",group.var="meta.cluster")
{
	out.tb <- as.data.table(ldply(c("P","N","T"),function(x){
					x.tb <- plotDistFromCellInfoTable(tb[loc==x,], plot.type="none",
									  cmp.var="cancerType",min.NTotal=30,
									  group.var=group.var,donor.var="patient.uid")
					x.tb[,loc:=x]
					x.tb[,stype:=stype]
				  }))
	if(type.return=="tb"){
		return(out.tb)
	}else if(type.return=="mtx"){
		d.tb <- out.tb
		ht.tb <- dcast(d.tb,group.var~loc+donor.var,value.var="freq",fill=0)
		ht.mtx <- as.matrix(ht.tb[,-1])
		rownames(ht.mtx) <- ht.tb[[1]]
		print(ht.mtx[,1:3])
		return(ht.mtx)
	}
}

if(!exists("../data/metaInfo/panC.freq.all.ht.tb.rds"))
{
    freq.CD8.ht.tb <- getFreqTable(meta.tb[usedForFreq=="Y" & treatment=="baseline" & stype=="CD8",],
                                   "CD8",type.return="tb")
    freq.CD4.ht.tb <- getFreqTable(meta.tb[usedForFreq=="Y" & treatment=="baseline" & stype=="CD4",],
                                   "CD4",type.return="tb")


    #### 
    freq.all.ht.tb <- rbind(freq.CD8.ht.tb,freq.CD4.ht.tb)
    saveRDS(freq.all.ht.tb,file=sprintf("%s.freq.all.ht.tb.rds",out.prefix))
    saveRDS(freq.CD8.ht.tb,file=sprintf("%s.freq.CD8.ht.tb.rds",out.prefix))
    saveRDS(freq.CD4.ht.tb,file=sprintf("%s.freq.CD4.ht.tb.rds",out.prefix))
    freq.all.ht.tb <- readRDS(sprintf("%s.freq.all.ht.tb.rds",out.prefix))
    freq.CD8.ht.tb <- readRDS(sprintf("%s.freq.CD8.ht.tb.rds",out.prefix))
    freq.CD4.ht.tb <- readRDS(sprintf("%s.freq.CD4.ht.tb.rds",out.prefix))
}else{
    ##### previously saved
    freq.all.ht.tb <- readRDS("../data/metaInfo/panC.freq.all.ht.tb.rds")
    freq.CD8.ht.tb <- readRDS("../data/metaInfo/panC.freq.CD8.ht.tb.rds")
    freq.CD4.ht.tb <- readRDS("../data/metaInfo/panC.freq.CD4.ht.tb.rds")
}

#### NImpactT (fig. S29B)
{
    make.cor.NImpactT.ind.plot <- function(freq.tb,mcls.oi,out.prefix)
    {
        freq.tb[,meta.cluster:=group.var]
        freq.tb[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
        ret.list <- llply(mcls.oi,function(mcls){
            ##mcls <- "CD8.c16.MAIT.SLC4A10"
            mcls.freq.N <- freq.tb[meta.cluster==mcls & loc %in% c("N"),
                       c("donor.var","cmp.var","meta.cluster","cluster.name","N","NTotal","freq","loc"),with=F]
            mcls.freq.T <- freq.tb[meta.cluster==mcls & loc %in% c("T"),
                       c("donor.var","cmp.var","meta.cluster","cluster.name","N","NTotal","freq","loc"),with=F]
            mcls.freq.cmp <- merge(mcls.freq.N,mcls.freq.T,by=c("donor.var","cmp.var","meta.cluster","cluster.name"))
            mcls.freq.cmp[,size:=(log10(N.x+1)+log10(N.y+1))/2]
            #mcls.freq.cmp[,isLIHC:=cmp.var %in% c("HCC","CHOL","STAD","UCEC","CRC")]
            p <- ggscatter(mcls.freq.cmp,x="freq.x",y="freq.y",size="size",color="cmp.var") +
                labs(x="Frequency In Normal",y="Frequency In Tumor",title=mcls.freq.cmp$cluster.name[1]) +
                geom_smooth(method='lm') +
                guides(size = guide_legend(order=2,direction="horizontal",
                               title.position = "top", label.position = "bottom"),
                   color=guide_legend(order=1,override.aes=list(size=5),ncol=2),
                   legend.direction="horizontal") +
                  geom_smooth(method='lm') +
                  coord_trans(clip="off") +
                  scale_size(breaks=log10(c(10,30,50,100)),range=c(1,10),labels=c(10,30,50,100),limits=c(0,10)) +
                  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
                  ##stat_cor(aes(label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")),) +
                  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
                  ###stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),vjust=2.5) +
                  stat_regline_equation(aes(label =  paste(..eq.label..)),vjust=2.5) +
                  ##stat_cor() +
                  #facet_wrap(~isLIHC,ncol=2,scales="free_x") +
                  scale_color_manual(values=colSet[["cancerType"]])
            pp <- plot_grid(p+theme(legend.position = "none"),get_legend(p),ncol=1,rel_heights=c(0.5,0.45),align="v")
            ggsave(sprintf("%s.tissueImpact.%s.pdf",out.prefix,mcls),width=4.5,height=7,useDingbats=F)
            ##ggsave(sprintf("%s.tissueImpact.%s.pdf",out.prefix,mcls),width=7.5,height=4,useDingbats=F)
            #mcls.freq.cmp$isLIHC %>% table
            return(p)
        },.parallel=T)
        names(ret.list) <- mcls.oi
        return(ret.list)
    }

    freq.all.ht.tb <- readRDS("../data/metaInfo/panC.freq.all.ht.tb.rds")
    freq.CD8.ht.tb <- readRDS("../data/metaInfo/panC.freq.CD8.ht.tb.rds")
    freq.CD4.ht.tb <- readRDS("../data/metaInfo/panC.freq.CD4.ht.tb.rds")
    p.list.cor.NImpactT.ind.CD8 <- make.cor.NImpactT.ind.plot(freq.tb=freq.CD8.ht.tb,
                                  mcls.oi=c("CD8.c09.Tk.KIR2DL4","CD8.c07.Temra.CX3CR1",
                                        "CD8.c16.MAIT.SLC4A10","CD8.c14.Tex.TCF7"),
                                  out.prefix=out.prefix)
    p.list.cor.NImpactT.ind.CD4 <- make.cor.NImpactT.ind.plot(freq.tb=freq.CD4.ht.tb,
                                  mcls.oi=c("CD4.c16.Tfh.CXCR5","CD4.c08.Tm.CREM",
                                        "CD4.c19.Treg.S1PR1", "CD4.c21.Treg.OAS1",
                                        "CD4.c18.Treg.RTKN2"),
                                  out.prefix=out.prefix)

    p.list.cor.NImpactT.merge <- c(p.list.cor.NImpactT.ind.CD8[c("CD8.c16.MAIT.SLC4A10")],
                       p.list.cor.NImpactT.ind.CD4[c("CD4.c16.Tfh.CXCR5")],
                       p.list.cor.NImpactT.ind.CD8[c("CD8.c09.Tk.KIR2DL4","CD8.c07.Temra.CX3CR1")])
    pp <- plot_grid(plot_grid(plotlist=llply(p.list.cor.NImpactT.merge,
                                             function(x){ x + theme(legend.position="none") }),
                    align="hv",nrow=2),
                    get_legend(p.list.cor.NImpactT.ind.CD8[["CD8.c16.MAIT.SLC4A10"]]),
                    align="v",ncol=1,rel_heights=c(0.6,0.4))
    ggsave(file=sprintf("%s.tissueImpact.cor.merged.v2.pdf",out.prefix),
           width=7,height=10,useDingbats=F)

}

#### clustering samples by T cell compositions (fig. S33A, Fig. 5A, fig. S37, fig. S33B)
#### prepare data for R^2 barplot (cancerType, immuneType)
{

	mapping.tb <- unique(freq.all.ht.tb[,c("donor.var","cmp.var"),with=F])
	mapping.vec <- structure(mapping.tb$cmp.var,names=mapping.tb$donor.var)

	freq.CD8.ht.tb[,meta.cluster:=group.var]
	freq.CD8.ht.tb[,group.var:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
	freq.CD4.ht.tb[,meta.cluster:=group.var]
	freq.CD4.ht.tb[,group.var:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]

	##### immune type .vs. cancer type
	do.ass.immuneType.cancerType <- function(meta.tb,cluster.sample.T.tb,out.prefix)
	{

	    dir.create(dirname(out.prefix),F,T)
	    meta.tb[1,]
	    donor.tb <-meta.tb[loc=="T",.(NCells=.N),by=c("cancerType","dataset","tech","patient.uid")]
	    donor.tb[,patient.uid:=gsub("_","\\.",patient.uid,perl=T)]
	    donor.immuneType.tb <- merge(donor.tb,cluster.sample.T.tb,by.x="patient.uid",by.y="donorID")
	    donor.immuneType.tb[,immuneType:=sprintf("C%02d",displayCluster)]
	    donor.immuneType.tb[,table(immuneType)]
	    donor.immuneType.tb[,table(cancerType)]
	    dist.immuneType.cancerType <- unclass(donor.immuneType.tb[,table(cancerType,immuneType)])

	    dist.immuneType.cancerType.melt <- test.dist.table(count.dist=dist.immuneType.cancerType,min.rowSum=5)
	    
	    dist.immuneType.cancerType.melt[adj.p.value < 0.1,]
	    dist.immuneType.cancerType.melt[,.(N=sum(count)),by=c("rid")][order(-N),]
	    dist.immuneType.cancerType.melt.sig <- dist.immuneType.cancerType.melt[p.value < 0.05 & OR > 0,
										   ][order(cid,p.value),]
	    write.table(dist.immuneType.cancerType.melt,row.names=F,
			file=sprintf("%s.dist.immuneType.cancerType.melt.txt",out.prefix),
			sep="\t",quote=F)
	    write.table(dist.immuneType.cancerType.melt.sig,row.names=F,
			file=sprintf("%s.dist.immuneType.cancerType.melt.sig.txt",out.prefix),
			sep="\t",quote=F)

	    library("webr")
	   
	    #out.prefix.plot <- sprintf("%s/immuneType.cancerType/%s",dirname(out.prefix),basename(out.prefix))
	    #dir.create(dirname(out.prefix.plot),F,T)

	    bcolor.tb <- unique(donor.immuneType.tb[,c("bcolor","immuneType"),with=F])
	    bcolor.vec <- structure(bcolor.tb$bcolor,names=bcolor.tb$immuneType)

	    cancerType.plot <- unique(dist.immuneType.cancerType.melt.sig[order(cid,p.value),][["rid"]])
	    p.list <- llply(cancerType.plot,function(cancerType.i){

		dat.plot.tmp <- dist.immuneType.cancerType.melt[rid==cancerType.i,][count > 0,]
		dat.plot.tmp.sig <- copy(dat.plot.tmp)
		dat.plot.tmp.sig <- dat.plot.tmp.sig[OR > 0 & p.value < 0.05,]
		dat.plot.tmp[cid %in% dat.plot.tmp.sig$cid, cid:=sprintf("%s*",cid)]
		dat.plot.tmp[,cid:=factor(cid,levels=sort(unique(cid)))]

		bcolor.vec.tmp <- bcolor.vec
		f.sig <- names(bcolor.vec.tmp) %in% dat.plot.tmp.sig$cid
		names(bcolor.vec.tmp)[f.sig] <- sprintf("%s*",names(bcolor.vec.tmp)[f.sig])

		p <- PieDonut(dat.plot.tmp,aes(cid,count=count),r0=0,
			      explode=which(levels(dat.plot.tmp$cid)==dat.plot.tmp[ OR >0 & p.value < 0.05,][["cid"]]),
			      showPieName=F,maxx=1.2,r1=1,r2=1,pieLabelSize=4.5) +
			annotate("text",x=0,y=1.2,label=cancerType.i,size=5) +
			scale_fill_manual(values=bcolor.vec.tmp) +
			theme(panel.border = element_blank())
#		ggsave(sprintf("%s.immuneType.cancerType.pieDonut.%s.pdf",out.prefix.plot,cancerType.i),
#		       width=5,height=5)
		return(p)

			})
	    names(p.list) <- cancerType.plot
	    pp <- plot_grid(plotlist=p.list,align="hv",ncol=3)
	    ggsave(sprintf("%s.immuneType.cancerType.pieDonut.merge.pdf",out.prefix),
		   width=12,
		   height=12)

	}

	#### barplot (immune types (grouping by frequency of meta-clusters))
    do.gen.RSq.by.immuneType <- function(freq.all.ht.tb,cluster.sample.T.mapping,out.prefix)
	{

	    mcls.freq.T <- freq.all.ht.tb[loc %in% c("T"),]
	    mcls.freq.T[,group.tumor:=cluster.sample.T.mapping[donor.var]]
	    mcls.freq.T <- mcls.freq.T[!is.na(group.tumor)]
	    mcls.freq.T$group.tumor %>% table(useNA="always")

	    mcls.freq.groupTumor.lm.df <- as.data.table(ldply(unique(mcls.freq.T$group.var),function(x){
								  res.lm <- mcls.freq.T[group.var==x,
											][,lm(freq~group.tumor)]
								  res.lm.s <- summary(res.lm)
								  nRich <- mcls.freq.T[group.var==x,
										       ][freq>0.1,nrow(.SD)]
								  data.table(meta.cluster=x,
									     cate="groupTumor",
									     r.squared=res.lm.s$r.squared,
									     adj.r.squared=res.lm.s$adj.r.squared,
									     p.value=anova(res.lm)$'Pr(>F)'[1],
									     nRich=nRich)
								  ##p.value=res.lm.s$coefficients["freq.x","Pr(>|t|)"])
					       }))
	    ###
	    mcls.freq.groupTumor.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
	    mcls.freq.groupTumor.lm.df[,adj.p.value:=p.adjust(p.value,"BH")]
	    mcls.freq.groupTumor.lm.df[adj.p.value < 0.2,]
	    ##write.table(mcls.freq.groupTumor.lm.df,file=sprintf("%s.data.lm.groupTumor.txt",out.prefix),
	    write.table(mcls.freq.groupTumor.lm.df,file=sprintf("%s.txt",out.prefix),
			row.names=F,sep="\t",quote=F)
	}

	#### immuen types, TH.cor: 0.35 (fig. S33A, Fig. 5A, fig. S37)
    #### use euclidean distance for samples
	{

	    cluster.sample.T <- do.plot.freq.heatmap(dat.plot.a=freq.CD8.ht.tb[loc=="T",],
						     mapping.vec=mapping.vec,
						     colSet=colSet,
						     group.var="cluster.name.full",
						     dat.plot.b=freq.CD4.ht.tb[loc=="T",],
						     #z.lo=0,z.hi=0.3,,z.len=10,
						     #z.lo=-1.5,z.hi=1.5,z.len=10,
						     out.prefix=sprintf("%s.CD4CD8.ht.T",out.prefix),
						     k.corMat=4,
						     k.column=8,k.row=2,pdf.width=16,pdf.height=10)
	    saveRDS(cluster.sample.T,file=sprintf("%s.cluster.sample.T.rds",out.prefix))
	    cluster.sample.T <- readRDS(file=sprintf("%s.cluster.sample.T.rds",out.prefix))

	    cluster.sample.T.tb <- data.table(donorID=names(cluster.sample.T$flt$obj.hclust.col$cluster),
					      donorCluster=cluster.sample.T$flt$obj.hclust.col$cluster,
					      order.dend=order.dendrogram(cluster.sample.T$flt$obj.hclust.col$branch))
	    cluster.sample.T.tb <- cluster.sample.T.tb[cluster.sample.T.tb$order.dend,]
	    cluster.sample.T.tb[,displayCluster:=as.integer(factor(donorCluster,levels=unique(donorCluster)))]
	    cluster.sample.T.tb$bcolor <- dendextend::get_leaves_branches_col(cluster.sample.T$flt$obj.hclust.col$branch)
	    cluster.sample.T.tb[displayCluster==7,]
	    write.table(cluster.sample.T.tb,file=sprintf("%s.sample.T.donorCluster.txt",out.prefix),
				    row.names=F,quote=F,sep="\t")

        #### pie plots ( fig. S37)
	    do.ass.immuneType.cancerType(meta.tb,cluster.sample.T.tb,
					 out.prefix=sprintf("%s/immuneType.cancerType.Th.cor.0.35/%s",
							    dirname(out.prefix),basename(out.prefix)))


	}
	

	### pair wise correlation (fig. S33B)
	{
	    #freq.TcTh.ht.tb[,cluster.name.full:=fetchMetaClusterID2CusterFullName()[as.character(group.var)]]
	    freq.all.ht.tb[,cluster.name.full:=fetchMetaClusterID2CusterFullName()[as.character(group.var)]]
	    makeFig.freqCor.scatter <- function(dat.tb,mcls,out.prefix,cex.point=0.5)
	    {
            require("ggpubr")
            require("ggplot2")

            dat.tb.a <- dat.tb[group.var==mcls[1],-c("stype")]
            dat.tb.b <- dat.tb[group.var==mcls[2],-c("stype")]
            #dat.plot <- merge(dat.tb.a,dat.tb.b,by=c("donor.var","cmp.var","loc"))
            dat.plot <- merge(dat.tb.a,dat.tb.b,by=c("donor.var","cmp.var"))
            dat.plot[,size:=(log10(N.x+1)+log10(N.y+1))*cex.point/2]
            dat.plot.debug <<- dat.plot
            print(summary(dat.plot$size))
            print(summary(dat.plot$freq.x))
            print(summary(dat.plot$freq.y))
            xvar.name <- dat.tb.a$cluster.name.full[1]
            yvar.name <- dat.tb.b$cluster.name.full[1]
            p <- ggscatter(dat.plot,x="freq.x",y="freq.y",size="size",color="cmp.var",alpha=0.75) +
                    labs(x=sprintf("%s",xvar.name),
                         y=sprintf("%s",yvar.name)) +
                    geom_smooth(method='lm') +
                    guides(size = guide_legend(order=2,direction="horizontal",
                                   title.position = "top", label.position = "bottom"),
                           color=guide_legend(order=1,override.aes=list(size=5),ncol=2),
                           legend.direction="horizontal") +
                    geom_smooth(method='lm') +
                    ###coord_cartesian(clip="off",ylim=c(0,rev(pretty(dat.plot$freq.y))[1])) + 
                    coord_cartesian(clip="off") + 
                    scale_size(breaks=log10(c(10,30,50,100))*cex.point,range=c(0.2,6),labels=c(10,30,50,100),
                           limits=c(0,10)*cex.point) +
                    theme(legend.position = "right", 
                          plot.title = element_text(hjust = 0.5)) +
                    ##stat_cor(aes(label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")),) +
                    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
                    stat_regline_equation(aes(label =  paste(..eq.label..)),vjust=2.5) +
                    #stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label..,
                    #					 sep = "~~~~")),vjust=2.5) +
                    ##stat_cor() +
                    #facet_wrap(~isLIHC,ncol=2,scales="free_x") +
                    scale_color_manual(values=colSet[["cancerType"]])
            pp <- plot_grid(p+theme(legend.position = "none"),get_legend(p),ncol=1,rel_heights=c(0.40,0.6),align="v")
            ggsave(sprintf("%s.pdf",out.prefix),width=3.5,height=7,useDingbats=F)
            return(p)
	    }
	   

	    ########## batch
	    cmp.freqCor.list <- list("CD8.Trm.Tex"=c("CD8.c10.Trm.ZNF683","CD8.c12.Tex.CXCL13"),
				     "CD8.Tem.Tex"=c("CD8.c06.Tem.GZMK","CD8.c12.Tex.CXCL13"),
				     "CD8.Two.Tex"=c("CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13"),
				     #"CD4.Two.Tfh"=c("CD4.c16.Tfh.CXCR5","CD4.c17.TfhTh1.CXCL13"),
				     "CD4.Tfh.1"=c("CD4.c02.Tn.PASK","CD4.c16.Tfh.CXCR5"),
				     "CD4.Tfh.2"=c("CD4.c02.Tn.PASK","CD4.c17.TfhTh1.CXCL13"),
				     "CD4.Tfh.3"=c("CD4.c16.Tfh.CXCR5","CD4.c17.TfhTh1.CXCL13"),
				     "all.Two.Th1"=c("CD4.c17.TfhTh1.CXCL13","CD8.c12.Tex.CXCL13"),
				     "all.Tex.Treg"=c("CD4.c20.Treg.TNFRSF9","CD8.c12.Tex.CXCL13"),
				     "all.MAIT.Th17"=c("CD8.c16.MAIT.SLC4A10","CD4.c15.Th17.IL23R"),
				     "all.Two.Temra"=c("CD8.c07.Temra.CX3CR1","CD4.c13.Temra.CX3CR1"),
				     "all.Two.Tn"=c("CD8.c01.Tn.MAL","CD4.c01.Tn.TCF7"))
	    p.list.freqCor.scatter <- llply(names(cmp.freqCor.list), function(x){
						makeFig.freqCor.scatter(dat.tb=freq.all.ht.tb[loc=="T" &
									donor.var %in% unique(cluster.sample.T.tb$donorID) &
									group.var %in% cmp.freqCor.list[[x]],],
								    mcls=cmp.freqCor.list[[x]],
								    out.prefix=sprintf("%s.fig.freqCor.Scatter.%s",
										       out.prefix,x))
							     },.parallel=T)
	    names(p.list.freqCor.scatter) <- names(cmp.freqCor.list)
	    saveRDS(p.list.freqCor.scatter,file=sprintf("%s.fig.freqCor.Scatter.merge.rds",out.prefix))

	    pp <- plot_grid(plot_grid(plotlist=llply(p.list.freqCor.scatter[c("CD8.Trm.Tex","CD8.Tem.Tex","CD8.Two.Tex","all.Two.Th1",
                                                                          "all.Tex.Treg","all.MAIT.Th17","all.Two.Temra","all.Two.Tn")],
                                                 function(x){ x + theme(legend.position="none") }),
							      align="hv",ncol=4),
					    get_legend(p.list.freqCor.scatter[["CD8.Trm.Tex"]]),
					    align="hv",ncol=1,rel_heights=c(0.8,0.6))
	    ggsave(file=sprintf("%s.fig.freqCor.Scatter.merge.pdf",out.prefix),width=13,height=9,useDingbats=F)

	}

	###### prepare data for R^2 barplot
	{

	    {
		cluster.sample.T.tb <- fread(sprintf("%s.sample.T.donorCluster.txt",out.prefix))
		cluster.sample.T.mapping <- structure(sprintf("G%02d",cluster.sample.T.tb$displayCluster),
										      names=cluster.sample.T.tb$donorID)

		#### barplot (immune types (grouping by frequency of meta-clusters))
		do.gen.RSq.by.immuneType(freq.all.ht.tb,cluster.sample.T.mapping,
					 out.prefix=sprintf("%s.data.lm.groupTumor",out.prefix))
	    }
	    
	    #### barplot (cancer type)
	    {
		mcls.freq.T <- freq.all.ht.tb[loc %in% c("T"),]
		mcls.freq.T$cmp.var %>% table(useNA="always")

		mcls.freq.cancerType.lm.df <- as.data.table(ldply(unique(mcls.freq.T$group.var),function(x){
								      res.lm <- mcls.freq.T[group.var==x,
											    ][,lm(freq~cmp.var)]
								      res.lm.s <- summary(res.lm)
								      nRich <- mcls.freq.T[group.var==x,
											   ][freq>0.1,nrow(.SD)]
								      data.table(meta.cluster=x,
										 cate="cancerType",
										 r.squared=res.lm.s$r.squared,
										 adj.r.squared=res.lm.s$adj.r.squared,
										 p.value=anova(res.lm)$'Pr(>F)'[1],
										 nRich=nRich)
								      ##p.value=res.lm.s$coefficients["freq.x","Pr(>|t|)"])
						   }))
		###
		mcls.freq.cancerType.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
		mcls.freq.cancerType.lm.df[,adj.p.value:=p.adjust(p.value,"BH")]
		mcls.freq.cancerType.lm.df[adj.p.value < 0.2,]

		write.table(mcls.freq.cancerType.lm.df,file=sprintf("%s.data.lm.cancerType.txt",out.prefix),
			    row.names=F,sep="\t",quote=F)
	    }

	}

}

############# R^2 barplots (fig. S28, fig. S36) ################
{
    #### prepare data for R^2 (tissue impact)
    {
        mcls.freq.N <- freq.all.ht.tb[loc %in% c("N"),]
        mcls.freq.T <- freq.all.ht.tb[loc %in% c("T"),]
        mcls.freq.cmp <- merge(mcls.freq.N,mcls.freq.T,by=c("donor.var","cmp.var","group.var","cluster.name.full","stype"))

        mcls.freq.tissueImpact.lm.df <- as.data.table(ldply(unique(mcls.freq.cmp$group.var),function(x){
                                    res.lm <- mcls.freq.cmp[group.var==x,][,lm(freq.y~freq.x)]
                                    res.lm.s <- summary(res.lm)
                                    nRich <- mcls.freq.cmp[group.var==x,][freq.x > 0.1,nrow(.SD)]
                                    data.table(meta.cluster=x,
                                           cate="tissueImpact",
                                           r.squared=res.lm.s$r.squared,
                                           adj.r.squared=res.lm.s$adj.r.squared,
                                           p.value=anova(res.lm)$'Pr(>F)'[1],
                                           nRich=nRich)
                                    ##p.value=res.lm.s$coefficients["freq.x","Pr(>|t|)"])
                           }))
        ###
        mcls.freq.tissueImpact.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name.full")[as.character(meta.cluster)]]
        mcls.freq.tissueImpact.lm.df[,adj.p.value:=p.adjust(p.value,"BH")]
        mcls.freq.tissueImpact.lm.df[adj.p.value < 0.2,]

        write.table(mcls.freq.tissueImpact.lm.df,file=sprintf("%s.data.lm.tissueImpact.txt",out.prefix),
                    row.names=F,sep="\t",quote=F)
    }

    #### prepare data for R^2 (age,gender,stage,BMI)
    {

        #clinical.info.tb <- fread("list/TableS1.clinicalInfo.20200609.txt")
        clinical.info.tb <- fread("../data/metaInfo/TableS1.clinicalInfo.P316.txt")
        colnames(clinical.info.tb)[1] <- "donor.var"
        colnames(clinical.info.tb)[5] <- "cmp.var"
        colnames(clinical.info.tb)[9] <- "stage"
        
        mcls.freq.clin.tb <- merge(freq.all.ht.tb,
                       clinical.info.tb[,c("donor.var","cmp.var","Age","Gender","stage","BMI"),with=T],
                       by=c("donor.var"))

        #### age and gender
        mcls.freq.AgeGender.lm.list <- list()
        for(aloc in unique(mcls.freq.clin.tb$loc))
        {
            mcls.freq.AgeGender.lm.list[[aloc]] <- as.data.table(ldply(unique(mcls.freq.clin.tb$group.var),
                                                                       function(x){
                                        dat.used.tb <- mcls.freq.clin.tb[group.var==x &
                                                    loc==aloc &
                                                    !(is.na(Age) & is.na(Gender)),]
                                        res.lm <- dat.used.tb[,lm(freq~Age+Gender)]
                                        res.lm.s <- summary(res.lm)
                                        nRich <- dat.used.tb[freq > 0.1,nrow(.SD)]
                                        data.table(meta.cluster=x,
                                               cate="AgeGender",
                                               r.squared=res.lm.s$r.squared,
                                               adj.r.squared=res.lm.s$adj.r.squared,
                                               p.value=anova(res.lm)$'Pr(>F)'[1],
                                               nRich=nRich)
                                        ##p.value=res.lm.s$coefficients["freq.x","Pr(>|t|)"])
                               }))
            mcls.freq.AgeGender.lm.list[[aloc]][,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name.full")[as.character(meta.cluster)]]
            mcls.freq.AgeGender.lm.list[[aloc]][,adj.p.value:=p.adjust(p.value,"BH")]
            write.table(mcls.freq.AgeGender.lm.list[[aloc]],file=sprintf("%s.data.lm.AgeGender.%s.txt",out.prefix,aloc),
                        row.names=F,sep="\t",quote=F)
        }

        ### stage
        mcls.freq.stage.lm.list <- list()
        for(aloc in unique(mcls.freq.clin.tb$loc))
        {
            mcls.freq.stage.lm.list[[aloc]] <- as.data.table(ldply(unique(mcls.freq.clin.tb$group.var),
                                                                   function(x){
                                        dat.used.tb <- mcls.freq.clin.tb[group.var==x &
                                                    loc==aloc &
                                                    !is.na(stage) &
                                                    stage %in% c("I","II","III","IV"),]
                                        res.lm <- dat.used.tb[,lm(freq~stage)]
                                        res.lm.s <- summary(res.lm)
                                        nRich <- dat.used.tb[freq > 0.1,nrow(.SD)]
                                        data.table(meta.cluster=x,
                                               cate="stage",
                                               r.squared=res.lm.s$r.squared,
                                               adj.r.squared=res.lm.s$adj.r.squared,
                                               p.value=anova(res.lm)$'Pr(>F)'[1],
                                               nRich=nRich)
                                        ##p.value=res.lm.s$coefficients["freq.x","Pr(>|t|)"])
                               }))
            mcls.freq.stage.lm.list[[aloc]][,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name.full")[as.character(meta.cluster)]]
            mcls.freq.stage.lm.list[[aloc]][,adj.p.value:=p.adjust(p.value,"BH")]
            write.table(mcls.freq.stage.lm.list[[aloc]],file=sprintf("%s.data.lm.stage.%s.txt",out.prefix,aloc),
                        row.names=F,sep="\t",quote=F)
        }

        ### BMI
        mcls.freq.BMI.lm.list <- list()
        for(aloc in unique(mcls.freq.clin.tb$loc))
        {
            mcls.freq.BMI.lm.list[[aloc]] <- as.data.table(ldply(unique(mcls.freq.clin.tb$group.var),
                                                                 function(x){
                                        dat.used.tb <- mcls.freq.clin.tb[group.var==x &
                                                    loc==aloc &
                                                    !is.na(BMI) ,]
                                        res.lm <- dat.used.tb[,lm(freq~BMI)]
                                        res.lm.s <- summary(res.lm)
                                        nRich <- dat.used.tb[freq > 0.1,nrow(.SD)]
                                        data.table(meta.cluster=x,
                                               cate="BMI",
                                               r.squared=res.lm.s$r.squared,
                                               adj.r.squared=res.lm.s$adj.r.squared,
                                               p.value=anova(res.lm)$'Pr(>F)'[1],
                                               nRich=nRich)
                               }))
            mcls.freq.BMI.lm.list[[aloc]][,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name.full")[as.character(meta.cluster)]]
            mcls.freq.BMI.lm.list[[aloc]][,adj.p.value:=p.adjust(p.value,"BH")]
            write.table(mcls.freq.BMI.lm.list[[aloc]],file=sprintf("%s.data.lm.BMI.%s.txt",out.prefix,aloc),
                        row.names=F,sep="\t",quote=F)
        }

        ## example
        makePlot.BMI.Freq <- function(out.prefix,in.tb,mcls,aloc,my.ylim=NULL,
                      par.stat_cor=list(),
                      par.stat_regline_equation=list(vjust=2.5),...)
        {
            mcls.freq.cmp <- in.tb[loc==aloc & group.var==mcls,]
            mcls.freq.cmp[,size:=(log10(N+1))]
            #mcls.freq.cmp[,isLIHC:=cmp.var %in% c("HCC","CHOL","STAD","UCEC","CRC")]
            p <- ggscatter(mcls.freq.cmp,x="BMI",y="freq",size="size",color="cmp.var.x",...) +
                labs(x="BMI",y="Frequency In Tumor",title=mcls.freq.cmp$group.var[1]) +
                geom_smooth(method='lm') +
                guides(size = guide_legend(order=2,direction="horizontal",
                               title.position = "top", label.position = "bottom"),
                       color=guide_legend(order=1,override.aes=list(size=5),ncol=2),
                       legend.direction="horizontal") +
                  geom_smooth(method='lm') +
                  coord_cartesian(clip="off",ylim=my.ylim) +
                  scale_size(breaks=log10(c(10,30,50,100)),range=c(1,10),
                         labels=c(10,30,50,100),limits=c(0,10)) +
                  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
            p <- p +
                  ##stat_cor(aes(label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")),) +
                  do.call(stat_cor,c(list(mapping=aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))),
                             par.stat_cor)) +
                  ###stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),vjust=2.5) +
                  do.call(stat_regline_equation,c(list(mapping=aes(label =  paste(..eq.label..))),
                                  par.stat_regline_equation)) +
                  ##stat_cor() +
                  #facet_wrap(~isLIHC,ncol=2,scales="free_x") +
                  scale_color_manual(values=colSet[["cancerType"]])
            pp <- plot_grid(p+theme(legend.position = "none"),get_legend(p),
                    ncol=1,rel_heights=c(0.5,0.45),align="v")
            ggsave(sprintf("%s.BMI.%s.pdf",out.prefix,mcls),width=4.5,height=7,useDingbats=F)
            ##ggsave(sprintf("%s.tissueImpact.%s.pdf",out.prefix,mcls),width=7.5,height=4,useDingbats=F)
            return(p)
        }

        ### panC
        p.list <- llply(c("CD4.c16.Tfh.CXCR5", "CD4.c17.TfhTh1.CXCL13",
                          "CD8.c12.Tex.CXCL13"),
                        function(x){
                            makePlot.BMI.Freq(out.prefix,mcls.freq.clin.tb,x,"T")
                        })
        p.tmp <- makePlot.BMI.Freq(out.prefix,mcls.freq.clin.tb,"CD8.c09.Tk.KIR2DL4","T",my.ylim=c(0,0.1),
                       par.stat_cor=list(label.y=0.10),
                       par.stat_regline_equation=list(label.y=0.09))
        p.list <- c(p.list,p.tmp)
        names(p.list) <- c("CD4.c16.Tfh.CXCR5", "CD4.c17.TfhTh1.CXCL13", "CD8.c12.Tex.CXCL13","CD8.c09.Tk.KIR2DL4")
        saveRDS(p.list,file=sprintf("%s.BMI.panC.p.list.rds",out.prefix))

        ### THCA
    ######    p.list <- llply(c("CD4.c16.Tfh.CXCR5", "CD4.c17.TfhTh1.CXCL13", "CD8.c12.Tex.CXCL13"),
    ######		    function(x){
    ######		    makePlot.BMI.Freq(out.prefix=sprintf("%s.THCA",out.prefix),
    ######				      mcls.freq.clin.tb[cmp.var.x=="THCA",],x,"T")
    ######		    })
    ######    p.tmp <- makePlot.BMI.Freq(out.prefix=sprintf("%s.THCA",out.prefix),
    ######			       mcls.freq.clin.tb[cmp.var.x=="THCA",],"CD8.c09.Tk.KIR2DL4","T",my.ylim=c(0,0.1),
    ######			       par.stat_cor=list(label.y=0.10),
    ######			       par.stat_regline_equation=list(label.y=0.09))
    ######    p.list <- c(p.list,p.tmp)
    ######    names(p.list) <- c("CD4.c16.Tfh.CXCR5", "CD4.c17.TfhTh1.CXCL13", "CD8.c12.Tex.CXCL13","CD8.c09.Tk.KIR2DL4")
    ######    saveRDS(p.list,file=sprintf("%s.BMI.THCA.p.list.rds",out.prefix))

    }

    #### prepare data for R^2 (combine all)
    {
        mcls.freq.cancerType.lm.df <- fread(sprintf("%s.data.lm.cancerType.txt",out.prefix))
        mcls.freq.tissueImpact.lm.df <- fread(sprintf("%s.data.lm.tissueImpact.txt",out.prefix))
        mcls.freq.groupTumor.lm.df <- fread(sprintf("%s.data.lm.groupTumor.txt",out.prefix))
        #mcls.freq.groupTumor.cor0.30.lm.df <- fread(sprintf("%s.data.lm.groupTumor.TH.cor.0.30.txt",out.prefix))
        #mcls.freq.groupTumor.cor0.30.lm.df[,cate:=sprintf("%s.0.30",cate)]
        mcls.freq.AgeGender.lm.df <- fread(sprintf("%s.data.lm.AgeGender.T.txt",out.prefix))
        mcls.freq.stage.lm.df <- fread(sprintf("%s.data.lm.stage.T.txt",out.prefix))
        mcls.freq.BMI.lm.df <- fread(sprintf("%s.data.lm.BMI.T.txt",out.prefix))
        mcls.freq.cancerType.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[meta.cluster] ]
        mcls.freq.tissueImpact.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[meta.cluster] ]
        mcls.freq.AgeGender.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[meta.cluster] ]
        mcls.freq.stage.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[meta.cluster] ]
        mcls.freq.BMI.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[meta.cluster] ]

        mcls.freq.TMB.lm.df <- fread("../data/metaInfo/TMB_lm_noRegCancerType.txt")
        mcls.freq.TMB.lm.df[,meta.cluster:=Cluster]
        mcls.freq.TMB.lm.df[,cate:="TMB"]
        mcls.freq.TMB.lm.df[,r.squared:=R2]
        mcls.freq.TMB.lm.df[,adj.r.squared:=adj.R2]
        mcls.freq.TMB.lm.df[,p.value:=P.val]
        mcls.freq.TMB.lm.df[,nRich:=100]
        mcls.freq.TMB.lm.df[,cluster.name:=fetchMetaClusterID2CusterFullName()[Cluster] ]
        mcls.freq.TMB.lm.df[,adj.p.value:=adj.P.val]
        mcls.freq.TMB.lm.df <- mcls.freq.TMB.lm.df[,colnames(mcls.freq.cancerType.lm.df),with=F]

        mcls.freq.tissueImpact.lm.df[adj.p.value < 0.2,]
        mcls.freq.cancerType.lm.df[adj.p.value < 0.2,]
        mcls.freq.groupTumor.lm.df[adj.p.value < 0.2,]

        mcls.freq.lm.df <- rbind(
                     mcls.freq.AgeGender.lm.df,
                     mcls.freq.stage.lm.df,
                     mcls.freq.BMI.lm.df,
                     mcls.freq.tissueImpact.lm.df,
                     mcls.freq.cancerType.lm.df,
                     mcls.freq.TMB.lm.df,
                     #mcls.freq.groupTumor.cor0.30.lm.df,
                     mcls.freq.groupTumor.lm.df)
        saveRDS(mcls.freq.lm.df,file=sprintf("%s.mcls.freq.lm.df.rds",out.prefix))
    }

    #### make the R^2 barplots (fig. S28, fig. S36)
    {

        mcls.freq.lm.df <- readRDS(file=sprintf("%s.mcls.freq.lm.df.rds",out.prefix))
        mcls.freq.lm.df[,cate:=factor(cate,levels=(c("AgeGender","stage","BMI",
                                                     "tissueImpact","cancerType","TMB", "groupTumor")))]
        mcls.freq.lm.df[,meta.cluster:=factor(as.character(meta.cluster),
                          levels=rev(as.character(mcls.freq.lm.df[cate=="groupTumor",
                                      ][order(-adj.r.squared),][["meta.cluster"]])))]
        mcls.freq.lm.df[,cluster.name:=factor(as.character(cluster.name),
                          levels=rev(as.character(mcls.freq.lm.df[cate=="groupTumor",
                                      ][order(-adj.r.squared),][["cluster.name"]])))]
        mcls.freq.lm.df[,pchar:=""]
        mcls.freq.lm.df[adj.p.value<0.1,pchar:="+"]
        mcls.freq.lm.df[adj.p.value<0.05,pchar:="*"]
        mcls.freq.lm.df[adj.p.value<0.01,pchar:="**"]
        mcls.freq.lm.df[adj.p.value<0.001,pchar:="***"]
        mcls.freq.lm.df[,pchar:=factor(pchar,levels=c("***","**","*","+",""))]
        mcls.freq.lm.df[,hjust:=ifelse(adj.r.squared>0, -0.1, 1.1)]
        mcls.freq.lm.df[,vjust:=ifelse(pchar=="+", 0.5, 0.8)]
        mcls.freq.lm.df[,size:=ifelse(pchar=="+", 8, 8)]
        mcls.freq.lm.df[,note.freq:=""]
        mcls.freq.lm.df[p.value<0.05 & nRich<2,note.freq:="(low freq)"]

        p.list <- llply(unique(mcls.freq.lm.df$cate),function(x){
            dat.plot <- mcls.freq.lm.df[cate==x,]
            dat.plot[,meta.cluster:=factor(as.character(meta.cluster),
                               levels=rev(as.character(dat.plot[order(-adj.r.squared),][["meta.cluster"]])))]
            dat.plot[,cluster.name:=factor(as.character(cluster.name),
                               levels=rev(as.character(dat.plot[order(-adj.r.squared),][["cluster.name"]])))]
            p <- ggbarplot(dat.plot,x="cluster.name",y="adj.r.squared",
                           color=NA,
                           fill="p.value") +
                labs(y=expression("Adjusted "*R^2),x="",title=x)+
                #coord_flip(ylim=c(0,if(x=="groupTumor") 0.8 else if(x=="cancerType") 0.4 else 0.5)) +
                coord_flip(ylim=c(0,0.7)) +
                geom_hline(yintercept=c(0.2,0.4,0.6),linetype="dashed",color="lightgray",alpha=0.8)+
                scale_fill_distiller(palette = "Reds",breaks=c(0,0.01,0.02,0.03,0.04,0.05),
                             limits=c(0,0.05),na.value="lightgray",
                             guide = guide_colorbar(label.theme = element_text(angle = 45,hjust = 1))) +
                #guides(size = NULL, fill=guide_legend(order=1,nrow=3,byrow=T,label.position="bottom")) +
                guides(size = FALSE) +
                geom_text(aes(label=pchar,hjust=hjust,vjust=vjust,size=size,group=meta.cluster)) +
                geom_text(aes(label=note.freq,y=Inf,hjust=1,vjust=0.5,size=size,group=meta.cluster)) +
                theme(axis.text.x=element_text(angle=60,hjust=1),axis.title.x=element_text(size=14),
                      legend.key.width=unit(0.8,"cm"))
            ggsave(file=sprintf("%s.mcls.freq.lm.%s.pdf",out.prefix,x),width=5,height=9)
            return(p)
                          })
        names(p.list) <- unique(mcls.freq.lm.df$cate)

        pp <- plot_grid(plot_grid(plotlist=llply(p.list[c("AgeGender","stage","tissueImpact",
                                  "cancerType","TMB","groupTumor")],
                             function(x){ x + theme(legend.position="none") }),
                      #labels=names(p.list),hjust=-1,
                      align="hv",ncol=3),
                get_legend(p.list[[1]]),
                ncol=1,
                rel_heights=c(0.9,0.1))
        ggsave(file=sprintf("%s.mcls.freq.lm.merge.pdf",out.prefix),width=15,height=18)

        pp <- plot_grid(plot_grid(plotlist=llply(p.list[c("AgeGender","stage","BMI",
                                  "tissueImpact","cancerType","TMB")],
                             function(x){ x + theme(legend.position="none") }),
                      #labels=names(p.list),hjust=-1,
                      align="hv",ncol=3),
                get_legend(p.list[[1]]),
                ncol=1,
                rel_heights=c(0.9,0.1))
        ggsave(file=sprintf("%s.mcls.freq.lm.merge.00.pdf",out.prefix),width=15,height=18)

    }

}

#### therapy (fig. S38B, Fig. 5C)
{

    tb.list <- list()

    do.therapy.freq.cmp <- function(ana.tb,out.prefix.therapy,donor.var="patient",min.ncell=50)
    {
        #### one patient may have multiple pre-treatment sample; it is reasonable to merge 
        #### those pre-treatment samples for baseline comparision.
        #### one patient may have multiple samples with different response

        #### filtering ana.tb
        blacklist.sample.tb <- ana.tb[,.N,by="sampleID"][N<min.ncell,]
        blacklist.sample <- blacklist.sample.tb[["sampleID"]]
        ana.tb <- ana.tb[!(sampleID %in% blacklist.sample),]
        ana.tb[,meta.cluster:=as.character(meta.cluster)]
        if(nrow(blacklist.sample.tb)>0){
            cat(sprintf("remove %d samples, %d cells\n",
                        length(blacklist.sample),blacklist.sample.tb[,sum(N)]))
        }

        ana.freq.tb <- sscVis::plotDistFromCellInfoTable(ana.tb[treatment=="baseline",],
                                 out.prefix=NULL,plot.type="none",test.method="wilcox.test",
                                 cmp.var="res",group.var="meta.cluster",donor.var=donor.var)
        ana.freq.ext.tb <- as.data.table(ldply(sort(unique(ana.freq.tb$group.var)),function(x){
                           dat.block <- ana.freq.tb[group.var==x,]
                           dat.block.mean <- dat.block[,.(freq.mean=mean(freq)),by=cmp.var]
                           res.wilcox <- wilcox.test(freq~cmp.var,dat.block)
                           data.table(group.var=x,
                                  freqDiff=dat.block.mean[cmp.var=="R", ][["freq.mean"]]-
                                       dat.block.mean[cmp.var=="NR", ][["freq.mean"]],
                                  p.value=res.wilcox$p.value)
                        }))
        ana.freq.ext.tb[,FDR:=p.adjust(p.value,"BH")]
        ana.freq.ext.tb <- ana.freq.ext.tb[order(FDR,p.value),]
        write.table(ana.freq.ext.tb,file=sprintf("%s.txt",out.prefix.therapy),row.names=F,sep="\t",quote=F)
        ana.freq.ext.sig.tb <- ana.freq.ext.tb[p.value<0.05,]
        #print("ana.freq.tb[,summary(NTotal)]")
        #print(ana.freq.tb[,summary(NTotal)])
        #print(ana.freq.tb[,summary(NTotal>=min.ncell)])
        ana.freq.tb <- ana.freq.tb[NTotal>=min.ncell,]

        l_ply(c(sort(unique(ana.freq.ext.sig.tb$group.var)),"CD8.c14.Tex.TCF7"),function(x){
                  p <- ggboxplot(ana.freq.tb[group.var==x,],x="cmp.var",y="freq",
                                 color="cmp.var",palette="npg",legend="none",
                                 add = "jitter", outlier.shape = NA,
                                 add.params=list(size=3,shape=16),
                                 title=x,xlab="",ylab="Frequency") +
                        stat_compare_means(label="p.format") +
                        theme(
                              #axis.text.x = element_text(angle = 60,hjust = 1),
                              plot.title=element_text(size=9)) +
                        coord_cartesian(clip = "off")
                  ggsave(sprintf("%s.therapy.%s.t.pdf",out.prefix.therapy,x),
                         width=2,height=3.5,useDingbats=F)
                                               })
        return(NULL)

        p <- sscVis::plotDistFromCellInfoTable(ana.tb[treatment=="baseline",],
                               out.prefix=NULL,
                               plot.type="boxplot2",test.method="wilcox.test",
                               cmp.var="res",group.var="meta.cluster",donor.var=donor.var)
        ggsave(file=sprintf("%s.therapy.pdf",out.prefix.therapy),width=6,height=25,useDingbats=F)

        ####
        ##cdata <- dcast(ana.tb,meta.cluster+res+patient~treatment,value.var="cellID")
        sample2treatment.tb <- ana.tb[,.N,by=c("sampleID","patient","treatment")
                                      ][order(patient,sampleID),]
        sample2treatment.tb[,sampleID:=gsub("_",".",sampleID)]
        sample2treatment.vec <- structure(sample2treatment.tb$treatment,
                                          names=sample2treatment.tb$sampleID)
        ana.freq.treatment.tb <- sscVis::plotDistFromCellInfoTable(ana.tb,
                                       out.prefix=NULL,
                                       plot.type="none",test.method="wilcox.test",
                                       cmp.var="res",group.var="meta.cluster",donor.var="sampleID")
        ana.freq.treatment.tb[,res:=factor(sample2treatment.vec[donor.var],
                                           levels=c("baseline","post.treatment"))]
        #print("ana.freq.treatment.tb[,summary(NTotal)]")
        #print(ana.freq.treatment.tb[,summary(NTotal)])
        #print(ana.freq.treatment.tb[,summary(NTotal>=min.ncell)])
        ana.freq.treatment.tb <- ana.freq.treatment.tb[NTotal>=min.ncell,]

        ana.freq.treatment.ext.tb <- ana.freq.treatment.tb[,.(p.value=wilcox.test(freq~res,.SD)$p.value,
                                      freqDiff=mean(.SD[res=="post.treatment",][["freq"]])-
                                      mean(.SD[res=="baseline",][["freq"]]) ),
                                    by=c("group.var","cmp.var")]
        ana.freq.treatment.ext.tb[,FDR:=p.adjust(p.value,"BH")]
        ana.freq.treatment.ext.tb <- ana.freq.treatment.ext.tb[order(FDR,-abs(freqDiff)),]
        write.table(ana.freq.treatment.ext.tb,
                file=sprintf("%s.treatment.txt",out.prefix.therapy),
                row.names=F,sep="\t",quote=F)
        ana.freq.treatment.ext.sig.tb <- ana.freq.treatment.ext.tb[p.value<0.05,]

        ana.freq.treatment.tb[,res:=factor(substring(as.character(res),1,4),levels=c("base","post"))]
        l_ply(seq_len(nrow(ana.freq.treatment.ext.sig.tb)),function(i){
              p <- ggboxplot(ana.freq.treatment.tb[group.var==ana.freq.treatment.ext.sig.tb$group.var[i] &
                             cmp.var==ana.freq.treatment.ext.sig.tb$cmp.var[i],],
                         x="res",y="freq",
                         color="res",palette="Paired",legend="none",
                         add = "jitter", outlier.shape = NA,
                         add.params=list(size=3,shape=16),
                         title=sprintf("%s(%s)",ana.freq.treatment.ext.sig.tb$group.var[i],
                                       ana.freq.treatment.ext.sig.tb$cmp.var[i]),
                             xlab="",ylab="Frequency") +
                    stat_compare_means(label="p.format") +
                    theme(
                          #axis.text.x = element_text(angle = 60,hjust = 1),
                          plot.title=element_text(size=9)) +
                    coord_cartesian(clip = "off")
              ggsave(sprintf("%s.therapy.treatment.%s.%s.pdf",out.prefix.therapy,
                             ana.freq.treatment.ext.sig.tb$group.var[i],
                             ana.freq.treatment.ext.sig.tb$cmp.var[i]),
                     width=2,height=3.5,useDingbats=F)
                                               })

        p <- ggboxplot(ana.freq.treatment.tb,x="res",y="freq",
                             color="res",palette="Paired",legend="none",
                             add = "jitter", outlier.shape = NA,
                             add.params=list(size=3,shape=16),
                             xlab="",ylab="Frequency") +
                    facet_wrap(group.var~cmp.var,ncol=8,scales="free_y") +
                    stat_compare_means(label="p.format") +
                    theme(axis.text.x = element_text(angle = 60,hjust = 1),
                          plot.title=element_text(size=9)) +
                    coord_cartesian(clip = "off")
        ggsave(file=sprintf("%s.therapy.treatment.pdf",out.prefix.therapy),useDingbats=F,
               width=14,height=20)



    }

    #meta.tb[,.(N=.N),by=c("dataset","treatment")][treatment!="baseline",]

    #### MELA.MosheSade-Feldman2018
    out.prefix.therapy <- sprintf("%s/therapy/%s.MELA.MosheSade-Feldman2018",dirname(out.prefix),basename(out.prefix))
    dir.create(dirname(out.prefix.therapy),F,T)
    MELA.therapy.tb <- fread("../data/metaInfo/MELA.MosheSade-Feldman2018.therapy.txt")
    MELA.therapy.vec <- structure(MELA.therapy.tb$res,names=MELA.therapy.tb$sampleID)
    ana.tb <- meta.tb[dataset=="MELA.MosheSade-Feldman2018" & loc=="T",]
    ana.tb$res <- MELA.therapy.vec[ana.tb$sampleID]
    #### by freq of meta-cluster
    do.therapy.freq.cmp(ana.tb,out.prefix.therapy)
    tb.list[["MELA.MosheSade-Feldman2018"]] <- ana.tb

    #### by immune type
    if(F)
    {
        cluster.sample.T.tb <- fread(file=sprintf("%s.sample.T.donorCluster.txt",out.prefix))
        cluster.sample.T.tb[1,]
        dat.for.fisher.tb <- ana.tb[dataset=="MELA.MosheSade-Feldman2018" & loc=="T",.N,
                                    by=c("cancerType","dataset",
                                        #"sampleID","treatment",
                                        "patient.uid","res"
                                        )]
        dat.for.fisher.tb <- merge(dat.for.fisher.tb,cluster.sample.T.tb[,
                                   c("donorID","displayCluster")],
                                by.x="patient.uid", by.y="donorID")
        f.dup <- dat.for.fisher.tb[duplicated(patient.uid),][["patient.uid"]]
        dat.for.fisher.tb[patient.uid %in% f.dup,]
        dat.for.fisher.tb <- dat.for.fisher.tb[!patient.uid %in% f.dup,]
        dat.for.fisher.tb[displayCluster %in% c(1,2),itype:="C0102"]
        dat.for.fisher.tb[!displayCluster %in% c(1,2),itype:="C03_08"]

        dat.for.fisher.tb[,table(displayCluster,res)]
        test.aa <- dat.for.fisher.tb[,table(itype,res)]
        test.res <- fisher.test(test.aa)
        test.aa.tb <- as.data.table(test.aa)

        p <- ggbarplot(test.aa.tb,x="res",y="N",group="itype",width=1,
               fill="itype",palette="npg") +
            labs(subtitle=sprintf("p:%4.2f, OR: %4.2f",test.res$p.value,test.res$estimate)) +
            #coord_cartesian(ylim=c(0,10)) +
            scale_y_continuous(name="N",breaks=seq(0,10,2),limits=c(0,10))
        ggsave(sprintf("%s.fisher.barplot.pdf",out.prefix.therapy),width=2,height=4)
    }

    #### BCC.KathrynEYost2019 & SCC.KathrynEYost2019 (no significant :-( )
#    out.prefix.therapy <- sprintf("%s/therapy/%s.BCC.SCC.KathrynEYost2019",dirname(out.prefix),basename(out.prefix))
#    dir.create(dirname(out.prefix.therapy),F,T)
#    ana.tb <- meta.tb[dataset %in% c("BCC.KathrynEYost2019","SCC.KathrynEYost2019") & loc=="T",]
#    ana.tb <- ana.tb[!(sampleID %in% c("scc.su010.post.cd39")),]
#    ana.tb$res <- "NR"
#    ana.tb[(patient %in% c("su001","su002","su003","su004","su009","su012","su011"))|
#	       (patient %in% c("su010") & cancerType=="SCC" ),res:="R"]
#    #ana.tb$res <- MELA.therapy.vec[ana.tb$sampleID]
#    do.therapy.freq.cmp(ana.tb,out.prefix.therapy)
#    tb.list[["BCC.SCC.KathrynEYost2019"]] <- ana.tb

    ### multi- datassets
#    ana.pan.tb <- do.call(rbind,tb.list)
#    out.prefix.therapy <- sprintf("%s/therapy/%s.panC",dirname(out.prefix),basename(out.prefix))
#    do.therapy.freq.cmp(ana.pan.tb,out.prefix.therapy)

}

