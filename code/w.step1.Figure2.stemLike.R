#!/usr/bin/env Rscript

library("data.table")
library("sscVis")
library("plyr")
library("ggplot2")
library("ggpubr")
library("ggplotify")
library("cowplot")
source("func.R")
#source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/inte.comb.miniClust.lib.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 8)

out.prefix <- "OUT_Fig2/stemLike/CD8.stemLike"
dir.create(dirname(out.prefix),F,T)

###

colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
#colSet$cluster.name <- colSet$meta.cluster
#names(colSet$cluster.name) <- fetchMetaClusterID2CusterFullName()[names(colSet$meta.cluster)]
#f.na <- is.na(names(colSet$cluster.name))
#colSet$cluster.name <- colSet$cluster.name[!f.na]

meta.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")
mapping.tb <- unique(meta.tb[,c("cancerType","dataset.old")])
datasetOld2cancerType <- c(structure(mapping.tb$cancerType,names=mapping.tb$dataset.old),
						   "PanC"="PanC")

######## gene expression data ####
sce.list.file <- "list/sce.CD8.fullPath.list"
sce.list.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.file),head=T)
###colnames(sce.list.tb) <- c("data.id", "measurement", "platform", "seufile", "scefile")

sce.list <- llply(seq_len(nrow(sce.list.tb)),function(i){ readRDS(sce.list.tb$scefile[i]) },.parallel=T)
names(sce.list) <- sce.list.tb$data.id

##############

calBinExp <- function(sce.list,
		       gene.to.test=c("FOXP3","CXCL13"),
		       dataset.x="OV.zhangLab5P",mcls="CD8.c12.Tex.CXCL13",
		       th.exprs=0.3, min.ncell=50,val.padding=0)
{

    sce <- sce.list[[dataset.x]]
    
    #### filter out dataset with cell number in Tex less than 50
    if(sum(sce[["meta.cluster"]]==mcls) < min.ncell){
	    return(NULL)
    }
    ####
    cat(sprintf("dataset (%s)\n",dataset.x))

    sce.plot <- ssc.scale(sce,gene.symbol=gene.to.test,adjB="batchV",do.scale=T)
    sce.plot <- sce.plot[,sce.plot$meta.cluster==mcls]

    dat.plot <- as.data.frame(t(as.matrix(assay(sce.plot,"norm_exprs.scale"))))
    colnames(dat.plot) <- rowData(sce.plot)$display.name
    col.dropout <- setdiff(gene.to.test,colnames(dat.plot))
    if(length(col.dropout)>0){
	dropout.mtx <- matrix(rep(val.padding,nrow(dat.plot)*length(col.dropout)),ncol=length(col.dropout))
	colnames(dropout.mtx) <- col.dropout
	dat.plot <- cbind(dat.plot,dropout.mtx)
    }
    dat.plot <- dat.plot[,gene.to.test]

    f.na <- is.na(dat.plot)
    dat.plot[f.na] <- val.padding

    for(gg in colnames(dat.plot)){
	    sce.plot[[sprintf("bin.%s",gg)]] <- as.integer(dat.plot[[gg]] > th.exprs)
    }

    ret.tb <- as.data.table(colData(sce.plot))[,c("cellID","cancerType","dataset","patient",
						  "cellID.uniq","miniCluster","meta.cluster",
						  sprintf("bin.%s",colnames(dat.plot)) ),with=F]

    #cat(sprintf("successful (%s)\n",dataset.x))

    return(ret.tb)

}

gene.in.box <- c("CXCL13","HAVCR2","FOXP3","TOX","TCF7","CXCR5",
		 "TIGIT","CTLA4","PDCD1","TNFRSF9","LAG3",
		 "CCR7","SELL","GZMK")
dataset.ana.vec <- names(sce.list)

### 
dat.bin.exp.other.tb <- as.data.table(ldply(c("CD8.c05.Tem.CXCR5","CD8.c14.Tex.TCF7"),function(aid){
    as.data.table(ldply(dataset.ana.vec,function(x){
				      calBinExp(sce.list,
						 gene.to.test=gene.in.box,
						 dataset.x=x,mcls=aid)
		    }))
}))


##############
run.it <- function(sce.list,out.prefix,
		   gene.to.test=c("FOXP3","CXCL13"), dataset.x="OV.zhangLab5P",mcls="CD8.c12.Tex.CXCL13",
		   th.exprs=0.3, th.ncell=2,
           ###show.all.expanded.clone=F,
		   note.txt="TCRSharing",min.ncell=50,
		   pdf.width=5,pdf.height=5)
{

    dir.create(dirname(out.prefix),F,T)

    sce <- sce.list[[dataset.x]]
    
    #### filter out dataset with cell number in Tex less than 50
    if(sum(sce[["meta.cluster"]]==mcls) < min.ncell){
	    return(NULL)
    }
    ####
    cat(sprintf("dataset (%s, %s)\n",dataset.x,note.txt))

    sce.plot <- ssc.scale(sce,gene.symbol=gene.to.test,adjB="batchV",do.scale=T)
    sce.plot <- sce.plot[,sce.plot$meta.cluster==mcls]
    sce.plot <- sce.plot[match(gene.to.test,rowData(sce.plot)$display.name),]

    dat.plot <- as.data.frame(t(as.matrix(assay(sce.plot,"norm_exprs.scale"))))
    colnames(dat.plot) <- rowData(sce.plot)$display.name
    dat.plot <- dat.plot[,gene.to.test]
    #dat.plot.melt <- melt(as.matrix(dat.plot))
    #colnames(dat.plot.melt) <- c("cell","gene","exprs.z")

    if(all(is.na(dat.plot[,1])) || all(is.na(dat.plot[,2]))){
	    return(NULL)
    }

    for(gg in colnames(dat.plot)){
	    sce.plot[[sprintf("bin.%s",gg)]] <- as.integer(dat.plot[[gg]] > th.exprs)
    }

    gene.1 <- colnames(dat.plot)[1]
    gene.2 <- colnames(dat.plot)[2]
    sce.plot$xtype <- "other"
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==1 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==0] <- sprintf("SP.%s",gene.1)
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==0 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==1] <- sprintf("SP.%s",gene.2)
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==1 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==1] <- sprintf("DP")
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==0 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==0] <- sprintf("DN")

    #print(table(sce.plot$xtype))


    #f.cell <- intersect(colnames(sce.plot),dat.tcr[["Cell_Name"]])
    #dat.tcr.use <- dat.tcr[colnames(sce.plot),]
    #sce.plot$cloneID <- dat.tcr.use$cloneID

    ### make a plot
    {

        unclass(table(sce.plot$xtype))
        dat.fisher <- c(sum(sce.plot$xtype=="DP"),
                sum(sce.plot$xtype==sprintf("SP.%s",gene.2)),
                sum(sce.plot$xtype==sprintf("SP.%s",gene.1)),
                sum(sce.plot$xtype=="DN"))
        res.fisher <- fisher.test(matrix(dat.fisher,ncol=2))
        dat.fisher.tb <- as.data.table(matrix(dat.fisher,nrow=1))
        colnames(dat.fisher.tb) <- c("DP",sprintf("SP.%s",gene.2),sprintf("SP.%s",gene.1),"DN")
        dat.fisher.tb <- cbind(data.table(dataset=dataset.x),
                       dat.fisher.tb,
                       data.table(OR=res.fisher$estimate,
                          p.value=res.fisher$p.value))
        if(F)
        {
            #sce.plot.debug <<- sce.plot
            #res.fisher.debug <<- res.fisher
            tryCatch({
                ssc.plot.GeneDensity(sce.plot,
                             out.prefix=sprintf("%s.clone.xtype.%s.cor",out.prefix,note.txt),
                             gene.symbol=rowData(sce.plot)$display.name,
                             assay.name="norm_exprs.scale",
                             expT = c(th.exprs, th.exprs),
                             ##expT = c(0.5, 0.5),
                             my.title=sprintf("%s\n(OR:%4.4f, p: %4.4f)",dataset.x,
                                      res.fisher$estimate,res.fisher$p.value))
            },error=function(cond) { print(cond) })
        }
        sce.plot.ret <- sce.plot
        #### 
        #sce.plot$xtype <- factor(sce.plot$xtype,levels=(c("DP", sprintf("SP.%s",gene.1), sprintf("SP.%s",gene.2), "DN")))

    }
    return(list("res.fisher"=res.fisher,
		"dat.fisher.tb"=dat.fisher.tb,
		"sce.mcls"=sce.plot.ret))
}


#####  show gene expression using 2D density plot (fig. S19A)
{
    gene.pairs.list <- list("CD8.c05.Tem.CXCR5"=c("CXCR5","TCF7"),
                            "CD8.c14.Tex.TCF7"=c("CXCR5","TCF7"))

    res.all <- llply(names(gene.pairs.list),function(aid){
        res.full <-llply(dataset.ana.vec,function(x){
                 run.it(sce.list,out.prefix=sprintf("%s/%s/%s.%s",
                                    dirname(out.prefix),
                                    sprintf("%s",aid),
                                    basename(out.prefix),x),
                    gene.to.test=gene.pairs.list[[aid]],
                    dataset.x=x,
                    mcls=aid,
                    pdf.width=3,
                    ###show.all.expanded.clone=T,
                    ##note.txt=sprintf("ShareWith%s.all",aid)
                    note.txt=sprintf("%s",aid))
            },.parallel=T)
        names(res.full) <- dataset.ana.vec
        cor.full.tb <- as.data.table(ldply(dataset.ana.vec,function(x){
                           if(!is.null(res.full[[x]])){
                               res.full[[x]]$dat.fisher.tb
                           }
                    }))

        #####
        sce.fisher.list <- llply(dataset.ana.vec,function(x){
                     if(!is.null(res.full[[x]])){
                         sce.ret <- res.full[[x]]$sce.mcls
                         dat.tmp <- assay(sce.ret,"norm_exprs.scale")
                         rownames(dat.tmp) <- unname(rowData(sce.ret)$display.name)
                         colnames(dat.tmp) <- sprintf("%s.%s",x,colnames(dat.tmp))
                         sce.ret <- ssc.build(dat.tmp,assay.name="norm_exprs.scale")
                         return(sce.ret)
            } })
        names(sce.fisher.list) <- dataset.ana.vec
        f.dataset <- sapply(dataset.ana.vec,function(x){ !is.null(res.full[[x]]) })
        sce.tmp.list <- sce.fisher.list[ dataset.ana.vec[f.dataset] ]
        sce.fisher <- do.call(cbind,sce.tmp.list)

        fisher.mtx <- matrix(c(sum(cor.full.tb[["DP"]]),
                   sum(cor.full.tb[[ sprintf("SP.%s",gene.pairs.list[[aid]][1]) ]]),
                   sum(cor.full.tb[[ sprintf("SP.%s",gene.pairs.list[[aid]][2]) ]]),
                   sum(cor.full.tb[["DN"]])),ncol=2)
        res.fisher <- fisher.test(fisher.mtx)

        ssc.plot.GeneDensity(sce.fisher,
                 out.prefix=sprintf("%s/%s/%s.%s",
                            dirname(out.prefix),
                            sprintf("%s",aid),
                            basename(out.prefix),"panC"),
                 gene.symbol=rowData(sce.fisher)$display.name,
                 assay.name="norm_exprs.scale",
                 expT = c(0.3,0.3),
                 my.title=sprintf("%s(OR:%4.4f, p: %s)","panC",
                          res.fisher$estimate,res.fisher$p.value))

        return(list("res.full"=res.full,"cor.full.tb"=cor.full.tb))
    })
    names(res.all) <- names(gene.pairs.list)

    #saveRDS(res.all,sprintf("%s.res.all.rds",out.prefix))

}


dat.bin.exp.stem.tb <- as.data.table(ldply(names(res.all),function(aid){
    dataset.use <- names(res.all[[aid]]$res.full)
    dat.aid <- as.data.table(ldply(dataset.use,function(x){
			    sce.a <- res.all[[aid]]$res.full[[x]]$sce.mcls
			    if(!is.null(sce.a)){
				as.data.table(colData(sce.a)[,c("cellID","cancerType","dataset","patient","cellID.uniq",
						  "miniCluster","meta.cluster",
						  "bin.CXCR5","bin.TCF7","xtype")])
			    }
    }))
    dat.aid
}))

dat.bin.exp.stem.tb[,.(N=.N),by=c("meta.cluster","xtype")]
dat.bin.exp.stem.tb[1:2,]

dat.bin.exp.other.tb[1:2,]

dat.bin.exp.tb <- merge(dat.bin.exp.other.tb,dat.bin.exp.stem.tb[,c("meta.cluster","cellID.uniq","xtype"),with=F],all=T)
dat.bin.exp.tb[,cancerType:=datasetOld2cancerType[dataset]]
dat.bin.exp.tb[dataset=="CHOL.YaoHe10X",cancerType:="CHOL"]
dat.bin.exp.tb[1:2,]
dim(dat.bin.exp.tb)
dim(dat.bin.exp.other.tb)
dim(dat.bin.exp.stem.tb)
dat.bin.exp.stem.tb[,.N,by=c("meta.cluster")]
dat.bin.exp.other.tb[,.N,by=c("meta.cluster")]
dat.bin.exp.tb[,.N,by=c("meta.cluster")]
dat.bin.exp.tb[is.na(xtype),.N,by=c("cancerType","dataset","meta.cluster")]

test.list <- list("union"=c("DP","SP.CXCR5","SP.TCF7"),
		  "intersection"=c("DP"),
		  "CXCR5"=c("DP","SP.CXCR5"),
		  "TCF7"=c("DP","SP.TCF7"))

#gene.in.box <- c("CXCL13","HAVCR2","FOXP3","TOX","TCF7","CXCR5","TIGIT","CTLA4","PDCD1")

makePairdPlot <- function(dat.bin.exp.tb,aid,out.prefix,gene.in.box,
			  use.facet=F,facet.nrow=2,pdf.width=8*0.8,pdf.height=8.5)
{
    dir.create(dirname(out.prefix),F,T)
    if(!use.facet)
    {
	plot.list <- llply(gene.in.box,function(a.gene){
	    dat.plot.freq <- dat.bin.exp.tb[xtype %in% test.list[[aid]],][,.(N=.N,freq=sum(.SD[[sprintf("bin.%s",a.gene)]])/.N),
							  by=c("cancerType","dataset","meta.cluster")
							  ][order(meta.cluster,freq),]

	    dat.plot.freq.pair <- dcast(dat.plot.freq,cancerType+dataset~meta.cluster,
					value.var="freq")[!is.na(CD8.c05.Tem.CXCR5) & !is.na(CD8.c14.Tex.TCF7),]
	    dat.plot.freq.pair.melt <- melt(dat.plot.freq.pair,id.vars=c("cancerType","dataset"),
					    variable.name = "meta.cluster", value.name = "Freq")

	    p <- ggpaired(dat.plot.freq.pair.melt,x="meta.cluster",y="Freq",
			   add = "jitter",outlier.shape=NA,
			   line.color = "gray", line.size=0.4,
			   #color="meta.cluster",
			   fill="meta.cluster") +
		    scale_fill_manual(values=colSet$meta.cluster) +
		    #scale_color_manual(values=colSet$meta.cluster) +
		    labs(x="",y="Freq",title=sprintf("%s",a.gene)) +
		    stat_compare_means(label="p.format",paired=T) +
		    theme(axis.text.x=element_text(angle=60,hjust=1),
			  plot.title=element_text(hjust=0.5),
			  legend.position="none")
	    ggsave(sprintf("%s.%s.box.pair.%s.pdf",out.prefix,aid,a.gene),width=2.0,height=4.5,useDingbats=F)
	    return(p)
	})
	names(plot.list) <- gene.in.box
	return(plot.list)
    }else{

	dat.plot.freq.pair.melt <- as.data.table(ldply(gene.in.box,function(a.gene){
				   ret.tb <- dat.bin.exp.tb[xtype %in% test.list[[aid]],
						  ][,.(N=.N,freq=sum(.SD[[sprintf("bin.%s",a.gene)]])/.N),
						      by=c("cancerType","dataset","meta.cluster")
						      ][order(meta.cluster,freq),]

				    ret.tb.pair <- dcast(ret.tb,cancerType+dataset~meta.cluster,
							 value.var="freq")[!is.na(CD8.c05.Tem.CXCR5) &
									    !is.na(CD8.c14.Tex.TCF7),]
				    ret.tb.pair.melt <- melt(ret.tb.pair,id.vars=c("cancerType","dataset"),
							     variable.name = "meta.cluster", value.name = "Freq")
				    ret.tb.pair.melt[,gene:=a.gene]
				    return(ret.tb.pair.melt)
	}))
	dat.plot.freq.pair.melt[,gene:=factor(gene,levels=gene.in.box)]

	#dat.plot.freq.pair.melt[,meta.cluster:=as.character(meta.cluster)]
	dat.plot.freq.pair.melt <- dat.plot.freq.pair.melt[order(meta.cluster,gene,cancerType,dataset),]
	dat.plot.freq.pair.melt[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
	#print(str(dat.plot.freq.pair.melt))

	#dat.plot.freq.pair.melt.debug <<- dat.plot.freq.pair.melt
	p <- ggpaired(dat.plot.freq.pair.melt,
		      #x="meta.cluster",y="Freq",
		      x="cluster.name",y="Freq",
		       add = "jitter",outlier.shape=NA,
		       line.color = "gray", line.size=0.4,
		       fill="cluster.name") +
		scale_fill_manual(values=colSet$cluster.name) +
		labs(x="",y="Freq") +
		stat_compare_means(label="p.format",paired=T) +
		facet_wrap(~gene,nrow=facet.nrow) +
		theme(axis.text.x=element_text(angle=60,hjust=1),
		      plot.title=element_text(hjust=0.5),
		      strip.text=element_text(size=14),
		      #strip.background=element_blank(),
		      legend.position="none")
	ggsave(sprintf("%s.%s.box.pair.pdf",out.prefix,aid),
	       #width=facet.ncol*2*0.8,
	       width=pdf.width,
	       height=pdf.height,useDingbats=F)
	return(p)

    }
}

#####  (fig. S19B)
for(aid in "union")
{
    makePairdPlot(dat.bin.exp.tb,aid,
		  sprintf("%s.%s.allInOne.v2",out.prefix,aid),
		  gene.in.box=c("TOX","PDCD1","TIGIT","CTLA4","TNFRSF9"),
		  facet.nrow=1,pdf.width=2*5*0.7,pdf.height=5,
		  use.facet=T)
}

#saveRDS(dat.bin.exp.tb,file=sprintf("%s.dat.bin.exp.tb.rds",out.prefix))


