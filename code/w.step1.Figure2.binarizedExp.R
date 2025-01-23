#!/usr/bin/env Rscript

library("data.table")
library("sscVis")
library("plyr")
library("ggplot2")
library("ggpubr")
library("ggplotify")
library("cowplot")
library("tictoc")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)
#source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/inte.comb.miniClust.lib.R")

out.prefix <- "OUT_Fig2/binExp.CD8.Tex/binExp.CD8"
dir.create(dirname(out.prefix),F,T)

###
{

    colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    colSet <- c(colSet,list("tech"=c("10X"="#7CAE00","SmartSeq2"="#00BFC4")))
    ##colSet <- c(colSet,list("tech"=c("10X"="#F8766D","SmartSeq2"="#00BFC4")))
    ###
    cellInfo.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")
    meta.tb <- cellInfo.tb
    #geneTablongLong <- readRDS("OUT.core.signature/CD8.multiAsTwo.minCell50/panC.core.signature.CD8.multiAsTwo.minCell50.geneTableLong.collapsed.rds")
    gene.desc.tb <- readRDS("../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds")
    #meta.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")
    mapping.tb <- unique(meta.tb[,c("cancerType","dataset.old")])
    datasetOld2cancerType <- c(structure(mapping.tb$cancerType,names=mapping.tb$dataset.old),
                               "PanC"="PanC")

    gene.desc.Tex.tb <- gene.desc.tb[meta.cluster=="CD8.c12.Tex.CXCL13",]
    ##setkey(gene.desc.Tex.tb,"geneID")
    gene.desc.Tex.tb[1,]
    #geneTablongLong.Tex <- geneTablongLong[meta.cluster=="CD8.c12.Tex.CXCL13",]
    ###

}

######## gene expression data ####
sce.list.file <- "list/sce.CD8.fullPath.list"
sce.list.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.file),head=T)
###colnames(sce.list.tb) <- c("data.id", "measurement", "platform", "seufile", "scefile")

sce.list <- llply(seq_len(nrow(sce.list.tb)),function(i){ readRDS(sce.list.tb$scefile[i]) },.parallel=T)
names(sce.list) <- sce.list.tb$data.id

##############

calBinExp <- function(sce.list,
		      cellIDToRID.tb,
		       gene.to.test=c("FOXP3","CXCL13"),locations=c("P","N","T"),
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
    #if(dataset.x=="CHOL.YaoHe10X"){
    #    m.tb <- cellIDToRID.tb[cancerType=="CHOL" & dataset=="HCC.QimingZhang2019.10X",]
    #}else
    {
        ##m.tb <- cellIDToRID.tb[dataset.old==dataset.x,]
        m.tb <- cellIDToRID.tb[dataset==dataset.x,]
    }
    m.vec <- structure(m.tb$rid,names=m.tb$cellID)
    sampleInfo.x.tb <- m.tb[,.N,by=setdiff(colnames(m.tb),c("cellID","cellID.uniq"))]

    sce.plot <- ssc.scale(sce,gene.symbol=gene.to.test,adjB="batchV",do.scale=T,add.padding=T)
    sce.plot <- sce.plot[,sce.plot$meta.cluster==mcls & sce.plot$loc %in% locations]
    assay(sce.plot,"binExp") <- assay(sce.plot,"norm_exprs.scale") > th.exprs

    ### only after subseting sce.plot, assign rid
    sce.plot$rid <- m.vec[colnames(sce.plot)]
    print(table(sce.plot$rid,useNA="always"))
    
    ##f.na.lib <- is.na(sce.plot$libraryID)
    ##sce.plot$libraryID[f.na.lib] <- "NA"
    ##col.avg <- c("patient","cancerType","dataset","batchV",
    ##		 "meta.cluster","loc")
    col.avg.n.tb <- as.data.table(colData(sce.plot))[,.(NCells.mcls=.N),by="rid"]
    sce.plot.avg <- ssc.average.cell(sce.plot,assay.name="binExp",
				     column="rid",
				     ret.type="sce",do.parallel=T)
    z.dat.mtx <- ssc.average.cell(sce.plot,assay.name="norm_exprs.scale",
				     column="rid",
				     ret.type="data.mtx",do.parallel=T)
    cat(sprintf("check colnames and rownames of sce.plot.avg and z.dat.mtx:\n"))
    print(all(colnames(sce.plot.avg)==colnames(z.dat.mtx)))
    print(all(rownames(sce.plot.avg)==rownames(z.dat.mtx)))
    assay(sce.plot.avg,"z.score") <- z.dat.mtx
    col.avg.n.df <- as.data.frame(merge(sampleInfo.x.tb,
					merge(as.data.table(colData(sce.plot.avg)),col.avg.n.tb),
					by="rid"))
    rownames(col.avg.n.df) <- col.avg.n.df$rid
    col.avg.n.df <- col.avg.n.df[colnames(sce.plot.avg),]
    colData(sce.plot.avg) <- DataFrame(col.avg.n.df)
    rownames(sce.plot.avg) <- rowData(sce.plot.avg)$display.name
    sce.plot.avg <- sce.plot.avg[gene.to.test,]
    sce.plot.avg$cid <- NULL

    return(sce.plot.avg)

}

gene.in.box <- gene.desc.tb[meta.cluster=="CD8.c12.Tex.CXCL13",][["geneID"]]

dataset.ana.vec <- names(sce.list)

col.group <- c("patient","cancerType","dataset",
				  ####"sampleID","batchV",
				  "treatment",
				  "patient.uid",
				  "cancerType.old","dataset.old",
				  "usedForFreq","dataSource","tech","tech.cate","pub",
				  "dataset.tech","loc")
sampleInfo.tb <- cellInfo.tb[stype=="CD8" & loc=="T",
			     .(NCells=.N),
			     by=col.group ]
sampleInfo.tb$rid <- sprintf("S%04d",seq_len(nrow(sampleInfo.tb)))
sampleInfo.tb[1:2,]

cellIDToRID.tb <- merge(cellInfo.tb[stype=="CD8" & loc=="T",c("cellID","cellID.uniq",col.group),with=F],sampleInfo.tb,all=T)
cellIDToRID.tb[1:2,]
cellIDToRID.tb[is.na(rid),]
##### change dataset names
cellIDToRID.tb[,dataset:=gsub("BRCA.ElhamAzizi2018.InDrop","BRCA.ElhamAzizi2018_Indrop",dataset)]
cellIDToRID.tb[,dataset:=gsub("HCC.QimingZhang2019.SS2","HCC.QimingZhang2019_SS2",dataset)]
cellIDToRID.tb[,dataset:=gsub("HCC.QimingZhang2019.10X","HCC.QimingZhang2019_10X",dataset)]
cellIDToRID.tb[,dataset:=gsub("NSCLC.DietherLambrechts2018","LC.DietherLambrechts2018",dataset)]
cellIDToRID.tb[,dataset:=gsub("NSCLC.RapolasZilionis2019","LC.RapolasZilionis2019",dataset)]
cellIDToRID.tb[,dataset:=gsub("BRCA.ElhamAzizi2018.10X","BRCA.ElhamAzizi2018_10X",dataset)]
cellIDToRID.tb[,dataset:=gsub("CRC.LeiZhang2020.10X","CRC.LeiZhang2020_10X",dataset)]
cellIDToRID.tb[,dataset:=gsub("NSCLC.QianqianSong2019","LC.QianqianSong2019",dataset)]
cellIDToRID.tb[,dataset:=gsub("NSCLC.XinyiGuo2018","LC.XinyiGuo2018",dataset)]
cellIDToRID.tb[,dataset:=gsub("STAD.BoxiKang2021","STAD.BoxiKang2019",dataset)]

## check
#sampleInfo.tb[,.N,by=c("cancerType","dataset","patient")]
#sampleInfo.tb[,.N,by=c("cancerType","dataset","patient")][N>1,]
#sampleInfo.tb[patient=="CRC.P1228",]
#sampleInfo.tb[dataset=="BRCA.PeterSavas2018" & patient=="S01",]
#sampleInfo.tb[patient=="Mel129",]
#sampleInfo.tb[dataset=="MELA.MosheSade-Feldman2018" & patient=="P3",]
#sampleInfo.tb[dataset=="AML.PeterVanGalen2019" & patient=="AML556",]
#sampleInfo.tb[dataset=="AML.PeterVanGalen2019" & patient=="AML870",]

### 
tic("calBinExp")
sce.binExp.list <- llply(dataset.ana.vec,function(x){
			     calBinExp(sce.list,cellIDToRID.tb, gene.to.test=gene.in.box,
				       locations="T",
				       dataset.x=x,mcls="CD8.c12.Tex.CXCL13")
		  },.parallel=F)
names(sce.binExp.list) <- dataset.ana.vec
toc()

sce.binExp.all <- NULL
for(x in names(sce.binExp.list)){
    if(!is.null(sce.binExp.list[[x]])){
        sce.tmp <- sce.binExp.list[[x]]
        if(is.null(sce.binExp.all)){
            sce.binExp.all <- sce.tmp
        }else{
            metadata(sce.tmp)$ssc <- NULL
            rowData(sce.tmp)$display.name <- NULL
            sce.binExp.all <- cbind(sce.binExp.all,sce.tmp)
        }
    }
}
as.data.table(colData(sce.binExp.all))[,all(NCells==N)]
sce.binExp.all$N <- NULL

saveRDS(sce.binExp.all,file=sprintf("%s.sce.binExp.rds",out.prefix))
#### this run
#sce.binExp.all <- readRDS(file=sprintf("%s.sce.binExp.rds",out.prefix))
#### original result
#sce.binExp.all <- readRDS(file="../data/expression/CD8/integration/binExp.CD8.sce.binExp.rds")

as.data.table(colData(sce.binExp.all))[tech=="10X" & NCells.mcls >= 10,.N,by=c("cancerType","dataset","tech")]
as.data.table(colData(sce.binExp.all))[tech=="10X" & NCells.mcls >= 10,.N,by=c("cancerType","dataset","tech")][,sum(N)]
as.data.table(colData(sce.binExp.all))[NCells.mcls >= 10,.N,by=c("cancerType","dataset","tech")][,sum(N)]

gene.list.plot <- c(
		    "TCF7",
		    "STAT2","IFIT1","IFI44L","ISG20","STAT1",
		    "FOXP3","TXK","KIR2DL3","IL1RL1",
		    "IL26","IL17A","IL17F","RORC","RORA","IL23R","KLRB1","CTSH", "IL17RE", "IL18RAP", "IFNGR1",
		    "CCR6", "TMIGD2", "IL4I1", "ZBTB16", "NCR3", "CCL20", "COLQ","SLC4A10",
		    "PDCD1","CD274","PDCD1LG2","CTLA4","CD80","CD86",
		    "ZNF683","GZMK","NTRK1","NTRK2","HIVEP1","TNFRSF4","TNFRSF9",
		    "CXCL13","HAVCR2","ENTPD1","LAYN","CTLA4",
		    "TOX","TOX2","ZBED2","RBPJ","SOX4",
		    "MKI67","STMN1",
		    "NFATC1","NFATC2","NR5A2","ARID5B","ETV1"
		    )

### boxplot
make.boxplot.cmpCancerType <- function(sce.binExp.all,gene.plot,out.prefix,min.NCell.mcls=10,
				       sort.by="rank.mean",
				       assay.name="binExp",
				       pdf.width=5,pdf.height=2.5,
				       split.by=NULL)
{
    p.list <- llply(gene.plot,function(a.gene){
	if(!a.gene %in% rownames(sce.binExp.all)){ 
	    print(sprintf("gene not found: %s",a.gene))
	    return(NULL) 
	}
	print(sprintf("plot gene: %s",a.gene))
	test.tb <- data.table(cancerType=as.character(sce.binExp.all$cancerType),
			      geneE=assay(sce.binExp.all,assay.name)[a.gene,],
			      treatment=sce.binExp.all$treatment,
			      NCells=sce.binExp.all$NCells,
			      patient=sce.binExp.all$patient,
			      tech=sce.binExp.all$tech,
			      NCells.mcls=sce.binExp.all$NCells.mcls)
	test.tb <- test.tb[NCells.mcls>=min.NCell.mcls & treatment=="baseline",]
	#summary(aov(geneE~cancerType,test.tb))
	#kruskal.test(data=test.tb,geneE~cancerType)

	if(is.null(split.by))
	{
	    test.tb <- test.tb[,geneE.rnk:=rank(geneE)]
	    if(sort.by=="rank.mean"){
		cancerType.order.tb <- test.tb[, .(medV=mean(geneE.rnk)),
					       by="cancerType"][order(medV),]
	    }else{
		cancerType.order.tb <- test.tb[, .(medV=median(geneE),meanV=mean(geneE)),
					       by="cancerType"][order(medV,meanV),]
	    }
	    test.tb[,cancerType:=factor(as.character(cancerType),
				    levels=cancerType.order.tb$cancerType)]
	    test.tb[,groupV:=cancerType]
	}else{
	    test.tb[["V.split.by"]] <- test.tb[[split.by]]
	    test.tb$groupV <- sprintf("%s.%s",test.tb$V.split.by,test.tb$cancerType)
	    cancerType.order.tb <- test.tb[, .(medV=median(geneE),meanV=mean(geneE)),
					   by=c("V.split.by","cancerType")][order(V.split.by,medV,meanV),]
	    cancerType.order.tb[,groupV:=sprintf("%s.%s",V.split.by,cancerType)]
	    test.tb[,groupV:=factor(as.character(groupV),
				    levels=cancerType.order.tb$groupV)]
	}
	test.debug <<- test.tb

	p <- ggboxplot(test.tb,x="groupV",y="geneE",group="groupV",
		       fill="cancerType", xlab="",ylab="Frequency Of Positive Cells",
		       title=a.gene)
	if(!is.null(split.by)){
	    p <- p + facet_grid(~V.split.by, scales="free_x", space="free")
	    ktest.out <- test.tb[,{ res <- kruskal.test(geneE~cancerType,data=.SD)
				    .(p.value=res$p.value)
				    } ,by="V.split.by"]
	    y.pretty <- pretty(test.tb$geneE)
	    p <- p + geom_text(data = ktest.out,
			       aes(x =  2, y = y.pretty[length(y.pretty)-1],
				   label = sprintf("p=%4.4g",p.value)))
	}else{
	    p <- p + stat_compare_means(label="p.format")
	}
	p <- p + scale_x_discrete(labels = function(x){ gsub("^(.+?)\\.", "", x,perl=T) })
	p <- p + scale_fill_manual(values=colSet$cancerType) +
		geom_hline(yintercept=0.05,color="lightgray",linetype="dashed",alpha=0.8) +
		#facet_wrap(~cluster.name,ncol=1, scales="free_y") +
		theme(axis.text.x=element_text(angle=45,hjust=1),
		      strip.text=element_text(size=12),
		      plot.title=element_text(hjust=0.5),
		      ####aspect.ratio=my.aspect.ratio,
		      legend.position="none")
	###ggsave(sprintf("%s.gene.bin.exp.tb.%s.pdf",out.prefix,a.gene),width=5,height=3,useDingbats=F)
	ggsave(sprintf("%s.gene.bin.exp.tb.%s.pdf",out.prefix,a.gene),width=pdf.width,height=pdf.height,useDingbats=F)
	return(p)
    })
    names(p.list) <- gene.plot
    return(p.list)
}

#### boxplot (Fig. 2G, fig. S21D, fig. S22A, fig. S24C ???)
{

    out.prefix.10X <- sprintf("%s/10X.nonpara/%s",dirname(out.prefix),basename(out.prefix))
    dir.create(dirname(out.prefix.10X),F,T)
    p.10X.list <- make.boxplot.cmpCancerType(sce.binExp.all[,sce.binExp.all$tech=="10X"],
			       gene.list.plot,
			       #"SLC4A10",
			       sprintf("%s.byMedian.10X",out.prefix.10X),
			       min.NCell.mcls=10,sort.by="median")
    saveRDS(p.10X.list,file=sprintf("%s.byMedian.boxplot.10X.rds",out.prefix.10X))

#    p.SmartSeq2.list <- make.boxplot.cmpCancerType(sce.binExp.all[,sce.binExp.all$tech=="SmartSeq2"],
#			       gene.list.plot,
#			       sprintf("%s.byMedian.SmartSeq2",out.prefix.10X),
#			       min.NCell.mcls=10,sort.by="median",pdf.width=2.0)
#    saveRDS(p.SmartSeq2.list,file=sprintf("%s.byMedian.boxplot.SmartSeq2.rds",out.prefix.10X))

#    p.Tech2.list <- make.boxplot.cmpCancerType(sce.binExp.all[,sce.binExp.all$tech %in% c("10X","SmartSeq2")],
#			       gene.list.plot,
#			       ##"IL26", 
#			       split.by="tech",
#			       ###assay.name="z.score",
#			       sprintf("%s.byMedian.Tech2.test",out.prefix.10X),
#			       min.NCell.mcls=10,sort.by="median",pdf.width=6.0,pdf.height=3)
#    saveRDS(p.Tech2.list,file=sprintf("%s.byMedian.boxplot.Tech2.rds",out.prefix.10X))

}


##############################################

