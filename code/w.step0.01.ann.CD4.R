#!/usr/bin/env Rscript

library("Seurat")
library("sscVis")
library("sscClust")
library("scPip")
library("ggplot2")
library("ggpubr")
library("data.table")
library("tictoc")
library("plyr")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 8)

exp.list.file <- "list/sce.CD4.fullPath.list"
sce <- readRDS("OUT.int.CD4/int.CD4.sce.merged.rds")
seu <- readRDS("OUT.int.CD4/int.CD4.seu.merged.rds")
meta.tb <- readRDS("OUT.int.CD4/int.CD4.meta.tb.rds")

out.prefix <- "OUT.int.CD4/int.CD4"

############ information including cancerType,... are required ##########
all(colnames(sce)==colnames(seu))
sce$cancerType <- unname(sapply(strsplit(sce$dataset,"\\.",perl=T),"[",1))
seu$cancerType <- sce$cancerType
mapping.mini2CancerType <- structure(sce$cancerType,names=colnames(sce))
#########################################################################

g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

#g.colSet <- list()
#cc.name <- c("dataset","cancerType")
#for(cc in cc.name){
#    cc.values <- sort(unique(colData(sce)[,cc]))
#    g.colSet[[cc]] <- structure(sscVis:::auto.colSet(length(cc.values),name = "Paired"), names = as.character(cc.values))
#}

p <- ssc.plot.tsne(sce, columns = "cancerType", 
                  reduced.name = "harmony.umap",
                  colSet=g.colSet,size=0.01,label=NULL,
                  vector.friendly=T,
                  par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                  par.geom_point = list(scale=0.8),
                  par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                  base_aspect_ratio = 1.35)
ggsave(sprintf("%s.harmony.umap.groupBy.%s.pdf",out.prefix,"cancerType"),width=6.0,height=4)

p <- ssc.plot.tsne(sce, columns = "dataset", 
                  reduced.name = "harmony.umap",
                  colSet=g.colSet,size=0.01,label=NULL,
                  vector.friendly=T,
                  par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                  par.geom_point = list(scale=0.8),
                  par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                  base_aspect_ratio = 1.35)
ggsave(sprintf("%s.harmony.umap.groupBy.%s.pdf",out.prefix,"dataset"),width=8.5,height=4)

p <- ssc.plot.tsne(sce, columns = "dataset", 
                   splitBy="dataset",
                  reduced.name = "harmony.umap",
                  colSet=g.colSet,size=0.01,label=NULL,
                  vector.friendly=T,
                  par.geom_point = list(scale=1),
                  par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                  base_aspect_ratio = 1.35)
ggsave(sprintf("%s.harmony.umap.groupBy.%s.split.pdf",out.prefix,"dataset"),width=14,height=8)


opt.res <- "RNA_snn_res.2.6"

sce$ClusterID <- sprintf("C%02d",as.integer(as.character(colData(sce)[[opt.res]])))
sce$SubClusterID <- sce$ClusterID


######

as.data.table(colData(sce))[,.N,by=c("ClusterID","SubClusterID")][order(ClusterID),]

{
    sce$meta.cluster <- ""

    ######### manually annotate
    sce$ClusterID <- as.character(sce$ClusterID)
    sce$meta.cluster <- (sce$ClusterID)

	sce$meta.cluster[sce$ClusterID %in% c("C23","C05","C07","C29","C09","C17")] <- "CD4.c01.Tn.TCF7"
	sce$meta.cluster[sce$ClusterID %in% c("C30")] <- "CD4.c02.Tn.PASK"
	sce$meta.cluster[sce$ClusterID %in% c("C32")] <- "CD4.c03.Tn.ADSL"
	sce$meta.cluster[sce$ClusterID %in% c("C33")] <- "CD4.c04.Tn.il7r"
	sce$meta.cluster[sce$ClusterID %in% c("C35")] <- "CD4.c05.Tm.TNF"
	sce$meta.cluster[sce$ClusterID %in% c("C02","C00","C19")] <- "CD4.c06.Tm.ANXA1"
	sce$meta.cluster[sce$ClusterID %in% c("C11")] <- "CD4.c07.Tm.ANXA2"
	sce$meta.cluster[sce$ClusterID %in% c("C25")] <- "CD4.c08.Tm.CREM"
	sce$meta.cluster[sce$ClusterID %in% c("C21")] <- "CD4.c09.Tm.CCL5"
	sce$meta.cluster[sce$ClusterID %in% c("C08")] <- "CD4.c10.Tm.CAPG"
	sce$meta.cluster[sce$ClusterID %in% c("C22")] <- "CD4.c11.Tm.GZMA"
	sce$meta.cluster[sce$ClusterID %in% c("C13")] <- "CD4.c12.Tem.GZMK"
	sce$meta.cluster[sce$ClusterID %in% c("C12")] <- "CD4.c13.Temra.CX3CR1"
	### Th17
	sce$meta.cluster[sce$ClusterID %in% c("C18")] <- "CD4.c14.Th17.SLC4A10"
	sce$meta.cluster[sce$ClusterID %in% c("C15")] <- "CD4.c15.Th17.IL23R"
	### Tfh like
	sce$meta.cluster[sce$ClusterID %in% c("C24","C06","C34")] <- "CD4.c16.Tfh.CXCR5"
	sce$meta.cluster[sce$ClusterID %in% c("C14")] <- "CD4.c17.TfhTh1.CXCL13"
	### Treg 
	sce$meta.cluster[sce$ClusterID %in% c("C16")] <- "CD4.c18.Treg.RTKN2"
	sce$meta.cluster[sce$ClusterID %in% c("C20")] <- "CD4.c19.Treg.S1PR1"
	sce$meta.cluster[sce$ClusterID %in% c("C10","C01","C03","C04")] <- "CD4.c20.Treg.TNFRSF9"
	sce$meta.cluster[sce$ClusterID %in% c("C31")] <- "CD4.c21.Treg.OAS1"
	sce$meta.cluster[sce$ClusterID %in% c("C27")] <- "CD4.c22.ISG.IFIT1"
	sce$meta.cluster[sce$ClusterID %in% c("C28")] <- "CD4.c23.Mix.NME1"
	sce$meta.cluster[sce$ClusterID %in% c("C26")] <- "CD4.c24.Mix.NME2"
	###

    sce$meta.cluster <- factor(sce$meta.cluster,
								      levels=unique(sort(sce$meta.cluster)))
    cluster.set <- unique(as.data.table(colData(sce)[,c("meta.cluster","ClusterID","SubClusterID")
					     ]))[order(meta.cluster,ClusterID),]
    sce$ClusterID <- factor(sce$ClusterID, levels=cluster.set$ClusterID)
    ####
    all(colnames(seu)==colnames(sce))
    seu$meta.cluster <- sce$meta.cluster
    seu$ClusterID <- sce$ClusterID

}

################ cluster annotation ###################
{

    #cc.name <- c("meta.cluster")
    #for(cc in cc.name){
    #    cc.values <- sort(unique(colData(sce)[,cc]))
    #    g.colSet[[cc]] <- structure(sscVis:::auto.colSet(length(cc.values),name = "Paired"), names = as.character(cc.values))
    #}

    p <- ssc.plot.tsne(sce, columns = "ClusterID", 
                      reduced.name = "harmony.umap",
                      colSet=list(),size=0.5,label=2,
                      vector.friendly=T,
                      par.geom_point = list(scale=0.8),
                      par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                      base_aspect_ratio = 1.35)
    ggsave(sprintf("%s.harmony.umap.groupBy.%s.pdf",out.prefix,"ClusterID"),width=6,height=4)

    p <- ssc.plot.tsne(sce, columns = "meta.cluster", 
                      reduced.name = "harmony.umap",
                      colSet=g.colSet,size=0.5,label=2,
                      par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                      vector.friendly=T,
                      #
                      fun.extra=function(p){ p + guides(color=guide_legend(ncol=2, override.aes=list(size=4))) },
                      par.geom_point = list(scale=0.8),
                      par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                      base_aspect_ratio = 1.35)
    ggsave(sprintf("%s.harmony.umap.groupBy.%s.pdf",out.prefix,"meta.cluster"),width=8.0,height=4)

}


################ save ###################
{

    all(colnames(sce)==colnames(seu))
    seu$ClusterID <- sce$ClusterID
    ##seu$SubClusterID <- sce$SubClusterID
    seu$meta.cluster <- sce$meta.cluster

    mapping.mini2majorCluster <- structure(sce$meta.cluster,names=colnames(sce))
    meta.tb[,meta.cluster:=mapping.mini2majorCluster[miniCluster]]
    meta.tb[,cancerType:=mapping.mini2CancerType[miniCluster]]
    idx.to.add <- meta.tb[,.(DIG.Score1=mean(DIG.Score1),
                             score.MALAT1=mean(score.MALAT1),
                             percent.mito=mean(percent.mito)),
                            by="miniCluster"]
    setkey(idx.to.add,"miniCluster")
    colData(sce) <- cbind(colData(sce),as.data.frame(idx.to.add[colnames(sce),-1]))

    saveRDS(sce,file=sprintf("%s.sce.merged.rds",out.prefix))
    saveRDS(seu,file=sprintf("%s.seu.merged.rds",out.prefix))
    saveRDS(meta.tb,file=sprintf("%s.meta.tb.rds",out.prefix))
    saveRDS(g.colSet,file=sprintf("%s.colSet.rds",out.prefix))
    #sce <- readRDS(file=sprintf("%s.sce.merged.rds",out.prefix))
    #seu <- readRDS(file=sprintf("%s.seu.merged.rds",out.prefix))
    #meta.tb <- readRDS(file=sprintf("%s.meta.tb.rds",out.prefix))
    meta.tb[,.N,by=c("meta.cluster","dataset")][order(meta.cluster),]
    meta.tb[,.N,by=c("cancerType","dataset")][order(cancerType),]

    exp.list.table <- fread(exp.list.file)

    loginfo("save sce objects by dataset ...")
    cn.tb <- scPip:::saveSCEPerDataSet(exp.list.table,meta.tb,out.prefix)
    loginfo("done.")

}






