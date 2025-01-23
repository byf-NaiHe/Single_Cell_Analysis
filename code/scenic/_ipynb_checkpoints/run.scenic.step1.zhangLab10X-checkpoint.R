#!/usr/bin/env Rscript

library("sscClust")
library("data.table")
library("plyr")
library("tictoc")
library("SCENIC")
library("SCopeLoomR")
library("ggplot2")
library("R.utils")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 20)

out.prefix <- "OUT.scenic/sce/scenic.zhangLab10X"
dir.create(dirname(out.prefix),F,T)

gene.exclude.file <- "../../data/external/exclude.gene.misc.misc.v3.withoutStress.RData"
env.misc <- loadToEnv(gene.exclude.file)

###
scenic.db.dir <- "../../data/external/scenic"
meta.tb.file <- "../../data/metaInfo/panC.freq.all.meta.tb.rds"
sce.xx.CD8 <- readRDS("../../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds")
sce.xx.CD4 <- readRDS("../../data/expression/CD4/integration/int.CD4.S35.sce.merged.rds")

sce.list <- list(
                 ### combined data of T cells from patients with solid tumors
                 "S"=readRDS("../../data/expression/zhangLab/S.10X.sce.rds"),
                 ### data of T cells from patients with BCL and MM
				 "BCL"=readRDS("../../data/expression/zhangLab/BCL.10X.sce.rds"),
				 "MM"=readRDS("../../data/expression/zhangLab/MM.10X.sce.rds"))

{
    x.col <- NULL
    rid <- NULL
    for(x in names(sce.list)){
        if(is.null(x.col)){
            x.col <- colnames(colData(sce.list[[x]]))
        }else{
            x.col <- intersect(x.col,colnames(colData(sce.list[[x]])))
        }
        if(is.null(rid)){
            rid <- rownames(sce.list[[x]])
        }else{
            rid <- intersect(rid,rownames(sce.list[[x]]))
        }
    }
    x.col <- x.col[!grepl("^SCT_snn",x.col)]

    for(x in names(sce.list)){
        colData(sce.list[[x]]) <- colData(sce.list[[x]])[,x.col]
        sce.list[[x]] <- sce.list[[x]][rid,]
    }

    sce <- do.call(cbind,sce.list)
}

meta.tb <- readRDS(meta.tb.file)
meta.tb <- meta.tb[dataset.tech=="zhangLab5P" ,]

all(meta.tb[["cellID"]] %in% colnames(sce))
sce <- sce[,meta.tb[["cellID"]]]
all(meta.tb[["cellID"]]==colnames(sce))
colData(sce) <- DataFrame(meta.tb)
colnames(sce) <- sce$cellID

miniClusterSize <- meta.tb[,.(miniCluster.size=.N),by=c("cancerType","dataset","dataset.tech",
									"ClusterID.harmony","meta.cluster","meta.cluster.coarse",
									"cancerType.old","dataset.old","stype",
									"dataSource","tech","tech.cate","pub",
									"miniCluster")]
##miniClusterSize[,cid:=sprintf("%s_%s_%s",dataset,stype,miniCluster)]
miniClusterSize[,cid:=sprintf("%s_%s",stype,miniCluster)]
setkey(miniClusterSize,"cid")

tic("ssc.average.cell")
sce.mini <- ssc.average.cell(sce,assay.name="norm_exprs",
							 ###column=c("dataset","stype","miniCluster"),
							 column=c("stype","miniCluster"),
							 ret.type="sce",do.parallel=T)
toc()

f.cell <- intersect(miniClusterSize$cid,colnames(sce.mini))
sce.mini <- sce.mini[,f.cell]
miniClusterSize <- miniClusterSize[f.cell,]
all(colnames(sce.mini)==miniClusterSize$cid)
colData(sce.mini) <- DataFrame(miniClusterSize,row.names=miniClusterSize$cid)

all(rownames(sce.mini)==rownames(sce))
rowData(sce.mini) <- rowData(sce)

saveRDS(sce.mini,file=sprintf("%s.sce.mini.rds",out.prefix))
#sce.mini <- readRDS(file=sprintf("%s.sce.mini.rds",out.prefix))

scenicOptions <- initializeScenic(org="hgnc",
								  dbDir=scenic.db.dir,
								  nCores=10)

#exprMat <- assay(sce,"counts")
#rownames(exprMat) <- rowData(sce)$display.name

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

sceToLoom <- function(obj,loom.file,rd.map,...)
{
	exprMat <- assay(obj)
	rownames(exprMat) <- rowData(obj)$display.name
	colnames(exprMat) <- gsub("^CD[48]_","",colnames(exprMat))
	genesKept <- geneFiltering(exprMat, scenicOptions,...)

	### filter out genes in blacklist
    f.gene.blackList.bool <- (rownames(exprMat) %in% env.misc$all.gene.ignore.df$geneSymbol) |
							grepl("^RP[LS]",rownames(exprMat),perl=T) |
							rownames(exprMat)=="MALAT1"
	f.gene.blackList <- unname(rownames(exprMat)[f.gene.blackList.bool])
	genesKept <- setdiff(genesKept,f.gene.blackList)
	cat(sprintf("final genes kept (remove those in black list): %d\n",length(genesKept)))

	##### check auc threshold
	nGenesDetectedPerCell <- colSums(exprMat[genesKept,] > 0)
	th.nGenesDetectedPerCell <- quantile(nGenesDetectedPerCell,c(0.01,0.05,0.1,0.5,1))

	p <- ggplot(data.table(nGenesDetectedPerCell=nGenesDetectedPerCell),aes(x=nGenesDetectedPerCell)) +
			geom_histogram(aes(y=..density..), colour=NA, fill="steelblue")+
			geom_density(color="orange") +
			geom_vline(xintercept=th.nGenesDetectedPerCell,
					   color="red", linetype="dashed", size=1) +
			theme_bw()
	ggsave(sprintf("%s.nGenesDetectedPerCell.pdf",gsub("loom","",loom.file)),width=5,height=3.5)
	#####
	##return(NULL)

	loom <- build_loom(loom.file, dgem=exprMat[genesKept,],
					   default.embedding = rd.map, default.embedding.name = "harmony.umap")
	for(x in setdiff(colnames(colData(obj)),c("miniCluster","cid"))){
		if(is.character(obj[[x]]) || is.factor(obj[[x]])){
			add_col_attr(loom=loom, key = x, value=as.character(obj[[x]]), as.annotation=T)
		}else{
			add_col_attr(loom=loom, key = x, value=obj[[x]], as.metric=T)
		}
	}
	close_loom(loom)
}

#### CD8
sce.dump <- sce.mini[,sce.mini$stype=="CD8"]
rd.map <- reducedDim(sce.xx.CD8,"harmony.umap")
all(sce.dump$miniCluster %in% rownames(rd.map))
rd.map <- rd.map[sce.dump$miniCluster,]
all(sce.dump$miniCluster==rownames(rd.map))

sceToLoom(sce.dump,
		  sprintf("%s.mini.CD8.loom",out.prefix),
		  rd.map,
		  minCountsPerGene=0.5,minSamples=3)

#### CD4
sce.dump <- sce.mini[,sce.mini$stype=="CD4"]
rd.map <- reducedDim(sce.xx.CD4,"harmony.umap")
all(sce.dump$miniCluster %in% rownames(rd.map))
rd.map <- rd.map[sce.dump$miniCluster,]
all(sce.dump$miniCluster==rownames(rd.map))

sceToLoom(sce.dump,
		  sprintf("%s.mini.CD4.loom",out.prefix),
		  rd.map,
		  minCountsPerGene=0.5,minSamples=3)





##saveRDS(sce.mini[,sce.mini$stype=="CD8"],file=sprintf("%s.sce.mini.CD8.rds",out.prefix))
##saveRDS(sce.mini[,sce.mini$stype=="CD4"],file=sprintf("%s.sce.mini.CD4.rds",out.prefix))


