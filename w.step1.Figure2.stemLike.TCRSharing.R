#.!/usr/bin/env Rscript

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
source("./func.R")
#source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/inte.comb.miniClust.lib.R")
#source("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/PanC.T/lib/TCR.ana.R")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
startrac.file <- "../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds"
out.prefix <- "OUT_Fig2/stemLike/CD8.stemLike.TCRSharing"

dir.create(dirname(out.prefix),F,T)

###  tcr data
{
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

    out <- readRDS(startrac.file)
    pIndex.sig.tran <- as.data.table(out@pIndex.sig.tran)
}



### sce
sce.list.file <- "list/sce.CD8.fullPath.list"
sce.list.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.file),head=T)
sce.list <- llply(seq_len(nrow(sce.list.tb)),function(i){ readRDS(sce.list.tb$scefile[i]) },.parallel=T)
names(sce.list) <- sce.list.tb$data.id

### probability matrix
{
    ds.tb <- in.dat[,.N,by=c("dataset","dataset.old")]
    prob.svm.mat.list <- llply(ds.tb$dataset.old,function(x){
	  prob.svm.mat <- readRDS(sprintf("../data/expression/CD8/integration/reClassify.v1/T.CD8.%s/T.CD8.prob.mat.svm.%s.rds",x,x))
	})
    names(prob.svm.mat.list) <- ds.tb$dataset
}


### clonotype level pIndex.tran
{

    obj.star <- new("Startrac",in.dat[stype=="CD8",],aid="panC")

    {

	### 
	{
	    object <- obj.star

	    mcls.1 <- "CD8.c14.Tex.TCF7"
	    mcls.2 <- c("CD8.c17.Tm.NME1")
	    ### pTrans per clonotype
	    clone.hi.pIndex.full.tb <- as.data.table(ldply(mcls.2,function(x){
            dat.block <- object@clonotype.dist.cluster[,c(x,mcls.1)]
            dat.block.clone.index <- Startrac:::mrow.entropy(dat.block)
            dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
            ##ff <- (dat.block.clone.index > 0.1)
            out.tb <- data.table(meta.cluster=x,
                         cloneID=names(dat.block.clone.index),
                         pIndex.tran=dat.block.clone.index,
                         ncells=rowSums(dat.block),
                         nTex=object@clonotype.dist.cluster[,mcls.1],
                         nO=object@clonotype.dist.cluster[,x],
                         nAll=rowSums(object@clonotype.dist.cluster))
            return(out.tb)
	    }))
	    ###clone.hi.pIndex.tb <- clone.hi.pIndex.full.tb[pIndex.tran > 0.1]
	    
	    clone.info.patient.tb <- in.dat[,.N,by=c("cancerType","dataset","dataset.old","patient","cloneID")]

	    ### LLR of clonotypes
	    {
            clone.hi.pIndex.tb <- clone.hi.pIndex.full.tb[pIndex.tran > 0]
            clone.hi.pIndex.tb <- clone.hi.pIndex.tb[order(-pIndex.tran),]
            clone.hi.pIndex.tb <- merge(clone.hi.pIndex.tb,clone.info.patient.tb,by="cloneID")
            clone.LLR.tb <- ldply(seq_len(nrow(clone.hi.pIndex.tb)),function(i){
                        calCloneLLR(in.dat[majorCluster %in% c("CD8.c14.Tex.TCF7","CD8.c17.Tm.NME1"),],
                        prob.mat=prob.svm.mat.list[[ clone.hi.pIndex.tb$dataset.old[i] ]],
                        cloneID.x=clone.hi.pIndex.tb$cloneID[i],verbose=F)
                        })
            clone.hi.pIndex.tb <- merge(clone.hi.pIndex.tb,clone.LLR.tb,by="cloneID")
            clone.hi.pIndex.tb <- clone.hi.pIndex.tb[order(-LLR,pIndex.tran),]
            clone.hi.pIndex.flt.tb <- clone.hi.pIndex.tb[LLR > 1 & G.best==G.obs,]
            saveRDS(clone.hi.pIndex.flt.tb,sprintf("%s.clone.hi.pIndex.flt.tb.rds", out.prefix))
            #clone.hi.pIndex.flt.tb <- readRDS(sprintf("%s.clone.hi.pIndex.flt.tb.rds", out.prefix))
            
            cinfo.clone.tb <- in.dat[stype=="CD8",][cloneID %in% clone.hi.pIndex.flt.tb$cloneID &
                                meta.cluster %in% c(mcls.1,mcls.2),]

	    }

	}

	### bar plots (fig. S18B ??)
    {
        tmp.colSet <- g.colSet["meta.cluster"]
        tmp.colSet$meta.cluster <- tmp.colSet$meta.cluster[c("CD8.c17.Tm.NME1","CD8.c14.Tex.TCF7")]
        out.Path.TCF7.NME1 <- ana.clonotypeAcrossMcls.moreThanTwo(obj.star,
                in.dat,
                out.prefix=sprintf("%s/Link.TCF7.NME1/%s",dirname(out.prefix),basename(out.prefix)),
                aid="Link.TCF7.NME1",
                par.barplot=list(pdf.width.byPatientF = 7, pdf.height.byPatientF = 5.0,
                         pdf.width.byPatientT = 7, pdf.height.byPatientT = 5.5),
                l.colSet=tmp.colSet,
                mcls.moi=c("CD8.c17.Tm.NME1","CD8.c14.Tex.TCF7"))

        tmp.colSet <- g.colSet["meta.cluster"]
        tmp.colSet$meta.cluster <- tmp.colSet$meta.cluster[c("CD8.c14.Tex.TCF7","CD8.c11.Tex.PDCD1")]
        out.Path.TCF7.NME1 <- ana.clonotypeAcrossMcls.moreThanTwo(obj.star,
                in.dat,
                out.prefix=sprintf("%s/Link.TCF7.GZMKTex/%s",dirname(out.prefix),basename(out.prefix)),
                aid="Link.TCF7.GZMKTex",
                par.barplot=list(pdf.width.byPatientF = 7, pdf.height.byPatientF = 5.0,
                         pdf.width.byPatientT = 7, pdf.height.byPatientT = 5.5),
                l.colSet=tmp.colSet,
                mcls.moi=c("CD8.c14.Tex.TCF7","CD8.c11.Tex.PDCD1"))
    }

	#### specific clonotype (fig. S18C ??)
	clone.example.tb <- clone.hi.pIndex.flt.tb[cloneID %in% c("UCEC-P20190305_C002457:164"),]
	cinfo.clone.tb <- in.dat[stype=="CD8",][cloneID %in% clone.example.tb$cloneID,]
	l_ply(seq_len(nrow(clone.example.tb)),function(i){
	    dataset.id <- clone.example.tb[i,][["dataset"]][1]
	    vis.clonotype(out.prefix=sprintf("%s.manualPick.%s.v4",
					     out.prefix,
					     clone.example.tb$meta.cluster[i]),
			  clone.x=clone.example.tb$cloneID[i],
			  sce.list=sce.list,cinfo.clone.tb=cinfo.clone.tb,
              gene.sig.list=list(),
			  gene.use=c("TCF7","CD27","CD28"),
			  prob.mat=prob.svm.mat.list[[dataset.id]],sortByProb=T,
			  pdf.width=18, pdf.height=3.5,
			  sig.pretty=seq(-1,1,0.5),
			  mcls.plot=c("CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13","CD8.c14.Tex.TCF7", "CD8.c17.Tm.NME1"))
				 })


    }
}



