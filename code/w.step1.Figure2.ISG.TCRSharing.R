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
source("./func.R")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

out.prefix <- "OUT_Fig2/ISG.TCRSharing/ISG.TCRSharing"
dir.create(dirname(out.prefix),F,T)

##### load data
{

    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
    startrac.file <- "../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds"

    in.dat <- readRDS(tcr.file)
    in.dat$majorCluster <- as.character(in.dat$meta.cluster)
    in.dat$clone.id <- in.dat$cloneID
    #### filter out patient.cluster with number of cell < 10
    ncell.patient.cluster <- sort(unclass(in.dat[,table(sprintf("%s.%s",patient,majorCluster))]))
    ##ncell.patient.cluster <- sort(unclass(in.dat[stype=="CD8", table(sprintf("%s.%s",patient,majorCluster))]))
    in.dat <- in.dat[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
    ####

    out <- readRDS(startrac.file)
    pIndex.sig.tran <- as.data.table(out@pIndex.sig.tran)

}

### pattern sorting
mat.patterSort <- function(dat.plot.mtx,dataset.x="panC")
{
    old.colname <- colnames(dat.plot.mtx)
    if(ncol(dat.plot.mtx)==3){
	colnames(dat.plot.mtx) <- c("DP","SP.geneA","SP.geneB")
	dat.plot.mtx.pattern <- as.data.table(dat.plot.mtx > 0)
	dat.plot.mtx.pattern$rid <- rownames(dat.plot.mtx)
	dat.plot.mtx.pattern$N.DP <- dat.plot.mtx[,"DP"]
	dat.plot.mtx.pattern$N.SP.geneA <- dat.plot.mtx[,"SP.geneA"]
	dat.plot.mtx.pattern$N.SP.geneB <- dat.plot.mtx[,"SP.geneB"]
	dat.plot.mtx.pattern <- dat.plot.mtx.pattern[order(-DP,-SP.geneA,-SP.geneB,-N.DP,-N.SP.geneA,-N.SP.geneB),]
	dat.plot.mtx.pattern[,Group:=sprintf("%s%s%s",as.integer(DP),as.integer(SP.geneA),as.integer(SP.geneB))]
	dat.plot.mtx <- dat.plot.mtx[dat.plot.mtx.pattern$rid,,drop=F]
	####
    }else if(ncol(dat.plot.mtx)==4){
	colnames(dat.plot.mtx) <- c("DP","SP.geneA","SP.geneB","DN")
	dat.plot.mtx.pattern <- as.data.table(dat.plot.mtx > 0)
	dat.plot.mtx.pattern$rid <- rownames(dat.plot.mtx)
	dat.plot.mtx.pattern$N.DP <- dat.plot.mtx[,"DP"]
	dat.plot.mtx.pattern$N.SP.geneA <- dat.plot.mtx[,"SP.geneA"]
	dat.plot.mtx.pattern$N.SP.geneB <- dat.plot.mtx[,"SP.geneB"]
	dat.plot.mtx.pattern$N.DN <- dat.plot.mtx[,"DN"]
	dat.plot.mtx.pattern <- dat.plot.mtx.pattern[order(-DP,-SP.geneA,-SP.geneB,-DN,
							   -N.DP,-N.SP.geneA,-N.SP.geneB,-N.DN),]
	dat.plot.mtx.pattern[,Group:=sprintf("%s%s%s%s",as.integer(DP),
					     as.integer(SP.geneA),as.integer(SP.geneB),
					     as.integer(DN))]
	dat.plot.mtx <- dat.plot.mtx[dat.plot.mtx.pattern$rid,,drop=F]
    }
    colnames(dat.plot.mtx) <- old.colname
    return(dat.plot.mtx)
}

### pattern: mutual exclusive or not ? (fig. S16B)
{
    
    obj.star <- new("Startrac",in.dat[stype=="CD8",],aid="panC")

    ana.pattern <- function(aid,obj.star,mcls.vec.id,out.prefix,gene.sig.slim.list,gene.use,
			    pdf.width.pattern=4.0,flt.LLR=F,plot.RepreClone=T,
			    min.cloneSize=2,bg.expandToOtherState=F)
    {

        dir.create(dirname(out.prefix),F,T)

        mcls.vec <- fetchMetaClusterID2CusterFullName()[mcls.vec.id]
        dat.block <- obj.star@clonotype.dist.cluster[,mcls.vec.id,drop=F]
        colnames(dat.block) <- fetchMetaClusterID2CusterFullName()[colnames(dat.block)]
        f.exp <- dat.block[,mcls.vec[1] ] >= min.cloneSize
        dat.block.exp <- dat.block[f.exp,,drop=F]
        if(bg.expandToOtherState){
            dat.block.exp.ext <- obj.star@clonotype.dist.cluster[rownames(dat.block.exp),,drop=F]
            f.crossState <- rowSums(dat.block.exp.ext >= 1) >= 2
            dat.block.crossState <- dat.block.exp.ext[f.crossState,,drop=F]
            dat.block.exp <- dat.block.exp[rownames(dat.block.crossState),,drop=F]
        }

        ### filtering by LLR
        if(flt.LLR){
            clone.info.tb <- in.dat[,.N,by=c("clone.id","cancerType","patient","dataset","dataset.old")]
            clone.LLR.tb <- ldply(seq_len(nrow(dat.block.exp)),function(i){
                    clone.id.i <- rownames(dat.block.exp)[i]
                    i.dataset.old <- clone.info.tb[clone.id==clone.id.i,][["dataset.old"]]
                    #dat.block.exp[clone.id.i,]
                    if(length(i.dataset.old)> 0){
                        my.calCloneLLR(cellData=in.dat[majorCluster %in% mcls.vec.id,],
                            prob.mat=prob.svm.mat.list[[ i.dataset.old ]],
                            cloneID.x=clone.id.i,verbose=F)
                    }
                    },.parallel=T)
            clone.LLR.tb <- as.data.table(clone.LLR.tb)
            clone.LLR.flt.tb <- clone.LLR.tb[(LLR > 1 & G.best==G.obs) | nStates==1,]
            clone.LLR.flt.tb <- clone.LLR.flt.tb[order(-LLR,-LL1),]
            dat.block.exp <- dat.block.exp[ clone.LLR.flt.tb$cloneID, ]
        }
        #####

        dat.block.exp <- mat.patterSort(dat.block.exp)
        dat.block.exp.clamp <- dat.block.exp
        dat.block.exp.clamp[ dat.block.exp.clamp > 9 ] <- 9

        ### numbers
        f.shared <- (dat.block.exp[,mcls.vec[1]] > 0) &
                (dat.block.exp[,mcls.vec[2]] > 0) &
                (dat.block.exp[,mcls.vec[3]] > 0)
        f.only1 <-  (dat.block.exp[,mcls.vec[1]] > 0) &
                (dat.block.exp[,mcls.vec[2]] > 0) &
                (dat.block.exp[,mcls.vec[3]]==0)
        f.only2 <-  (dat.block.exp[,mcls.vec[1]] > 0) &
                (dat.block.exp[,mcls.vec[2]] == 0) &
                (dat.block.exp[,mcls.vec[3]] > 0)
        f.null <-   (dat.block.exp[,mcls.vec[1]] > 0) &
                (dat.block.exp[,mcls.vec[2]] == 0) &
                (dat.block.exp[,mcls.vec[3]]==0)
        nclone.tb <- data.table(dataset="panC",
                    aid=aid,
                    nclone.shared=sum(f.shared),
                    nclone.only1=sum(f.only1),
                    nclone.only2=sum(f.only2),
                    nclone.null=sum(f.null))
        nclone.tb <- cbind(nclone.tb,as.data.table(ldply(seq_len(nrow(nclone.tb)),function(i){
                        fisher.out <- fisher.test(matrix(c(nclone.tb$nclone.shared[i],
                                           nclone.tb$nclone.only1[i],
                                           nclone.tb$nclone.only2[i],
                                           nclone.tb$nclone.null[i]),ncol=2))
                        out.tb <- data.table(OR=fisher.out$estimate,
                                 p.value=fisher.out$p.value)
                    })))
        write.table(nclone.tb,file=sprintf("%s.pattern.me.00.txt",out.prefix),
                sep="\t",quote=F,row.names=F)

        ann.row.vec <- rep(NA,nrow(dat.block.exp))
        ann.row.vec[f.shared] <- "Both"
        ann.row.vec[f.only1] <- mcls.vec[2]
        ann.row.vec[f.only2] <- mcls.vec[3]
        ann.row.vec[f.null] <- "NULL"
        ann.row <- rowAnnotation(cate=ann.row.vec,
                      col=list(cate=structure(RColorBrewer::brewer.pal(4,"Set2"),
                                       names=c("Both",mcls.vec[2],mcls.vec[3],"NULL"))),
                      annotation_legend_param=list(),
                      border=T)
        sscVis::plotMatrix.simple(dat.block.exp.clamp,mytitle="TCR sharing",
                        out.prefix=sprintf("%s.pattern.01",out.prefix),
                        col.ht=rev(structure(c("lightgray",
                                       sscVis:::getColorPaletteFromNameContinuous("YlOrRd")),
                                     names=0:9)),
                        #show.number=dat.plot.mtx,
                        par.heatmap = list(cex.row=0, left_annotation = ann.row),
                        my.cell_fun=function(j, i, x, y, w, h, col){
                            nn <- dat.block.exp[i,j]
                            if(nn>25){
                                n.fontsize <- 8
                            }else{
                                n.fontsize <- 10
                            }
                            if(nrow(dat.block.exp)<=10 && nn>0){
                                if(dat.block.exp[i,j]>=7){
                                grid::grid.text(dat.block.exp[i, j], x, y,gp=grid::gpar(col="white",
                                                               fontsize=n.fontsize))
                                }else{
                                grid::grid.text(dat.block.exp[i, j], x, y,gp=grid::gpar(fontsize=n.fontsize))
                                }
                            } },
                        exp.name="CellNumber",returnHT=F,
                        pdf.width=pdf.width.pattern,pdf.height=8)


        ### representative clones
        if(plot.RepreClone)
        {
            f.rep <- rowSums(dat.block.exp > 1)==3
            dat.block.exp[f.rep,]
            clone.rep.vec <- head(rownames(dat.block.exp[f.rep,]),n=10)

            sce.list <- list()
            ##OV−P20190304_C000163:228
            ##THCA−P20181226_C001600:2693
            ##UCEC−P20190213_C002471:107
            ##UCEC−P20190305_C002457:164
            ##UCEC−P20190312_C002534:59
            ##ESCA−P20190613_C000174:45
            ### STAD.P181019_C000323:71
            ### UCEC−P20190312_C002965:79
            l_ply(seq_along(clone.rep.vec),function(i){
                          vis.clonotype(out.prefix=sprintf("%s.%s", out.prefix, "rep"),
                                clone.x=clone.rep.vec[i],
                                sce.list=sce.list,
                                cinfo.clone.tb=in.dat,
                                gene.sig.list=gene.sig.slim.list,
                                ##sig.pretty=seq(-1,1,0.5),
                                gene.use=gene.use,
                                mcls.plot=mcls.vec.id)
                            })
        }
    }

    #### used in figure
    ana.pattern(aid="CD8.c02.Tm.IL7R",
		obj.star=obj.star,
		mcls.vec.id=c("CD8.c15.ISG.IFIT1","CD8.c02.Tm.IL7R", "CD8.c12.Tex.CXCL13"),
		gene.sig.slim.list=list("CD8.c02(IL7R+ Tm)"="IL7R","ZFP36L2","CXCR4","ZFP36","ANXA1","GPR183",
				     "CD8.c12(terminal Tex)"=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4"),
				     "CD8.c15(ISG+ CD8+ T)"="IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		gene.use=c("IL7R","ZFP36L2","CXCR4","ZFP36","ANXA1","GPR183",
			   "CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
			   "IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		flt.LLR=F,
		plot.RepreClone=F,
		min.cloneSize=2,
		bg.expandToOtherState=T,
		out.prefix=sprintf("%s/pattern.CD8.c02.Tm.IL7R.minCloneSize.expandToOtherState/%s",
                           dirname(out.prefix),basename(out.prefix)))

    ### used in figure
    ana.pattern(aid="CD8.c10.Trm.ZNF683",
		obj.star=obj.star,
		mcls.vec.id=c("CD8.c15.ISG.IFIT1","CD8.c10.Trm.ZNF683", "CD8.c12.Tex.CXCL13"),
		gene.sig.slim.list=list("CD8.c10(ZNF683+CXCR6+ Trm)"=c("ZNF683","HOPX","CAPG"),
				     "CD8.c12(terminal Tex)"=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4"),
				     "CD8.c15(ISG+ CD8+ T)"="IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		gene.use=c("ZNF683","HOPX","CAPG",
			   "CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
			   "IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		pdf.width.pattern=4.6,
		flt.LLR=F,
		plot.RepreClone=F,
		min.cloneSize=2,
		bg.expandToOtherState=T,
		out.prefix=sprintf("%s/pattern.CD8.c10.Trm.ZNF683.minCloneSize.expandToOtherState/%s",
                           dirname(out.prefix),basename(out.prefix)))

    #### used in figure
    ana.pattern(aid="CD8.c06.Tem.GZMK",
		obj.star=obj.star,
		mcls.vec.id=c("CD8.c15.ISG.IFIT1","CD8.c06.Tem.GZMK", "CD8.c12.Tex.CXCL13"),
		gene.sig.slim.list=list("CD8.c06(GZMK+ Tem)"=c("GZMK", "CD74", "CST7"),
				     "CD8.c12(terminal Tex)"=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4"),
				     "CD8.c15(ISG+ CD8+ T)"="IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		gene.use=c("GZMK", "CD74", "CST7",
			   "CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
			   "IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		pdf.width.pattern=4.1,
		flt.LLR=F,
		plot.RepreClone=F,
		min.cloneSize=2,
		bg.expandToOtherState=T,
		out.prefix=sprintf("%s/pattern.CD8.c06.Tem.GZMK.minCloneSize2.expandToOtherState/%s",
                           dirname(out.prefix),basename(out.prefix)))

    ### used in figure
    ana.pattern(aid="CD8.c11.Tex.PDCD1",
		obj.star=obj.star,
		mcls.vec.id=c("CD8.c15.ISG.IFIT1","CD8.c11.Tex.PDCD1", "CD8.c12.Tex.CXCL13"),
		gene.sig.slim.list=list("CD8.c11(GZMK+ Tex)"=c("GZMK", "CD74", "CST7"),
				     "CD8.c12(terminal Tex)"=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4"),
				     "CD8.c15(ISG+ CD8+ T)"="IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		gene.use=c("GZMK", "CD74", "CST7",
			   "CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
			   "IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
		pdf.width.pattern=4.1,
		flt.LLR=F,
		plot.RepreClone=F,
		min.cloneSize=2,
		bg.expandToOtherState=T,
		out.prefix=sprintf("%s/pattern.CD8.c11.Tex.PDCD1.minCloneSize.expandToOtherState/%s",
                           dirname(out.prefix),basename(out.prefix)))

}

################################


