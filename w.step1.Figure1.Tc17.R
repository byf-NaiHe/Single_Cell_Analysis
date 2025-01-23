#!/usr/bin/env Rscript

library("data.table")
library("sscVis")
library("plyr")
library("ggplot2")
library("ggpubr")
library("ggplotify")
library("cowplot")
library("tictoc")
source("./func.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 8)

out.prefix <- "OUT_Fig1/Tc17/Tc17"
dir.create(dirname(out.prefix),F,T)

### load data
{

    dat.tcr <- readRDS("../data/tcr/byCell/tcr.zhangLab.comb.flt.rds")
    colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    cellInfo.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")
    ######## gene expression data ####
    sce.list.file <- "list/sce.CD8.fullPath.list"
    sce.list.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.file),head=T)
    sce.list <- llply(seq_len(nrow(sce.list.tb)),function(i){ readRDS(sce.list.tb$scefile[i]) },.parallel=T)
    names(sce.list) <- sce.list.tb$data.id

}

##############

calZExp <- function(sce.list,
		      #cellIDToRID.tb,
		       gene.to.test=c("FOXP3","CXCL13"),locations=c("P","N","T"),
		       dataset.x="OV.zhangLab5P",
		       mcls="CD8.c12.Tex.CXCL13",mcls.focus="CD8.c16.MAIT.SLC4A10",
		       th.exprs=0.3, min.ncell=50,val.padding=0)
{

    sce <- sce.list[[dataset.x]]
    
    if(mcls=="all"){
        if(sum(sce[["meta.cluster"]]==mcls.focus) < min.ncell){
            return(NULL)
        }
    }else{
        #### filter out dataset with cell number in Tex less than 50
        if(sum(sce[["meta.cluster"]]==mcls) < min.ncell){
            return(NULL)
        }
    }

    exp.mat <- assay(sce,"norm_exprs")
    rnames <- rowData(sce)$display.name
    rownames(exp.mat) <- unname(rnames)
    colnames(exp.mat) <- sce$cellID.uniq
    sce.plot <- ssc.build(exp.mat,assay.name="norm_exprs")
    cdata <- colData(sce)[,c("patient","cellID","libraryID","cancerType","loc",
				"batchV","dataset","cellID.uniq","miniCluster","meta.cluster")]
    rownames(cdata) <- sce$cellID.uniq
    colData(sce.plot) <- cdata
    sce.plot <- ssc.scale(sce.plot,gene.symbol=gene.to.test,adjB="batchV",do.scale=T,add.padding=T)
    if(mcls=="all"){
        sce.plot <- sce.plot[,sce.plot$loc %in% locations]
    }else{
        sce.plot <- sce.plot[,sce.plot$meta.cluster==mcls & sce.plot$loc %in% locations]
    }
    assay(sce.plot,"binExp") <- assay(sce.plot,"norm_exprs.scale") > th.exprs
    sce.plot$meta.cluster <- as.character(sce.plot$meta.cluster)
    sce.plot <- sce.plot[intersect(gene.to.test,rownames(sce.plot)),]
    cat(sprintf("dataset done (%s)\n",dataset.x))
    return(sce.plot)
}

gene.in.box <- c("IL26","IL17A","IL17F","RORC","IL23R","KLRB1","CTSH","CCR6","TMIGD2","IL4I1","ZBTB16","SLC4A10",
		 "NCR3", "LTK", "CA2", "ME1", "AQP3", "ADAM12",
		 "CXCR6", "RUNX2", "ID2", "CA10", "CAPG", "PTPN13", "MGAT4A", "ERN1", "ANKRD28")

dataset.ana.vec <- names(sce.list)

### 
tic("calZ")
sce.ZExp.list <- llply(dataset.ana.vec,function(x){
			     calZExp(sce.list,gene.to.test=gene.in.box,
				     mcls="all",
				     dataset.x=x)
		  })
names(sce.ZExp.list) <- dataset.ana.vec
toc()

sce.ZExp.all <- NULL
for(x in names(sce.ZExp.list)){
    if(!is.null(sce.ZExp.list[[x]])){
        sce.tmp <- sce.ZExp.list[[x]]
        metadata(sce.tmp)$ssc <- NULL
        if(is.null(sce.ZExp.all)){
            sce.ZExp.all <- sce.tmp
        }else{
            sce.ZExp.all <- cbind(sce.ZExp.all,sce.tmp)
        }
        cat(sprintf("dataset done (%s)\n",x))
    }
}

dat.tcr[,cellID.uniq:=sprintf("%s.%s",dataset.old,Cell_Name)]
f.cell <- intersect(dat.tcr$cellID.uniq,colnames(sce.ZExp.all))
setkey(dat.tcr,"cellID.uniq")
sce.ZExp.all <- sce.ZExp.all[,f.cell]
colData(sce.ZExp.all) <- cbind(colData(sce.ZExp.all),dat.tcr[f.cell,c("cloneID","cloneSize","cluster.name","invariantA")])

sce.ZExp.all$invariantA.old <- sce.ZExp.all$invariantA

f.mait <- sce.ZExp.all$invariantA.old=="MAIT" & sce.ZExp.all$meta.cluster=="CD8.c16.MAIT.SLC4A10"
sce.ZExp.all$invariantA[f.mait] <- "MAIT"
f.Tc17 <- sce.ZExp.all$invariantA.old=="" & sce.ZExp.all$meta.cluster=="CD8.c16.MAIT.SLC4A10"
sce.ZExp.all$invariantA[f.Tc17] <- "Tc17"
f.oo <- sce.ZExp.all$meta.cluster!="CD8.c16.MAIT.SLC4A10"
sce.ZExp.all$invariantA[f.oo] <- "other"
table(sce.ZExp.all$invariantA)

# iNKT  MAIT other  Tc17
#    1  1411 90763  1579

saveRDS(sce.ZExp.all,file=sprintf("%s.sce.rds",out.prefix))
#sce.ZExp.all <- readRDS(file=sprintf("%s.sce.rds",out.prefix))

sce.ZExp.miniC <- ssc.average.cell(sce.ZExp.all,assay.name="norm_exprs.scale",
                                     column=c("invariantA","cancerType","miniCluster"),
                                     ret.type="sce",do.parallel=T)
saveRDS(sce.ZExp.miniC,file=sprintf("%s.sce.miniC.rds",out.prefix))
#sce.ZExp.miniC <- readRDS(file=sprintf("%s.sce.miniC.rds",out.prefix))
table(sce.ZExp.miniC$invariantA)
# iNKT  MAIT other  Tc17
#    1   145  4982   142

dat.plot.tb <- as.data.table(colData(sce.ZExp.all))
col.cinfo <- colnames(dat.plot.tb)
dat.plot.tb <- cbind(dat.plot.tb,t(assay(sce.ZExp.all,"norm_exprs.scale")))
dat.plot.tb[,.N,by=c("invariantA","miniCluster")][order(-N),]
N.tb <- dcast(dat.plot.tb[invariantA!="iNKT",.N,by=c("invariantA","cancerType")][order(-N),],
	      cancerType~invariantA,value.var="N")

#    cancerType MAIT Tc17 other
# 1:         BC   60   49  3245
# 2:        BCL  101   48  2592
# 3:        CRC   46   78  3038
# 4:       ESCA   82  132  8356
# 5:        HCC  153   21  1135
# 6:       LUNG   76   41  3376
# 7:         MM  137   21  7146
# 8:         OV    6   71  2925
# 9:       PACA   91  251  4385
#10:         RC   87   86 12927
#11:       STAD   13   43  1775
#12:       THCA  294  385 24582
#13:       UCEC  265  353 15281


a.clamp <- c(-0.5,1.5)

### expression of marker genes (fig. S6B)
{

    b.clamp <- c(-0.5,4)

    sce.plot.tmp <- sce.ZExp.miniC[,sce.ZExp.miniC$invariantA!="iNKT"]
    sce.plot.tmp$invariantA[sce.plot.tmp$invariantA=="other"] <- "Other"
    sce.plot.tmp$invariantA <- factor(sce.plot.tmp$invariantA,levels=c("MAIT","Tc17","Other"))

    p <- ssc.plot.violin(sce.plot.tmp,
			 assay.name="norm_exprs.scale",
			 #splitBy="cancerType",
			 #gene=gene.in.box,
			 gene=c("RORC","KLRB1","IL23R","CCR6","ZBTB16","ME1"),
			 group.var="invariantA",
			 p.ncol=6,
			 add="boxplot",
			 clamp=b.clamp) +
		    xlab("") + ylab("z-score") +
		    coord_cartesian(ylim=b.clamp) +
		    geom_hline(yintercept=0.3,linetype="dashed") +
		    #theme_pubr() +
		    theme(
			  strip.text=element_text(size=12),
			  axis.text=element_text(size=12),
			  legend.position="right",
			  panel.grid.minor=element_blank())
    ggsave(sprintf("%s.violin.example.splitCancertype.FALSE.miniC.fig.panC.pdf",out.prefix),width=13,height=2.5)

}

#################


#################

