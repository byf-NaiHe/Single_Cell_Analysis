#!/usr/bin/env Rscript

library("scPip")
library("sscClust")
library("data.table")
library("tictoc")
library("Seurat")
library("R.utils")
library("plyr")
library("dplyr")
library("tibble")
library("doParallel")
library("Matrix")
library("gplots")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("limma")
library("reticulate")
source("./func.R")

dat.ext.dir <- system.file("extdata",package="scPip")
gene.exclude.file <- sprintf("%s/exclude.gene.misc.misc.v3.RData",dat.ext.dir)

env.misc <- loadToEnv(gene.exclude.file)
env.misc$all.gene.ignore.df %>% head

#### run Seurat for OV data
{
    sce <- readRDS("../data/expression/CD8/byDataset/OV.thisStudy.sce.rds")
    out.prefix <- "./OUT_Fig2/OV/CD8.OV"
    dir.create(dirname(out.prefix),F,T)

    sce.full <- sce
    sce <- sce[,sce$meta.cluster=="CD8.c12.Tex.CXCL13"]

    tic("run.Seurat3")
    obj.list <- run.Seurat3(seu=NULL,sce=sce,out.prefix,
                            gene.exclude.df=env.misc$all.gene.ignore.df,n.top=1500,
                            measurement="counts",platform="10X",use.sctransform=T,
                            opt.res=0.8, opt.npc=15,ncores=12,run.stage=100)
    toc()
}

#### original result
sce.plot <- readRDS("../data/expression/CD8/integration/CD8.OV.sce.rds")
#### from this run
#### sce.plot <- readRDS(sprintf("%s.sce.rds",out.prefix))

### some genes (OV data) (Fig. 2E ???)
{
	gene.to.plot <- c("HAVCR2","CTLA4","CXCL13", "FOXP3","TXK","KIR2DL3")
	sce.tmp <- ssc.scale(sce.full,gene.symbol=gene.to.plot, adjB=NULL,do.scale=T)
	sce.tmp <- sce.tmp[,colnames(sce.plot)]
	reducedDim(sce.tmp,"seurat.umap") <- reducedDim(sce.plot,"seurat.umap")
	p <- ssc.plot.tsne(sce.tmp,assay.name="norm_exprs.scale",
						   gene=gene.to.plot,
						   reduced.name="seurat.umap", 
						   vector.friendly=T,size=3,
						   #par.geom_point=list(raster.width=6,raster.height=6),
						   clamp=c(-0.5,1.5),
						   p.ncol=3,
						   palette.name = "YlOrRd",
						   theme.use=theme_void,
						   ##par.geneOnTSNE=list(scales="fixed",pt.order="random",pt.alpha = 0.5))
						   par.geneOnTSNE=list(scales="fixed",pt.order="value",pt.alpha = 1))
	ggsave(sprintf("%s.umap.gene.Treg.man.01.pdf",out.prefix),width=7.5,height=4.5)
	#saveRDS(sce.tmp,file=sprintf("%s.sce.for.umap.gene.Treg.man.01.rds",out.prefix))

}


sce.merged.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds"),
                        "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.sce.merged.rds"))
g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

### some genes in panC data, used as a reference
{

    sce.panC <- sce.merged.list$CD8
	p <- ssc.plot.tsne(sce.panC[,sce.panC$meta.cluster %in% c("CD8.c12.Tex.CXCL13","CD8.c09.Tk.KIR2DL4")],
					   assay.name="exprs",
					   gene=gene.to.plot,
					   reduced.name="harmony.umap", 
					   vector.friendly=T,size=0.4,
					   #par.geom_point=list(raster.width=6,raster.height=6),
					   clamp=c(-0.5,1.5),
					   p.ncol=3,
					   palette.name = "YlOrRd",
					   theme.use=theme_void,
					   par.geneOnTSNE=list(scales="fixed",pt.order="random",pt.alpha = 0.5))
	ggsave(sprintf("%s.umap.gene.Treg.man.panC.gene.pdf",out.prefix),width=7.5,height=4.5)

    l.colSet <- g.colSet
    l.colSet$meta.cluster <- l.colSet$meta.cluster[c("CD8.c12.Tex.CXCL13","CD8.c09.Tk.KIR2DL4")]
	p <- ssc.plot.tsne(sce.panC[,sce.panC$meta.cluster %in% c("CD8.c12.Tex.CXCL13","CD8.c09.Tk.KIR2DL4")],
					   assay.name="exprs",
					   columns="meta.cluster",
					   reduced.name="harmony.umap", 
					   vector.friendly=T, size=0.1,
					   #par.geom_point=list(raster.width=6,raster.height=6),
					   #par.geom_point=list(scale=1),
					   #clamp=c(-0.5,1.5),
					   p.ncol=1,
					   colSet=l.colSet,label=3,
					   #palette.name = "YlOrRd",
					   theme.use=theme_void,
					   par.geneOnTSNE=list(scales="fixed",pt.order="random",pt.alpha = 0.5))
		#theme(legend.position="none")
	ggsave(sprintf("%s.umap.gene.Treg.man.panC.meta.cluster.pdf",out.prefix),width=4.0,height=2.5,useDingbats=F)

}

### ectopic expression
{
    sce.CD8 <- sce.merged.list$CD8
    sce.CD8$cluster.name <- fetchMetaClusterID2CusterFullName()[as.character(sce.CD8$meta.cluster)]
    sce.CD8$cluster.id <- gsub("^CD8\\.","",fetchMetaClusterID2CusterFullName("cluster.id")[as.character(sce.CD8$meta.cluster)])
    sce.CD8$cluster.name.short <- fetchMetaClusterID2CusterFullName("cluster.name")[as.character(sce.CD8$meta.cluster)]

    sce.CD4 <- sce.merged.list$CD4
    sce.CD4$cluster.name <- fetchMetaClusterID2CusterFullName()[as.character(sce.CD4$meta.cluster)]
    sce.CD4$cluster.id <- gsub("^CD4\\.","",fetchMetaClusterID2CusterFullName("cluster.id")[as.character(sce.CD4$meta.cluster)])
    sce.CD4$cluster.name.short <- fetchMetaClusterID2CusterFullName("cluster.name")[as.character(sce.CD4$meta.cluster)]

    ####KIR (TXK, KIR2DL3) (fig. S21A ??)
    {
        sce.tmp <- sce.CD8
        sce.tmp$group <- sce.tmp$cluster.name.short
        sce.tmp$group[!(sce.tmp$cluster.name.short %in% c("terminal Tex","KIR+TXK+ NK-like"))] <- "Other"
        sce.tmp$group <- factor(sce.tmp$group,levels=c("terminal Tex","KIR+TXK+ NK-like","Other"))
        table(sce.tmp$group)

        b.clamp <- c(-0.5,2.5)
        gene.plot <- c("HAVCR2","CTLA4","CXCL13","TXK","KIR2DL3")
        p <- ssc.plot.violin(sce.tmp,gene=gene.plot,group.var="group",
                 clamp=b.clamp,
                 add="boxplot") +
            xlab("") + ylab("z-score") +
            coord_cartesian(ylim=b.clamp) +
            geom_hline(yintercept=0.3,linetype="dashed") +
            #theme_pubr() +
            theme(
              #axis.ticks.x=element_blank(),
              axis.text=element_text(size=12),
              strip.text=element_text(size=12),
              legend.position="right",
              panel.grid.minor=element_blank())
        ggsave(sprintf("%s.CD8.KIR.TXK.NKLike.00.pdf",out.prefix),width=3,height=7)

    }

    #### Treg(FOXP3) (fig. S21A ??)
    {

        sce.tmp.1 <- sce.CD8
        sce.tmp.2 <- sce.CD4[,sce.CD4$meta.cluster %in% c("CD4.c18.Treg.RTKN2","CD4.c19.Treg.S1PR1",
                                                          "CD4.c20.Treg.TNFRSF9","CD4.c21.Treg.OAS1")]
        rowData(sce.tmp.1) <- NULL
        rowData(sce.tmp.2) <- NULL
        sce.tmp <- cbind(sce.tmp.1,sce.tmp.2)
        all(rownames(sce.tmp)==rownames(sce.CD8))
        rowData(sce.tmp) <- rowData(sce.CD8)
        
        sce.tmp$meta.cluster <- as.character(sce.tmp$meta.cluster)
        sce.tmp$group <- "Other"
        sce.tmp$group[grepl("^CD4",sce.tmp$meta.cluster,perl=T)] <- "CD4+ Treg"
        sce.tmp$group[sce.tmp$meta.cluster=="CD8.c12.Tex.CXCL13"] <- "terminal Tex"
        sce.tmp$group <- factor(sce.tmp$group,levels=c("terminal Tex","CD4+ Treg","Other"))
        table(sce.tmp$group)

        b.clamp <- c(-0.5,2.5)
        gene.plot <- c("HAVCR2","CTLA4","CXCL13","FOXP3")
        p <- ssc.plot.violin(sce.tmp,gene=gene.plot,group.var="group",
                 clamp=b.clamp,
                 add="boxplot") +
            xlab("") + ylab("z-score") +
            coord_cartesian(ylim=b.clamp) +
            geom_hline(yintercept=0.3,linetype="dashed") +
            #theme_pubr() +
            theme(
              #axis.ticks.x=element_blank(),
              axis.text=element_text(size=12),
              strip.text=element_text(size=12),
              legend.position="right",
              panel.grid.minor=element_blank())
        ggsave(sprintf("%s.CD8.Treg.00.pdf",out.prefix),width=3,height=5.5)

    }

    #### Tc17/MAIT (Fig. 2F ???)
    {

        ##
        l_ply(c("RORC","IL23R","IL17A","IL17F","SLC4A10","IL26"),function(a.gene) {

            p.list <- list()
            sce.tmp <- sce.CD8[rowData(sce.CD8)$display.name==a.gene,]
            rowData(sce.tmp)$display.name <- "PanC"
            rownames(sce.tmp) <- "PanC"
            p.list[[sprintf("%s.overall",a.gene)]] <- ssc.plot.violin(sce.tmp,
                                    #gene=a.gene,group.var="cluster.id",
                                    gene="PanC",group.var="cluster.id",
                                    clamp=c(-0.5,1.5),
                                    add="boxplot") +
                xlab("") +
                coord_cartesian(ylim=c(-0.5,1.5)) +
                #theme_pubr() +
                theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
                  strip.text=element_text(size=12),
                  legend.position="right",
                  panel.grid.minor=element_blank())
            for(a.cancerType in c("ESCA","STAD","SCC"))
            {
                sce.tmp <- sce.CD8[rowData(sce.CD8)$display.name==a.gene,sce.CD8$cancerType==a.cancerType]
                rowData(sce.tmp)$display.name <- a.cancerType
                rownames(sce.tmp) <- a.cancerType
                p.list[[sprintf("%s.%s",a.gene,a.cancerType)]] <- ssc.plot.violin(sce.tmp,
                                     #gene=a.gene,group.var="cluster.id",
                                     gene=a.cancerType,group.var="cluster.id",
                                     ###clamp=c(-1.5,3), 
                                     clamp=c(-0.5,1.5), 
                                     add="boxplot") +
                xlab("") +
                coord_cartesian(ylim=c(-0.5,1.5)) +
                #theme_pubr() +
                theme(strip.text=element_text(size=12),
                      legend.position="right",
                      panel.grid.minor=element_blank())
            }

            p <- plot_grid(plotlist=p.list,align="hv",ncol=1)
            ##ggsave(sprintf("%s.ectopic.IL26.ESCA.01.pdf",out.prefix),width=6.20,height=3.5)
            ggsave(sprintf("%s.ectopic.%s.3cancer.01.pdf",out.prefix,a.gene),width=6.20,height=7)

        },.parallel=T)

    }



}





