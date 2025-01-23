#!/usr/bin/env Rscript

library("ggpubr")
library("ggplot2")
library("data.table")
library("sscClust")
library("fmsb")
library("grid")
library("ComplexHeatmap")
library("plyr")
source("./func.R")

g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

cluster.name.tb <- fread("../data/metaInfo/name.conversion.txt",head=T)
mcls2Name <- structure(factor(cluster.name.tb$cluster.name.full,
                              levels=cluster.name.tb$cluster.name.full),
                       names=cluster.name.tb$meta.cluster)


##### REACTOME (fig. S08)
{

    th.FDR <- 0.05

    gsea.tb.file <- "./gsea/gsea.tb.C2.CP.txt"
    out.prefix <- "./OUT_Fig1/gsea/gsea.REACTOME"

    dir.create(dirname(out.prefix),F,T)

    gsea.tb <- fread(gsea.tb.file)
    colnames(gsea.tb) <- make.names(colnames(gsea.tb))

    f.gset <- unique(gsea.tb[FDR.q.val < 0.05,][["NAME"]])

    dat.plot.NES <- dcast(gsea.tb,NAME~meta.cluster,value.var="NES")
    dat.plot.FDR <- dcast(gsea.tb,NAME~meta.cluster,value.var="FDR.q.val")
    
    dat.plot.NES.mtx <- as.matrix(dat.plot.NES[,-1])
    rownames(dat.plot.NES.mtx) <- dat.plot.NES[[1]]
    dat.plot.NES.mtx <- dat.plot.NES.mtx[f.gset,]
    dat.plot.NES.mtx <- dat.plot.NES.mtx[grepl("^REACTOME_",rownames(dat.plot.NES.mtx),perl=T),]
    
    dat.plot.FDR.mtx <- as.matrix(dat.plot.FDR[,-1])
    rownames(dat.plot.FDR.mtx) <- dat.plot.FDR[[1]]
    dat.plot.FDR.mtx <- dat.plot.FDR.mtx[f.gset,]
    dat.plot.FDR.mtx <- dat.plot.FDR.mtx[grepl("^REACTOME_",rownames(dat.plot.FDR.mtx),perl=T),]


    #### OR
    {
        ### output from w.step1.Figure1.Barplot.Boxplot.OR.Index.R
        OR.df <- read.csv("../data/metaInfo/Fig1.OR.all.baseline.csv")
        OR.df[ OR.df > 3 ] <- 3
        print(head(OR.df))

        OR.hclust.row <- run.cutree(as.matrix(OR.df),k=4,method.distance="cosine",method.hclust="ward.D2")
        OR.hclust.row$branch <- dendextend::set(OR.hclust.row$branch,"branches_lwd", 2)
       
        col.OR <- circlize::colorRamp2(c(0, 1, 3), viridis::viridis(3))
        ca <- ComplexHeatmap::HeatmapAnnotation(df = OR.df, col=list("P"=col.OR,"N"=col.OR,"T"=col.OR),
                            show_legend=T,annotation_height=unit(0.5,"cm"),
                            annotation_legend_param=list(color_bar="continuous", at=seq(0,3,1),
                                         grid_width = unit(0.4,"cm"),
                                         grid_height = unit(0.4, "cm"),
                                         legend_width = unit(6, "cm"),
                                         #legend_height = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 10),
                                         direction="horizontal",
                                         title_gp = gpar(fontsize = 12, fontface = "bold")))
    }

    f.TCR.signaling <- c("REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
			 "REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE",
			 "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
			 "REACTOME_DOWNSTREAM_TCR_SIGNALING",
			 "REACTOME_TCR_SIGNALING")
    
    #### clustering using reactome pathways of high sd
    {

        dat.plot.NES.mtx[1:4,1:3]
        dat.plot.FDR.mtx[1:4,1:3]
        row.var <- rowSds(dat.plot.NES.mtx)
        #f.goi <- row.var >= quantile(row.var,0.75)
        f.goi <- row.var >= quantile(row.var,0.8)

        dat.plot.NES.var.mtx <- dat.plot.NES.mtx[f.goi,]
        dat.plot.FDR.var.mtx <- dat.plot.FDR.mtx[f.goi,]

        colnames(dat.plot.NES.var.mtx) <- fetchMetaClusterID2CusterFullName()[colnames(dat.plot.NES.var.mtx)]
        rownames(dat.plot.NES.var.mtx) <- gsub("^REACTOME_","",rownames(dat.plot.NES.var.mtx))
        colnames(dat.plot.FDR.var.mtx) <- fetchMetaClusterID2CusterFullName()[colnames(dat.plot.FDR.var.mtx)]
        rownames(dat.plot.FDR.var.mtx) <- gsub("^REACTOME_","",rownames(dat.plot.FDR.var.mtx))
        dat.plot.NES.var.mtx <- dat.plot.NES.var.mtx[,rownames(OR.df)]
        dat.plot.FDR.var.mtx <- dat.plot.FDR.var.mtx[,rownames(OR.df)]

        print(all(rownames(dat.plot.NES.var.mtx)==rownames(dat.plot.FDR.var.mtx)))
        print(all(colnames(dat.plot.NES.var.mtx)==colnames(dat.plot.FDR.var.mtx)))

        zz <- 3.0
        dat.plot.NES.var.mtx.z <- dat.plot.NES.var.mtx
        dat.plot.NES.var.mtx.z[ dat.plot.NES.var.mtx > zz] <- zz
        dat.plot.NES.var.mtx.z[ dat.plot.NES.var.mtx < -zz] <- -zz

        obj.hclust.col <- run.cutree(t(dat.plot.NES.var.mtx),k=4,method.distance="cosine")
        obj.hclust.col$branch <- dendextend::set(obj.hclust.col$branch,"branches_lwd", 1)
        obj.hclust.row <- run.cutree(dat.plot.NES.var.mtx,k=1,method.distance="cosine")
        obj.hclust.row$branch <- dendextend::set(obj.hclust.row$branch,"branches_lwd", 1)

        ### whited out
        dat.plot.NES.var.mtx.z[ dat.plot.FDR.var.mtx > th.FDR ] <- 0


        f.TCR.signaling
            gene.highlight <- gsub("^REACTOME_","",f.TCR.signaling)
            gene.highlight.at <- match(gene.highlight,rownames(dat.plot.NES.var.mtx.z))
        ha <- rowAnnotation(geneH = anno_mark(at = gene.highlight.at,
                              labels_gp=gpar(fontsize=10),
                              labels = rownames(dat.plot.NES.var.mtx.z)[gene.highlight.at]))

        #### slim (show pathways)
        do.it <- function(show.highlight=F,pdf.width=10,pdf.height=14)
        {
            sscVis:::plotMatrix.simple(dat.plot.NES.var.mtx.z,
                       show.dendrogram=T, z.lo=-zz,z.hi=zz,
                       par.legend=list(at=seq(-zz,zz,1),
                               direction="horizontal", legend_width = unit(6, "cm"),
                               labels_gp = gpar(fontsize = 10),
                               title_gp = gpar(fontsize = 12, fontface = "bold")),
                       clust.row=obj.hclust.row$branch,
                       ##clust.column=obj.hclust.col$branch,
                       ##clust.column=OR.hclust.row.pre$branch,
                       #clust.column=OR.hclust.row$branch,
                       clust.column=obj.hclust.col$branch,
                       palatte=rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")),
                       top_annotation=ca,
                       par.heatmap=list(
                                column_names_gp=gpar(fontsize=10,
                                         col=g.colSet$cluster.name[colnames(dat.plot.NES.var.mtx.z)]),
                                row_names_gp=gpar(fontsize=0),
                                ##cex.column=1.25, cex.row=1.2,
                                ###border=T,
                                row_gap=unit(0.0,"cm"),
                                column_gap=unit(0.0,"cm"),
                                column_dend_height = unit(3, "cm"),
                                row_dend_width = unit(3, "cm"),
                                #row_split=k.row,
                                #row_split=gset.info.tb$Group,
                                #column_split=k.column,
                                right_annotation = if(show.highlight) ha else NULL,
                                use_raster = F,raster_quality=5),
                       out.prefix=sprintf("%s.REACTOME.NES.slim.var.hightligh.%s",out.prefix,show.highlight),
                       mytitle="REACTOME",show.number=F,
                       #heatmap_legend_side = "left", annotation_legend_side = "left",
                       heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                       exp.name="NES",pdf.width=pdf.width,pdf.height=pdf.height)
        }

        ###do.it(show.highlight=F,pdf.width=10,pdf.height=14)
        do.it(show.highlight=T,pdf.width=14,pdf.height=12)

    }

    ### focus on TCR signaling
    {
        dat.plot.NES.mtx[1:4,1:3]
        dat.plot.FDR.mtx[1:4,1:3]
        f.goi <- c("REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
               "REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE",
               "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
               "REACTOME_DOWNSTREAM_TCR_SIGNALING",
               "REACTOME_TCR_SIGNALING")

        dat.plot.NES.TCR.mtx <- dat.plot.NES.mtx[f.goi,]
        dat.plot.FDR.TCR.mtx <- dat.plot.FDR.mtx[f.goi,]

        colnames(dat.plot.NES.TCR.mtx) <- fetchMetaClusterID2CusterFullName()[colnames(dat.plot.NES.TCR.mtx)]
        rownames(dat.plot.NES.TCR.mtx) <- gsub("^REACTOME_","",rownames(dat.plot.NES.TCR.mtx))
        colnames(dat.plot.FDR.TCR.mtx) <- fetchMetaClusterID2CusterFullName()[colnames(dat.plot.FDR.TCR.mtx)]
        rownames(dat.plot.FDR.TCR.mtx) <- gsub("^REACTOME_","",rownames(dat.plot.FDR.TCR.mtx))
        dat.plot.NES.TCR.mtx <- dat.plot.NES.TCR.mtx[,rownames(OR.df)]
        dat.plot.FDR.TCR.mtx <- dat.plot.FDR.TCR.mtx[,rownames(OR.df)]

        print(all(rownames(dat.plot.NES.TCR.mtx)==rownames(dat.plot.FDR.TCR.mtx)))
        print(all(colnames(dat.plot.NES.TCR.mtx)==colnames(dat.plot.FDR.TCR.mtx)))

        zz <- 3.0
        dat.plot.NES.TCR.mtx.z <- dat.plot.NES.TCR.mtx
        dat.plot.NES.TCR.mtx.z[ dat.plot.NES.TCR.mtx > zz] <- zz
        dat.plot.NES.TCR.mtx.z[ dat.plot.NES.TCR.mtx < -zz] <- -zz

        ## use that from var pathways
        #obj.hclust.col <- run.cutree(t(dat.plot.NES.TCR.mtx),k=2,method.distance="cosine")
        #obj.hclust.col$branch <- dendextend::set(obj.hclust.col$branch,"branches_lwd", 4)
        #obj.hclust.row <- run.cutree(dat.plot.NES.TCR.mtx,k=1,method.distance="cosine")
        #obj.hclust.row$branch <- dendextend::set(obj.hclust.row$branch,"branches_lwd", 2)

        ### whited out
        if(T){
            dat.plot.NES.TCR.mtx.z[ dat.plot.FDR.TCR.mtx > th.FDR ] <- 0
        }
        
        #### slim (show a few pathways)
        sscVis:::plotMatrix.simple(dat.plot.NES.TCR.mtx.z,
                       show.dendrogram=T, z.lo=-zz,z.hi=zz,
                       par.legend=list(at=seq(-zz,zz,1),
                               direction="horizontal", legend_width = unit(6, "cm"),
                               labels_gp = gpar(fontsize = 10),
                               title_gp = gpar(fontsize = 12, fontface = "bold")),
                       ###clust.row=obj.hclust.row$branch,
                       ##clust.column=obj.hclust.col$branch,
                       ##clust.column=OR.hclust.row.pre$branch,
                       #clust.column=OR.hclust.row$branch,
                       clust.column=obj.hclust.col$branch,
                       palatte=rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")),
                       top_annotation=ca,
                       par.heatmap=list(
                            column_names_gp=gpar(fontsize=10,
                                         col=g.colSet$cluster.name[colnames(dat.plot.NES.TCR.mtx.z)]),
                            row_names_gp=gpar(fontsize=10),
                            ##cex.column=1.25, cex.row=1.2,
                            ###border=T,
                            row_gap=unit(0.0,"cm"),
                            column_gap=unit(0.0,"cm"),
                            column_dend_height = unit(2, "cm"),
                            row_dend_width = unit(3, "cm"),
                            #row_split=k.row,
                            #row_split=gset.info.tb$Group,
                            #column_split=k.column,
                            #right_annotation =  ha,
                            use_raster = F,raster_quality=5),
                       out.prefix=sprintf("%s.REACTOME.NES.slim.TCR.hclustSD",out.prefix),
                       mytitle="REACTOME",show.number=F,
                       #heatmap_legend_side = "left", annotation_legend_side = "left",
                       heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                       exp.name="NES",pdf.width=12,pdf.height=6.8)

    }

}

##### KEGG.flt (fig. S12)
{

    th.FDR <- 0.05
    th.nSig <- 1
    gsea.tb.file <- "./gsea/gsea.tb.KEGG.flt.txt"
    out.prefix <- "./OUT_Fig1/gsea/gsea.KEGG"

    dir.create(dirname(out.prefix),F,T)

    gsea.tb <- fread(gsea.tb.file)
    colnames(gsea.tb) <- make.names(colnames(gsea.tb))
    #### examples
    gsea.plot.tb <- gsea.tb
    gsea.plot.tb[,full.name:=as.character(mcls2Name[meta.cluster])]
    gsea.plot.tb[,pchar:=""]
    gsea.plot.tb[NOM.p.val<0.05 & FDR.q.val < 0.1,pchar:="*"]
    gsea.plot.tb[,vjust:=ifelse(NES>0, 0.5, 1.1)]

    #### Tex
    {
        pathway.use <- c(
                 "KEGG_OXIDATIVE_PHOSPHORYLATION",
                 "KEGG_CITRATE_CYCLE_TCA_CYCLE", 
                 "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                 "KEGG_FATTY_ACID_METABOLISM")
        gsea.plot.tmp <- gsea.plot.tb[NAME %in% pathway.use &
                          meta.cluster %in% c("CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13",
                                  "CD8.c13.Tex.myl12a","CD8.c14.Tex.TCF7"),]
        gsea.plot.tmp[,cluster.name.short:=fetchMetaClusterID2CusterFullName("cluster.name")[as.character(meta.cluster)]]
        gsea.plot.tmp[,NAME:=gsub("^KEGG_","",NAME)]
        gsea.plot.tmp[NOM.p.val<0.01 & FDR.q.val < 0.01,pchar:="*"]
        pathway.use.s <- gsub("^KEGG_","",pathway.use)
        p.list <- llply(pathway.use.s,function(x){
            d.tb <- gsea.plot.tmp[NAME==x,]
            nes.pretty <- pretty(d.tb$NES)
            nes.pretty.step <- nes.pretty[2]-nes.pretty[1]
            p <- ggbarplot(d.tb,
                   x="cluster.name.short",y="NES",
                   ###x="NAME",y="NES",
                   color=NA,fill="cluster.name.short",
                   #palette="npg",
                   ###position=position_dodge2(),
                   xlab="",title=x) +
                #facet_wrap(~NAME,ncol=2,scales="free_y") +
                scale_fill_manual(values=g.colSet$cluster.name.short.CD8) +
                geom_text(aes(label=pchar,vjust=vjust),size=10) +
                coord_cartesian(clip="off",
                        ylim=c(nes.pretty[1]-nes.pretty.step,
                           nes.pretty[length(nes.pretty)]+nes.pretty.step)) +
                theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
                  strip.background=element_blank(),
                  #strip.text=element_text(size=10),
                  legend.position="none")
                                  },.parallel=T)
        names(p.list) <- pathway.use.s
        p <- cowplot::plot_grid(plotlist=p.list,nrow=1,align="hv")
        ggsave(sprintf("%s.barplot.example.Tex.00.pdf",out.prefix),width=6.5,height=4)
    }

}

##############################


