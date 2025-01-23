#!/usr/bin/env Rscript

library("Startrac")
library("scPip")
library("sscVis")
library("data.table")
library("plyr")
library("R.utils")
library("ggpubr")
library("ggplot2")
source("./func.R")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

out.prefix <- "OUT_Fig2/regenUMAP/regenUMAP"
dir.create(dirname(out.prefix),F,T)

### data
{
    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    sce.Path.file <- "../data/expression/CD8/integration/int.CD8.S35.HVG.continumOnly.v1.sce.Path.rds" 
    startrac.file <- "../data/tcr/startrac/CD8/CD8.out.startrac.nperm1000.rds"
}

if(!file.exists(sce.Path.file))
{
    ### expression data
    {
        sce.file <- "../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds"
        sce <- readRDS(sce.file)

        exp.list.file <- "list/obj.inte.all.CD8.post.limma.continumOnly.v1.list"
        exp.list.table <- fread(cmd=sprintf("awk '!/^#/' %s",exp.list.file),header=F)
        colnames(exp.list.table) <- c("stype","data.id","measurement","platform","efile","dfile")

        dat.ext.dir <- system.file("extdata",package="scPip")
        gene.exclude.file <- sprintf("%s/exclude.gene.misc.misc.v3.RData",dat.ext.dir)
        env.misc <- loadToEnv(gene.exclude.file)

    }

    ############### HVG #####
    {
        gene.de.list <- list()
        for(i in seq_len(nrow(exp.list.table))){
            id.d <- exp.list.table$data.id[i]
            ifile <- exp.list.table$dfile[i]
            de.out <- readRDS(ifile)
            gene.de.list[[id.d]] <- de.out$all
            gene.de.list[[id.d]]$geneID <- gene.de.list[[id.d]]$geneSymbol
        }
        names(gene.de.list) <- exp.list.table$data.id

        gene.common <- rownames(sce)

        gene.rank.tb <- as.data.table(ldply(names(gene.de.list),function(x){
                            ret.tb <- unique(gene.de.list[[x]][,c("geneID","F.rank")])
                            ret.tb$dataset.id <- x
                            ret.tb <- ret.tb[geneID %in% gene.common,]
                            return(ret.tb) }))

        gene.rank.tb <- dcast(gene.rank.tb,geneID~dataset.id,value.var="F.rank",fill=1)
        gene.rank.tb$median.F.rank <- rowMedians(as.matrix(gene.rank.tb[,-c("geneID"),with=F]))
        gene.rank.tb <- gene.rank.tb[order(median.F.rank),]
        ##gene.rank.tb <- gene.rank.tb[geneID %in% gene.common,]
        #rowData(sce.pb)$median.F.rank <- gene.rank.tb[["median.F.rank"]][match(rownames(sce.pb),gene.rank.tb$geneID)]

        #### ASH1L MALAT1  CLDN15  ZNF117 NF1 ZNF546 HLA-A HLA-B HLA-C HLA-E PDIA3 IL2RG PTPRC B2M 
        f.gene.blackList <- (gene.rank.tb$geneID %in% env.misc$all.gene.ignore.df$geneSymbol) |
        grepl("^RP[LS]",gene.rank.tb$geneID,perl=T) |
        gene.rank.tb$geneID=="MALAT1"
            
        #### select genes
        #gene.de.common <- head(gene.rank.tb[!f.gene.blackList,][["geneID"]],n=1500)
        #saveRDS(gene.rank.tb[!f.gene.blackList,],sprintf("%s.gene.rank.tb.flt.rds",out.prefix))
        gene.de.common.tmp.tb <- HVG.From.GeneRankTb(gene.rank.tb[!f.gene.blackList,],
                             n.common=1500,n.specific=0000,th.rank=0.1)
        gene.de.common.tmp.tb <- gene.de.common.tmp.tb[!is.na(geneID),]

        write.table(gene.de.common.tmp.tb,
            file=sprintf("%s.gene.de.common.1500.0000.tb",out.prefix),
            row.names=F,sep="\t",quote=F)
        gene.de.common.tmp.tb <- fread(sprintf("%s.gene.de.common.1500.0000.tb",out.prefix))

    }

    ############### UMAP ####
    ### re-calculate the UMAP
    {
        f.cell <- sce$meta.cluster %in% c("CD8.c01.Tn.MAL","CD8.c02.Tm.IL7R",
                          "CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK",
                          "CD8.c10.Trm.ZNF683",
                          "CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13")
        sce.Path <- sce[,f.cell]
        exp.dat <- assay(sce.Path,"exprs")[gene.de.common.tmp.tb$geneID,]
        mdata.df <- as.data.frame(colData(sce.Path))
        mdata.df[1:2,]
        seu <- CreateSeuratObject(exp.dat,project="panC", meta.data=mdata.df)
        seu <- SetAssayData(seu,"scale.data",exp.dat)

        npc <- 15
        seu <- RunPCA(seu,features = rownames(seu), npcs = npc, verbose = FALSE)
        seu <- RunUMAP(seu,reduction="pca",dims=1:npc)
        seu <- RunHarmony(seu, c("dataset"))
        seu <- RunUMAP(seu,reduction = "harmony",reduction.name = "harmony.umap", dims = 1:npc)

        for(x in c("pca","umap","harmony","harmony.umap"))
        {
            rd.x <- Embeddings(object = seu, reduction = x)
            print(all(colnames(sce.Path)==rownames(rd.x)))
            reducedDim(sce.Path,sprintf("recal.%s",x)) <- rd.x
        }

        saveRDS(sce.Path,file=sce.Path.file)
    }
}else{ 
    sce.Path <- readRDS(file=sce.Path.file)
}


{

    ### 
    {
        #### cluster in the continumn
        mcls.vec <- c("CD8.c01.Tn.MAL","CD8.c02.Tm.IL7R",
                  "CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK",
                  "CD8.c10.Trm.ZNF683",
                  "CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13")
        colSet <- g.colSet
        colSet$meta.cluster <- colSet$meta.cluster[mcls.vec]

        #### pTrans network
        {
            ### pairwise transition network
            out <- readRDS(startrac.file)
            dat.export.tb <- as.data.table(out@pIndex.tran)[aid=="panC",]
            dat.export.mx <- as.matrix(dat.export.tb[,-c("aid","majorCluster"),with=F])
            rownames(dat.export.mx) <- dat.export.tb[["majorCluster"]]
            dat.export.mx[lower.tri(dat.export.mx)] <- NA
            dat.export.flt.tb <- as.data.frame(dat.export.mx)
            setDT(dat.export.flt.tb,T)
            dat.export.flt.tb <- melt(dat.export.flt.tb,id.var="rn")
            dat.export.flt.tb <- dat.export.flt.tb[!is.na(value),]
            colnames(dat.export.flt.tb) <- c("mcls.1","mcls.2","pTran")
            pTrans.tb <- dat.export.flt.tb[pTran > 0.1,]
            
            pTrans.tb <- pTrans.tb[mcls.1 %in% mcls.vec & mcls.2 %in% mcls.vec,]
            pTrans.tb[,linetype:="solid"]
            pTrans.patch.tb <- rbind(pTrans.tb,data.table(mcls.1=c("CD8.c01.Tn.MAL","CD8.c01.Tn.MAL"),
                                      mcls.2=c("CD8.c02.Tm.IL7R","CD8.c12.Tex.CXCL13"),
                                      pTran=0,
                                      linetype=c("dashed","direct")))
        }

        #### umap coordinate
        dat.map <- reducedDim(sce.Path,"recal.harmony.umap")
        dat.map.tb <- cbind(data.table(cellID=colnames(sce.Path),
                           meta.cluster=as.character(sce.Path$meta.cluster)),
                    dat.map)
        dat.map.center.tb <- dat.map.tb[,.(Dim1=median(harmony.umap_1),
                           Dim2=median(harmony.umap_2)
                           ),by="meta.cluster"][order(meta.cluster),]
        #### distances between cluster centers
        dat.map.center.dist <- as.matrix(dist(dat.map.center.tb[,-c("meta.cluster")],diag=T,upper=T))
        dat.map.center.dist.tb <- cbind(data.table(mcls.1=rownames(dat.map.center.dist)),as.data.table(dat.map.center.dist))
        dat.map.center.dist.melt.tb <- melt(dat.map.center.dist.tb,id.vars="mcls.1",variable.name="mcls.2",value.name="gep.dist")
        dat.map.center.dist.melt.tb[,mcls.1:=dat.map.center.tb$meta.cluster[as.integer(as.character(mcls.1))] ]
        dat.map.center.dist.melt.tb[,mcls.2:=dat.map.center.tb$meta.cluster[as.integer(as.character(mcls.2))] ]
        dat.map.center.dist.melt.tb[mcls.1=="CD8.c01.Tn.MAL",]
        ##            mcls.1             mcls.2   gep.sim
        ## 1: CD8.c01.Tn.MAL     CD8.c01.Tn.MAL  0.000000
        ## 2: CD8.c01.Tn.MAL    CD8.c02.Tm.IL7R  3.723553
        ## 3: CD8.c01.Tn.MAL  CD8.c05.Tem.CXCR5  5.402684
        ## 4: CD8.c01.Tn.MAL   CD8.c06.Tem.GZMK  8.676101
        ## 5: CD8.c01.Tn.MAL CD8.c10.Trm.ZNF683  8.480233
        ## 6: CD8.c01.Tn.MAL  CD8.c11.Tex.PDCD1 11.309992
        ## 7: CD8.c01.Tn.MAL CD8.c12.Tex.CXCL13 13.111907

        #### graph data
        {
            dat.graph.tb <- merge(pTrans.patch.tb,dat.map.center.dist.melt.tb)
            dat.graph.tb <- merge(dat.graph.tb,dat.map.center.tb,by.x="mcls.1",by.y="meta.cluster")
            dat.graph.tb <- merge(dat.graph.tb,dat.map.center.tb,by.x="mcls.2",by.y="meta.cluster")
            dat.graph.tb[,dis.pTran:=1-pTran]
            #### STARTRAC-based
            dat.graph.startrac.tb <- copy(dat.graph.tb[linetype!="direct",])
            dat.graph.startrac.tb <- dat.graph.startrac.tb[c(1,2,4,5,9,11,12),]

            #dat.graph.Dijkstra.tb <- copy(dat.graph.tb[linetype!="direct",])
            #dat.graph.Dijkstra.tb <- dat.graph.Dijkstra.tb[c(1,3,9,12),]

            dat.graph.direct.tb <- dat.graph.tb[linetype=="direct",]

            dat.graph.PAGA.tb <- dat.graph.tb[c(1,2,4,5,7,9,13),]
            
            dat.graph.Slingshot.tb <- dat.graph.tb[c(1,4,5,7,9,13),]
        }

        ##### evaluation of methods using ratio of SS (fig. S14E)
        {
            ### v, w, p should be data.table/list with columns "x" and "y"
            distSqr <- function(v, w) { return((v$x - w$x)^2 + (v$y - w$y)^2) }
            distToSegmentSquared <- function(p, v, w) {
                l2 <- distSqr(v, w)
                if (l2 == 0) return(distSqr(p, v))
                tt <- as.numeric(c(p$x-v$x,p$y-v$y) %*% c(w$x-v$x,w$y-v$y) / l2)
                tt = max(0, min(1, tt))
                return(distSqr(p, data.table(x=v$x + tt * (w$x - v$x),
                                 y=v$y + tt * (w$y - v$y) )))
            }
            distToSegment <- function(p, v, w) { return(sqrt(distToSegmentSquared(p, v, w))) }

            #distToSegmentSquared(p=data.table(x=3,y=2),v=data.table(x=0,y=0),w=data.table(x=0,y=3))

            datToGraphSS <- function(out.prefix,dat.map.tb,dat.graph,dat.map.center.tb,sce.plot)
            {
                SS <- as.data.table(ldply(seq_len(nrow(dat.map.tb)),function(i){
                      p.tb <- dat.map.tb[i,]
                      p.tb <- data.table(x=p.tb[["harmony.umap_1"]],y=p.tb[["harmony.umap_2"]])
                      data.table(distS=min(sapply(seq_len(nrow(dat.graph)),function(j){
                             d.tb <- dat.graph[j,]
                             v.tb <- data.table(x=d.tb[["Dim1.x"]],y=d.tb[["Dim2.x"]])
                             w.tb <- data.table(x=d.tb[["Dim1.y"]],y=d.tb[["Dim2.y"]])
                             #print(str(v.tb))
                             #print(str(w.tb))
                             distToSegmentSquared(p.tb,v.tb,w.tb)
                               })))
                               },.parallel=T))
                ret.list <- list(dat.map=cbind(dat.map.tb,SS), sumSS=sum(SS$distS))
                if(!is.null(sce.plot)){
                    m.tb <- ret.list$dat.map
                    setkey(m.tb,"cellID")
                    sce.plot$distS <- m.tb[colnames(sce.plot),][["distS"]]

                    p.multi <- ssc.plot.tsne(sce.plot,columns = c("distS"),reduced.name = "recal.harmony.umap",
                      vector.friendly=T,legend.w=1.2,
                      theme.use=theme_pubr,verbose=T,
                      size=0.6,
                      palette.name="Oranges",
                      clamp=c(0,10),
                      par.legend=list(breaks=seq(0,10,2),limits=c(0,10)),
                      par.geom_point=list(raster.width=7,raster.height=6),
                      par.geneOnTSNE = list(pt.order = "random"),
                      colSet=colSet)
                    p <- p.multi$list[[1]] + theme(legend.position="right")
                    p <- p + geom_point(data=dat.map.center.tb,color="black") +
                        geom_segment(data=dat.graph,
                                 aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y),
                                 colour = "black") +
                        theme_bw(base_size=12) +
                        theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) +
                        labs(x="UMAP1",y="UMAP2") +
                        coord_fixed(ratio=1, xlim = c(-10, 10), ylim = c(-10, 10), expand = FALSE)
                    ggsave(sprintf("%s.vecFri.graph.00.pdf",out.prefix),
                       width=5,height=3.5,useDingbats=FALSE)

                }
                return(ret.list)
            }

            res.SS.ourInference <- datToGraphSS(out.prefix=sprintf("%s.distS.ourInference",out.prefix),
                         dat.map.tb,dat.graph.startrac.tb,dat.map.center.tb,sce.Path)

#            res.SS.Dijkstra <- datToGraphSS(out.prefix=sprintf("%s.distS.Dijkstra",out.prefix),
#                         dat.map.tb,dat.graph.Dijkstra.tb,dat.map.center.tb,sce.Path)

            res.SS.direct <- datToGraphSS(out.prefix=sprintf("%s.distS.direct",out.prefix),
                         dat.map.tb,dat.graph.direct.tb,dat.map.center.tb,sce.Path)
            
            res.SS.PAGA <- datToGraphSS(out.prefix=sprintf("%s.distS.PAGA",out.prefix),
                         dat.map.tb,dat.graph.PAGA.tb,dat.map.center.tb,sce.Path)
            
            res.SS.Slingshot <- datToGraphSS(out.prefix=sprintf("%s.distS.Slingshot",out.prefix),
                         dat.map.tb,dat.graph.Slingshot.tb,dat.map.center.tb,sce.Path)

            dat.plot.sumSS <- data.table(method=c("Our Inference","PAGA","Slingshot","Direct Link"),
                         sumSS=c(res.SS.ourInference$sumSS/res.SS.ourInference$sumSS,
                             res.SS.PAGA$sumSS/res.SS.ourInference$sumSS,
                             res.SS.Slingshot$sumSS/res.SS.ourInference$sumSS,
                             res.SS.direct$sumSS/res.SS.ourInference$sumSS))
            dat.plot.sumSS[,sumSS:=1/sumSS]
            dat.plot.sumSS[,method:=factor(method,levels=dat.plot.sumSS$method)]

            p <- ggbarplot(dat.plot.sumSS,x="method",y="sumSS",color=NA,
                   fill=RColorBrewer::brewer.pal(9,"Blues")[7]) +
                labs(y="Ratio Of SS",x="") +
                theme(axis.text.x=element_text(angle=60,hjust=1))
            ggsave(sprintf("%s.distS.sumSS.barplot.00.pdf",out.prefix),width=3,height=3.5)

        }

    }

    ###### connection with pTrans
    {
        p.multi <- ssc.plot.tsne(sce.Path,columns = c("meta.cluster"),reduced.name = "recal.harmony.umap",
                  vector.friendly=T,legend.w=1.2,
                  theme.use=theme_pubr,verbose=T,
                  size=0.6,
                  fun.extra=function(p){
                p + theme_bw(base_size=12) +
                    theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=14)) +
                    labs(x="UMAP1",y="UMAP2") +
                    coord_fixed(ratio=1,
                        xlim = c(-10, 10), ylim = c(-10, 10),
                        expand = FALSE)
                  },
                  ####xlim=c(-10,10), ylim=c(-10,10),
                  #par.geom_point=list(raster.width=7,raster.height=6),
                  par.geom_point=list(scale=0.4),
                  par.geneOnTSNE = list(pt.order = "random"),
                  colSet=colSet)
        p <- p.multi$list[[1]] + theme(legend.position="right")
        ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)
        p <- p + geom_point(data=dat.map.center.tb,color="black") + 
            geom_segment(data=dat.graph.tb,
                 aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y,
                     #size=pTran,
                     linetype=linetype),
                 colour = "black") +
            scale_linetype_manual(values=c("solid"="solid","dashed"="dashed"))
        ###ggsave(sprintf("%s.meta.cluster.CD8.vecFri.graph.lim.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)
        ggsave(sprintf("%s.meta.cluster.CD8.vecFri.graph.lim.v1.pdf",out.prefix),width=5,height=3.5,useDingbats=FALSE)

    }

    ##### show the graphs inferred by eahc method (fig. S14D)
    {
        ###### startrac-based graph
        {
            p.multi <- ssc.plot.tsne(sce.Path,columns = c("meta.cluster"),reduced.name = "recal.harmony.umap",
                      vector.friendly=T,legend.w=1.2,
                      theme.use=theme_pubr,verbose=T,
                      size=0.6,
                      ##xlim=c(-10,10), ylim=c(-10,10),
                      #par.geom_point=list(raster.width=7,raster.height=6),
                      par.geom_point=list(scale=0.4),
                      par.geneOnTSNE = list(pt.order = "random"),
                      colSet=colSet)
            p <- p.multi$list[[1]] + theme(legend.position="right")
            #ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)
            p <- p + geom_point(data=dat.map.center.tb,color="black") + 
                geom_segment(data=dat.graph.startrac.tb,
                     aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y
                         #size=pTran,
                         #linetype=linetype
                         ),
                     colour = "black") +
                theme_bw(base_size=12) +
                theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14)) +
                labs(x="UMAP1",y="UMAP2") +
                #scale_linetype_manual(values=c("solid"="solid","dashed"="dashed")) +
                coord_fixed(ratio=1,
                    xlim = c(-10, 10),
                    ylim = c(-10, 10),
                    expand = FALSE)
            ggsave(sprintf("%s.meta.cluster.CD8.vecFri.graph.lim.startrac.01.pdf",out.prefix),width=5,height=3.5,useDingbats=FALSE)

        }

        ##### PAGA
        {
            p.multi <- ssc.plot.tsne(sce.Path,columns = c("meta.cluster"),reduced.name = "recal.harmony.umap",
                      vector.friendly=T,legend.w=1.2,
                      theme.use=theme_pubr,verbose=T,
                      size=0.6,
                      ##xlim=c(-10,10), ylim=c(-10,10),
                      #par.geom_point=list(raster.width=7,raster.height=6),
                      par.geom_point=list(scale=0.4),
                      par.geneOnTSNE = list(pt.order = "random"),
                      colSet=colSet)
            p <- p.multi$list[[1]] + theme(legend.position="right")
            #ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)
            p <- p + geom_point(data=dat.map.center.tb,color="black") + 
                geom_segment(data=dat.graph.PAGA.tb,
                     aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y
                         #size=pTran,
                         #linetype=linetype
                         ),
                     colour = "black") +
                theme_bw(base_size=12) +
                theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14)) +
                labs(x="UMAP1",y="UMAP2") +
                #scale_linetype_manual(values=c("solid"="solid","dashed"="dashed")) +
                coord_fixed(ratio=1,
                    xlim = c(-10, 10),
                    ylim = c(-10, 10),
                    expand = FALSE)
            ggsave(sprintf("%s.meta.cluster.CD8.vecFri.graph.lim.PAGA.01.pdf",out.prefix),width=5,height=3.5,useDingbats=FALSE)

        }

        ##### Slingshot
        {
            p.multi <- ssc.plot.tsne(sce.Path,columns = c("meta.cluster"),reduced.name = "recal.harmony.umap",
                      vector.friendly=T,legend.w=1.2,
                      theme.use=theme_pubr,verbose=T,
                      size=0.6,
                      ##xlim=c(-10,10), ylim=c(-10,10),
                      #par.geom_point=list(raster.width=7,raster.height=6),
                      par.geom_point=list(scale=0.4),
                      par.geneOnTSNE = list(pt.order = "random"),
                      colSet=colSet)
            p <- p.multi$list[[1]] + theme(legend.position="right")
            #ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)
            p <- p + geom_point(data=dat.map.center.tb,color="black") + 
                geom_segment(data=dat.graph.Slingshot.tb,
                     aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y
                         #size=pTran,
                         #linetype=linetype
                         ),
                     colour = "black") +
                theme_bw(base_size=12) +
                theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14)) +
                labs(x="UMAP1",y="UMAP2") +
                #scale_linetype_manual(values=c("solid"="solid","dashed"="dashed")) +
                coord_fixed(ratio=1,
                    xlim = c(-10, 10),
                    ylim = c(-10, 10),
                    expand = FALSE)
            ggsave(sprintf("%s.meta.cluster.CD8.vecFri.graph.lim.Slingshot.01.pdf",out.prefix),width=5,height=3.5,useDingbats=FALSE)

        }

        ###### Direct Link
        {
            p.multi <- ssc.plot.tsne(sce.Path,columns = c("meta.cluster"),reduced.name = "recal.harmony.umap",
                      vector.friendly=T,legend.w=1.2,
                      theme.use=theme_pubr,verbose=T,
                      size=0.6,
                      ##xlim=c(-10,10), ylim=c(-10,10),
                      #par.geom_point=list(raster.width=7,raster.height=6),
                      par.geom_point=list(scale=0.4),
                      par.geneOnTSNE = list(pt.order = "random"),
                      colSet=colSet)
            p <- p.multi$list[[1]] + theme(legend.position="right")
            #ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)
            p <- p + geom_point(data=dat.map.center.tb,color="black") + 
                geom_segment(data=dat.graph.direct.tb,
                     aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y
                         #size=pTran,
                         #linetype=linetype
                         ),
                     colour = "black") +
                theme_bw(base_size=12) +
                theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14)) +
                labs(x="UMAP1",y="UMAP2") +
                #scale_linetype_manual(values=c("solid"="solid","dashed"="dashed")) +
                coord_fixed(ratio=1,
                    xlim = c(-10, 10),
                    ylim = c(-10, 10),
                    expand = FALSE)
            ggsave(sprintf("%s.meta.cluster.CD8.vecFri.graph.lim.direct.01.pdf",out.prefix),width=5,height=3.5,useDingbats=FALSE)

        }


    }

    ###### a few genes (fig. S14A)
    {

        p <- ssc.plot.tsne(sce.Path,gene=c("CXCL13","IL7R","GZMK","ZNF683"),
                   reduced.name="recal.harmony.umap",
                   vector.friendly=T,clamp=c(-0.5,1.5),
                   p.ncol=4,
                   fun.extra=function(p){
                       p + theme_bw(base_size=12) +
                       theme(axis.text=element_text(size=12),
                         plot.title = element_text(hjust=0.5)) +
                       labs(x="UMAP1",y="UMAP2") +
                       coord_fixed(ratio=1, xlim = c(-10, 10), ylim = c(-10, 10), expand = FALSE)
                   },
                   theme.use=theme_void,size=0.55,
                   par.geneOnTSNE=list(scales="fixed",
                               pt.order="random",pt.alpha = 0.5))
        ##ggsave(sprintf("%s.geneExample.CD8.vecFri.pdf",out.prefix), width=6.5,height=4.5,useDingbats=FALSE)
        ##ggsave(sprintf("%s.geneExample.CD8.vecFri.test.pdf",out.prefix), width=10,height=8,useDingbats=FALSE)
        ggsave(sprintf("%s.geneExample.CD8.vecFri.pdf",out.prefix), width=15,height=3.5,useDingbats=FALSE)

        p <- ssc.plot.tsne(sce.Path,gene=c("IL7R","GZMK","ZNF683","CXCL13"),
                   reduced.name="recal.harmony.umap",
                   vector.friendly=T,clamp=c(-0.5,1.5),
                   p.ncol=2,
                   fun.extra=function(p){
                       p + theme_bw(base_size=12) +
                       theme(axis.text=element_text(size=12),
                         plot.title = element_text(hjust=0.5)) +
                       labs(x="UMAP1",y="UMAP2") +
                       coord_fixed(ratio=1, xlim = c(-10, 10), ylim = c(-6, 6), expand = FALSE)
                   },
                   theme.use=theme_void,size=0.55,
                   par.geneOnTSNE=list(scales="fixed",
                               pt.order="random",pt.alpha = 0.5))
        ##ggsave(sprintf("%s.geneExample.CD8.vecFri.pdf",out.prefix), width=6.5,height=4.5,useDingbats=FALSE)
        ##ggsave(sprintf("%s.geneExample.CD8.vecFri.test.pdf",out.prefix), width=10,height=8,useDingbats=FALSE)
        ggsave(sprintf("%s.geneExample.CD8.vecFri.v2.pdf",out.prefix), width=6.5,height=4,useDingbats=FALSE)

        p <- ssc.plot.tsne(sce.Path,gene=c("IL7R","GZMK","ZNF683","CXCL13"),
                   reduced.name="recal.harmony.umap",
                   vector.friendly=T,clamp=c(-0.5,1.5),
                   p.ncol=2,
                   fun.extra=function(p){
                       p + theme_void(base_size=12) +
                       theme(
                             #axis.text=element_text(size=12),
                             axis.text=element_blank(),
                         plot.title = element_text(hjust=0.5)) +
                       labs(x="UMAP1",y="UMAP2") +
                       coord_fixed(ratio=1, xlim = c(-10, 10), ylim = c(-6, 6), expand = FALSE)
                   },
                   theme.use=theme_void,size=0.55,
                   par.geneOnTSNE=list(scales="fixed",
                               pt.order="random",pt.alpha = 0.5))
        ##ggsave(sprintf("%s.geneExample.CD8.vecFri.pdf",out.prefix), width=6.5,height=4.5,useDingbats=FALSE)
        ##ggsave(sprintf("%s.geneExample.CD8.vecFri.test.pdf",out.prefix), width=10,height=8,useDingbats=FALSE)
        ggsave(sprintf("%s.geneExample.CD8.vecFri.v3.pdf",out.prefix), width=6.5,height=4,useDingbats=FALSE)

    }

}




