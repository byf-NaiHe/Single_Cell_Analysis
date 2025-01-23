#!/usr/bin/env Rscript


# import libraries

library("magrittr")
library("sscVis")
library("data.table")
library("R.utils")
library("ggpubr")
library("ggplot2")
library("plyr")
library("grid")
library("cowplot")
library("clusterProfiler")
library("ComplexHeatmap")
source("./func.R")
ncores <- 12
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = ncores)

out.prefix <- "OUT_Fig3/Treg/Fig3.Treg"
dir.create(dirname(out.prefix),F,T)

{

    colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    gene.desc.top.slim <- readRDS("../data/expression/CD4/integration/int.CD4.S35.gene.tb.rds")
    gene.desc.top.slim.sig <- gene.desc.top.slim[sig==T,]
    gene.desc.top.slim.sig[,table(meta.cluster)]
    gene.desc.top.slim.sig.Treg <- gene.desc.top.slim.sig[meta.cluster=="CD4.c20.Treg.TNFRSF9",]

    gene.file.foreign <- list(
                  "Plitas2016.mmc2"="../data/external/Plitas2016.mmc2.UP.txt",
                  "Plitas2016.mmc3"="../data/external/Plitas2016.mmc3.UP.txt",
                  "DeSimone2016"="../data/external/DeSimone_2016.Treg.txt",
                  "Tirosh2016"="../data/external/Tirosh2016.Treg.txt",
                  "Zheng2017"="../data/external/Zheng2017.Treg.txt",
                  "Guo2018"="../data/external/Guo2018.Treg.txt"
                  )

    gene.tb.foreign <- as.data.table(ldply(names(gene.file.foreign),function(x){
                           a.tb <- fread(gene.file.foreign[[x]])
                           r.tb <- data.table(study=x,
                                      geneSymbol=a.tb[[2]])
                           return(r.tb)
                        }))
    gene.tb.this <- data.table(study="this study",
                   geneSymbol=gene.desc.top.slim.sig.Treg[["geneSymbol"]])

}

### pattern sorting
{
    dat.plot.tb <- rbind(gene.tb.this,gene.tb.foreign)
    dat.plot.tb[,sig:=1]
    dat.plot.tb[,study:=factor(study,levels=c("this study","Plitas2016.mmc2","Plitas2016.mmc3","DeSimone2016",
					      "Zheng2017","Guo2018","Tirosh2016"))]
    dat.plot.tb <- dcast(dat.plot.tb,geneSymbol~study,value.var="sig",fill=0)
    dat.plot.tb <- merge(dat.plot.tb,gene.desc.top.slim[meta.cluster=="CD4.c20.Treg.TNFRSF9",c("geneSymbol","comb.ES"),with=F])
    dat.plot.tb.sort <- setorderv(dat.plot.tb,colnames(dat.plot.tb)[-1],order=-1)
    dat.plot.tb.sort.mtx <- as.matrix(dat.plot.tb.sort[,-c("geneSymbol","comb.ES")])
    rownames(dat.plot.tb.sort.mtx) <- dat.plot.tb.sort[["geneSymbol"]]

}

#### select genes to highlight
{

    f.gene.core <- which(rowSums(dat.plot.tb.sort[,-c("geneSymbol","comb.ES")]==0)==0)
    dat.plot.tb.sort[f.gene.core,]

    f.gene <- which(rowSums(dat.plot.tb.sort[,3:8])==0 & dat.plot.tb[[2]]>0)
    dat.plot.tb.sort[f.gene,]
    gene.desc.top.slim.sig.Treg$novelty <- gene.desc.top.slim.sig.Treg$geneSymbol %in% dat.plot.tb.sort[f.gene,][["geneSymbol"]]
    write.table(gene.desc.top.slim.sig.Treg,file=sprintf("%s.sig.withNovelty.txt",out.prefix),
		row.names=F,sep="\t",quote=F)

    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.25,
				][,1:8][,-c(1,3)] %>% as.data.frame
    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.20 & geneSet.TF==T,
				][,1:8][,-c(1,3)]
    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.20 & geneSet.kinase==T & geneSet.membrane==T,
				][,1:8][,-c(1,3)]
    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.20 & (geneSet.drug.target.FDA==T | geneSet.drug.target.potential==T),
				][,c(1:8,30)][,-c(1,3)]
    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.15 & (geneSet.cytokine==T |geneSet.cytokineReceptor==T),
				][,c(1:8)][,-c(1,3)]
    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.15 & (geneSet.chemokine==T|geneSet.chemokineReceptor==T),
				][,c(1:8)][,-c(1,3)]
    gene.desc.top.slim.sig.Treg[novelty==T,
				][comb.ES>0.15 & (geneSet.CD.molecular==T),
				][,c(1:8)][,-c(1,3)]

    gene.to.highlight.list <- list("core"=dat.plot.tb.sort[f.gene.core,][["geneSymbol"]],
				   "TF"=c("BCL3", "HIVEP1","TGIF2"),
				   "membrane.kinase"=c("CAMK1","IGF2R"),
				   "drug.target"=c("IFNAR2","TOP1"), ###  "CD247"
				   "cytokineReceptor"=c("IL15RA"), ### "CCR6",
				   "CD.molecule"=c("SELL","CD320"),
				   "other"=c("MIR181A1HG","PHLDA2")
				   )

}

#### "venn"-like heatmpa showing the overlap of signature genes (fig. 27A)
{

    f.gene.specific <- which(rowSums(dat.plot.tb.sort[,3:8])==0 & dat.plot.tb.sort[,2]>0)
    f.gene.common <- which(dat.plot.tb.sort[,2]==1 & !(rowSums(dat.plot.tb.sort[,3:8])==0 & dat.plot.tb.sort[,2]>0) )
    gene.category <- rep("other",nrow(dat.plot.tb.sort))
    gene.category[f.gene.specific] <- "specific"
    gene.category[f.gene.common] <- "common"
    gene.category <- factor(gene.category,levels=c("common","specific","other"))

    ES.range <- c(-0.3,0.6)
    ES.step <- 0.15
    a.palette <- sscVis:::getColorPaletteFromNameContinuous("RdYlBu")
    z.len <- length(a.palette)
    ann.dataset <- HeatmapAnnotation(dataset=c("scRNAseq",
						    "bulkRNAseq","bulkRNAseq","bulkRNAseq",
						    "scRNAseq","scRNAseq","scRNAseq"),
					  col = list(dataset = c("scRNAseq"=RColorBrewer::brewer.pal(9,"Pastel1")[2],
								 "bulkRNAseq"=RColorBrewer::brewer.pal(9,"Pastel1")[3])))
    ann.gene <- rowAnnotation(comb.ES=dat.plot.tb.sort$comb.ES,
			      gene.category=gene.category,
			      col=list(comb.ES=circlize::colorRamp2(seq(ES.range[1],ES.range[2], length = z.len),
							  colorRampPalette(a.palette)(z.len)),
				       gene.category=structure(RColorBrewer::brewer.pal(3,"Greys"),
							       names=c("other","specific","common"))),
			      annotation_legend_param=list(comb.ES=list(at = seq(ES.range[1],ES.range[2],ES.step),
								   grid_width = unit(0.4,"cm"),
								   grid_height = unit(0.4, "cm"),
								   legend_height = unit(6, "cm"))),
			      border=T)
    gene.highlight.at <- match(unlist(gene.to.highlight.list),rownames(dat.plot.tb.sort.mtx))
    gene.highlight.at <- gene.highlight.at[!is.na(gene.highlight.at)]
    ha <- rowAnnotation(geneH = anno_mark(at = gene.highlight.at,labels_gp = gpar(fontsize=10),
    				      labels = rownames(dat.plot.tb.sort.mtx)[gene.highlight.at]))
    sscVis::plotMatrix.simple(dat.plot.tb.sort.mtx,
			      out.prefix=sprintf("%s.patternSort.00",out.prefix),
			      mytitle=sprintf("%s","Treg"),
			      col.ht=c("0"="lightgray","1"="steelblue"),
			      #row.split=gene.spe.tb$ES.cancerType,
			      #row_gap = unit(0, "mm"),
			      #row_title_gp = gpar(fontsize = 10),
			      clust.row=F,
			      z.lo=0,z.hi=1,
			      #par.legend=list(at = seq(0,1,0.5)),
			      pdf.width=6.5,
			      pdf.height=10,
			      top_annotation=ann.dataset,
			      par.heatmap=list(cex.column=1.2,cex.row=0.5,
					       row_title_rot=0,border=T,
					       raster_device="png",
					       raster_quality = 5,
					       right_annotation =  ha,
					       left_annotation = ann.gene
							   ),
			      exp.name="significance")

}

#### novel sig genes heatmpa of effect size (Fig. 3D)
{

    gene.to.highlight.list.novel <- gene.to.highlight.list[setdiff(names(gene.to.highlight.list),"core")]
    gene.plot.tb <- gene.desc.top.slim[geneID %in% unname(unlist(gene.to.highlight.list.novel)),]

    gene.plot.vec <- unique(gene.plot.tb$geneID)
    dat.plot.tb.sort[geneSymbol %in% gene.plot.vec,]
    gene.plot.vec <- dat.plot.tb.sort[geneSymbol %in% gene.plot.vec,][["geneSymbol"]]

    gene.plot.tb[,cluster.name:=fetchMetaClusterID2CusterFullName("cluster.name")[as.character(meta.cluster)]]
    ##gene.plot.tb[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
    cluster.name.tb <- gene.plot.tb[,.N,by=c("meta.cluster","cluster.name")][order(meta.cluster),]
    mcls.vec <- cluster.name.tb$cluster.name
    gene.plot.tb[,y:=factor(geneID,levels=rev(gene.plot.vec))]
    gene.plot.tb[,x:=factor(cluster.name,levels=mcls.vec)]
    gene.plot.tb[,ES:=comb.ES]
    ES.pretty <- c(-0.15,0,0.15,0.3)
    #ES.pretty <- c(-0.3,0,0.3,0.6)
    gene.plot.tb[comb.ES > ES.pretty[length(ES.pretty)],ES:=ES.pretty[length(ES.pretty)]]
    gene.plot.tb[comb.ES < ES.pretty[1],ES:=ES.pretty[1]]
    gene.plot.tb[,Group:="novel"]

    p <- ggplot(gene.plot.tb,aes(x,y)) +
		    geom_point(aes(size=ES,color=ES),shape=16) +
		    facet_grid(Group ~ ., scales = "free", space = "free") +
		    scale_size(breaks=ES.pretty) +
		    scale_colour_distiller(palette = "RdYlBu",
					   guide = guide_colourbar(draw.ulim = T, draw.llim = T),
					   breaks = ES.pretty) +
		    labs(x="",y="") +
		    theme_pubr() + 
		    theme(strip.text.y = element_blank(),
			      axis.line.x=element_blank(),
			      axis.line.y=element_blank(),
			      panel.background = element_rect(colour = "black", fill = "white"),
			      #panel.grid = element_line(colour = "grey", linetype = "dashed"),
			      #panel.grid.major = element_line( colour = "grey", linetype = "dashed", size = 0.2),
			      axis.text.y = element_text(size=10),
			      axis.text.x = element_text(angle = 45,size=10, hjust = 1))
    ##ggsave(sprintf("%s.sigGene.novel.example.01.pdf",out.prefix),width=7,height=6,useDingbats=F)
    ggsave(sprintf("%s.sigGene.novel.example.01.pdf",out.prefix),width=7,height=5,useDingbats=F)
    #saveRDS(p,file=sprintf("%s.sigGene.novel.example.00.rds",out.prefix))

}

#### GO analysis (fig. 27B)
{

    my.parallel <- F

    f.gene <- which(rowSums(dat.plot.tb.sort[,3:8])==0 & dat.plot.tb.sort[,2]>0)
    f.gene.common <- which(dat.plot.tb.sort[,2]==1 & !(rowSums(dat.plot.tb.sort[,3:8])==0 & dat.plot.tb.sort[,2]>0) )
    dat.plot.tb.sort[f.gene,]
    dat.plot.tb.sort[f.gene.common,]
    gene.list.bg <- unique(gene.desc.top.slim[meta.cluster=="CD4.c20.Treg.TNFRSF9",][["geneSymbol"]])

    gene.list.test <- list("Treg.all"=structure(gene.desc.top.slim.sig.Treg[["comb.ES"]],
						names=gene.desc.top.slim.sig.Treg[["geneSymbol"]]),
			   "Treg.common"=structure(dat.plot.tb.sort[f.gene.common,][["comb.ES"]],
						  names=dat.plot.tb.sort[f.gene.common,][["geneSymbol"]]),
			   "Treg.novel"=structure(dat.plot.tb.sort[f.gene,][["comb.ES"]],
						  names=dat.plot.tb.sort[f.gene,][["geneSymbol"]])
			   )
    gene.list.test.gname <- llply(gene.list.test,names)
    names(gene.list.test.gname) <- names(gene.list.test)
    print(str(gene.list.test))

    hsGO <- llply(c("CC","MF","BP"),function(x){ GOSemSim::godata('org.Hs.eg.db', ont=x) },
			      .parallel=my.parallel,
			      .paropts=list(.packages=c("GOSemSim","clusterProfiler","AnnotationDbi"),
							    .export="org.Hs.eg.db"))
    names(hsGO) <- c("CC","MF","BP")

    out.prefix.plot <- sprintf("%s/clusterProfiler/%s",dirname(out.prefix),basename(out.prefix))
    dir.create(dirname(out.prefix.plot),F,T)
   
    th.qvalue <- 0.1
    th.showCategory <- 30
    font.size.cmp <- 10
    th.geneRatio <- c(0,0.1)
    pdf.width.dotplot <- NULL

    out.enrich.go <- list()
    for(gene.cate in names(gene.list.test))
    {
        out.enrich.go[[gene.cate]] <- llply(c("CC","MF","BP"),function(x){
                es.x <- gene.list.test[[gene.cate]]
                es.x[es.x > 2 ] <- 2
                run.clusterProfiler(sprintf("%s.%s.GO.%s",out.prefix.plot,gene.cate,x),
                            names(gene.list.test[[gene.cate]]),
                            gene.list.bg,
                            my.title=sprintf("%s (%s)",x,gene.cate),
                            db.name="GO",qvalueCutoff=0.1,
                            OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont=x,semData=hsGO[[x]],
                            es=es.x,
                            verbose=F,
                            pdf.width=c(9,7,7,18))
                    },.parallel=my.parallel,
                .paropts=list(.packages=c("clusterProfiler","org.Hs.eg.db","GOSemSim","AnnotationDbi")))
        names(out.enrich.go[[gene.cate]]) <- c("CC","MF","BP")
    }

    saveRDS(out.enrich.go,file=sprintf("%s.out.enrich.GO.rds",out.prefix.plot))
    #out.enrich.go <- readRDS(file=sprintf("%s.out.enrich.GO.rds",out.prefix.plot))
    cat(sprintf("enrich.GO done.\n"))

    #### GO comparison ####
    {
        out.cmp.go <- llply(c("CC","MF","BP"),function(x){
            ck <- compareCluster(geneCluster = gene.list.test.gname, fun = "enrichGO",
                         keyType="SYMBOL",ont=x,universe=gene.list.bg,
                         qvalueCutoff=th.qvalue, OrgDb='org.Hs.eg.db')
            ck.slim <- clusterProfiler::simplify(gofilter(ck,c(4,5)),cutoff=0.6, by="pvalue", select_fun=min,
                                measure="Jiang",semData=hsGO[[x]])
            ck.slim.fortify <- fortify(ck.slim, showCategory=30, by="geneRatio",
                                       includeAll=TRUE, split=NULL)
            ###write.table(ck.slim@compareClusterResult,file=sprintf("%s.cmp.GO.%s.txt",out.prefix.plot,x),
            ###			row.names=F,sep="\t",quote=F)
            #####
            ck.tmp <- compareCluster(geneCluster = gene.list.test.gname, fun = "enrichGO",
                         keyType="SYMBOL",ont=x,universe=gene.list.bg,
                         qvalueCutoff=Inf,pvalueCutoff=Inf, OrgDb='org.Hs.eg.db')
            enrich.ref <- clusterProfiler::enrichGO(gene.list.test.gname[[1]],OrgDb=org.Hs.eg.db,
                                universe=gene.list.bg,keyType="SYMBOL",ont=x,
                                pvalueCutoff=Inf,qvalueCutoff=Inf)
            enrich.ref.tmp <- clusterProfiler.dplyr::filter(enrich.ref,
                                    as.character(Description) %in%
                                        as.character(ck.slim.fortify$Description))
            enrich.ref.tmp <- clusterProfiler::simplify(gofilter(enrich.ref.tmp,c(4,5)),
                                    cutoff=0.6,by="pvalue",select_fun=min,
                                    measure="Jiang",semData=hsGO[[x]])
            ck.slim.plot <- clusterProfiler.dplyr::filter(ck.tmp,as.character(Description) %in%
                                    as.character(enrich.ref.tmp@result$Description) &
                                      p.adjust < 0.05 & qvalue<th.qvalue)
            write.table(ck.slim.plot@compareClusterResult,file=sprintf("%s.cmp.GO.%s.txt",out.prefix.plot,x),
                        row.names=F,sep="\t",quote=F)
            p1 <- dotplot(ck.slim.plot, showCategory=th.showCategory,font.size=font.size.cmp) +
                ggtitle(x) +
                scale_size(breaks=c(0.025,0.05,0.075),range=c(0.2,6),labels=c("2.5%","5%","7.5%"),
                           limits=th.geneRatio) +
                scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
            ggsave(sprintf("%s.cmp.GO.%s.dotplot.pdf",out.prefix.plot,x),
                   width=if(!is.null(pdf.width.dotplot)) pdf.width.dotplot else { if(length(gene.list.test)>2) 9.5 else 8 },
                   height=8,useDingbats=F)
            ####
            #p1 <- dotplot(ck.slim, showCategory=30,font.size=8) +
            #	ggtitle(x) +
            #	scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
            #ggsave(sprintf("%s.cmp.GO.%s.dotplot.pdf",out.prefix.plot,x),width=9,height=7)

            p2 <- NULL
            if(F){
                p2 <- emapplot(ck.slim,showCategory=th.showCategory,pie="count",
                               ##### version difference
                               #pie_scale=1.0,
                               layout="kk")
                ggsave(sprintf("%s.cmp.GO.%s.emaplot.pdf",out.prefix.plot,x),width=9,height=7)
            }

            return(list("plot"=list("p1"=p1,"p2"=p2),
                    "cmpResult"=ck.slim.plot))
        },.parallel=my.parallel)
        names(out.cmp.go) <- c("CC","MF","BP")
        saveRDS(out.cmp.go,file=sprintf("%s.out.cmp.GO.rds",out.prefix.plot))
        
        ggsave(sprintf("%s.cmp.GO.BP.dotplot.man.00.pdf",out.prefix.plot),
               plot=out.cmp.go$BP$plot$p1,
               width=9.0,
               height=9,useDingbats=F)

#        pp <- plot_grid(out.cmp.go$BP$plot$p1+theme(legend.position = "none"),
#                out.cmp.go$MF$plot$p1+theme(legend.position = "none"),
#                out.cmp.go$CC$plot$p1+theme(legend.position = "none"),
#                get_legend(out.cmp.go$BP$plot$p1),
#                #rel_widths=c(1,1,1,0.3),
#                axis = "l",
#                align="hv",nrow=2)
#        ggsave(sprintf("%s.cmp.GO.merge.dotplot.pdf",out.prefix.plot),width=15,height=15,useDingbats=F)

#        pp <- plot_grid(out.cmp.go$MF$plot$p1+theme(legend.position = "none"),
#                out.cmp.go$CC$plot$p1+theme(legend.position = "none"),
#                get_legend(out.cmp.go$BP$plot$p1),
#                rel_widths=c(1,1,0.3),
#                axis = "l",
#                align="hv",nrow=1)
#        ggsave(sprintf("%s.cmp.GO.merge.n2.dotplot.pdf",out.prefix.plot),width=15,height=7.5,useDingbats=F)
        ############################

    }

    ####

}


###### examples of regulons (generate files for Cytoscape, Fig. 3E)
{

    TF.list <- c("HIVEP1(+)")
    network.file <- "../data/expression/CD4/integration/scenic/comb.CD4.pyScenic.reg.tb.txt"
    importance.file <- "../data/expression/CD4/integration/scenic/comb.CD4.adj.rds"
    ### edge
    importance.tb <- readRDS(importance.file)
    network.tb <- fread(network.file)
    network.tb[is.na(w.zhangLab10X.CD4),w.zhangLab10X.CD4:=0]
    network.tb[is.na(w.zhangLabSS2.CD4),w.zhangLabSS2.CD4:=0]
    network.plot.tb <- network.tb[TF %in% TF.list,]
    network.plot.tb[,TF:=gsub("\\(\\+\\)","",TF)]
    network.plot.tb$hl <- "gray"
    network.plot.tb[target %in% c("TNFRSF4","TNFRSF9","VDR","IL21R","GATA3","CTNNB1","BCL3","ID3"),hl:="green"]
    network.plot.tb <- merge(network.plot.tb,importance.tb)
    ### node 
    g.info.tb <- gene.desc.top.slim[meta.cluster=="CD4.c20.Treg.TNFRSF9",]
    setkey(g.info.tb,"geneID")
    node.tb <- unique(network.plot.tb[,c("target","hl"),with=F][target!="HIVEP1",])
    colnames(node.tb)[1] <- "geneID"
    node.tb <- rbind(node.tb,data.table(geneID="HIVEP1",hl="red"))
    node.tb <- merge(node.tb,g.info.tb[node.tb$geneID,])
    node.tb[sig==T,]

    write.table(network.plot.tb,file=sprintf("%s.network.HIVEP1.txt",out.prefix),
		row.names=F,sep="\t",quote=F)
    write.table(node.tb,file=sprintf("%s.network.HIVEP1.node.txt",out.prefix),
		row.names=F,sep="\t",quote=F)

}



###################### immune checkpoints (Treg)
{
####    #gene.long.tb <- readRDS(sprintf("%s.geneTableLong.rds",in.prefix))
####    ### CD4 
####    gene.long.collapsed.CD4 <- readRDS(sprintf("%s.geneTableLong.collapsed.rds",in.prefix))
####    gene.desc.top.slim.CD4 <- readRDS(sprintf("%s.gene.desc.tb.slim.rds",in.prefix))
####
####    ### CD8
####    in.prefix.CD8 <- "OUT.core.signature/CD8.multiAsTwo.minCell50/panC.core.signature.CD8.multiAsTwo.minCell50"
####    gene.long.collapsed.CD8 <- readRDS(sprintf("%s.geneTableLong.collapsed.rds",in.prefix.CD8))
####    gene.desc.top.slim.CD8 <- readRDS(sprintf("%s.gene.desc.tb.slim.rds",in.prefix.CD8))
####
####    gene.long.collapsed <- rbind(gene.long.collapsed.CD4,gene.long.collapsed.CD8)
####    gene.desc.top.slim <- rbind(gene.desc.top.slim.CD4,gene.desc.top.slim.CD8)
####
####    doit <- function(mcls,gene.long.collapsed.tb,gene.desc.top.slim,aid=NULL,
####             pdf.width=8,pdf.height=7,
####             my.fontsize=12,my.k=2,palette.name="RdYlBu")
####    {
####
####        ES.range <- c(0,1)
####        ES.step <- 0.2
####        #palette.name <- "RdYlBu"
####        #pdf.width <- 8
####        #pdf.height <- 7
####        ##pdf.height <- 10
####
####        if(is.null(aid)){ aid <- mcls }
####        gene.IC.file <- "/lustre1/zeminz_pkuhpc/zhenglt/proj.database/database/DeSimone_2016.Treg.TableS5.IC.only.txt"
####        gene.IC.tb <- fread(gene.IC.file)
####        gene.IC.tb[1:2,]
####        gene.long.collapsed.tb[1,]
####        dat.plot.IC.tb <- gene.long.collapsed.tb[geneID %in% gene.IC.tb[["symbol"]] & meta.cluster %in% mcls,]
####        dat.plot.IC.tb[,sig.bool:=as.logical(sig)]
####        print(setdiff(gene.IC.tb[["symbol"]],dat.plot.IC.tb[["geneID"]]))
####        ## "BTNL2"   "VSIR"    "HHLA2"   "IDO2"    "TNFSF15" "TNFSF18" "VTCN1"
####
####        ## heatmap plot
####        {
####
####        dat.list <- llply(mcls,function(x) {
####            dat.plot.IC.mtx.tb <- dcast(dat.plot.IC.tb[meta.cluster==x,],geneID~cancerType,value.var="dprime")
####            dat.plot.IC.mtx.tb <- merge(dat.plot.IC.mtx.tb,
####                        gene.desc.top.slim[meta.cluster==x,c("geneID","comb.ES","sig")])
####            setkey(dat.plot.IC.mtx.tb,"geneID")
####
####            dat.plot.IC.mtx <- as.matrix(dat.plot.IC.mtx.tb[,-c("geneID","comb.ES","sig")])
####            rownames(dat.plot.IC.mtx) <- dat.plot.IC.mtx.tb[["geneID"]]
####
####            dat.plot.IC.mtx.sig.tb <- dcast(dat.plot.IC.tb[meta.cluster==x,],geneID~cancerType,value.var="sig.bool")
####            dat.plot.IC.mtx.sig <- as.matrix(dat.plot.IC.mtx.sig.tb[,-c("geneID")])
####            rownames(dat.plot.IC.mtx.sig) <- dat.plot.IC.mtx.sig.tb[["geneID"]]
####
####            #f.gene <- rowSums(dat.plot.IC.mtx.sig)==0
####            #f.gene <- rownames(dat.plot.IC.mtx.sig)[!f.gene]
####            #dat.plot.IC.mtx.sig <- dat.plot.IC.mtx.sig[f.gene,]
####            #dat.plot.IC.mtx <- dat.plot.IC.mtx[f.gene,]
####
####            ES.hclust.col <- run.cutree(t(dat.plot.IC.mtx),k=my.k,method.distance="cosine",method.hclust="ward.D2")
####            dat.plot.IC.mtx <- dat.plot.IC.mtx[,ES.hclust.col$hclust$order]
####            dat.plot.IC.mtx.sig <- dat.plot.IC.mtx.sig[,colnames(dat.plot.IC.mtx)]
####            
####            ## white
####            print("all(rownames(dat.plot.IC.mtx.sig)==rownames(dat.plot.IC.mtx))")
####            print(all(rownames(dat.plot.IC.mtx.sig)==rownames(dat.plot.IC.mtx)))
####            print("all(colnames(dat.plot.IC.mtx.sig)==colnames(dat.plot.IC.mtx))")
####            print(all(colnames(dat.plot.IC.mtx.sig)==colnames(dat.plot.IC.mtx)))
####            dat.plot.IC.mtx.white <- dat.plot.IC.mtx
####            #dat.plot.IC.mtx.white[dat.plot.IC.mtx < 0.15 ] <- 0
####            dat.plot.IC.mtx.white[!dat.plot.IC.mtx.sig] <- 0
####            
####            dat.plot.ann <- structure(dat.plot.IC.mtx.tb[["comb.ES"]],
####                          names=dat.plot.IC.mtx.tb[["geneID"]])
####            dat.plot.ann[ dat.plot.IC.mtx.tb[["sig"]]==F ] <- 0
####
####            return(list(dat.plot.IC.mtx=dat.plot.IC.mtx,
####                dat.plot.IC.mtx.sig=dat.plot.IC.mtx.sig,
####                dat.plot.IC.mtx.white=dat.plot.IC.mtx.white,
####                dat.plot.ann=dat.plot.ann))
####        })
####        names(dat.list) <- mcls
####
####        #### check rownames
####        print("dat.list[[x]]$dat.plot.IC.mtx have the same rownames:")
####        all(apply(t(ldply(names(dat.list),function(x){
####              (rownames(dat.list[[x]]$dat.plot.IC.mtx))
####        })),1,function(xx){ all(xx==xx[1]) }))
####
####        dat.for.clust <- do.call(cbind,llply(names(dat.list),function(x){ dat.list[[x]]$dat.plot.IC.mtx } ))
####        dat.for.plot <- do.call(cbind,llply(names(dat.list),function(x){ dat.list[[x]]$dat.plot.IC.mtx.white } ))
####
####        f.zero <- rowSums(dat.for.plot>0)==0
####        dat.for.clust <- dat.for.clust[!f.zero,]
####        dat.for.plot <- dat.for.plot[!f.zero,]
####
####        ES.hclust.row <- run.cutree(dat.for.clust,k=2,method.distance="cosine",method.hclust="ward.D2")
####        gene.order <- rownames(dat.for.clust)[ES.hclust.row$hclust$order]
####        gene.use <- rownames(dat.for.clust)
####        #ES.hclust.col <- run.cutree(t(dat.plot.IC.mtx),k=2,method.distance="cosine",method.hclust="ward.D2")
####
####        a.palette <- sscVis:::getColorPaletteFromNameContinuous(palette.name)
####        z.len <- length(a.palette)
####
####        ht.list <- llply(seq_along(dat.list),function(i){
####            x <- names(dat.list)[i]
####            ann.gene.ES <- rowAnnotation(panC=dat.list[[x]]$dat.plot.ann[gene.use],
####                          show_legend=F,
####                          col=list(panC=circlize::colorRamp2(seq(ES.range[1],ES.range[2], length = z.len),
####                                      colorRampPalette(a.palette)(z.len))),
####                          annotation_legend_param=list(at = seq(ES.range[1],ES.range[2],ES.step),
####                                       grid_width = unit(0.4,"cm"),
####                                       grid_height = unit(0.4, "cm"),
####                                       legend_height = unit(6, "cm")),
####                          border=T)
####            ### whiten
####            ht <- sscVis::plotMatrix.simple(dat.list[[x]]$dat.plot.IC.mtx.white[gene.use,],
####                          returnHT=T,
####                          out.prefix=NULL,
####                          #out.prefix=sprintf("%s.IC.%s.%s.00",out.prefix,mcls,palette.name),
####                          mytitle=sprintf("%s",x),
####                          #row.split=gene.spe.tb$ES.cancerType,
####                          #row_gap = unit(0, "mm"),
####                          #row_title_gp = gpar(fontsize = 10),
####                          clust.row=ES.hclust.row$branch,
####                          #clust.column=ES.hclust.col$branch,
####                          clust.column=F,
####                          #par.warterfall=list(method.distance="cosine",k=2),
####                          #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
####                          palatte=a.palette,
####                          #z.lo=-0.3,z.hi=0.6,
####                          #par.legend=list(at = seq(-0.3,0.6,0.3)),
####                          z.lo=ES.range[1],z.hi=ES.range[2],
####                          par.legend=list(at = seq(ES.range[1],ES.range[2],ES.step)),
####                          pdf.width=pdf.width,
####                          pdf.height=pdf.height,
####                          #show_row_dend=T,
####                          show.dendrogram=if(i==1) T else F,
####                          row_dend_width = unit(4.0, "cm"),
####                          par.heatmap=list(
####                                   column_title=fetchMetaClusterID2CusterFullName()[x],
####                                   show_heatmap_legend = if(i==1) T else F,
####                                   column_names_gp=gpar(fontsize=my.fontsize),
####                                   row_names_gp=gpar(fontsize=my.fontsize),
####                                   #cex.column=1.2,cex.row=1.5,
####                                   #right_annotation =  ha,
####                                   #left_annotation = ann.gene,
####                                   #show_column_dend=T,
####                                   column_dend_height = unit(2.0, "cm"),
####                                   left_annotation = ann.gene.ES,
####                                   row_title_rot=0,border=T
####                                           ),
####                          exp.name="Effect Size")
####            return(ht)
####
####        })
####        names(ht.list) <- names(dat.list)
####
####        ht.merge <- NULL
####        for(x in names(ht.list)){
####            ht.merge <- ht.merge + ht.list[[x]]
####        }
####
####        pdf(sprintf("%s.IC.%s.%s.00.pdf",out.prefix,aid,palette.name),width=pdf.width,height=pdf.height)
####        opar <- par(mar = c(4, 2, 2, 4))
####        plot.new()
####            title(main = "", cex.main = 2)
####            vps <- gridBase::baseViewports()
####            pushViewport(vps$inner, vps$figure, vps$plot)
####        ComplexHeatmap::draw(ht.merge, newpage = FALSE, merge_legends = TRUE)
####        dev.off()
####
####        }
####
####    }
####
####    doit("CD4.c20.Treg.TNFRSF9",gene.long.collapsed,gene.desc.top.slim,pdf.width=7.5,pdf.height=8.5)
####    doit("CD4.c17.TfhTh1.CXCL13",gene.long.collapsed,gene.desc.top.slim,pdf.width=7.5,pdf.height=9.0)
####    doit("CD8.c12.Tex.CXCL13",gene.long.collapsed,gene.desc.top.slim,pdf.width=7.5,pdf.height=9.0)
####
####    doit(mcls=c("CD4.c20.Treg.TNFRSF9","CD4.c17.TfhTh1.CXCL13","CD8.c12.Tex.CXCL13"),
####         gene.long.collapsed,gene.desc.top.slim,
####         aid="merge3",pdf.width=12,pdf.height=8,my.fontsize=10,palette.name="plasma")
####
####
####    ########### binarized expression frequency of IC
####    ######## gene expression data ####
####    sce.list.CD8.file <- "list/obj.inte.all.CD8.post.limma.list"
####    sce.list.CD8.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.CD8.file),head=F)
####    colnames(sce.list.CD8.tb) <- c("stype","aid", "measurement", "platform", "scefile", "defile")
####
####    sce.list.CD4.file <- "list/obj.inte.all.CD4.post.limma.list"
####    sce.list.CD4.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.CD4.file),head=F)
####    colnames(sce.list.CD4.tb) <- c("stype","aid", "measurement", "platform", "scefile", "defile")
####
####    RhpcBLASctl::omp_set_num_threads(1)
####    doParallel::registerDoParallel(cores = 8)
####
####    sce.CD8.list <- llply(seq_len(nrow(sce.list.CD8.tb)),function(i){ readRDS(sce.list.CD8.tb$scefile[i]) },.parallel=T)
####    names(sce.CD8.list) <- sce.list.CD8.tb$aid
####
####    sce.CD4.list <- llply(seq_len(nrow(sce.list.CD4.tb)),function(i){ readRDS(sce.list.CD4.tb$scefile[i]) },.parallel=T)
####    names(sce.CD4.list) <- sce.list.CD4.tb$aid
####
####
####    calBinExp <- function(sce.list,
####                  gene.to.test=c("FOXP3","CXCL13"),
####                  dataset.x="OV.zhangLab5P",
####                  mcls="CD8.c12.Tex.CXCL13",
####                  th.exprs=0.3, min.ncell=50,val.padding=0)
####    {
####
####        sce <- sce.list[[dataset.x]]
####        
####        #### filter out dataset with cell number in Tex less than 50
####        if( (is.null(mcls) && ncol(sce) < min.ncell  ) ||
####            (!is.null(mcls) && sum(sce[["meta.cluster"]]==mcls) < min.ncell) )
####        {
####            return(NULL)
####        }
####        ####
####        cat(sprintf("dataset (%s)\n",dataset.x))
####
####        sce.plot <- ssc.scale(sce,gene.symbol=gene.to.test,adjB="batchV",do.scale=T)
####        if(!is.null(mcls)){
####        sce.plot <- sce.plot[,sce.plot$meta.cluster==mcls]
####        }
####
####        dat.plot <- as.data.frame(t(as.matrix(assay(sce.plot,"norm_exprs.scale"))))
####        colnames(dat.plot) <- rowData(sce.plot)$display.name
####        col.dropout <- setdiff(gene.to.test,colnames(dat.plot))
####        if(length(col.dropout)>0){
####        dropout.mtx <- matrix(rep(val.padding,nrow(dat.plot)*length(col.dropout)),ncol=length(col.dropout))
####        colnames(dropout.mtx) <- col.dropout
####        dat.plot <- cbind(dat.plot,dropout.mtx)
####        }
####        dat.plot <- dat.plot[,gene.to.test]
####
####        f.na <- is.na(dat.plot)
####        dat.plot[f.na] <- val.padding
####
####        for(gg in colnames(dat.plot)){
####        sce.plot[[sprintf("bin.%s",gg)]] <- as.integer(dat.plot[[gg]] > th.exprs)
####        sce.plot[[sprintf("z.%s",gg)]] <- dat.plot[[gg]]
####        }
####
####        ret.tb <- as.data.table(colData(sce.plot))[,c("cellID","cancerType","dataset","patient",
####                              "cellID.uniq","miniCluster","meta.cluster",
####                              "loc","libraryID","batchV",
####                              sprintf("bin.%s",colnames(dat.plot)),
####                              sprintf("z.%s",colnames(dat.plot))
####                              ),with=F]
####
####        #cat(sprintf("successful (%s)\n",dataset.x))
####
####        return(ret.tb)
####
####    }
####
####
####    gene.in.box <- c("CXCL13","HAVCR2","FOXP3",
####             "NTRK1",
####             "CD28","CD40","CD44","IDO1",
####             "CTLA4","CD80","CD86",
####             "PDCD1","PDCD1LG2","CD274",
####             "CD276","TNFRSF14","ADORA2A")
####    dataset.ana.vec <- names(sce.CD4.list)
####
####    ####### check distribution and get threshold
####    {
####        dat.bin.exp.check.CD4.tb <- as.data.table(ldply(dataset.ana.vec,function(x){
####                              calBinExp(sce.CD4.list,
####                                 gene.to.test=gene.in.box,
####                                 dataset.x=x,mcls=NULL)
####                    }))
####        dat.bin.exp.check.CD4.tb <- changeSomeNames(dat.bin.exp.check.CD4.tb)
####        dat.bin.exp.check.CD4.tb <- correctCellInfo(dat.bin.exp.check.CD4.tb)
####
####        dat.bin.exp.check.CD8.tb <- as.data.table(ldply(dataset.ana.vec,function(x){
####                              calBinExp(sce.CD8.list,
####                                 gene.to.test=gene.in.box,
####                                 dataset.x=x,mcls=NULL)
####                    }))
####        dat.bin.exp.check.CD8.tb <- changeSomeNames(dat.bin.exp.check.CD8.tb)
####        dat.bin.exp.check.CD8.tb <- correctCellInfo(dat.bin.exp.check.CD8.tb)
####
####        make.zscore.densityPlot <- function(aid,dat.bin.exp.check.tb,out.prefix,gene.in.box,th.tb=NULL)
####        {
####            #a.gene <- "FOXP3"
####            l_ply(gene.in.box,function(a.gene){
####                    p <- ggplot(dat.bin.exp.check.tb,aes_string(sprintf("z.%s",a.gene))) +
####                        geom_density(aes(color=dataset,group=dataset,fill=dataset),alpha=0.5,size=0.5) +
####                        scale_fill_manual(values=colSet$dataset) +
####                        scale_color_manual(values=colSet$dataset) +
####                        geom_vline(xintercept=c(0.3),linetype="dashed",color="lightgray",alpha=0.5)
####                    if(is.null(th.tb)){
####                    p <- p + geom_vline(xintercept=c(0.3),linetype="dashed",color="red")
####                    }else{
####                    p <- p + geom_vline(aes(xintercept=TH.z),
####                            data=th.tb[gene==a.gene & stype==aid,],
####                            linetype="dashed",color="red")
####                    }
####                    p <- p + xlab("z-score") +
####                        #facet_wrap(~dataset,ncol=5,scales="free_y") +
####                        facet_wrap(~dataset,ncol=5,scales="free") +
####                        theme_pubr() +
####                        theme(legend.position="none")
####                    ggsave(sprintf("%s.IC.gene.zscore.%s.%s.pdf",out.prefix,a.gene,aid),width=10,height=12)
####                 },.parallel=T)
####
####        }
####
####        out.prefix.bin.exp <- sprintf("%s/bin.exp/%s",dirname(out.prefix),basename(out.prefix))
####        dir.create(dirname(out.prefix.bin.exp),F,T)
####
####        thres.tb <- as.data.table(expand.grid(dataset=sort(unique(dat.bin.exp.check.CD4.tb$dataset)),
####                              gene=gene.in.box))
####        thres.tb$stype <- "CD4"
####        thres.tb$TH.z <- 0.3
####        thres.tb[dataset %in% c("PACA.thisStudy",
####                    "THCA.thisStudy","UCEC.thisStudy", "OV.thisStudy",
####                    "CRC.LeiZhang2018","NSCLC.XinyiGuo2018","HCC.QimingZhang2019.SS2",
####                    "BCC.KathrynEYost2019","BRCA.PeterSavas2018",
####                    "SCC.KathrynEYost2019","HNSCC.SidharthVPuram2017","NPC.YangLiu2020") & gene=="FOXP3",TH.z:=0]
####        thres.tb[dataset %in% c("FTC.thisStudy") & gene=="FOXP3",TH.z:=0.5]
####        thres.tb[dataset %in% c("ESCA.thisStudy","BRCA.PeterSavas2018") & gene=="FOXP3",TH.z:=-0.3]
####        thres.tb[gene=="FOXP3",]
####        sort(unique(dat.bin.exp.check.CD4.tb$dataset))
####
####        thres.CD8.tb <- as.data.table(expand.grid(dataset=sort(unique(dat.bin.exp.check.CD8.tb$dataset)),
####                              gene=gene.in.box))
####        thres.CD8.tb$stype <- "CD8"
####        thres.CD8.tb$TH.z <- 0.3
####
####        make.zscore.densityPlot(aid="CD4",dat.bin.exp.check.CD4.tb,out.prefix.bin.exp,gene.in.box,th.tb=thres.tb)
####
####        make.zscore.densityPlot(aid="CD8",dat.bin.exp.check.CD8.tb,out.prefix.bin.exp,gene.in.box,th.tb=thres.CD8.tb)
####
####
####        ########  classsification
####        dat.bin.exp.CD4.tb <- as.data.table(ldply(c("CD4.c20.Treg.TNFRSF9","CD4.c17.TfhTh1.CXCL13"),function(aid){
####            as.data.table(ldply(dataset.ana.vec,function(x){
####                              calBinExp(sce.CD4.list,
####                                 gene.to.test=gene.in.box,
####                                 dataset.x=x,mcls=aid)
####                    }))
####        }))
####
####        dat.bin.exp.CD8.tb <- as.data.table(ldply(c("CD8.c12.Tex.CXCL13"),function(aid){
####            as.data.table(ldply(dataset.ana.vec,function(x){
####                              calBinExp(sce.CD8.list,
####                                 gene.to.test=gene.in.box,
####                                 dataset.x=x,mcls=aid)
####                    }))
####        }))
####
####        dat.bin.exp.tb <- rbind(dat.bin.exp.CD4.tb,dat.bin.exp.CD8.tb)
####        dat.bin.exp.tb[1,]
####        dat.bin.exp.tb[,.N,by=c("meta.cluster")]
####        dat.bin.exp.tb <- changeSomeNames(dat.bin.exp.tb)
####        dat.bin.exp.tb <- correctCellInfo(dat.bin.exp.tb)
####
####        saveRDS(dat.bin.exp.tb,file=sprintf("%s.IC.gene.bin.exp.tb.rds",out.prefix))
####
####
####        make.bin.exp.boxplot <- function(dat.bin.exp.tb,out.prefix.bin.exp,pmode=NULL,use.onlyT=T)
####        {
####            dat.bin.exp.long.tb <- as.data.table(ldply(gene.in.box,function(a.gene){
####                    #a.gene <- "PDCD1LG2"
####                    dat.plot <- dat.bin.exp.tb[loc %in% (if(use.onlyT) c("T") else c("P","N","T","L")),
####                                   .(NPos=sum(.SD[[sprintf("bin.%s",a.gene)]]),
####                                     NTotal=.N),
####                                   by=c("cancerType","dataset","meta.cluster","patient")]
####                    dat.plot[,FreqPos:=NPos/NTotal]
####                    dat.plot[,gene:=a.gene]
####                    dat.plot <- dat.plot[ NTotal >= 10,]
####                    cancerType.order.tb <- dat.plot[meta.cluster=="CD4.c20.Treg.TNFRSF9",
####                                    .(medV=median(FreqPos)),
####                                    by="cancerType"][order(medV),]
####                    dat.plot[,cancerType:=factor(cancerType,levels=cancerType.order.tb$cancerType)]
####                    dat.plot[,cluster.name:=fetchMetaClusterID2CusterFullName()[as.character(meta.cluster)]]
####
####                    if(is.null(pmode))
####                    {
####                        p <- ggboxplot(dat.plot,
####                               x="cancerType",y="FreqPos",fill="cancerType",
####                               xlab="",ylab="Frequency Of Positive Cells",title=a.gene) +
####                            scale_fill_manual(values=colSet$cancerType) +
####                            geom_hline(yintercept=0.05,color="lightgray",linetype="dashed",alpha=0.8) +
####                            stat_compare_means(label="p.format") +
####                            facet_wrap(~cluster.name,ncol=1,
####                                   scales="free_y") +
####                            theme(axis.text.x=element_text(angle=60,hjust=1),
####                              strip.text=element_text(size=12),
####                              plot.title=element_text(hjust=0.5),
####                              legend.position="none")
####                        ggsave(sprintf("%s.IC.gene.bin.exp.tb.%s.pdf",out.prefix.bin.exp,a.gene),width=5,height=6,useDingbats=F)
####                    }else{
####                        mcls.vec <- unique(dat.plot$meta.cluster)
####                        p.list <- llply(mcls.vec,function(x){
####                                a.title <- fetchMetaClusterID2CusterFullName()[as.character(x)]
####                                dat.x <- dat.plot[meta.cluster==x,]
####                                cancerType.order.tb <- dat.x[,
####                                                .(medV=median(FreqPos)),
####                                                by="cancerType"][order(medV),]
####                                dat.x[,cancerType:=factor(as.character(cancerType),
####                                              levels=cancerType.order.tb$cancerType)]
####                                p <- ggboxplot(dat.x,
####                                           x="cancerType",y="FreqPos",fill="cancerType",
####                                           xlab="",ylab="Frequency Of Positive Cells",title=a.title) +
####                                    scale_fill_manual(values=colSet$cancerType) +
####                                    geom_hline(yintercept=0.05,color="lightgray",linetype="dashed",
####                                           alpha=0.8) +
####                                    stat_compare_means(label="p.format") +
####                                    theme(axis.text.x=element_text(angle=60,hjust=1),
####                                          strip.text=element_text(size=12),
####                                          plot.title=element_text(hjust=0.5),
####                                          legend.position="none")
####                                return(p)
####                                })
####                        pp <- cowplot::plot_grid(plotlist=p.list,align="hv",
####                                     labels=a.gene,hjust=0,
####                                     nrow=1)
####                        ggsave(sprintf("%s.IC.gene.bin.exp.tb.%s.pdf",out.prefix.bin.exp,a.gene),
####                           width=13,height=3.2,useDingbats=F)
####                    }
####
####                    return(dat.plot)
####                     },.parallel=T))
####        }
####
####        make.bin.exp.boxplot(dat.bin.exp.tb,out.prefix.bin.exp,pmode="freeX")
####
####    }
####
}



