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
source("./func.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)

out.prefix <- "OUT_Fig1/Fig1"
dir.create(dirname(out.prefix),F,T)

# load data
{
    dir.int.list <- list("CD8"="../data/expression/CD8/integration",
                         "CD4"="../data/expression/CD4/integration")
    sce.merged.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds"),
                            "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.sce.merged.rds"))
    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    colSet.CD8 <- g.colSet["meta.cluster"]
    colSet.CD8$meta.cluster <- colSet.CD8$meta.cluster[ grepl("^CD8\\.c",names(colSet.CD8$meta.cluster))]
    colSet.CD4 <- g.colSet["meta.cluster"]
    colSet.CD4$meta.cluster <- colSet.CD4$meta.cluster[ grepl("^CD4\\.c",names(colSet.CD4$meta.cluster))]
    colSet.list <- list("CD8"=colSet.CD8,"CD4"=colSet.CD4)

    gene.desc.tb.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds"),
                              "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.gene.tb.rds"))
}


###### UMAP
# Fig. 1 B and C
{

    p.multi <- ssc.plot.tsne(sce.merged.list[["CD8"]],columns = c("meta.cluster"),reduced.name = "harmony.umap",
                              vector.friendly=T,legend.w=1.2,
                              theme.use=theme_pubr,verbose=T,
                              size=1.0,
                              label = 2,
                              par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                              par.geom_point=list(scale=0.6),
                              par.geneOnTSNE = list(pt.order = "random"),
                              colSet=colSet.CD8)
    p <- p.multi$list[[1]] + theme(legend.position="right")
    #print(p)
    #ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=9,height=5,useDingbats=FALSE)
    ggsave(sprintf("%s.meta.cluster.CD8.vecFri.pdf",out.prefix),width=6.5,height=4.25,useDingbats=FALSE)

}

{
    p.multi <- ssc.plot.tsne(sce.merged.list[["CD4"]],columns = c("meta.cluster"),reduced.name = "harmony.umap",
                              vector.friendly=T,legend.w=1.2,
                              theme.use=theme_pubr,verbose=T,
                              size=1.0,
                              label = 2,
                              par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                              par.geom_point=list(scale=0.6),
                              par.geneOnTSNE = list(pt.order = "random"),
                              colSet=colSet.CD4)
    p <- p.multi$list[[1]] + theme(legend.position="right")
    #print(p)
    #ggsave(sprintf("%s.meta.cluster.CD4.vecFri.pdf",out.prefix),width=10.5,height=5,useDingbats=FALSE)
    ggsave(sprintf("%s.meta.cluster.CD4.vecFri.pdf",out.prefix),width=8.0,height=4.25,useDingbats=FALSE)

}


#### compare steps to evaluate integration
# UMAP (fig.S02 C to D ??)
{
    ###sce.merged.CD8 <- changeSomeNames(sce.merged.CD8)

    lim.list <- list("umap"=list("xlim"=c(-5,6.5),
				 "ylim"=c(-6,6)),
		     "harmony.umap"=list("xlim"=c(-7,7),
				 "ylim"=c(-7.5,7.5)))

    sce.merged <- sce.merged.list[["CD8"]]
    
    p <- ssc.plot.tsne(sce.merged,columns="dataset",reduced.name="umap",
                       colSet=g.colSet,verbose=T)
    pp <- p$list[[1]]
    pp <- pp + guides(colour=guide_legend(override.aes = list(size=4),
                                          label.theme = element_text(size=8),ncol=4))
    p.legend <- get_legend(pp)
    pdf(sprintf("%s.dataset.CD8.legend.v1.pdf",out.prefix),width=8,height=3,useDingbats = F)
    grid.newpage()
    grid.draw(p.legend)
    dev.off()


    ##### cancer type examples
    {

        cancerType.multDatasets.vec <- c("BRCA","LC","MELA","CRC","PACA","RC")
        ctype.name <- cancerType.multDatasets.vec
        l_ply(c("umap","harmony.umap"),function(rd){
              plist.byCancerType <- llply(ctype.name,function(ctype){
                      p <- ssc.plot.tsne(sce.merged[,sce.merged$cancerType==ctype],columns="dataset",
                             reduced.name=rd,colSet=g.colSet,
                             vector.friendly=T,legend.w=0,
                             xlim=lim.list[[rd]][["xlim"]],ylim=lim.list[[rd]][["ylim"]],
                             size=2.5,
                             label=2.0,
                             par.geom_point=list(scale=0.2),
                             fun.extra=function(p){ p + labs(title = sprintf("datasets(%s)",ctype),
                                             x=sprintf("%s1",toupper(rd)),
                                             y=sprintf("%s2",toupper(rd))) +
                                        theme(axis.title=element_text(size=12),
                                              axis.text=element_text(size=10)) },
                             #theme.use=theme_void,
                             theme.use=theme_classic,
                             par.geneOnTSNE = list(pt.order="random"))
                      #ggsave(filename=sprintf("%s.%s.figS.byCancerType.miniC.%s.v0.pdf",out.prefix,rd,ctype),
                      #       width=4,height=4.2)
                      return(p)
                     },.parallel=T)
            names(plist.byCancerType) <- ctype.name
            saveRDS(plist.byCancerType,file=sprintf("%s.%s.figS.byCancerType.miniC.plist.CD8.rds",out.prefix,rd))
        })


        ### version2
        l_ply(c("umap","harmony.umap"),function(rd){
            plist.byCancerType <- readRDS(sprintf("%s.%s.figS.byCancerType.miniC.plist.CD8.rds",out.prefix,rd))
            p.list <- list()
            p.list[["dataset"]] <- ssc.plot.tsne(sce.merged,columns="dataset",
                         reduced.name=rd,colSet=colSet.list,
                         vector.friendly=T,legend.w=0,
                         xlim=lim.list[[rd]][["xlim"]],ylim=lim.list[[rd]][["ylim"]],
                         ##label=2.0,
                         par.geom_point=list(scale=0.25),
                         fun.extra=function(p){ p + labs(title = sprintf("datasets(all)"),
                                         x=sprintf("%s1",toupper(rd)),
                                         y=sprintf("%s2",toupper(rd))) +
                                        theme(axis.title=element_text(size=12),
                                          axis.text=element_text(size=10)) },
                         #theme.use=theme_void,
                         theme.use=theme_classic,
                         par.geneOnTSNE = list(pt.order="random"))
            p.list[["dataset(MELA)"]] <- plist.byCancerType[["MELA"]]
            p.list[["gene"]] <- ssc.plot.tsne(sce.merged,gene=c("TCF7","CX3CR1","CXCL13","TYROBP","SLC4A10"),
                              reduced.name=rd,
                              vector.friendly=T,clamp=c(-0.5,1.5),p.ncol=5,
                              size=1.5,
                              par.geom_point=list(scale=0.25),
                              xlim=lim.list[[rd]][["xlim"]],ylim=lim.list[[rd]][["ylim"]],
                              fun.extra=function(p){ p + labs(
                                         x=sprintf("%s1",toupper(rd)),
                                         y=sprintf("%s2",toupper(rd))) +
                                        theme(axis.title=element_text(size=12),
                                          axis.text=element_text(size=10)) },
                              #theme.use=theme_void,
                              theme.use=theme_classic,
                              par.geneOnTSNE=list(scales="fixed",pt.order="random",pt.alpha = 0.5))
            pp <- cowplot::plot_grid(plotlist=p.list,nrow = 1,align = "hv",rel_widths=c(2.8,2.8,16))
            ggsave(filename=sprintf("%s.%s.figS.CD8.v2.pdf",out.prefix,rd),width=16,height=2.2)
        },.parallel=T)



    }

}

# UMAP (fig.S03 C to D ??)
{

    lim.list <- list("umap"=list("xlim"=c(-7.5,6),
                                 "ylim"=c(-6,6)),
                     "harmony.umap"=list("xlim"=c(-8,8),
                                 "ylim"=c(-6,6)))

    sce.merged <- sce.merged.list[["CD4"]]
    
    p <- ssc.plot.tsne(sce.merged,columns="dataset",reduced.name="umap",
                       colSet=g.colSet,verbose=T)
    pp <- p$list[[1]]
    pp <- pp + guides(colour=guide_legend(override.aes = list(size=4),
                                          label.theme = element_text(size=8),ncol=4))
    p.legend <- get_legend(pp)
    pdf(sprintf("%s.dataset.CD4.legend.v1.pdf",out.prefix),width=8,height=3,useDingbats = F)
    grid.newpage()
    grid.draw(p.legend)
    dev.off()


    ##### cancer type examples
    {

      cancerType.multDatasets.vec <- c("BRCA","LC","MELA","CRC","PACA","RC")
	  ctype.name <- cancerType.multDatasets.vec
	  l_ply(c("umap","harmony.umap"),function(rd){
	      plist.byCancerType <- llply(ctype.name,function(ctype){
		      p <- ssc.plot.tsne(sce.merged[,sce.merged$cancerType==ctype],columns="dataset",
				     reduced.name=rd,colSet=g.colSet,
				     vector.friendly=T,legend.w=0,
				     xlim=lim.list[[rd]][["xlim"]],ylim=lim.list[[rd]][["ylim"]],
				     size=2.5,
				     label=2.0,
				     par.geom_point=list(scale=0.2),
				     fun.extra=function(p){ p + labs(title = sprintf("datasets(%s)",ctype),
								     x=sprintf("%s1",toupper(rd)),
								     y=sprintf("%s2",toupper(rd))) +
								theme(axis.title=element_text(size=12),
								      axis.text=element_text(size=10)) },
				     #theme.use=theme_void,
				     theme.use=theme_classic,
				     par.geneOnTSNE = list(pt.order="random"))
		      #ggsave(filename=sprintf("%s.%s.figS.byCancerType.miniC.%s.v0.pdf",out.prefix,rd,ctype),
		      #       width=4,height=4.2)
		      return(p)
				 },.parallel=T)
	    names(plist.byCancerType) <- ctype.name
	    saveRDS(plist.byCancerType,file=sprintf("%s.%s.figS.byCancerType.miniC.plist.CD4.rds",out.prefix,rd))
	  })

	  ### version2
	  l_ply(c("umap","harmony.umap"),function(rd){
	    plist.byCancerType <- readRDS(sprintf("%s.%s.figS.byCancerType.miniC.plist.CD4.rds",out.prefix,rd))
	    p.list <- list()
	    p.list[["dataset"]] <- ssc.plot.tsne(sce.merged,columns="dataset",
					 reduced.name=rd,colSet=colSet.list,
					 vector.friendly=T,legend.w=0,
					 xlim=lim.list[[rd]][["xlim"]],ylim=lim.list[[rd]][["ylim"]],
					 ##label=2.0,
				     par.geom_point=list(scale=0.25),
					 fun.extra=function(p){ p + labs(title = sprintf("datasets(all)"),
									 x=sprintf("%s1",toupper(rd)),
									 y=sprintf("%s2",toupper(rd))) +
								    theme(axis.title=element_text(size=12),
									  axis.text=element_text(size=10)) },
					 #theme.use=theme_void,
					 theme.use=theme_classic,
					 par.geneOnTSNE = list(pt.order="random"))
	    p.list[["dataset(MELA)"]] <- plist.byCancerType[["MELA"]]
	    p.list[["gene"]] <- ssc.plot.tsne(sce.merged,gene=c("TCF7","GZMK","CX3CR1","CXCL13","FOXP3"),
					      reduced.name=rd,
					      vector.friendly=T,clamp=c(-0.5,1.5),p.ncol=5,
					      size=1.5,
				          par.geom_point=list(scale=0.25),
					      xlim=lim.list[[rd]][["xlim"]],ylim=lim.list[[rd]][["ylim"]],
					      fun.extra=function(p){ p + labs(
									 x=sprintf("%s1",toupper(rd)),
									 y=sprintf("%s2",toupper(rd))) +
								    theme(axis.title=element_text(size=12),
									  axis.text=element_text(size=10)) },
					      #theme.use=theme_void,
					      theme.use=theme_classic,
					      par.geneOnTSNE=list(scales="fixed",pt.order="random",pt.alpha = 0.5))
	    pp <- cowplot::plot_grid(plotlist=p.list,nrow = 1,align = "hv",rel_widths=c(2.8,2.8,16))
	    ggsave(filename=sprintf("%s.%s.figS.CD4.v2.pdf",out.prefix,rd),width=16,height=2.2)
			     },.parallel=T)

    }
}

# LISI index
{

    #### calculate lisi index
    for(stype in c("CD8","CD4"))
    {    
        if(!file.exists(sprintf("%s/int.%s.S35.lisi.RData",dir.int.list[[stype]],stype)))
        {
            sce.merged <- sce.merged.list[[stype]]
            lisi.pca <- lisi::compute_lisi(reducedDim(sce.merged,"pca"),
                                           colData(sce.merged),c("dataset","dataset.tech","ClusterID.pca"))
            lisi.harmony <- lisi::compute_lisi(reducedDim(sce.merged,"harmony"),
                                               colData(sce.merged),c("dataset","dataset.tech","ClusterID.harmony"))

            save(lisi.pca,lisi.harmony,file=sprintf("%s/int.%s.S35.lisi.RData",dir.int.list[[stype]],stype))
        }
    }

    ############
    for(stype in c("CD8","CD4"))
    {

        lisi.sc.raw <- loadToEnv(sprintf("../data/expression/%s/integration/sc/int.harmony.%s.raw.S35.lisi.RData",stype,stype))
        lisi.sc <- loadToEnv(sprintf("../data/expression/%s/integration/sc/int.harmony.%s.S35.lisi.RData",stype,stype))
        load(sprintf("%s/int.%s.S35.lisi.RData",dir.int.list[[stype]],stype))

        lisi.pca.sc.raw.tb <- cbind(data.table(cellID=rownames(lisi.sc.raw$lisi.pca),
                           rd="PCA.sc.raw",step="A"),
                 lisi.sc.raw$lisi.pca[,c("dataset"),drop=F])
        lisi.pca.sc.tb <- cbind(data.table(cellID=rownames(lisi.sc$lisi.pca),
                           rd="PCA.sc",step="B"),
                 lisi.sc$lisi.pca[,c("dataset"),drop=F])
        lisi.pca.tb <- cbind(data.table(cellID=rownames(lisi.pca),
                        rd="PCA",step="C"),
                 lisi.pca[,c("dataset"),drop=F])
        lisi.harmony.tb <- cbind(data.table(cellID=rownames(lisi.harmony),
                        rd="Harmony",step="D"),
                 lisi.harmony[,c("dataset"),drop=F])
        lisi.merge.tb <- rbind(lisi.pca.sc.raw.tb,lisi.pca.sc.tb,lisi.pca.tb,lisi.harmony.tb)
        lisi.merge.tb[,.(mean(dataset)),by="step"]
        lisi.merge.tb[,.(median(dataset)),by="step"]

        p <- ggboxplot(lisi.merge.tb,x="step",y="dataset",
               fill="step",alpha=0.8) +
            stat_compare_means(comparisons=list(c("A","B"),c("B","C"),c("C","D"),c("B","D"))) +
            ylab("LISI") +
            theme(legend.position="right")
        ggsave(sprintf("%s.LISI.dataset.merge.%s.01.pdf",out.prefix,stype),width=3.2,height=4)

    }

}

#######################

# dot plot showing effect size of signature genes
{

    #### all CD8 clusters
    {

        gene.desc.top <- gene.desc.tb.list[["CD8"]]
        gene.plot.example.tb <- fread("../data/expression/CD8/integration/gene.plot.CD8.example.00.txt",head=T)
        gene2group <- structure(factor(gene.plot.example.tb$Group,levels=unique(gene.plot.example.tb$Group)),
                                names=gene.plot.example.tb$geneID)
        gene.plot.tb <- gene.desc.top[geneID %in% gene.plot.example.tb$geneID,]
        gene.plot.tb[,y:=factor(geneID,levels=rev(gene.plot.example.tb$geneID))]
        gene.plot.tb[,x:=cluster.name ]
        gene.plot.tb[,ES:=comb.ES]
        gene.plot.tb[comb.ES > 0.5,ES:=0.5]
        gene.plot.tb[comb.ES < -0.25,ES:=-0.25]
        gene.plot.tb[,Group:=gene2group[geneID]]

        p <- ggplot(gene.plot.tb,aes(x,y)) +
                geom_point(aes(size=ES,color=ES),shape=16) +
                facet_grid(Group ~ ., scales = "free", space = "free") +
                scale_colour_distiller(palette = "RdYlBu") +
                labs(x="",y="") +
                theme_pubr() + 
                theme(strip.text.y = element_blank(),
                      axis.line.x=element_blank(),
                      axis.line.y=element_blank(),
                      panel.background = element_rect(colour = "black", fill = "white"),
                      #panel.grid = element_line(colour = "grey", linetype = "dashed"),
                      #panel.grid.major = element_line( colour = "grey", linetype = "dashed", size = 0.2),
                      axis.text.y = element_text(size=10),
                      axis.text.x = element_text(angle = 60,size=10, hjust = 1))
        ggsave(sprintf("%s,sigGene.top.Fig1SXX.CD8.01.pdf",out.prefix),width=6,height=16,useDingbats=F)
        saveRDS(p,file=sprintf("%s,sigGene.top.Fig1SXX.CD8.01.rds",out.prefix))
    }

    #### Tex
    {
        gene.desc.top <- gene.desc.tb.list[["CD8"]]
        gene.plot.example.tb <- fread("../data/expression/CD8/integration/gene.plot.CD8.Tex.example.00.txt",head=T)
        gene2group <- structure(factor(gene.plot.example.tb$Group,levels=unique(gene.plot.example.tb$Group)),
                                names=gene.plot.example.tb$geneID)
        gene.plot.tb <- gene.desc.top[geneID %in% gene.plot.example.tb$geneID,]
        gene.plot.tb <- gene.plot.tb[meta.cluster %in% c("CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13",
                                                         "CD8.c13.Tex.myl12a","CD8.c14.Tex.TCF7"),]
        gene.plot.tb[,y:=factor(geneID,levels=rev(gene.plot.example.tb$geneID))]
        gene.plot.tb[,x:=cluster.name ]
        gene.plot.tb[,ES:=comb.ES]
        gene.plot.tb[comb.ES > 0.5,ES:=0.5]
        gene.plot.tb[comb.ES < -0.25,ES:=-0.25]
        gene.plot.tb[,Group:=gene2group[geneID]]

        p <- ggplot(gene.plot.tb,aes(x,y)) +
                geom_point(aes(size=ES,color=ES),shape=16) +
                facet_grid(Group ~ ., scales = "free", space = "free") +
                scale_colour_distiller(palette = "RdYlBu") +
                labs(x="",y="") +
                #theme_pubr() + 
                theme_void() + 
                guides(colour=guide_colorbar(barheight=10,
                                             title.theme=element_text(size=14),
                                             label.theme=element_text(size=12)),
                       size=guide_legend(title.theme=element_text(size=14),
                                         label.theme=element_text(size=12))) +
                theme(strip.text.y = element_blank(),
                      axis.line.x=element_blank(),
                      axis.line.y=element_blank(),
                      panel.background = element_rect(colour = "black", fill = "white"),
                      #panel.grid = element_line(colour = "grey", linetype = "dashed"),
                      #panel.grid.major = element_line( colour = "grey", linetype = "dashed", size = 0.2),
                      axis.text.y = element_text(size=12),
                      axis.text.x = element_text(angle = 60,size=10, vjust = 1,hjust=1))
        ##ggsave(sprintf("%s,sigGene.top.Fig1SXX.Tex.01.pdf",out.prefix),width=6,height=12,useDingbats=F)
        ggsave(sprintf("%s,sigGene.top.Fig1SXX.Tex.01.pdf",out.prefix),width=3.0,height=14,useDingbats=F)
        saveRDS(p,file=sprintf("%s,sigGene.top.Fig1SXX.Tex.01.rds",out.prefix))
    }

    #### all CD4 clusters
    {
        gene.desc.top <- gene.desc.tb.list[["CD4"]]
        gene.plot.example.tb <- fread("../data/expression/CD4/integration/gene.plot.CD4.example.00.txt",head=T)
        gene2group <- structure(factor(gene.plot.example.tb$Group,levels=c("Th","Treg0","Tn","Tm","Th17","Tfh","Treg","NME1+")),
                                names=gene.plot.example.tb$geneID)
        gene.plot.tb <- gene.desc.top[geneID %in% gene.plot.example.tb$geneID,]
        gene.plot.tb[,y:=factor(geneID,levels=rev(gene.plot.example.tb$geneID))]
        gene.plot.tb[,x:=cluster.name ]
        gene.plot.tb[,ES:=comb.ES]
        gene.plot.tb[comb.ES > 0.5,ES:=0.5]
        gene.plot.tb[comb.ES < -0.25,ES:=-0.25]
        gene.plot.tb[,Group:=gene2group[geneID]]

        p <- ggplot(gene.plot.tb,aes(x,y)) +
                geom_point(aes(size=ES,color=ES),shape=16) +
                facet_grid(Group ~ ., scales = "free", space = "free") +
                scale_colour_distiller(palette = "RdYlBu") +
                labs(x="",y="") +
                theme_pubr() + 
                theme(strip.text.y = element_blank(),
                      axis.line.x=element_blank(),
                      axis.line.y=element_blank(),
                      panel.background = element_rect(colour = "black", fill = "white"),
                      #panel.grid = element_line(colour = "grey", linetype = "dashed"),
                      #panel.grid.major = element_line( colour = "grey", linetype = "dashed", size = 0.2),
                      axis.text.y = element_text(size=10),
                      axis.text.x = element_text(angle = 60,size=10, hjust = 1))
        ggsave(sprintf("%s,sigGene.top.Fig1SXX.CD4.01.pdf",out.prefix),width=6,height=16,useDingbats=F)
        saveRDS(p,file=sprintf("%s,sigGene.top.Fig1SXX.CD4.01.rds",out.prefix))
    }


}



