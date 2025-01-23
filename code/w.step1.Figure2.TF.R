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
source("./func.R")
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)

out.prefix <- "OUT_Fig2/TF/Fig2.TF"
dir.create(dirname(out.prefix),F,T)

# load data
{
#    dir.int.list <- list("CD8"="../data/expression/CD8/integration",
#                         "CD4"="../data/expression/CD4/integration")
#    sce.merged.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds"),
#                            "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.sce.merged.rds"))
    g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")

    gene.desc.tb.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds"),
                              "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.gene.tb.rds"))
    gene.long.cancerType.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.geneLong.cancerType.tb.rds"),
                                      "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.geneLong.cancerType.tb.rds"))
    gene.long.dataset.list <- list("CD8"=readRDS("../data/expression/CD8/integration/int.CD8.S35.geneLong.dataset.tb.rds"),
                                   "CD4"=readRDS("../data/expression/CD4/integration/int.CD4.S35.geneLong.dataset.tb.rds"))
}


###### top universal transcription factors
{
    gene.desc.top.slim <- gene.desc.tb.list$CD8
    dat.plot <- head(gene.desc.top.slim[meta.cluster=="CD8.c12.Tex.CXCL13" & sig==T &
				   cancerType.sig.freq > 0.8 & geneSet.TF==T,][order(-comb.ES),],n=10)
    dat.plot[,geneSymbol:=factor(geneSymbol,levels=rev(geneSymbol))]
    dat.plot[1,]
    ggbarplot(dat.plot,x="geneSymbol",y="comb.ES",
	      color=NA,fill="steelblue") +
        labs(x="",y="Effect Size") +
        geom_hline(yintercept=0.15,linetype="dashed",color="lightgray",alpha=0.8) +
        coord_flip()
    ggsave(sprintf("%s.top.universal.TF.bar.00.pdf",out.prefix),width=2.5,height=3.6)

}


#### compare effect size of genes across cancer types
{

    colSet.cancerType <- g.colSet$cancerType
    colSet.cancerType["panC"] <- "#f2f3f4"

    gene.to.plot <- c("CXCL13","HAVCR2","ENTPD1","LAYN","PDCD1",
		      "RBPJ", "TOX", "ZBED2", "TOX2", "VDR", "BATF", "IKZF4", "PRDM1", "STAT3", "IFI16",
		      "KLRD1", "SOX4","TP73","IRF4", "IRF8", "ETV7","EGR2","TSHZ1",
		      "FOXP3","KLRC3","IL17RB","IL1RL1",
		      "IL26","IL17A","IL17F","RORA","RORC","ZBTB16","CCL20",
		      "NR4A2", "ZNF282",
		      "NFATC1","NFATC2","NR5A2","ARID5B","ETV1",
		      "CD274","PDCD1LG2","NTRK1",
		      "KLRC1","XCL1")
    mcls.plot <- "CD8.c12.Tex.CXCL13"

    out.prefix.plot <- sprintf("%s/example.mod.3/%s",dirname(out.prefix),basename(out.prefix))
    dir.create(dirname(out.prefix.plot),F,T)
    makeFig.ExampleGeneBarplot(sprintf("%s.CD8.Tex",out.prefix.plot),gene.to.plot,
						       gene.long.dataset.list[["CD8"]],
                               gene.long.cancerType.list[["CD8"]],
                               gene.desc.tb.list[["CD8"]],mcls.plot,mod.sort=3)
    p.list <- readRDS(sprintf("%s.CD8.Tex.cmp.gene.example.fig.rds",out.prefix.plot))

    ## w: 3.5 * ?
    ## h: 2.25 * ? + 2.25
    p <- plot_grid(plot_grid(plotlist=llply(p.list[c("CXCL13","TOX","NFATC1",
						     "HAVCR2","TOX2","NFATC2",
						     "ENTPD1","ZBED2","NR5A2",
						     "LAYN","RBPJ","ARID5B",
						     "PDCD1","SOX4","ETV1")],
					    function(x){ x+theme(legend.position = "none") }),
			     align="hv",ncol=3),
		   get_legend(p.list[[1]]),
		   ncol=1,
		   rel_heights=c(3,1))
    ggsave(sprintf("%s.merged.01.pdf",out.prefix.plot),width=10.5,height=13.5)

    p <- plot_grid(plot_grid(plotlist=llply(p.list[c("RBPJ","BATF","NFATC1",
						     "TOX","IKZF4","NFATC2",
						     "ZBED2","PRDM1","NR5A2",
						     "TOX2","STAT3","ARID5B",
						     "VDR","IFI16","ETV1")],
					    function(x){ x+theme(legend.position = "none") }),
			     align="hv",ncol=3),
		   get_legend(p.list[[1]]),
		   ncol=1,
		   rel_heights=c(3,1))
    ggsave(sprintf("%s.merged.02.pdf",out.prefix.plot),width=10.5,height=13.5)

    p <- plot_grid(plot_grid(plotlist=llply(p.list[c("RBPJ","BATF","NR5A2",
						     "TOX","IKZF4","ARID5B",
						     "ZBED2","PRDM1","ETV1",
						     "TOX2","STAT3","SOX4",
						     "VDR","IFI16","FOXP3")],
					    function(x){ x+theme(legend.position = "none") }),
			     align="hv",ncol=3),
		   get_legend(p.list[[1]]),
		   ncol=1,
		   rel_heights=c(3,1))
    ggsave(sprintf("%s.merged.03.pdf",out.prefix.plot),width=10.5,height=13.5)

}


