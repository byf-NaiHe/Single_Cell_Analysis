#!/usr/bin/env Rscript

library("tidyverse")
library("tibble")
library("data.table")
library("plyr")
library("ggpubr")
library("ggplot2")
source("./func.R")

out.prefix <- "./OUT_Fig4/sharedGenes/sharedGenes"
colSet.file <- "../data/metaInfo/panC.colSet.list.rds"
gene.core.CD8.file <- "../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds"
gene.core.CD4.file <- "../data/expression/CD4/integration/int.CD4.S35.gene.tb.rds"
dir.create(dirname(out.prefix),F,T)

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

colSet <- readRDS(colSet.file)
gene.core.CD8.tb <- readRDS(gene.core.CD8.file)
gene.core.CD4.tb <- readRDS(gene.core.CD4.file)
gene.core.tb <- rbind(gene.core.CD8.tb, gene.core.CD4.tb)
gene.bg <- unique(gene.core.CD8.tb$geneID)
###############################

############## 
cmp.list=list("Tex.CD4CD8"=c("CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13"),
      "CD8Tex.ActTreg"=c("CD8.c12.Tex.CXCL13","CD4.c20.Treg.TNFRSF9"),
      "TfhTh1.ActTreg"=c("CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9"),
      "TwoCD8Tex"=c("CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13"),
      "TwoTfh"=c("CD4.c16.Tfh.CXCR5","CD4.c17.TfhTh1.CXCL13"),
      "TwoCD8Ts"=c("CD8.c05.Tem.CXCR5","CD8.c14.Tex.TCF7"),
      "TwoTreg"=c("CD4.c18.Treg.RTKN2","CD4.c20.Treg.TNFRSF9"),
      "ISG.TcTh"=c("CD8.c15.ISG.IFIT1","CD4.c22.ISG.IFIT1"),
      "ISG.TcTreg"=c("CD8.c15.ISG.IFIT1","CD4.c21.Treg.OAS1"),
      "ISG.ThTreg"=c("CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1"),
      "MAIT.Th17"=c("CD8.c16.MAIT.SLC4A10","CD4.c15.Th17.IL23R"))

########## venn of CD8+ Tex .vs. TfhTh1 .vs. ActTreg (fig. S34A right, Fig. 4C)
{
    sig.list.all <- list("CD8+Tex"=gene.core.tb[meta.cluster=="CD8.c12.Tex.CXCL13" & sig==T,][["geneSymbol"]],
             "CD4+Tfh/Th1"=gene.core.tb[meta.cluster=="CD4.c17.TfhTh1.CXCL13" & sig==T,][["geneSymbol"]],
             "CD4+ActTreg"=gene.core.tb[meta.cluster=="CD4.c20.Treg.TNFRSF9" & sig==T,][["geneSymbol"]])
    sigGeneVennPlot(sig.list.all,
            background.list=unique(gene.core.tb$geneSymbol),
            col.venn=colSet$meta.cluster[c("CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9")],
            fill.venn=colSet$meta.cluster[c("CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9")],
            out.prefix=sprintf("%s.venn.%s",out.prefix,"CD8Tex.TfhTh1.ActTreg"))

    sig.list.TF <- list("CD8+Tex"=gene.core.tb[meta.cluster=="CD8.c12.Tex.CXCL13" & sig==T & geneSet.TF==T,][["geneSymbol"]],
            "CD4+Tfh/Th1"=gene.core.tb[meta.cluster=="CD4.c17.TfhTh1.CXCL13" & sig==T & geneSet.TF==T,
                           ][["geneSymbol"]],
            "CD4+ActTreg"=gene.core.tb[meta.cluster=="CD4.c20.Treg.TNFRSF9" & sig==T & geneSet.TF==T,
                           ][["geneSymbol"]])
    sigGeneVennPlot(sig.list.TF,
            background.list=unique(gene.core.tb[geneSet.TF==T,][["geneSymbol"]]),
            col.venn=colSet$meta.cluster[c("CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9")],
            fill.venn=colSet$meta.cluster[c("CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9")],
            out.prefix=sprintf("%s.venn.%s.TF",out.prefix,"CD8Tex.TfhTh1.ActTreg"))
    sig.TF.cmp3.tb <-  sigGeneVennTable(gene.core.tb,
                    cmp=c("CD8.c12.Tex.CXCL13",
                          "CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9"),
                    only.sig=T)
    conn <- gzfile(sprintf("%s.cmp.CD8Tex.CD4DualFun.CD4ActTreg.txt.gz",out.prefix),"w")
    write.table(sig.TF.cmp3.tb,file=conn,row.names=F,sep="\t",quote=F)
    close(conn)

}

########## venn of ISG (fig. S34 A and B, middle, ISG)
{
    sig.list.all <- list("Tc.ISG"=gene.core.tb[meta.cluster=="CD8.c15.ISG.IFIT1" & sig==T,][["geneSymbol"]],
                 "Th.ISG"=gene.core.tb[meta.cluster=="CD4.c22.ISG.IFIT1" & sig==T,][["geneSymbol"]],
                 "Tr.ISG"=gene.core.tb[meta.cluster=="CD4.c21.Treg.OAS1" & sig==T,][["geneSymbol"]])
    sigGeneVennPlot(sig.list.all,
            background.list=unique(gene.core.tb$geneSymbol),
            col.venn=colSet$meta.cluster[c("CD8.c15.ISG.IFIT1","CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1")],
            fill.venn=colSet$meta.cluster[c("CD8.c15.ISG.IFIT1","CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1")],
            out.prefix=sprintf("%s.venn.%s",out.prefix,"ISG.g3"))

    sig.list.TF <- list("Tc.ISG"=gene.core.tb[meta.cluster=="CD8.c15.ISG.IFIT1" & sig==T & geneSet.TF==T,
                ][["geneSymbol"]],
                "Th.ISG"=gene.core.tb[meta.cluster=="CD4.c22.ISG.IFIT1" & sig==T & geneSet.TF==T,
                          ][["geneSymbol"]],
                "Tr.ISG"=gene.core.tb[meta.cluster=="CD4.c21.Treg.OAS1" & sig==T & geneSet.TF==T,
                          ][["geneSymbol"]])
    sigGeneVennPlot(sig.list.TF,
            background.list=unique(gene.core.tb$geneSymbol),
            col.venn=colSet$meta.cluster[c("CD8.c15.ISG.IFIT1","CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1")],
            fill.venn=colSet$meta.cluster[c("CD8.c15.ISG.IFIT1","CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1")],
            out.prefix=sprintf("%s.venn.%s.TF",out.prefix,"ISG.g3"))
    sig.TF.cmp3.tb <-  sigGeneVennTable(gene.core.tb,
                        cmp=c("CD8.c15.ISG.IFIT1",
                          "CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1"),
                        only.sig=T)
    conn <- gzfile(sprintf("%s.cmp.ISG.g3.txt.gz",out.prefix),"w")
    write.table(sig.TF.cmp3.tb,file=conn,row.names=F,sep="\t",quote=F)
    close(conn)
}

######### nichenet of overlapped genes among 3 meta-clusters (three ISG meta-clusters)
{
    ### nichenet for ISG
    {
        gene.cmp.ISG.tb <- fread(sprintf("%s.cmp.ISG.g3.txt.gz",out.prefix))
        #gene.bg <- unique(gene.core.CD8.tb$geneID)
        #### comp3
        gene.oi <- gene.cmp.ISG.tb[sig_CD4.c21.Treg.OAS1==T & sig_CD4.c22.ISG.IFIT1==T & sig_CD8.c15.ISG.IFIT1==T,
                                  ][["geneID"]]
        es.tb <- dcast(gene.core.tb[meta.cluster %in% c("CD8.c15.ISG.IFIT1","CD4.c22.ISG.IFIT1","CD4.c21.Treg.OAS1"),
                     c("geneID","meta.cluster","comb.ES")],geneID~meta.cluster,value.var="comb.ES")
        es.avg.df <- tibble(gene=es.tb$geneID,avg.comb.ES=rowMeans(es.tb[,-("geneID")]))
        tic("run.nichenet (ISG.g3)")
        res.nichenet.ISG.g3 <- run.nichenet(gene.oi,gene.bg,out.prefix=sprintf("%s.nichenet.cmp.ISG.g3",out.prefix),
                            #n.top=20,comb.height=6,comb.rel.height=c(5,2),lr.height=4,
                            n.top=10,comb.height=4,comb.rel.height=c(5,2),lr.height=4,
                            ligands_all=c("IFNA1","IFNB1"),es.df=es.avg.df,do.eval=T)
        toc()
    }

}

########## venn of Tc17 and Th17 (fig. S34 A and B, left, T*17; fig. S35, T*17)
{
    res.nichenet.MAIT.Th17.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c16.MAIT.SLC4A10","CD4.c15.Th17.IL23R"),
                         sname=c("MAIT","Th17"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.rel.height=c(5,2),lr.height=4,lr.width=8,
                         #n.top=20,comb.height=6,comb.rel.height=c(5,2),lr.height=4,lr.width=8,
                         ligands_all="IL23A",
                         #es.df=es.avg.df,
                         do.eval=T)
}

########### venn and nichenet of other meta-cluster pairs (fig. S35)
{

    res.nichenet.TwoTemra.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c07.Temra.CX3CR1","CD4.c13.Temra.CX3CR1"),
                         sname=c("CD8.Temra","CD4.Temra"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.rel.height=c(5,2),lr.height=4,lr.width=10,
                         #n.top=20,comb.height=6,comb.rel.height=c(5,2),lr.height=4,lr.width=10,
                         ligands_all=c("IL12A", "AGT", "IL21", "IL2", "IL15"),
                         #es.df=es.avg.df,
                         do.eval=T)

    res.nichenet.CD8TexActTreg.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c12.Tex.CXCL13","CD4.c20.Treg.TNFRSF9"),
                         sname=c("CD8Tex","ActTreg"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=15,comb.rel.height=c(5,2),
                         comb.rel.width=c(0.12,1.00),
                         lr.height=4,lr.width=10,
                         ligands_all=c("IFNB1"),
                         #es.df=es.avg.df,
                         do.eval=F)

    res.nichenet.CD8TexCD4TfhTh1.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c12.Tex.CXCL13","CD4.c17.TfhTh1.CXCL13"),
                         sname=c("CD8Tex","CD4TfhTh1"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=14,comb.rel.height=c(5,2),
                         comb.rel.width=c(0.1,0.9),
                         lr.height=4,lr.width=9,
                         ligands_all=c("TGFB1"),
                         #es.df=es.avg.df,
                         do.eval=F)

    res.nichenet.CD4TfhTh1ActTreg.g2 <- run.venn.nicheNet.g2(mcls=c("CD4.c17.TfhTh1.CXCL13","CD4.c20.Treg.TNFRSF9"),
                         sname=c("CD4TfhTh1","ActTreg"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=14,comb.rel.height=c(5,2),
                         comb.rel.width=c(0.1,0.9),
                         lr.height=4,lr.width=9,
                         #ligands_all=c("TGFB1"),
                         ligands_all=c("IFNB1"),
                         #es.df=es.avg.df,
                         do.eval=F)

    res.nichenet.TwoTn.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c01.Tn.MAL","CD4.c01.Tn.TCF7"),
                         sname=c("CD8.Tn","CD4.Tn"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=7.5,comb.rel.height=c(5,2),
                         comb.rel.width=c(0.12,0.63),
                         lr.height=4,lr.width=6,
                         ligands_all=NULL,
                         #es.df=es.avg.df,
                         do.eval=F)

   res.nichenet.CD8TnCD8Tex.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c01.Tn.MAL","CD8.c12.Tex.CXCL13"),
                         sname=c("CD8.Tn","CD8.Tex"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=7.5,comb.rel.height=c(5,2),
                         comb.rel.width=c(0.12,0.63),
                         lr.height=4,lr.width=6,
                         ligands_all=NULL,
                         #es.df=es.avg.df,
                         do.eval=F)

   res.nichenet.CD8TemraCD8Tex.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c07.Temra.CX3CR1","CD8.c12.Tex.CXCL13"),
                         sname=c("CD8.Temra","CD8.Tex"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=7.5,comb.rel.height=c(5,2),
                         comb.rel.width=c(0.12,0.63),
                         lr.height=4,lr.width=6,
                         ligands_all=NULL,
                         #es.df=es.avg.df,
                         do.eval=F)

    res.nichenet.TwoCD8Ts.g2 <- run.venn.nicheNet.g2(mcls=c("CD8.c05.Tem.CXCR5","CD8.c14.Tex.TCF7"),
                         sname=c("CD8.Tem.CXCR5","CD8.Ts.TCF7"),out.prefix=out.prefix,gene.bg,
                         n.top=10,comb.height=4,comb.width=7.5,
                         comb.rel.height=c(5,2),comb.rel.width=c(0.12,0.63),
                         lr.height=4,lr.width=6,
                         ligands_all=NULL,
                         #es.df=es.avg.df,
                         do.eval=F)

}

{
    
    save(res.nichenet.ISG.g3,
         res.nichenet.MAIT.Th17.g2,res.nichenet.TwoTemra.g2,res.nichenet.TwoTn.g2,res.nichenet.TwoCD8Ts.g2,
         res.nichenet.CD8TexActTreg.g2,res.nichenet.CD8TexCD4TfhTh1.g2,res.nichenet.CD4TfhTh1ActTreg.g2,
	 file=sprintf("%s.res.nichenet.RData",out.prefix))
    
    lname <- load(sprintf("%s.res.nichenet.RData",out.prefix))

    #### print something
    l_ply(list(res.nichenet.ISG.g3,res.nichenet.CD8TexCD4TfhTh1.g2),function(x){
				    print(str(x$eval.list))
			     })
    ########################

}

