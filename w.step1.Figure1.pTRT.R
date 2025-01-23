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

g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
g.colSet$cancerType <- c(g.colSet$cancerType,c("panC"="#f2f3f4"))
###g.colSet$cluster.name <- c(g.colSet$cluster.name,c("CD8.c16(MAIT/Tc17)"=unname(g.colSet$cluster.name["CD8.c16(Tc17)"])))
tcr.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
pTRT.tb.file <- "../data/tcr/pTRT/TRT.signaling.tcr.ext.sc.tb.rds"
#startrac.file <- "./OUT.startrac.fltMinCell10/CD8/tcr.zhangLab.comb.CD8.fltMinCell10.out.startrac.nperm1000.rds"
####sample.usedForFreq.file <- "/lustre1/zeminz_pkuhpc/zhenglt/work/panC/ana/zhangLab.10X/inte.metaClust.20200111/OUT.freq/panC.freq.all.meta.tb.summary.count.usedForFreq.baseline.rds"
out.prefix <- "OUT_Fig1/pTRT/pTRT.signaling"

dir.create(dirname(out.prefix),F,T)

#### tcr data
{
    in.dat <- readRDS(tcr.file)
    in.dat$majorCluster <- as.character(in.dat$meta.cluster)
    in.dat$clone.id <- in.dat$cloneID
    #### filter out patient.cluster with number of cell < 10
    ncell.patient.cluster <- sort(unclass(in.dat[,table(sprintf("%s.%s",patient,majorCluster))]))
    ##ncell.patient.cluster <- sort(unclass(in.dat[stype=="CD8", table(sprintf("%s.%s",patient,majorCluster))]))
    in.dat <- in.dat[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
    #in.dat <- in.dat[stype=="CD8",]
    ####
    in.dat[dataset=="NSCLC.XinyiGuo2018",dataset:="LC.XinyiGuo2018"]
    in.dat[dataset=="STAD.BoxiKang2020",dataset:="STAD.BoxiKang2019"]


}


####### gene sets
{

    gset.prol.file <- "../data/external/gene.ProliferationScore.list"
    gset.file <- "../data/external/c2.cp.v7.0.symbols.gmt"
    gset.tb <- as.data.table(clusterProfiler::read.gmt(gset.file))
    gset.prol.tb <- fread(gset.prol.file)

    fname.TCR.signaling <- c("REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
			     "REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE",
			     "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
			     "REACTOME_DOWNSTREAM_TCR_SIGNALING",
			     "REACTOME_TCR_SIGNALING")

    gset.list <- llply(fname.TCR.signaling,function(x){
				gset.tb[term==x,][["gene"]]
			     })
    names(gset.list) <- fname.TCR.signaling
    gset.list[["proliferation"]] <- gset.prol.tb$geneSymbol

}

######### main analysis ##################
### pTRTs data
if(!file.exists(pTRT.tb.file))
{

    gsea.act.tb <- as.data.table(ldply(c("CD8","CD4"),function(stype){
        gsea.act.list <- llply(unique(in.dat$dataset),function(data.id){
            gsea.out.list <- readRDS(file=sprintf("%s.%s.gsea.out.list.zscore.mini.%s.rds",
                              out.prefix,stype,data.id))

            gsea.act.NES.tb <- do.call(cbind,llply(seq_len(length(gsea.out.list)),function(i){
                 out.tb <- gsea.out.list[[i]][,"NES"]
                 names(out.tb) <- rownames(gsea.out.list[[i]])
                 gset.padding <- setdiff(names(gset.list),names(out.tb))
                 out.tb <- c(out.tb,
                         structure(rep(0,length(gset.padding)),
                               names=gset.padding))[names(gset.list)]
                 return(out.tb)
                       }))
            colnames(gsea.act.NES.tb) <- names(gsea.out.list)
            gsea.act.NES.tb[,1:5]

            gsea.act.pvalue.tb <- do.call(cbind,llply(seq_len(length(gsea.out.list)),function(i){
                 out.tb <- gsea.out.list[[i]][,"pvalue"]
                 names(out.tb) <- rownames(gsea.out.list[[i]])
                 gset.padding <- setdiff(names(gset.list),names(out.tb))
                 out.tb <- c(out.tb,
                         structure(rep(1,length(gset.padding)),
                               names=gset.padding))[names(gset.list)]
                 return(out.tb)
                       }))
            colnames(gsea.act.pvalue.tb) <- names(gsea.out.list)
            gsea.act.pvalue.tb[,1:5]

            return(list("gsea.act.NES.tb"=gsea.act.NES.tb,"gsea.act.pvalue.tb"=gsea.act.pvalue.tb))
                     })
        names(gsea.act.list) <- unique(in.dat$dataset)

        ####
        out.tb <- as.data.table(ldply(names(gsea.act.list),function(x){
                dat.NES.mtx <- gsea.act.list[[x]][["gsea.act.NES.tb"]]
                dat.pvalue.mtx <- gsea.act.list[[x]][["gsea.act.pvalue.tb"]]
                print("all(colnames(dat.NES.mtx)==colnames(dat.pvalue.mtx))")
                print(all(colnames(dat.NES.mtx)==colnames(dat.pvalue.mtx)))
                data.table(dataset=x,
                       stype=stype,
                       miniCluster=colnames(dat.NES.mtx),
                       NES=dat.NES.mtx["REACTOME_TCR_SIGNALING",],
                       pvalue=dat.pvalue.mtx["REACTOME_TCR_SIGNALING",],
                       prol.NES=dat.NES.mtx["proliferation",],
                       prol.pvalue=dat.pvalue.mtx["proliferation",]
                       )
                     }))
        return(out.tb)

    }))

    tcr.ext.tb <- merge(in.dat[,.(miniCluster.size=.N),
			    by=c("cancerType","dataset","miniCluster",
				 "meta.cluster","cluster.name","stype")],
			gsea.act.tb,
			by.x=c("dataset","miniCluster","stype"),
			by.y=c("dataset","miniCluster","stype"))

    tcr.ext.summary.tb <- tcr.ext.tb[,.(NES.med=median(NES),
		  pvalue.med=median(pvalue),
		  nsig=sum(pvalue<0.05 & NES>0),
		  ntotal=.N,
		  prol.NES.med=median(prol.NES),
		  nsig.prol=sum(prol.NES > 0 & prol.pvalue < 0.05)
		  ),
	       by=c("meta.cluster")]
    tcr.ext.summary.tb[,freq.sig:=nsig/ntotal]
    tcr.ext.summary.tb[,prol.freq.sig:=nsig.prol/ntotal]
    tcr.ext.summary.tb[order(NES.med),]

    ### single-cell level
    {
        tcr.ext.sc.tb <- merge(in.dat,tcr.ext.tb[,c("stype","miniCluster","miniCluster.size",
                                "NES","pvalue","prol.NES","prol.pvalue"),with=F],
                       by.x=c("stype","miniCluster"),
                       by.y=c("stype","miniCluster"),all.x=T)

        ### frequency of proliferative cells. replace Fig.1G ? (hi: freq.prol > 0.15)
        tcr.ext.sc.prol.summary.tb <- tcr.ext.sc.tb[,.(n.total=.N,
                                   n.prol=sum(prol.NES>0 & prol.pvalue<0.05)),
                                by="cluster.name"]
        tcr.ext.sc.prol.summary.tb[,freq.prol:=n.prol/n.total]
        tcr.ext.sc.prol.summary.tb[order(freq.prol),]
    }

    #saveRDS(tcr.ext.tb,file=sprintf("%s.tcr.ext.tb.rds",out.prefix))
    #saveRDS(tcr.ext.sc.tb,file=sprintf("%s.tcr.ext.sc.tb.rds",out.prefix))
    saveRDS(tcr.ext.sc.tb,file=pTRT.tb.file)

}else{
    tcr.ext.sc.tb <- readRDS(file=pTRT.tb.file)
    #tcr.ext.sc.tb <- readRDS(file=sprintf("%s.tcr.ext.sc.tb.rds",out.prefix))
    tcr.ext.sc.tb[dataset=="NSCLC.XinyiGuo2018",dataset:="LC.XinyiGuo2018"]
    tcr.ext.sc.tb[dataset=="STAD.BoxiKang2020",dataset:="STAD.BoxiKang2019"]
}

###### make the plots (fig. S9)
{
    
    run.one <- function(in.dat,tcr.ext.sc.bear.TCR_hi.tb,
			out.prefix,
			par.loc="T")
    {

        dim(tcr.ext.sc.bear.TCR_hi.tb)
        col.use <- c("cancerType","dataset","Tech", "meta.cluster","cluster.name","patient")
        col.use.P <- c("cancerType","dataset","Tech", "patient")
        tcr.perPatient.onlyT.tb <- in.dat[loc==par.loc,.(n.total=.N), by=col.use]
        tcr.perPatient.TCR_hi.onlyT.tb <- tcr.ext.sc.bear.TCR_hi.tb[loc==par.loc,][,.(n.TCR_hi=.N), by=col.use]
        tcr.perPatient.ana.onlyT <- merge(tcr.perPatient.TCR_hi.onlyT.tb,
                          tcr.perPatient.onlyT.tb,
                          all.y=T)
        tcr.perPatient.ana.onlyT[is.na(n.TCR_hi),n.TCR_hi:=0]
        tcr.perPatient.ana.onlyT[,freq.TCR_hi:=n.TCR_hi/n.total]
        tcr.perPatient.ana.onlyT <- tcr.perPatient.ana.onlyT[,.(meta.cluster=.SD$meta.cluster,
                        cluster.name=.SD$cluster.name,
                        n.TCR_hi=.SD$n.TCR_hi,
                        n.total=.SD$n.total,
                        patient.ncells.TCR_hi=sum(.SD$n.TCR_hi),
                        patient.ncells=sum(.SD$n.total),
                        freq.TCR_hi=.SD$freq.TCR_hi,
                        dist.TCR_hi=if(sum(.SD$n.TCR_hi)==0) 0 else .SD$n.TCR_hi/sum(.SD$n.TCR_hi)),
                     by=col.use.P]
        tcr.perPatient.ana.onlyT.full <- tcr.perPatient.ana.onlyT
        tcr.perPatient.ana.onlyT.full[,summary(lm(n.TCR_hi~n.total+meta.cluster))]
    #	tcr.perPatient.ana.onlyT.panC <- copy(tcr.perPatient.ana.onlyT)
    #	tcr.perPatient.ana.onlyT.panC[,cancerType:="panC"]
    #	tcr.perPatient.ana.onlyT.panC[,dataset:="thisStudy"]
    #	tcr.perPatient.ana.onlyT.full <- rbind(tcr.perPatient.ana.onlyT,
    #					       tcr.perPatient.ana.onlyT.panC)

        tcr.perPatient.TCR_hi.onlyT.strict.tb <- tcr.ext.sc.TCR_hi.tb[loc==par.loc,][,.(n.TCR_hi=.N),by=col.use]
        tcr.perPatient.ana.onlyT.strict <- merge(tcr.perPatient.TCR_hi.onlyT.strict.tb,
                             tcr.perPatient.onlyT.tb,all.y=T)
        tcr.perPatient.ana.onlyT.strict[is.na(n.TCR_hi),n.TCR_hi:=0]
        tcr.perPatient.ana.onlyT.strict[,freq.TCR_hi:=n.TCR_hi/n.total]
        tcr.perPatient.ana.onlyT.strict <- tcr.perPatient.ana.onlyT.strict[,.(meta.cluster=.SD$meta.cluster,
                                              cluster.name=.SD$cluster.name,
                                              n.TCR_hi=.SD$n.TCR_hi,
                                              n.total=.SD$n.total,
                                              patient.ncells.TCR_hi=sum(.SD$n.TCR_hi),
                                              patient.ncells=sum(.SD$n.total),
                                              freq.TCR_hi=.SD$freq.TCR_hi,
                                              dist.TCR_hi=.SD$n.TCR_hi/sum(.SD$n.TCR_hi)),
                                            by=col.use.P]
        tcr.perPatient.ana.onlyT.full.strict <- tcr.perPatient.ana.onlyT.strict

        #########################################################################
        ## in each meta.cluster, what fraction of cells are potentail tumor-reactive T cells
        ## (sharing TCR with cells with high TCR signaling)
        {
            ### compare meta.cluster
            #out.prefix.cmpMetaCluster <- sprintf("%s/freq.cmpMetaCluster/%s",
            #						 dirname(out.prefix),basename(out.prefix))
            out.prefix.freqTRT <- sprintf("%s/freq.TRT/%s",
                             dirname(out.prefix),basename(out.prefix))
            dir.create(dirname(out.prefix.freqTRT),F,T)
            data.plot <- tcr.perPatient.ana.onlyT.full[cancerType!="panC",][n.total>=30,]
            data.plot.med <- data.plot[,.(freq.TCR_hi.med=median(freq.TCR_hi)),
                           by="cluster.name"][order(freq.TCR_hi.med),]
            data.plot[,cluster.name:=factor(cluster.name,levels=data.plot.med$cluster.name)]

            p <- ggboxplot(data.plot,x="cluster.name",y="freq.TCR_hi",
               fill = "cluster.name",
               color = "cluster.name",
               legend="none",title="",
               alpha=0.8, xlab="",ylab="Freq. Of pTRTs",
               add = "jitter",outlier.shape=NA) +
            scale_fill_manual(values=g.colSet$cluster.name) +
            scale_color_manual(values=g.colSet$cluster.name) +
            stat_compare_means(label="p.format") +
            coord_cartesian(clip="off",ylim=c(0,1)) +
            theme(plot.title = element_text(hjust = 0.5,size=14),
              axis.text.x = element_text(angle = 60, hjust = 1,vjust=1))
            ggsave(sprintf("%s.freq.TRT.%s.pdf",out.prefix.freqTRT, "meta.cluster"),
               width=6,height=5,useDingbats=F)

            ### strict version (only consider cells showing high TCR signaling, not include those)
            ### sharing TCR with them
            data.plot <- tcr.perPatient.ana.onlyT.full.strict[cancerType!="panC",][n.total>=30,]
            data.plot.med <- data.plot[,.(freq.TCR_hi.med=median(freq.TCR_hi)),by="cluster.name"][order(freq.TCR_hi.med),]
            data.plot[,cluster.name:=factor(cluster.name,levels=data.plot.med$cluster.name)]

            p <- ggboxplot(data.plot,x="cluster.name",y="freq.TCR_hi",
              fill = "cluster.name",
              color = "cluster.name",
              legend="none",title="",
              alpha=0.8, xlab="",ylab="NES",
              #add = "jitter",
              outlier.shape=NA) +
            scale_fill_manual(values=g.colSet$cluster.name) +
            scale_color_manual(values=g.colSet$cluster.name) +
            stat_compare_means(label="p.format",label.y=0.25) +
            stat_compare_means(label="p.signif",ref.group="CD8.c01(Tn)") +
            coord_cartesian(clip="off") +
            theme(plot.title = element_text(hjust = 0.5,size=14),
              axis.text.x = element_text(angle = 60, hjust = 1,vjust=1))
            ggsave(sprintf("%s.freq.TRT.strict.%s.pdf",out.prefix.freqTRT, "meta.cluster"),
               width=6,height=5,useDingbats=F)

            ### compare cancerType
            #out.prefix.cmpCancerType <- sprintf("%s/freq.cmpCancerType/%s",
            #				       dirname(out.prefix),basename(out.prefix))
            #dir.create(dirname(out.prefix.cmpCancerType),F,T)

            plot.list <- llply(unique(tcr.perPatient.ana.onlyT.full$cluster.name),function(mcls){
                data.plot <- tcr.perPatient.ana.onlyT.full[cluster.name==mcls,]
                data.plot.med <- data.plot[,.(freq.TCR_hi.med=median(freq.TCR_hi)),
                           by="cancerType"][order(freq.TCR_hi.med),]
                ##data.plot[,cancerType:=factor(cancerType, levels=c(setdiff(data.plot.med$cancerType,"panC"), "panC"))]
                data.plot[,cancerType:=factor(cancerType, levels=data.plot.med$cancerType)]
            p <- ggboxplot(data.plot,
                     x="cancerType",y="freq.TCR_hi",
                     fill = "cancerType", legend="none",title=mcls,
                     alpha=0.8, xlab="",ylab="Freq. Of pTRTs",
                     add = "jitter",outlier.shape=NA) +
                    scale_fill_manual(values=g.colSet$cancerType) +
                    stat_compare_means(label="p.format") +
                    coord_cartesian(clip="off",ylim=c(0,1)) +
                    theme(plot.title = element_text(hjust = 0.5,size=14),
                      axis.text.x = element_text(angle = 60, hjust = 1,vjust=1))
            ggsave(sprintf("%s.freq.TRT.%s.pdf",out.prefix.freqTRT,
                       as.character(data.plot$meta.cluster[1])),
                   width=5,height=3,useDingbats=F)
            return(p)
            })
            names(plot.list) <- unique(tcr.perPatient.ana.onlyT.full$cluster.name)
        }

        #########################################################################
        ## in each patient/tumor, in which meta.clusters the potentail tumor-reactive T cells
        ## (sharing TCR with cells with high TCR signaling) distribute ?
        {
            ### compare meta.cluster
            out.prefix.distTRT <- sprintf("%s/dist.TRT/%s",
                             dirname(out.prefix),basename(out.prefix))
            #out.prefix.cmpMetaCluster <- sprintf("%s/freq.cmpMetaCluster/%s",
            #					 dirname(out.prefix),basename(out.prefix))
            dir.create(dirname(out.prefix.distTRT),F,T)
            ##data.plot <- tcr.perPatient.ana.onlyT.full[cancerType!="panC",][patient.ncells>=30,]
            data.plot <- tcr.perPatient.ana.onlyT.full[cancerType!="panC",][patient.ncells.TCR_hi>=10,]
            data.plot.med <- data.plot[,.(dist.TCR_hi.med=median(dist.TCR_hi),
                          dist.TCR_hi.q75=quantile(dist.TCR_hi,0.75)),
                           by="cluster.name"]
            data.plot.med <- data.plot.med[order(dist.TCR_hi.med,dist.TCR_hi.q75,cluster.name),]
            data.plot[,cluster.name:=factor(cluster.name,levels=data.plot.med$cluster.name)]
            data.plot.debug <<- data.plot
            print(data.plot.med)

            p <- ggboxplot(data.plot,x="cluster.name",y="dist.TCR_hi",
               fill = "cluster.name",
               color = "cluster.name",
               legend="none",title="",
               alpha=0.8, xlab="",ylab="Freq. Of pTRTs",
               add.params=list(size=0.5),
               ###### strang ! when large number of 0s, jitter in y direction cause some negative frequencies (y values) !
               add = "jitter",
               outlier.shape=NA) +
            scale_fill_manual(values=g.colSet$cluster.name) +
            scale_color_manual(values=g.colSet$cluster.name) +
            stat_compare_means(label="p.format") +
            coord_cartesian(clip="off",ylim=c(0,1)) +
            theme(plot.title = element_text(hjust = 0.5,size=14),
              axis.text.x = element_text(angle = 60, hjust = 1,vjust=1))
            ggsave(sprintf("%s.dist.TRT.%s.test.pdf",out.prefix.distTRT, "meta.cluster"),
               width=6,height=5,useDingbats=F)

            write.table(data.plot,file=sprintf("%s.dist.TRT.%s.dump.txt",out.prefix.distTRT, "meta.cluster"),
                quote=F, row.names=F,sep="\t")
            saveRDS(p,sprintf("%s.dist.TRT.%s.rds",out.prefix.distTRT, "meta.cluster"))

            ### compare cancerType
            #out.prefix.freq <- sprintf("%s/freq.cmpCancerType/%s",
            #				       dirname(out.prefix),basename(out.prefix))
            #dir.create(dirname(out.prefix.freq),F,T)

            plot.list <- llply(unique(tcr.perPatient.ana.onlyT.full$cluster.name),function(mcls){
            ##data.plot <- tcr.perPatient.ana.onlyT.full[cluster.name==mcls,]
            data.plot <- tcr.perPatient.ana.onlyT.full[patient.ncells.TCR_hi>=10,][cluster.name==mcls,]
            data.plot.med <- data.plot[,.(dist.TCR_hi.med=median(dist.TCR_hi)),
                           by="cancerType"][order(dist.TCR_hi.med),]
            #data.plot[,cancerType:=factor(cancerType, levels=c(setdiff(data.plot.med$cancerType,"panC"), "panC"))]
            data.plot[,cancerType:=factor(cancerType, levels=data.plot.med$cancerType)]
            p <- ggboxplot(data.plot,
                     x="cancerType",y="dist.TCR_hi",
                     fill = "cancerType",
                     color = "cancerType",
                     legend="none",title=mcls,
                     alpha=0.8, xlab="",ylab="Freq. Of pTRTs",
                     add = "jitter",outlier.shape=NA) +
                    scale_fill_manual(values=g.colSet$cancerType) +
                    scale_color_manual(values=g.colSet$cancerType) +
                    stat_compare_means(label="p.format") +
                    coord_cartesian(clip="off") +
                    theme(plot.title = element_text(hjust = 0.5,size=14),
                      axis.text.x = element_text(angle = 60, hjust = 1,vjust=1))
            ggsave(sprintf("%s.dist.TRT.%s.pdf",out.prefix.distTRT,
                       as.character(data.plot$meta.cluster[1])),
                   width=5,height=3,useDingbats=F)
            return(p)
            })
            names(plot.list) <- unique(tcr.perPatient.ana.onlyT.full$cluster.name)
            saveRDS(plot.list,file=sprintf("%s.dist.TRT.cmpCancerType.plot.list.rds",out.prefix.distTRT))
        }

    }

    #### pTRTs by high level of TCR signaling | proliferation (fig. S9)
    {

        tcr.ext.sc.onlyT.tb <- tcr.ext.sc.tb[loc=="T",]
        tcr.ext.sc.onlyT.size.tb <- tcr.ext.sc.onlyT.tb[,.(cloneSizeInT=.N),
                                by=c("stype","patient","cloneID","cancerType","loc",
                                     "cloneSize","dataset","dataset.old","clone.id")]
        tcr.ext.sc.onlyT.size.tb[,cloneDist:=cloneSizeInT/cloneSize]
        tcr.ext.sc.onlyT.size.flt.tb <- tcr.ext.sc.onlyT.size.tb[cloneSizeInT >= 2 & cloneDist >= 0,]
        dim(tcr.ext.sc.onlyT.size.flt.tb)
        ### 6593 (cloneDist >= 1)
        ### 8376 (cloneDist >= 0.5)
        ### 8985 (cloneDist >= 0)

        tcr.ext.sc.TCR_hi.tb <- tcr.ext.sc.tb[ ( (NES > 0 & pvalue<0.05) |
                             (prol.NES > 0 & prol.pvalue<0.05)
                               ) &
                            clone.id %in% tcr.ext.sc.onlyT.size.flt.tb$clone.id &
                            loc=="T",]
        tcr.ext.sc.TCR_hi.tb[,table(loc)]
        tcr.ext.sc.TCR_hi.tb[,table(cloneSize>=2)]
        ### loc
        ###     T
        ### 13218 (cloneDist >= 1)
        ### 21381 (cloneDist >= 0.5)
        ### 22861 (cloneDist >= 0)
        tcr.ext.sc.bear.TCR_hi.tb <- tcr.ext.sc.tb[cloneID %in% tcr.ext.sc.TCR_hi.tb$cloneID,]
        dim(tcr.ext.sc.bear.TCR_hi.tb)
        ### [1] 21518   38 (cloneDist >= 1)
        ### [1] 43141   38 (cloneDist >= 0.5)
        ### [1] 59039   38 (cloneDist >= 0)

        l_ply(c("CD8","CD4"),function(x){
            run.one(in.dat[stype==x,],tcr.ext.sc.bear.TCR_hi.tb[stype==x,],
                out.prefix=sprintf("%s/TCR_hi.or.Prol_hi.Tumor.v2/TRT.signaling.%s",dirname(out.prefix),x), par.loc="T")

            run.one(in.dat[stype==x,],tcr.ext.sc.bear.TCR_hi.tb[stype==x,],
                out.prefix=sprintf("%s/TCR_hi.or.Prol_hi.Blood.v2/TRT.signaling.%s",dirname(out.prefix),x), par.loc="P")
            },.parallel=T)

        p1 <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Tumor.v2/dist.TRT/TRT.signaling.CD8.dist.TRT.meta.cluster.rds",
                              dirname(out.prefix)))
        p2 <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Tumor.v2/dist.TRT/TRT.signaling.CD4.dist.TRT.meta.cluster.rds",
                              dirname(out.prefix)))
        p3 <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Blood.v2/dist.TRT/TRT.signaling.CD8.dist.TRT.meta.cluster.rds",
                              dirname(out.prefix)))
        p4 <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Blood.v2/dist.TRT/TRT.signaling.CD4.dist.TRT.meta.cluster.rds",
                              dirname(out.prefix)))
        p <- cowplot::plot_grid(p1,p2,p3,p4,align="hv", rel_widths=c(0.85,1))
        ggsave(sprintf("%s/TCR_hi.or.Prol_hi.fig.00.cloneEnrichedInT.v2.pdf",dirname(out.prefix)),
               width=11,height=10,useDingbats=FALSE)

        #### cmpCancerType
        {

            plist.t <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Tumor.v2/dist.TRT/TRT.signaling.CD8.dist.TRT.cmpCancerType.plot.list.rds",dirname(out.prefix)))
            plist.b <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Blood.v2/dist.TRT/TRT.signaling.CD8.dist.TRT.cmpCancerType.plot.list.rds",dirname(out.prefix)))
            p <- cowplot::plot_grid(plotlist=c(plist.t[c("CD8.c10(ZNF683+CXCR6+ Trm)","CD8.c12(terminal Tex)")],
                               plist.b[c("CD8.c07(Temra)")]),
                        align="hv", rel_widths=c(1,1,0.5),nrow=1)
            ggsave(sprintf("%s/TCR_hi.or.Prol_hi.fig.CD8.01.cloneEnrichedInT.v2.pdf",dirname(out.prefix)),
               width=12,height=3.0,useDingbats=FALSE)

            plist.t <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Tumor.v2/dist.TRT/TRT.signaling.CD4.dist.TRT.cmpCancerType.plot.list.rds",dirname(out.prefix)))
            plist.b <- readRDS(sprintf("%s/TCR_hi.or.Prol_hi.Blood.v2/dist.TRT/TRT.signaling.CD4.dist.TRT.cmpCancerType.plot.list.rds",dirname(out.prefix)))
            p <- cowplot::plot_grid(plotlist=c(plist.t[c("CD4.c17(IFNG+ Tfh/Th1)","CD4.c20(TNFRSF9+ Treg)")],
                               plist.b[c("CD4.c13(Temra)")]),
                        align="hv", rel_widths=c(1,1,0.5),nrow=1)
            ggsave(sprintf("%s/TCR_hi.or.Prol_hi.fig.CD4.01.cloneEnrichedInT.v2.pdf",dirname(out.prefix)),
               width=12,height=3.0,useDingbats=FALSE)

        }

    }

}

###########################

