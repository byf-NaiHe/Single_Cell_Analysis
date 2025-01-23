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

dir.metaInfo <- "../data/metaInfo"
dir.create((dir.metaInfo),F,T)
meta.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")

##### calculate frequency using only baseline samples
{

    getFreqTable <- function(tb,stype,type.return="mtx",group.var="meta.cluster")
    {
        out.tb <- as.data.table(ldply(c("P","N","T"),function(x){
                        x.tb <- plotDistFromCellInfoTable(tb[loc==x,], plot.type="none",
                                          cmp.var="cancerType",min.NTotal=30,
                                          group.var=group.var,donor.var="patient.uid")
                        x.tb[,loc:=x]
                        x.tb[,stype:=stype]
                      }))
        if(type.return=="tb"){
            return(out.tb)
        }else if(type.return=="mtx"){
            d.tb <- out.tb
            ht.tb <- dcast(d.tb,group.var~loc+donor.var,value.var="freq",fill=0)
            ht.mtx <- as.matrix(ht.tb[,-1])
            rownames(ht.mtx) <- ht.tb[[1]]
            print(ht.mtx[,1:3])
            return(ht.mtx)
        }
    }


    freq.CD8.ht.tb <- getFreqTable(meta.tb[usedForFreq=="Y" & treatment=="baseline" & stype=="CD8",],
                                   "CD8",type.return="tb")
    freq.CD4.ht.tb <- getFreqTable(meta.tb[usedForFreq=="Y" & treatment=="baseline" & stype=="CD4",],
                                   "CD4",type.return="tb")

    freq.all.ht.tb <- rbind(freq.CD8.ht.tb,freq.CD4.ht.tb)
    saveRDS(freq.all.ht.tb,file=sprintf("%s/panC.freq.all.ht.tb.rds",dir.metaInfo))
    saveRDS(freq.CD8.ht.tb,file=sprintf("%s/panC.freq.CD8.ht.tb.rds",dir.metaInfo))
    saveRDS(freq.CD4.ht.tb,file=sprintf("%s/panC.freq.CD4.ht.tb.rds",dir.metaInfo))

}

