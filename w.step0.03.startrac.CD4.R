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

in.file <- "../data/tcr/byCell/tcr.zhangLab.comb.flt.rds"
dir.startrac <- "../data/tcr/startrac/CD4"
dir.create(dir.startrac,F,T)

g.colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
in.dat <- readRDS(in.file)

in.dat$majorCluster <- as.character(in.dat$meta.cluster)
in.dat$clone.id <- in.dat$cloneID

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

{

    #### filter out patient.cluster with number of cell < 10
    ncell.patient.cluster <- sort(unclass(in.dat[,table(sprintf("%s.%s",patient,majorCluster))]))
    #ncell.patient.cluster <- sort(unclass(in.dat[stype=="CD4",table(sprintf("%s.%s",patient,majorCluster))]))
    in.dat <- in.dat[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
    ####
    dim(in.dat)
    ### [1] 164271     33

    if(!file.exists(sprintf("%s/CD4.out.startrac.nperm1000.rds",dir.startrac))){
        tic("Startrac.run")
        out <- Startrac.run(in.dat[stype=="CD4",], proj="panC",verbose=F,cores=6,n.perm=1000)
        toc()
        saveRDS(out,sprintf("%s/CD4.out.startrac.nperm1000.rds",dir.startrac))
    }
    out <- readRDS(sprintf("%s/CD4.out.startrac.nperm1000.rds",dir.startrac))

}


for(a.loc in c("T","P","N"))
{
    if(!file.exists(sprintf("%s/CD4.out.only%s.startrac.nperm1000.rds",dir.startrac,a.loc))){
        tic("Startrac.run")
        in.dat.flt <- in.dat[stype=="CD4" & loc==a.loc,]
        ncell.patient.cluster <- sort(unclass(in.dat.flt[,table(sprintf("%s.%s",patient,majorCluster))]))
        in.dat.flt <- in.dat.flt[ncell.patient.cluster[sprintf("%s.%s",patient,majorCluster)]>=10,]
        out.onlyALoc <- Startrac.run(in.dat.flt, proj="panC",verbose=F,cores=10,n.perm=1000)
        toc()
        saveRDS(out.onlyALoc,sprintf("%s/CD4.out.only%s.startrac.nperm1000.rds",dir.startrac,a.loc))
    }
}



############
cancerType.vec <- in.dat[stype=="CD4",unique(cancerType)]

if(!file.exists(sprintf("%s/CD4.out.startrac.byCancerType.rds",dir.startrac)))
{
    res.byCancerType <- llply(cancerType.vec,function(x){
                      Startrac.run(in.dat[stype=="CD4" & cancerType==x,],
                           proj=x,verbose=F,cores=4,n.perm=NULL)
                            })
    names(res.byCancerType) <- cancerType.vec
    saveRDS(res.byCancerType,sprintf("%s/CD4.out.startrac.byCancerType.rds",dir.startrac))
}

if(!file.exists(sprintf("%s/CD4.out.startrac.byCancerType.onlyT.rds",dir.startrac)))
{
    res.byCancerType.onlyT <- llply(cancerType.vec,function(x){
                        Startrac.run(in.dat[stype=="CD4" & cancerType==x & loc=="T",],
                             proj=x,verbose=F,cores=4,n.perm=NULL)
                            })
    names(res.byCancerType.onlyT) <- cancerType.vec
    saveRDS(res.byCancerType.onlyT,sprintf("%s/CD4.out.startrac.byCancerType.onlyT.rds",dir.startrac))
}
############







