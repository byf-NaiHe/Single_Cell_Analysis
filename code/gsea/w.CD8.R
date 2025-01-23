#!/usr/bin/env Rscript

library("plyr")
library("doParallel")
library("data.table")
source("../func.R")
#####

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(10)

cur.dir <- system2("pwd",stdout=T)
gene.tb.file <- "../../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds"
out.prefix <- sprintf("%s/OUT.gsea/gsea.CD8",cur.dir)
dir.create(dirname(out.prefix),F,T)
sh.dir <- "./sh.gsea.CD8"
dir.create(sh.dir,F,T)

### scoring_scheme of GSEA: classic,weighted,weighted_p2,weighted_p1.5
args.scoring.scheme <- "classic"
args.max <- 500

###### 
gene.desc.top <- readRDS(gene.tb.file)

###################################
## C2.CP: curated gene sets, Canonical pathways
db.file.list <- list(
					 "C2.CP"=sprintf("%s/../../data/external/c2.cp.v7.0.symbols.gmt",cur.dir),
                     "KEGG.flt"=sprintf("%s/../../data/external/MSigDB.KEGG.flt.gmt",cur.dir)
					 )

l_ply(names(db.file.list),function(x){
		  gen.gsea.script(gene.desc.top,
						  sprintf("%s/%s",sh.dir,x),
						  sprintf("%s/%s/%s",dirname(out.prefix),x,basename(out.prefix)),
						  db.file.list[[x]])
					 })

