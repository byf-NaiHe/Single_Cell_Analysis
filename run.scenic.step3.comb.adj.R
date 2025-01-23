#!/usr/bin/env Rscript

library("data.table")
library("plyr")
library("matrixStats")

file.CD4.list <- list("zhangLab10X.CD4"="OUT.scenic/grn/zhangLab10X.CD4.adj.csv.gz",
					  "zhangLabSS2.CD4"="OUT.scenic/grn/zhangLabSS2.CD4.adj.csv.gz")
				
file.CD8.list <- list("zhangLab10X.CD8"="OUT.scenic/grn/zhangLab10X.CD8.adj.csv.gz",
					  "zhangLabSS2.CD8"="OUT.scenic/grn/zhangLabSS2.CD8.adj.csv.gz")
				
combAdj <- function(adj.file.list,out.prefix)
{
	adj.tb.list <- llply(names(adj.file.list),function(x){ fread(adj.file.list[[x]]) })
	names(adj.tb.list) <- names(adj.file.list)

	adj.comb.tb <- merge(adj.tb.list[[1]],adj.tb.list[[2]],by=c("TF","target"))
	adj.comb.tb[,importance.median:=rowMedians(as.matrix(.SD[,c("importance.x","importance.y"),with=F]),na.rm=T) ]
	adj.comb.tb[,importance.mins:=rowMins(as.matrix(.SD[,c("importance.x","importance.y"),with=F]),na.rm=T) ]
	adj.comb.tb <- adj.comb.tb[order(TF,-importance.median),]
	adj.comb.flt <- adj.comb.tb[importance.median > 1,]
	saveRDS(adj.comb.flt,file=sprintf("%s.rds",out.prefix))
	adj.comb.out <- adj.comb.flt[,c("TF","target","importance.median"),with=F]
	colnames(adj.comb.out)[3] <- "importance"
	conn <- gzfile(sprintf("%s.csv.gz",out.prefix),"w")
	write.table(adj.comb.out,file=conn,quote=F,sep=",",row.names=F)
	close(conn)
}

combAdj(file.CD4.list,out.prefix=sprintf("OUT.scenic/grn/comb.CD4.adj"))

combAdj(file.CD8.list,out.prefix=sprintf("OUT.scenic/grn/comb.CD8.adj"))


