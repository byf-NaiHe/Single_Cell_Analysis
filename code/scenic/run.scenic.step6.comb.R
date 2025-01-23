#!/usr/bin/env Rscript

library("sscClust")
library("SCopeLoomR")
library("SCENIC")
library("data.table")
library("tictoc")
library("plyr")
library("R.utils")
library("ggplot2")
library("ggrepel")
library("igraph")
library("ggraph")
library("ggrastr")
library("matrixStats")

source("../func.R")

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 16)

colSet <- readRDS("../../data/metaInfo/panC.colSet.list.rds")

geneSet.list <- fetechGeneSetList()

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

#out.prefix <- "OUT.scenic/grn/zhangLab10X.CD8"

plotRSS <- function(out.prefix,gene.desc,mcls,gene.highlight=NULL)
{
	set.seed(123456)
	dat.plot <- gene.desc[Group==mcls,][order(-rss.median),]
	dat.plot[,rrank:=seq_len(nrow(dat.plot))]
	if(is.null(gene.highlight)){
		dat.plot[,top:=rrank<=10]
	}else{
		dat.plot[,top:=(geneID %in% gene.highlight)]
	}
	p <- ggplot(dat.plot,aes(x=rrank,y=rss.median)) +
		geom_point_rast(aes(color=top),shape=16,show.legend=F,raster.width=2.5,raster.height=6) +
		scale_color_manual(values=c("TRUE"="#F8766D","FALSE"="#619CFF")) +
		geom_text_repel(data=dat.plot[top==T,],aes(label=geneID),color="#F8766D",
						nudge_x= 200, direction= "y",hjust=0,
						force=1) +
		xlab("Rank") + ylab("Regulon Specificity Score") +
		theme_pubr()
	ggsave(sprintf("%s.rss.rank.%s.pdf",out.prefix,mcls),width=2.2,height=3.5,useDingbats=F)
}


highlightRegulon <- function(out.prefix,gene.desc,reg.tb,motifs.tb,geneSig.tb,highlight.regulon=NULL,
							 par.node.text=list(repel=T,force=1,size=1.5),
							 show.motif=T,opt.layout="fr")
{
	if(is.null(highlight.regulon)){
		mcls <- sort(unique(gene.desc$Group))
		highlight.regulon <- llply(mcls,function(x){
					       gene.desc[Group==x,][order(-rss.median),head(geneID,n=10)]
								 })
		names(highlight.regulon) <- mcls
	}

	htmlOuput <- function(tableSubset, colsToShow,out.prefix,meta = NULL, cacheable = NA) {
		rmarkdown::render('./showMotifTb.rmd',
						  output_file = sprintf("%s.MotifTb.html",out.prefix))
	}
	print(str(highlight.regulon))
	l_ply(names(highlight.regulon),function(x){

			  reg.plot.tb <- reg.tb[TF %in% highlight.regulon[[x]],]
			  reg.plot.tb[,TF:=gsub("\\(\\+\\)","",TF)]

		      ###htmlOuput(motifs.tb[TF %in% gsub("\\(\\+\\)","",unname(unlist(highlight.regulon))),],
			  if(show.motif==T)
			  {
				  htmlOuput(motifs.tb[TF %in% reg.plot.tb$TF,],
							colsToShow=colnames(motifs.tb)[-c(4,10)],
							out.prefix=sprintf("%s.%s",out.prefix,x))
			  }

			  g <- graph_from_data_frame(reg.plot.tb, directed = TRUE, vertices = NULL)
			  g <- set.vertex.attribute(graph = g,name = 'isRegulator',
										value = vertex_attr(g,"name") %in% reg.plot.tb$TF)
			  g <- set.vertex.attribute(graph = g,name = 'isTF',
										value = vertex_attr(g,"name") %in% geneSet.list$TF )
			  Betweenness=igraph::betweenness(graph = g,directed=FALSE,normalized = TRUE)
			  g <- set.vertex.attribute(graph = g,name = 'Betweenness',
										value = Betweenness[vertex_attr(g,"name")])
			 
			  #opt.layout <- "dh"
			  #opt.layout <- "fr"
			  #opt.layout <- "graphopt"
			  #opt.layout <- "lgl"
			  set.seed(123456)
			  p <- ggraph(g,layout=opt.layout) + 
				  geom_edge_link(alpha=0.8,colour="lightgray",edge_width=0.1,
								 arrow = arrow(length = unit(1, 'mm'))) + 
				  geom_node_point(aes_string(
								shape="isRegulator",
								color="isTF",
								fill="isRegulator"),alpha=0.8,size=0.5) +
		          scale_shape_manual(values=c("TRUE"=21,"FALSE"=22)) +
		          #scale_shape_manual(values=c("TRUE"=16,"FALSE"=15)) +
				  scale_fill_manual(values=c("TRUE"="#F8766D","FALSE"="#619CFF")) +
				  scale_color_manual(values=c("TRUE"="#F8766D","FALSE"="#619CFF")) +
				  #guides(color=guide_legend(ncol=2))+
			      theme_void()
			  text.hi <- geneSig.tb[sig==T & meta.cluster==x,][ comb.ES>0.25 & rank(-comb.ES)<=50 ,]
			  p <- p + do.call(geom_node_text,c(list(mapping=aes_string(label="name",color="isRegulator"),
											   data=p$data[(p$data$Betweenness>0.000 & p$data$name %in% text.hi$geneSymbol) |
														   (p$data$isRegulator==T),],
											   show.legend=F), par.node.text))
			  ggsave(sprintf("%s.ggraph.%s.%s.pdf",out.prefix,x,opt.layout),
					 width=4.2, height=3.3)

								 },.parallel=F)

}


run.it <- function(out.prefix,mytitle,loom.file.list,sce.list,
				   rss.file.list,reg.file.list,geneSig.tb)
{
	#out.prefix <- "OUT.scenic/grn/zhangLab10X.CD8.pyScenic.activity.ht"
	dir.create(dirname(out.prefix),F,T)


	motifs.tb.list <- llply(names(reg.file.list),function(x){
				    motifs.tb <- data.table::fread(cmd=sprintf("gzip -cd %s | sed '1,3d' ",reg.file.list[[x]]),
								   header=F)
				    colnames(motifs.tb) <- c("TF","MotifID","AUC","Annotation","Context",
							     "MotifSimilarityQvalue","NES","OrthologousIdentity",
							     "RankAtMax","TargetGenes")
							return(motifs.tb)
						})
	names(motifs.tb.list) <- names(reg.file.list)

	mdata.tb.list <- llply(names(sce.list),function(x){
							   m.tb <- colData(sce.list[[x]])
							   rownames(m.tb) <- gsub("^CD[48]_","",rownames(m.tb))
							   m.tb$meta.cluster <- as.character(m.tb$meta.cluster)
							   return(m.tb)
				   })
	names(mdata.tb.list) <- names(sce.list)

	obj.list <- llply(names(loom.file.list),function(x){
						loom <- open_loom(loom.file.list[[x]], mode="r")
						#exprMat <- get_dgem(loom)
						#cellInfo <- get_cellAnnotation(loom)
						regulons_incidMat <- get_regulons(loom)
						regulons <- regulonsToGeneLists(regulons_incidMat)
						regulonsAUC <- get_regulonsAuc(loom)
						regulonsAucThresholds <- get_regulonThresholds(loom)
						close_loom(loom)

						mdata.tb <- mdata.tb.list[[x]]

						sce.auc <- ssc.build(assay(regulonsAUC))
						row.sd <- apply(assay(sce.auc),1,sd,na.rm=T)
						print(table(row.sd > 0))
						sce.auc <- sce.auc[row.sd>0,]
						f.cell <- intersect(colnames(sce.auc),rownames(mdata.tb))
						sce.auc <- sce.auc[,f.cell]
						mdata.tb <- mdata.tb[f.cell,]
						#mdata.tb <- mdata.tb[colnames(sce),]
						#all(colnames(sce.auc)==rownames(mdata.tb))
						colData(sce.auc) <- mdata.tb

						sce.auc <- ssc.scale(sce.auc,gene.id=rownames(sce.auc),assay.name="exprs",do.scale=T)

						gene.actFreq <- as.data.table(ldply(sort(unique(sce.auc$meta.cluster)),function(x){
								  obj.x <- sce.auc[,sce.auc$meta.cluster==x]
								  data.table(geneID=rownames(obj.x),
											 Group=x, 
											 actFreq=rowSums(assay(obj.x,"exprs") > 0)/ncol(obj.x))
													 }))

						####
						#regulons
						reg.tb <- as.data.table(ldply(names(regulons),function(x){ data.table(TF=x,target=regulons[[x]],w=1)  }))
						#write.table(reg.tb,file=sprintf("%s.reg.gene.txt",out.prefix),row.names=F,quote=F,sep="\t")
						
						return(list("reg.tb"=reg.tb,"sce.auc"=sce.auc,"gene.actFreq"=gene.actFreq))
				   })
	names(obj.list) <- names(loom.file.list)

	sce.auc.list <- llply(names(obj.list),function(x){ obj.list[[x]]$sce.auc })
	names(sce.auc.list) <- names(obj.list)

	### reg
	reg.tb <- merge(obj.list[[1]]$reg.tb,obj.list[[2]]$reg.tb,by=c("TF","target"),all=T,
					suffixes=sprintf(".%s",names(obj.list)[c(1,2)]))

	### actFreq
	actFreq.tb <- merge(obj.list[[1]]$gene.actFreq,obj.list[[2]]$gene.actFreq,by=c("geneID","Group"),
					suffixes=sprintf(".%s",names(obj.list)[c(1,2)]))
	actFreq.tb[,actFreq.median:=rowMedians(as.matrix(.SD)),.SDcol=grep("^actFreq\\.",colnames(actFreq.tb),value=T)]

	### rss
	rss.tb.list <- llply(names(rss.file.list),function(x){
						rss.tb <- fread(rss.file.list[[x]])
						cnames <- rss.tb[[1]]
						rss.mtx <- t(rss.tb[,-1])
						colnames(rss.mtx) <- cnames
						rss.mtx <- rss.mtx[,sort(colnames(rss.mtx))]
						f.gene <- apply(rss.mtx,1,function(x){ all(is.na(x)) })
						rss.mtx <- rss.mtx[!f.gene,]

						gene.rss <- as.data.table(ldply(colnames(rss.mtx),function(x){
										 #xx <- head(sort(rss.mtx[,x],decreasing=T),n=10)
										 xx <- rss.mtx[,x]
										 data.table(geneID=names(xx),
													Group=x,
													rss=xx,
													#rrank=rank(-xx)/length(xx))
													rrank=rank(-xx))
									}))
						return(gene.rss)
				   })
	names(rss.tb.list) <- names(rss.file.list)
	
	rss.tb <- merge(rss.tb.list[[1]],rss.tb.list[[2]],by=c("geneID","Group"),
					suffixes=sprintf(".%s",names(rss.tb.list)[c(1,2)]))
	rss.tb[,rss.median:=rowMedians(as.matrix(.SD)),.SDcol=grep("^rss\\.",colnames(rss.tb),value=T)]

	#rss.tb[Group=="CD8.c12.Tex.CXCL13",][order(rrank.zhangLab10X.CD8),]
	#rss.tb[Group=="CD8.c12.Tex.CXCL13",][order(-rss.median),]

	gene.desc <- merge(rss.tb,actFreq.tb,by=c("geneID","Group"))
	gene.desc <- gene.desc[order(Group,-rss.median),]
	#print(head(gene.desc[Group=="CD8.c12.Tex.CXCL13",],n=2))

	group.vec <- sort(unique(gene.desc$Group))
	gene.hl <- llply(group.vec,function(x){
			  ### top10
			  plotRSS(out.prefix,gene.desc,x,gene.highlight=NULL)
			  ### extrem value
			  res.outlier <- sscClust:::classify.outlier(gene.desc[Group==x,][["rss.median"]],
												 out.prefix=sprintf("%s.extremeV.%s",out.prefix,x))
			  th.x <- res.outlier$K$limit[["Right"]]
			  gene.hl.x <- gene.desc[Group==x & rss.median >= th.x,][["geneID"]]
			  if(length(gene.hl.x)>0){
				  plotRSS(out.prefix=sprintf("%s.extremeV",out.prefix),gene.desc,x,
						  gene.highlight=gene.hl.x)
			  }
			  ###
			  return(gene.hl.x)
					},.parallel=T)
	names(gene.hl) <- group.vec
	#
	print("XXXX")
	print(gene.hl)
	l_ply(names(gene.hl),function(x){
			  ### top 10
			  highlightRegulon(out.prefix=out.prefix,
							   gene.desc=gene.desc[Group==x,],
							   reg.tb=reg.tb,
							   motifs.tb=motifs.tb.list[[1]],
							   geneSig.tb=geneSig.tb,
							   highlight.regulon=NULL,
							   par.node.text=list(repel=T,force=1,size=1.5))
			  ### extremeV
			  print("Test:")
			  print(gene.hl[[x]])
			  if(length(gene.hl[[x]])>0){
				  highlightRegulon(out.prefix=sprintf("%s.extremeV",out.prefix),
								   gene.desc=gene.desc[Group==x,],
								   reg.tb=reg.tb,
								   motifs.tb=motifs.tb.list[[1]],
								   geneSig.tb=geneSig.tb,
								   highlight.regulon=gene.hl[x],
								   par.node.text=list(repel=T,force=1,size=1.5))
			  }
								 },.parallel=F)

	#####
	write.table(reg.tb,file=sprintf("%s.reg.tb.txt",out.prefix),quote=F,sep="\t",row.names=F)
	write.table(gene.desc,file=sprintf("%s.gene.desc.txt",out.prefix),quote=F,sep="\t",row.names=F)
	
	###### heatamp show regulons activities
	#gene.desc.flt <- gene.desc[order(Group,-rss,-actFreq),][actFreq>=0.2,][,head(.SD,n=10),by="Group"]
	#gene.desc.flt <- gene.desc.flt[!duplicated(geneID),]
	
	return(list("reg.tb"=reg.tb,
				"sce.auc.list"=sce.auc.list,
				"motifs.tb.list"=motifs.tb.list,
				#"gene.desc.flt"=gene.desc.flt,
				"gene.desc"=gene.desc))

}

sce.list <- list("zhangLab10X"=readRDS("./OUT.scenic/sce/scenic.zhangLab10X.sce.mini.rds"),
				 "zhangLabSS2"=readRDS("./OUT.scenic/sce/scenic.zhangLabSS2.sce.mini.rds")
				 )
sce.CD8.list <- llply(sce.list,function(x){ x[,x$stype=="CD8"] })
names(sce.CD8.list) <- c("zhangLab10X.CD8","zhangLabSS2.CD8")
sce.CD4.list <- llply(sce.list,function(x){ x[,x$stype=="CD4"] })
names(sce.CD4.list) <- c("zhangLab10X.CD4","zhangLabSS2.CD4")

geneSig.CD8.tb=readRDS("../../data/expression/CD8/integration/int.CD8.S35.gene.tb.rds")
geneSig.CD4.tb=readRDS("../../data/expression/CD4/integration/int.CD4.S35.gene.tb.rds")


res.CD8 <- run.it(out.prefix="OUT.scenic/grn.comb.pos/comb.CD8/comb.CD8.pyScenic",
					mytitle="CD8",
					loom.file.list=list("zhangLab10X.CD8"="OUT.scenic/grn.comb.pos/zhangLab10X.CD8.pyScenic.loom",
										"zhangLabSS2.CD8"="OUT.scenic/grn.comb.pos/zhangLabSS2.CD8.pyScenic.loom"),
					sce.list=sce.CD8.list,
					rss.file.list=list("zhangLab10X.CD8"="OUT.scenic/grn.comb.pos/zhangLab10X.CD8.pyScenic.RSS.meta.cluster.csv.gz",
								  "zhangLabSS2.CD8"="OUT.scenic/grn.comb.pos/zhangLabSS2.CD8.pyScenic.RSS.meta.cluster.csv.gz"),
					reg.file.list=list("zhangLab10X.CD8"="OUT.scenic/grn.comb.pos/zhangLab10X.CD8.reg.csv.gz",
								  "zhangLabSS2.CD8"="OUT.scenic/grn.comb.pos/zhangLabSS2.CD8.reg.csv.gz"),
				    geneSig.tb=geneSig.CD8.tb
					)
saveRDS(res.CD8,file="OUT.scenic/grn.comb.pos/comb.CD8/comb.CD8.pyScenic.rds")
#res.CD8 <- readRDS(file="OUT.scenic/grn.comb.pos/comb.CD8/comb.CD8.pyScenic.rds")

regulon.Tex.top <- res.CD8$gene.desc[Group=="CD8.c12.Tex.CXCL13",head(.SD$geneID,n=50)]
reg.Tex.top.tb <- res.CD8$reg.tb[TF %in% regulon.Tex.top & sprintf("%s(+)",target) %in% regulon.Tex.top,]
reg.Tex.top.tb[,TF:=gsub("[\\(\\+\\)]","",TF,perl=T)]

highlightRegulon(out.prefix="OUT.scenic/grn.comb.pos/comb.CD8/comb.CD8.pyScenic.Tex.top.network",
				 gene.desc=res.CD8$gene.desc,
				 reg.tb=reg.Tex.top.tb,
				 motifs.tb=NULL,
				 geneSig.tb=geneSig.CD8.tb,
				 #highlight.regulon=list("CD8.c12.Tex.CXCL13"=c("NR5A2(+)", "SOX4(+)", "ETV1(+)", "IRF5(+)")),
				 highlight.regulon=list("CD8.c12.Tex.CXCL13"=reg.Tex.top.tb$TF),opt.layout="graphopt",
				 par.node.text=list(repel=T,force=1,size=1.5),show.motif=F)

#aaa <- sscClust:::classify.outlier(res.CD8$gene.desc[Group=="CD8.c09.Tk.KIR2DL4",][["rss.median"]],
#								  out.prefix="OUT.scenic/grn.comb.pos/comb.CD8/comb.CD8.pyScenic.ttt")

#highlightRegulon(out.prefix="OUT.scenic/grn.comb.pos/comb.CD8/comb.CD8.pyScenic.t",
#				 gene.desc=res.CD8$gene.desc[Group=="CD8.c12.Tex.CXCL13",],
#				 reg.tb=res.CD8$reg.tb,
#				 motifs.tb=res.CD8$motifs.tb.list$zhangLab10X.CD8,
#				 geneSig.tb=geneSig.CD8.tb,
#				 highlight.regulon=NULL,
#				 par.node.text=list(repel=T,force=1,size=1.5))

res.CD4 <- run.it(out.prefix="OUT.scenic/grn.comb.pos/comb.CD4/comb.CD4.pyScenic",
					mytitle="CD4",
					loom.file.list=list("zhangLab10X.CD4"="OUT.scenic/grn.comb.pos/zhangLab10X.CD4.pyScenic.loom",
										"zhangLabSS2.CD4"="OUT.scenic/grn.comb.pos/zhangLabSS2.CD4.pyScenic.loom"),
					sce.list=sce.CD4.list,
					rss.file.list=list("zhangLab10X.CD4"="OUT.scenic/grn.comb.pos/zhangLab10X.CD4.pyScenic.RSS.meta.cluster.csv.gz",
								  "zhangLabSS2.CD4"="OUT.scenic/grn.comb.pos/zhangLabSS2.CD4.pyScenic.RSS.meta.cluster.csv.gz"),
					reg.file.list=list("zhangLab10X.CD4"="OUT.scenic/grn.comb.pos/zhangLab10X.CD4.reg.csv.gz",
								  "zhangLabSS2.CD4"="OUT.scenic/grn.comb.pos/zhangLabSS2.CD4.reg.csv.gz"),
				    geneSig.tb=geneSig.CD4.tb
					)
saveRDS(res.CD4,file="OUT.scenic/grn.comb.pos/comb.CD4/comb.CD4.pyScenic.rds")

###################################### edges table

make.reg.detail.table <- function(stype,trim.TF=T)
{
    #stype <- "CD8"
    out.prefix.tb <- sprintf("OUT.scenic/grn.comb.pos/comb.%s/comb.%s.pyScenic.reg.tb.detail",stype,stype)

    in.tb <- fread(sprintf("OUT.scenic/grn.comb.pos/comb.%s/comb.%s.pyScenic.reg.tb.txt",stype,stype))
    f.na <- is.na(in.tb[[sprintf("w.zhangLab10X.%s",stype)]])
    in.tb[[sprintf("w.zhangLab10X.%s",stype)]][f.na] <- 0
    f.na <- is.na(in.tb[[sprintf("w.zhangLabSS2.%s",stype)]])
    in.tb[[sprintf("w.zhangLabSS2.%s",stype)]][f.na] <- 0

    adj.tb <- readRDS(sprintf("OUT.scenic/grn/comb.%s.adj.rds",stype))
    colnames(adj.tb) <- c("TF","target","importance.zhangLab10X","importance.zhangLabSS2",
			  "importance.median","importance.mins")

    if(trim.TF){
	in.tb[,TF:=gsub("\\(\\+\\)","",TF,perl=T)]
    }

    out.tb <- merge(in.tb,adj.tb)
    conn <- gzfile(sprintf("%s.txt.gz",out.prefix.tb),"w")
    write.table(out.tb,conn,row.names=F,sep="\t",quote=F)
    close(conn)

}

make.reg.detail.table(stype="CD8",trim.TF=T)
make.reg.detail.table(stype="CD4",trim.TF=T)

######## correlation between TOX and others
####cor.dat.list <- readRDS("../OUT.GRN/CD8/panC.GRN.CD8.dat.cor.list.rds")
####
####all(colnames(cor.dat.list$dat.cor.p)==colnames(cor.dat.list$dat.cor.s))
####cor.TOX.tb <- data.table(gene.p=colnames(cor.dat.list$dat.cor.p),
####			 cor.pearson=cor.dat.list$dat.cor.p["TOX",],
####			 cor.spearman=cor.dat.list$dat.cor.s["TOX",])
####
####adj.TOX.a.tb <- adj.tb[TF=="TOX",]
####adj.TOX.a.tb <- merge(adj.TOX.a.tb,cor.TOX.tb,by.x="target",by.y="gene.p")
####adj.TOX.a.tb <- adj.TOX.a.tb[,c("TF","target",colnames(adj.TOX.a.tb)[3:8]),with=F][order(-importance.median),]
####
####adj.TOX.b.tb <- adj.tb[target=="TOX",]
####adj.TOX.b.tb <- merge(adj.TOX.b.tb,cor.TOX.tb,by.x="TF",by.y="gene.p")[order(-importance.median),]
####
####adj.TOX.tb <- rbind(adj.TOX.a.tb,adj.TOX.b.tb)
####adj.TOX.tb$gene.p <- adj.TOX.tb$TF
####adj.TOX.tb[TF=="TOX",gene.p:=target]
####
####geneSet.list <- fetechGeneSetList(1)
####
####for(x in names(geneSet.list)){
####    adj.TOX.tb[[sprintf("geneSet.%s",x)]] <- adj.TOX.tb[["gene.p"]] %in% geneSet.list[[x]]
####}
####
####write.table(adj.TOX.tb,sprintf("%s.TOX.adj.txt",out.prefix.tb),row.names=F,sep="\t",quote=F)
####




