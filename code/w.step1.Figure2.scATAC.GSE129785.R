#!/usr/bin/env Rscript

library("sscVis")
library("chromVARmotifs")
library("motifmatchr")
library("BSgenome.Hsapiens.UCSC.hg19")
library("data.table")
library("plyr")
library("Gviz")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("BiocParallel")
register(SerialParam())
set.seed(1)

genome <- BSgenome.Hsapiens.UCSC.hg19
data("human_pwms_v1")


out.prefix <- "./OUT_Fig2/scATAC.GSE129785/scATAC.GSE129785"
dir.create(dirname(out.prefix),F,T)

sc.TCells.file <- "../data/external/scATAC_TME_TCells_SummarizedExperiment.final.rds"
se <- readRDS(sc.TCells.file)
rownames(se) <- sprintf("%s.%s",rowData(se)[["SYMBOL"]],rowData(se)[["name"]])

uf <- readRDS("../data/external/scATAC.GSE129785.Unique_Peaks.rds")

########################## 
## plot function
plot.track.scATACSeq <- function(se.plot,
				 dat.scATAC.cls.tb,
				 out.prefix,
				 txdb,
				 ft.diff,
				 chr,
				 cluster.show=c("c12","c13","c14","c15","c16","c17"),
				 genome.mid=60031767,
				 genome.width=10000,
				 tf.name=NULL,
				 dat.out.tb=NULL,
				 lim.y=c(0,200),
				 cluster.id2name,cluster.id2col,highlight.tb,
				 pdf.width=8,pdf.height=7)
{
    {
	track.Tx <- GeneRegionTrack(txdb, chromosome = chr,
				    ##start = 60031767-1e6, end = 60031767+1e6,
				    start = genome.mid-genome.width, end = genome.mid+genome.width,
				    background.title="transparent",
				    size=1,
				    col.line="darkblue",col=NULL,fill="darkblue")

    }

    {
	#my.col.global <- "darkgray"
	my.col.global <- "black"

	#chr <- as.character(unique(seqnames(rowRanges(se.plot))))

	track.genome <- GenomeAxisTrack(cex=1,size=2.5)
	track.diffPeak <- AnnotationTrack(rowRanges(se.plot)[ft.diff,], name = "diffPeak",
				      stacking="squish",shape="box",
				      col.axis=my.col.global,
				      col.title=my.col.global,
				      background.title="transparent",
				      fill="#800000",
				      cex.title=1,
				      frame=T,
				      fontface.title=1,
				      rotation.title=0,
				      size=1,
				      col=NULL)

	##track.tf <- llply(tf.goi,function(x)
	##track.tf <- llply(c("ARID5B","IRF4","IRF8","ETV1","NR5A2","NFATC1","NFATC2"),function(x)
	##track.tf <- llply(c("ARID5B","IRF4","IRF8","ETV1","NR5A2"),function(x)

	track.tf <- NULL
	if(!is.null(tf.name) && !is.null(dat.out.tb)){
	    track.tf <- llply(tf.name,function(x)
			      {
					    dat.x <- dat.out.tb[TF==x,][!duplicated(name),]
					    dat.x[,chromosome:=seqnames]
					    AnnotationTrack(dat.x,name=x,shape="box",
								    col.axis=my.col.global,
								    col.title=my.col.global,
								    background.title="transparent",
								    fill=my.col.global,
								    cex.title=1,
								    frame=T,
								    size=1,
								    fontface.title=1,
								    rotation.title=0,
								    col=NULL)
				       })
	}


	track.list.data <- llply(cluster.show,function(cls.x){
					#cls.x <- "c17"
					dat.x <- dat.scATAC.cls.tb[T.Cluster==cls.x,]
					track.data <- DataTrack(data = dat.x[["y"]],
								start = dat.x[["start"]],
								end = dat.x[["end"]],
								chromosome = dat.x[["seqnames"]], genome = "hg19",
								col=NULL,
								col.axis=my.col.global,
								col.title=my.col.global,
								background.title="transparent",
								#col.border.title="black",
								col.histogram="transparent",
								frame=T,
								fontface.title=1,
								fontface=1,
								size=2.5,
								rotation.title=0,
								ylim=lim.y,
								cex.title=1,cex.axis=0.84,
								fill.histogram=cluster.id2col[cls.x],
								name = cluster.id2name[cls.x])
					return(track.data)
				   })

	track.ht <- HighlightTrack(trackList = c(list(track.diffPeak),track.tf,track.list.data),
			     start = highlight.tb$start,end=highlight.tb$end, chromosome = chr)


	pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
	plotTracks(c(track.genome,track.ht,track.Tx),type="histogram",
	###plotTracks(c(track.ht,track.Tx),type="histogram",
		   margin=42,
		   from=genome.mid-genome.width,to=genome.mid+genome.width)
	dev.off()

    }

}

## some settings
{
    cell2cluster <- structure(gsub("^Cluster","c",se$T_Cell_Cluster,perl=T),
			  names=colnames(se))
    cellInfo.tab <- data.table(cellID=colnames(se),
			   totalCountsPerCell=colSums(assay(se)),
			   T.Cluster=gsub("^Cluster","c",se$T_Cell_Cluster,perl=T))
    clsInfo.tab <- cellInfo.tab[,.(ncounts.tot=sum(totalCountsPerCell)),by=c("T.Cluster")]
    clsInfo.tab[,sf:=10e6/ncounts.tot]

    cluster.id2name <- c("c8"="Treg 1","c9"="Treg 2","c10"="Treg 3","c11"="Treg 4",
			 "c12"="Effector","c13"="Naive","c14"="Memory",
			 "c15"="Early Tex","c16"="Intermediate\nTex","c17"="Terminal\nTex")
    cluster.id2col <- c("c8"="#829ECD","c9"="#FFE50E","c10"="#263063","c11"="#007175",
			"c12"="#734995","c13"="#D365A1","c14"="#FF7209",
			 "c15"="#EFC1D7","c16"="#E6A35E","c17"="#72D8E4")

    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    tx.feature <- txdb %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique

}

####################### TOX (fig. S23B) ################################
{

    se.plot <- se[which(rowData(se)$SYMBOL=="TOX"),]

    uf$mat[grep("^TOX\\.",rownames(uf$mat),value=T),1:5]
    ft.TOX <- grep("^TOX\\.",rownames(uf$mat),value=T)

    sort(rowRanges(se.plot)[ft.TOX,])
    write.table(as.data.frame(sort(rowRanges(se.plot)[ft.TOX,]))[,c(1,2,3,6,4,5)],
		sprintf("%s.diffPeak.TOX.bed",out.prefix),
		row.names=F,sep="\t",quote=F,col.names=F)

    human_pwms_v1.tf <- sapply(seq_along(human_pwms_v1),function(i){ human_pwms_v1[[i]]@name })
    tf.goi <- c("ARID5B","CREM","ELF1","EOMES","ETV1", "ETV7","IKZF2", "IRF2","IRF4","IRF8",
		"NR5A2","PRDM1",
		"NFATC1","NFATC2")
    f.goi <- which(human_pwms_v1.tf %in% tf.goi)

    matches <- matchMotifs(human_pwms_v1[f.goi], rowRanges(se.plot), genome = "BSgenome.Hsapiens.UCSC.hg19",
			   p.cutoff=1e-04)

    ## peaks with matched motif
    dat.out.tb <- as.data.table(cbind(as.data.frame(rowRanges(matches)), as.matrix(assay(matches))))
    dat.out.tb <- melt(dat.out.tb,id.vars=c("seqnames","start","end","width","strand","name","SYMBOL","GC"))[value==T,]
    dat.out.tb[,variable:=as.character(variable)]
    dat.out.tb[,TF:=human_pwms_v1.tf[variable]]
    dat.out.tb[,start:=start-1]

    l_ply(unique(dat.out.tb$TF),function(x){
		write.table(dat.out.tb[TF==x,][,c(1,2,3,6,4,5)],
			    sprintf("%s.TOX.TF.%s.bed",out.prefix,x),
			    row.names=F,sep="\t",quote=F,col.names=F)
			   })

    ### peaks data
    dat.scATAC.tb <- as.data.table(cbind(as.data.frame(rowRanges(se.plot)), as.matrix(assay(se.plot))))
    dat.scATAC.tb <- melt(dat.scATAC.tb,id.vars=c("seqnames","start","end","width","strand",
						  "name","SYMBOL","GC"))
    dat.scATAC.tb[,variable:=as.character(variable)]

    dat.scATAC.tb[,T.Cluster:=cell2cluster[variable]]
    dat.scATAC.tb[1:2,]
    dat.scATAC.cls.tb <- dat.scATAC.tb[,.(ncounts=sum(value)),
				       by=c("seqnames","start","end","width","strand","name",
					    "SYMBOL","GC","T.Cluster")
				       ][order(T.Cluster,seqnames,start,end),]

    dat.scATAC.cls.tb <- merge(dat.scATAC.cls.tb,clsInfo.tab)
    dat.scATAC.cls.tb[,y:=ncounts*sf]
    dat.scATAC.cls.tb[,id:=sprintf("%s.%s",SYMBOL,name)]

    dat.scATAC.cls.tb[id %in% c("TOX.scATAC_274938"),]

    ############## plot
    {
	which(mcols(tx.feature)[,"tx_name"]=="uc003xtw.1") %>% tx.feature[.,]
	#GRanges object with 1 range and 2 metadata columns:
	#      seqnames    ranges strand |     tx_id     tx_name
	#         <Rle> <IRanges>  <Rle> | <integer> <character>
	#  [1]     chr8  60031767      - |     33508  uc003xtw.1
	#  -------
	#  seqinfo: 93 sequences (1 circular) from hg19 genome

	###
	seq.TOX <- DNAStringSet(genome$chr8[(60031767-10000):(60031767+10000)])
	names(seq.TOX) <- sprintf("chr8:%s-%s",60031767-10000,60031767+10000)
	writeXStringSet(seq.TOX,sprintf("%s.seq.TOX.fa",out.prefix))

	###

	as.data.table(as.data.frame(ranges(rowRanges(se.plot)[ft.TOX,])))[start > 60031767-10000 & end < 60031767+10000,]
	as.data.table(as.data.frame(ranges(rowRanges(se.plot)[ft.TOX,])))[start > 60031767-20000 & end < 60031767+20000,]

	highlight.tb <- rbind(dat.out.tb[TF=="ARID5B" & start > 60031767-10000 & end < 60031767+10000,][3,],
				  dat.out.tb[TF=="ETV1" & start > 60031767-10000 & end < 60031767+10000,],
				  dat.out.tb[TF=="NR5A2" & start > 60031767-10000 & end < 60031767+10000,][3,])
	highlight.tb <- data.table(seqnames="chr8",
				       start=c(60027545,60031546),
				       end=c(60028046,60032557))

	plot.track.scATACSeq(se.plot, dat.scATAC.cls.tb,
			     out.prefix=sprintf("%s.track.TOX.pdf",out.prefix),
			     txdb=txdb, ft.diff=ft.TOX,
			     chr="chr8", cluster.show=c("c12","c13","c14","c15","c16","c17"),
			     genome.mid=60031767, genome.width=10000,
			     #tf.name=c("ARID5B","IRF4","IRF8","ETV1","NR5A2"),
			     #genome.mid=60031767, genome.width=20000,
			     #tf.name=c("ARID5B","IRF4","IRF8","ETV1","ETV7","PRDM1","NR5A2"),
			     tf.name=c("ARID5B","ETV1","NR5A2"),pdf.width=8,pdf.height=6.5,
			     dat.out.tb=dat.out.tb,
			     cluster.id2name=cluster.id2name,cluster.id2col=cluster.id2col,
			     highlight.tb=highlight.tb)

    }

}


