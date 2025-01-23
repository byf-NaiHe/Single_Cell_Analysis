#!/usr/bin/env Rscript

library("data.table")
library("sscVis")
library("plyr")
library("ggplot2")
library("ggpubr")
library("ggplotify")
library("cowplot")
library("ComplexHeatmap")

out.prefix <- "OUT_Fig2/intraTex.TCRSharing/intraTex.TCRSharing"
dir.create(dirname(out.prefix),F,T)

###### load data
{
    colSet <- readRDS("../data/metaInfo/panC.colSet.list.rds")
    meta.tb <- readRDS("../data/metaInfo/panC.freq.all.meta.tb.rds")
    mapping.tb <- unique(meta.tb[,c("cancerType","dataset.old")])
    datasetOld2cancerType <- c(structure(mapping.tb$cancerType,names=mapping.tb$dataset.old),
                               "PanC"="PanC")

    ######## gene expression data ####
    sce.list.file <- "list/sce.CD8.fullPath.list"
    sce.list.tb <- fread(cmd=sprintf("awk '!/^#/' %s ",sce.list.file),head=T)

    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = 8)

    sce.list <- llply(seq_len(nrow(sce.list.tb)),function(i){ readRDS(sce.list.tb$scefile[i]) },.parallel=T)
    names(sce.list) <- sce.list.tb$data.id

    ######### TCR #######
    dat.tcr <- readRDS("../data/tcr/byCell/tcr.zhangLab.comb.flt.rds")
    setkey(dat.tcr,"Cell_Name")

    dat.tcr[dataset=="NSCLC.XinyiGuo2018",dataset:="LC.XinyiGuo2018"]
    dat.tcr[dataset=="STAD.BoxiKang2020",dataset:="STAD.BoxiKang2019"]
    dataset.ana.vec <- unique(dat.tcr$dataset)
    all(dataset.ana.vec %in% names(sce.list))
}
##############

##############
mat.patterSort <- function(dat.plot.mtx)
{
    old.colname <- colnames(dat.plot.mtx)
    if(ncol(dat.plot.mtx)==3){
	colnames(dat.plot.mtx) <- c("DP","SP.geneA","SP.geneB")
	dat.plot.mtx.pattern <- as.data.table(dat.plot.mtx > 0)
	dat.plot.mtx.pattern$rid <- rownames(dat.plot.mtx)
	dat.plot.mtx.pattern$N.DP <- dat.plot.mtx[,"DP"]
	dat.plot.mtx.pattern$N.SP.geneA <- dat.plot.mtx[,"SP.geneA"]
	dat.plot.mtx.pattern$N.SP.geneB <- dat.plot.mtx[,"SP.geneB"]
	dat.plot.mtx.pattern <- dat.plot.mtx.pattern[order(-DP,-SP.geneA,-SP.geneB,-N.DP,-N.SP.geneA,-N.SP.geneB),]
	dat.plot.mtx.pattern[,Group:=sprintf("%s%s%s",as.integer(DP),as.integer(SP.geneA),as.integer(SP.geneB))]
	dat.plot.mtx <- dat.plot.mtx[dat.plot.mtx.pattern$rid,,drop=F]
    }else if(ncol(dat.plot.mtx)==4){
	colnames(dat.plot.mtx) <- c("DP","SP.geneA","SP.geneB","DN")
	dat.plot.mtx.pattern <- as.data.table(dat.plot.mtx > 0)
	dat.plot.mtx.pattern$rid <- rownames(dat.plot.mtx)
	dat.plot.mtx.pattern$N.DP <- dat.plot.mtx[,"DP"]
	dat.plot.mtx.pattern$N.SP.geneA <- dat.plot.mtx[,"SP.geneA"]
	dat.plot.mtx.pattern$N.SP.geneB <- dat.plot.mtx[,"SP.geneB"]
	dat.plot.mtx.pattern$N.DN <- dat.plot.mtx[,"DN"]
	dat.plot.mtx.pattern <- dat.plot.mtx.pattern[order(-DP,-SP.geneA,-SP.geneB,-DN,
							   -N.DP,-N.SP.geneA,-N.SP.geneB,-N.DN),]
	dat.plot.mtx.pattern[,Group:=sprintf("%s%s%s%s",as.integer(DP),
					     as.integer(SP.geneA),as.integer(SP.geneB),
					     as.integer(DN))]
	dat.plot.mtx <- dat.plot.mtx[dat.plot.mtx.pattern$rid,,drop=F]
    }
    colnames(dat.plot.mtx) <- old.colname
    return(dat.plot.mtx)
}

run.it <- function(sce.list,dat.tcr,out.prefix,
		   gene.to.test=c("FOXP3","CXCL13"), dataset.x="OV.zhangLab5P",mcls="CD8.c12.Tex.CXCL13",
		   th.exprs=0.3, th.ncell=2, 
		   show.all.expanded.clone=F,
		   note.txt="TCRSharing",
		   min.ncell=10,loc.vec=c("P","N","T"),
		   pdf.width=5,pdf.height=5)
{

    dir.create(dirname(out.prefix),F,T)

    sce <- sce.list[[dataset.x]]
    
    #### filter out dataset with cell number in Tex less than 50
    if(sum(sce[["meta.cluster"]]==mcls & (sce[["loc"]] %in% loc.vec) ) < min.ncell){
	    return(NULL)
    }
    ####
    cat(sprintf("dataset (%s, %s)\n",dataset.x,note.txt))

    sce.plot <- ssc.scale(sce,gene.symbol=gene.to.test,adjB="batchV",do.scale=T)
    sce.plot <- sce.plot[,sce.plot$meta.cluster==mcls & sce.plot$loc %in% loc.vec]
    
    dat.plot <- as.data.frame(t(as.matrix(assay(sce.plot,"norm_exprs.scale"))))
    colnames(dat.plot) <- rowData(sce.plot)$display.name
    dat.plot <- dat.plot[,gene.to.test]
    #dat.plot.melt <- melt(as.matrix(dat.plot))
    #colnames(dat.plot.melt) <- c("cell","gene","exprs.z")

    if(all(is.na(dat.plot[,1])) || all(is.na(dat.plot[,2])) ||
       all(dat.plot[,1]==0) || all(dat.plot[,2]==0) ){
	    return(NULL)
    }

    for(gg in colnames(dat.plot)){
	sce.plot[[sprintf("bin.%s",gg)]] <- as.integer(dat.plot[[gg]] > th.exprs)
	if(all(sce.plot[[sprintf("bin.%s",gg)]]==0)) { return(NULL) }
    }


    gene.1 <- colnames(dat.plot)[1]
    gene.2 <- colnames(dat.plot)[2]
    sce.plot$xtype <- "other"
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==1 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==0] <- sprintf("SP.%s",gene.1)
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==0 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==1] <- sprintf("SP.%s",gene.2)
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==1 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==1] <- sprintf("DP")
    sce.plot$xtype[sce.plot[[sprintf("bin.%s",gene.1)]]==0 &
			       sce.plot[[sprintf("bin.%s",gene.2)]]==0] <- sprintf("DN")

    #print(table(sce.plot$xtype))


    #f.cell <- intersect(colnames(sce.plot),dat.tcr[["Cell_Name"]])
    dat.tcr.use <- dat.tcr[colnames(sce.plot),]
    sce.plot$cloneID <- dat.tcr.use$cloneID

    ### make a plot
    {

	unclass(table(sce.plot$xtype))
	dat.fisher <- c(sum(sce.plot$xtype=="DP"),
			sum(sce.plot$xtype==sprintf("SP.%s",gene.2)),
			sum(sce.plot$xtype==sprintf("SP.%s",gene.1)),
			sum(sce.plot$xtype=="DN"))
	res.fisher <- fisher.test(matrix(dat.fisher,ncol=2))
	#print(dat.fisher)
	#print(res.fisher)
	dat.fisher.tb <- as.data.table(matrix(dat.fisher,nrow=1))
	colnames(dat.fisher.tb) <- c("DP",sprintf("SP.%s",gene.2),sprintf("SP.%s",gene.1),"DN")
	dat.fisher.tb <- cbind(data.table(dataset=dataset.x),
			       dat.fisher.tb,
			       data.table(OR=res.fisher$estimate,
					  p.value=res.fisher$p.value))
	
	#cat(sprintf("???????????\n"))
	#sce.plot.debug <<- sce.plot
	ssc.plot.GeneDensity(sce.plot,
			     out.prefix=sprintf("%s.clone.xtype.%s.cor",out.prefix,note.txt),
			     gene.symbol=rowData(sce.plot)$display.name,
			     assay.name="norm_exprs.scale",
			     expT = c(th.exprs, th.exprs),
			     ##expT = c(0.5, 0.5),
			     my.title=sprintf("%s\n(OR:%4.4f, p: %4.4f)",dataset.x,
					      res.fisher$estimate,res.fisher$p.value))

	sce.plot.ret <- sce.plot
	#### sharing plot
	#sce.plot <- sce.plot[,sce.plot$xtype!="DN"]
	
	##### expression correlation
	###dat.test.norm <- as.data.frame(t(as.matrix(assay(sce.plot,"norm_exprs"))))
	###colnames(dat.test.norm) <- rowData(sce.plot)$display.name
	###res.cor <- cor.test(dat.test.norm[[1]],dat.test.norm[[2]])

	###
	sce.plot$xtype <- factor(sce.plot$xtype,levels=(c("DP", sprintf("SP.%s",gene.1), sprintf("SP.%s",gene.2), "DN")))
	clone.dist <- unclass(as.data.table(colData(sce.plot[,!is.na(sce.plot$cloneID)]))[,table(cloneID,xtype)])
	#clone.dist <- clone.dist[,levels(sce.plot$xtype)]

	#clone.dist <- clone.dist[,-which(colnames(clone.dist)=="DN"),drop=F]
	#clone.dist <- clone.dist[rowSums(clone.dist)>0,]

	### expanded clones with positive gene.2 cells
	f.clone <- (clone.dist[,sprintf("SP.%s",gene.2)] > 0 | clone.dist[,sprintf("DP")] > 0  ) &
			rowSums(clone.dist) >= th.ncell
#	if(show.all.expanded.clone){
#	    f.clone <- rowSums(clone.dist) >= th.ncell
#	}else{
#	    ### clones across multiple coluns or DP
#	    f.clone <- (rowSums(clone.dist > 0) > 1 )
#	    if("DP" %in% colnames(clone.dist)){
#		    f.clone <- f.clone | clone.dist[,"DP"] > 0
#	    }
#	}

	if(sum(f.clone)==0){
	    cat(sprintf("no clone found (%s, %s)\n",dataset.x,note.txt))
	    return(NULL)
	}
	dat.plot.mtx <- clone.dist[f.clone,,drop=F]

    
#	if((all(dat.plot.mtx[,"DP"]==0) && all(dat.plot.mtx[,sprintf("SP.%s",gene.1)]==0) ) ||
#	   (all(dat.plot.mtx[,"DP"]==0) && all(dat.plot.mtx[,sprintf("SP.%s",gene.2)]==0) )){
#	    return(NULL)
#	}
	
	dat.plot.mtx <- mat.patterSort(dat.plot.mtx)

	dat.plot.mtx.clamp <- dat.plot.mtx
	dat.plot.mtx.clamp[ dat.plot.mtx.clamp > 9 ] <- 9
	### numbers
	nclone.tb <- data.table(dataset=dataset.x,
				nclone.shared=sum((dat.plot.mtx[,"DP"] > 0) |
						  (dat.plot.mtx[,sprintf("SP.%s",gene.2)] > 0 &
						   dat.plot.mtx[,sprintf("SP.%s",gene.1)] > 0 ) ),
				nclone.gene1=sum((dat.plot.mtx[,"DP"] == 0) &
						 (dat.plot.mtx[,sprintf("SP.%s",gene.1)] > 0 &
						  dat.plot.mtx[,sprintf("SP.%s",gene.2)]==0)),
				nclone.gene2=sum((dat.plot.mtx[,"DP"] == 0) &
						 (dat.plot.mtx[,sprintf("SP.%s",gene.1)] == 0 &
						  dat.plot.mtx[,sprintf("SP.%s",gene.2)]>0)),
				nclone.null=sum((dat.plot.mtx[,"DP"] == 0) &
						(dat.plot.mtx[,sprintf("SP.%s",gene.1)] == 0 &
						 dat.plot.mtx[,sprintf("SP.%s",gene.2)]==0)))
	nclone.gene2 <- sum(dat.plot.mtx[,sprintf("SP.%s",gene.2)] > 0 | dat.plot.mtx[,"DP"] > 0)
	v.shared <- (dat.plot.mtx[,sprintf("SP.%s",gene.2)] > 0 & dat.plot.mtx[,sprintf("SP.%s",gene.1)] > 0 ) |
			     dat.plot.mtx[,"DP"] > 0
	nclone.shared <- sum(v.shared)
	####
	ann.row <- rowAnnotation(shared=v.shared,
			      col=list(shared=structure(RColorBrewer::brewer.pal(3,"Greys")[c(1,3)],
							       names=c("FALSE","TRUE"))),
			      annotation_legend_param=list(),
			      border=T)
	ht <- sscVis::plotMatrix.simple(dat.plot.mtx.clamp,mytitle=note.txt,
					out.prefix=NULL,
					col.ht=rev(structure(c("lightgray",sscVis:::getColorPaletteFromNameContinuous("YlOrRd")),
							     names=0:9)),
					show.number=dat.plot.mtx,
					par.heatmap = list(cex.row=0, left_annotation = ann.row),
					my.cell_fun=function(j, i, x, y, w, h, col){
					    nn <- dat.plot.mtx[i,j]
					    if(nn>25){
						n.fontsize <- 8
					    }else{
						n.fontsize <- 10
					    }
					    if(nrow(dat.plot.mtx)<=10 && nn>0){
						if(dat.plot.mtx[i,j]>=7){
						    grid::grid.text(dat.plot.mtx[i, j], x, y,gp=grid::gpar(col="white",
													   fontsize=n.fontsize))
						}else{
						    grid::grid.text(dat.plot.mtx[i, j], x, y,gp=grid::gpar(fontsize=n.fontsize))
						}
					    } },
					exp.name="CellNumber",returnHT=T,
					pdf.width=pdf.width,pdf.height=pdf.height)

    }
    return(list("ht"=ht,"dat.plot.mtx"=dat.plot.mtx,
			    "res.fisher"=res.fisher,"dat.fisher.tb"=dat.fisher.tb,
			    "sce.mcls"=sce.plot.ret,
			    "nclone.tb"=nclone.tb,
			    "nclone.show"=nrow(dat.plot.mtx),
			    "nclone.gene2"=nclone.gene2,
			    "nclone.shared"=nclone.shared))
}


##### 
gene.pairs.list <- list("Treg"=c("CXCL13","FOXP3"),
			"Tk"=c("CXCL13","KIR2DL3"),
			"Tk.2"=c("CXCL13","TXK"),
			"Test"=c("HAVCR2","CTLA4"),
			"Tc17.IL26"=c("CXCL13","IL26"),
			"Tc17.IL17A"=c("CXCL13","IL17A"),
			#"Tc17.IL17F"=c("CXCL13","IL17F"),
			#"Tc17.IL23R"=c("CXCL13","IL23R"),
			#"Tc17.SLC4A10"=c("CXCL13","SLC4A10"),
			"Tc17.RORC"=c("CXCL13","RORC"))

##### test it!
{
#####    tmp.aid <- "Tc17.IL26"
#####    tmp.x <- "BC.zhangLab5P"
#####    res.full <-l_ply(dataset.ana.vec,function(tmp.x){
#####	tmp.list <- run.it(sce.list,dat.tcr,out.prefix=sprintf("%s/%s.test/%s.%s", dirname(out.prefix),
#####										sprintf("ShareWith%s",tmp.aid),
#####										basename(out.prefix),tmp.x),
#####			   gene.to.test=gene.pairs.list[[tmp.aid]], dataset.x=tmp.x,
#####			   pdf.width=3, mcls="CD8.c12.Tex.CXCL13", 
#####			   #show.all.expanded.clone=T,
#####			   loc.vec=c("T"),
#####			   note.txt=sprintf("ShareWith%s.all",tmp.aid))
#####			})
}

##### run it! (fig. S20B, S21B )
{
    res.all <- llply(names(gene.pairs.list),function(aid){
        res.full <-llply(dataset.ana.vec,function(x){
                     run.it(sce.list,dat.tcr,out.prefix=sprintf("%s/%s/%s.%s",
                                        dirname(out.prefix),
                                        sprintf("ShareWith%s",aid),
                                        basename(out.prefix),x),
                        gene.to.test=gene.pairs.list[[aid]],
                        dataset.x=x,
                        pdf.width=3,
                        loc.vec=c("T"),
                        mcls="CD8.c12.Tex.CXCL13",
                        ####show.all.expanded.clone=T,
                        note.txt=sprintf("ShareWith%s.all",aid))
                },.parallel=F)
        names(res.full) <- dataset.ana.vec
        cor.full.tb <- as.data.table(ldply(dataset.ana.vec,function(x){
                               if(!is.null(res.full[[x]])){
                               res.full[[x]]$dat.fisher.tb
                               } }))
        nclone.all.tb <- as.data.table(ldply(dataset.ana.vec,function(x){
                             if(!is.null(res.full[[x]])){
                                 res.full[[x]]$nclone.tb
                             } }))
        tmp.tb <- as.data.table(matrix(colSums(nclone.all.tb[,-1]),nrow=1))
        colnames(tmp.tb) <- colnames(nclone.all.tb)[-1]
        nclone.all.tb <- nclone.all.tb[order(nclone.shared/(nclone.shared+nclone.gene2+nclone.gene1+nclone.null)),]
        nclone.all.tb <- rbind(nclone.all.tb,cbind(data.table(dataset="PanC"),tmp.tb))
        nclone.all.tb[,cancerType:=datasetOld2cancerType[dataset]]
        nclone.all.tb[,NClones:=rowSums(nclone.all.tb[,2:5])]
        colnames(nclone.all.tb)[2:5] <- c("Both",
                          sprintf("Only.%s",gene.pairs.list[[aid]][1]),
                          sprintf("Only.%s",gene.pairs.list[[aid]][2]),
                          "Neither")
            nclone.all.freq.mtx <- as.matrix(nclone.all.tb[,2:5])
            nclone.all.freq.mtx <- sweep(nclone.all.freq.mtx,1,rowSums(nclone.all.freq.mtx),"/")
            nclone.all.freq.tb <- cbind(nclone.all.tb[,"cancerType"],as.data.table(nclone.all.freq.mtx))
            nclone.all.melt.tb <- melt(nclone.all.freq.tb,id.vars="cancerType")
            nclone.all.melt.tb[,cancerType:=factor(cancerType,levels=nclone.all.tb$cancerType)]
            nclone.all.melt.tb[,variable:=factor(variable,levels=colnames(nclone.all.tb)[c(5,3,4,2)])]

        p3 <- ggbarplot(nclone.all.melt.tb,x="cancerType",y="value",group="variable",
                       #position=position_dodge2(width=0.8),
                       fill="variable",color=NA,
                       palette="npg") +
                geom_text(data=nclone.all.tb, aes(label=NClones),y=1,vjust=-0.1) +
                coord_cartesian(clip = "off") +
                labs(x="",y="Frequency")+
                theme(legend.position="right",
                      axis.text.x=element_text(angle=45,hjust=1))
        #ggsave(sprintf("%s.%s.bar.pdf",out.prefix,aid),plot=p3,width=6,height=2.0)

        #####
        sce.fisher.list <- llply(dataset.ana.vec,function(x){
                         if(!is.null(res.full[[x]])){
                         sce.ret <- res.full[[x]]$sce.mcls
                         dat.tmp <- assay(sce.ret,"norm_exprs.scale")
                         rownames(dat.tmp) <- unname(rowData(sce.ret)$display.name)
                         dat.tmp <- dat.tmp[gene.pairs.list[[aid]],]
                         sce.ret <- ssc.build(dat.tmp,assay.name="norm_exprs.scale")
                         return(sce.ret)
                      } })
            names(sce.fisher.list) <- dataset.ana.vec
            f.dataset <- sapply(dataset.ana.vec,function(x){ !is.null(res.full[[x]]) })
            sce.tmp.list <- sce.fisher.list[ dataset.ana.vec[f.dataset] ]
            sce.fisher <- do.call(cbind,sce.tmp.list)

            fisher.mtx <- matrix(c(sum(cor.full.tb[["DP"]]),
                       sum(cor.full.tb[[ sprintf("SP.%s",gene.pairs.list[[aid]][1]) ]]),
                       sum(cor.full.tb[[ sprintf("SP.%s",gene.pairs.list[[aid]][2]) ]]),
                       sum(cor.full.tb[["DN"]])),ncol=2)
            res.fisher <- fisher.test(fisher.mtx)

            ssc.plot.GeneDensity(sce.fisher,
                     out.prefix=sprintf("%s/%s/%s.%s",
                            dirname(out.prefix),
                            sprintf("ShareWith%s",aid),
                            basename(out.prefix),"panC"),
                     gene.symbol=rowData(sce.fisher)$display.name,
                     assay.name="norm_exprs.scale",
                     expT = c(0.3,0.3),
                     my.title=sprintf("%s(OR:%4.4f, p: %s)","panC",
                              res.fisher$estimate,res.fisher$p.value))

            ht.list <- llply(dataset.ana.vec,function(x){
                     if(!is.null(res.full[[x]]$ht)){
                     as.ggplot(res.full[[x]]$ht,scale=0.95)
                     } })
            names(ht.list) <- dataset.ana.vec
            ht.list <- ht.list[!sapply(ht.list,is.null)]
            ht.list <- ht.list[ order(sapply(names(ht.list),function(x){ res.full[[x]]$nclone.shared })) ]
            nclone.show <- sapply(names(ht.list),function(x){ res.full[[x]]$nclone.show })
            nclone.shared <- sapply(names(ht.list),function(x){ res.full[[x]]$nclone.shared })
            nclone.gene2 <- sapply(names(ht.list),function(x){ res.full[[x]]$nclone.gene2 })
        ht.list.nrow <- ceiling(length(ht.list)/4)
            p <- plot_grid(plotlist=ht.list,
                   hjust=-0.05,
                   ncol=4,
                   label_size=10,labels=sprintf("%s(%d,%d/%d)", names(ht.list),
                                nclone.show,nclone.shared,nclone.gene2) )
            ggsave(sprintf("%s/%s/%s.%s.ht.pdf",
                       dirname(out.prefix),
                       sprintf("ShareWith%s",aid),
                       basename(out.prefix),"panC"),
               width=10,height=3*ht.list.nrow)

            dat.plot.merged.mtx <- do.call(rbind,llply(names(res.full),function(x){ res.full[[x]]$dat.plot.mtx }))
        nclone.show.panC <- sum(nclone.show)
        nclone.shared.panC <- sum(nclone.shared)
        nclone.gene2.panC <- sum(nclone.gene2)
        #nclone.show.panC <- sum(sapply(names(res.full),function(x){ res.full[[x]]$nclone.show }))
        #nclone.shared.panC <- sum(sapply(names(res.full),function(x){ res.full[[x]]$nclone.shared }))
        #nclone.gene2.panC <- sum(sapply(names(res.full),function(x){ res.full[[x]]$nclone.gene2 }))
            dat.plot.merged.mtx <- mat.patterSort(dat.plot.merged.mtx)

            dat.plot.merged.mtx.clamp <- dat.plot.merged.mtx
            dat.plot.merged.mtx.clamp[ dat.plot.merged.mtx.clamp > 9 ] <- 9

        v.shared.panC <- dat.plot.merged.mtx[,1] > 0 | (dat.plot.merged.mtx[,2] > 0 & dat.plot.merged.mtx[,3] > 0)
        ann.row.panC <- rowAnnotation(shared=v.shared.panC,
                      col=list(shared=structure(RColorBrewer::brewer.pal(3,"Greys")[c(1,3)],
                                       names=c("FALSE","TRUE"))),
                      annotation_legend_param=list(),
                      border=T)
            ht.all <- sscVis::plotMatrix.simple(dat.plot.merged.mtx.clamp,
                            mytitle=sprintf("panC(%d,%d/%d)",nclone.show.panC,
                                    nclone.shared.panC,nclone.gene2.panC),
                        out.prefix=sprintf("%s/%s/%s.%s.ht.all.pdf",
                                   dirname(out.prefix),
                                   sprintf("ShareWith%s",aid),
                                   basename(out.prefix),"panC"),
                        col.ht=rev(structure(c("lightgray",sscVis:::getColorPaletteFromNameContinuous("YlOrRd")),
                                     names=0:9)),
                        show.number=dat.plot.merged.mtx,
                        par.heatmap = list(cex.row=0,left_annotation=ann.row.panC),
                        my.cell_fun=function(j, i, x, y, w, h, col){
                            nn <- dat.plot.merged.mtx[i,j]
                            if(nrow(dat.plot.merged.mtx)<35 && nn > 0){
                            if(nn>=7){
                                grid::grid.text(nn, x, y,
                                        gp=grid::gpar(col="white",fontsize=10))
                            }else{
                                grid::grid.text(nn, x, y,
                                        gp=grid::gpar(fontsize=10))
                            }
                            }
                        },
                        exp.name="CellNumber",returnHT=T,
                        pdf.width=3.5,pdf.height=8)


        return(list("res.full"=res.full,"cor.full.tb"=cor.full.tb,
                    "ht.list"=ht.list,
                    "nclone.shared"=nclone.shared,
                    "nclone.show"=nclone.show,
                    "nclone.gene2"=nclone.gene2,
                    "nclone.all.tb"=nclone.all.tb,
                    "ht.all"=ht.all))
    })
    names(res.all) <- names(gene.pairs.list)

}

######################################


