
changeSomeNames <- function(obj,col.mcls="meta.cluster",col.ctype="cancerType",col.dataset="dataset")
{
    #### change cluster name here
    if(!is.null(obj[[col.mcls]])){
#		obj[[sprintf("%s.old",col.mcls)]] <- obj[[sprintf("%s",col.mcls)]]
#		obj[[col.mcls]][ obj[[sprintf("%s.old",col.mcls)]]=="CD8.c05.Tem.TNF" ] <- "CD8.c05.Tact.TNF"
#		obj[[col.mcls]][ obj[[sprintf("%s.old",col.mcls)]]=="CD4.c09.Trm.CAPG" ] <- "CD4.c09.Tm.CAPG"
#		obj[[col.mcls]][ obj[[sprintf("%s.old",col.mcls)]]=="CD4.c04.Tm.TNF" ] <- "CD4.c04.Tact.TNF"
    }

    #### change cancerType here
    if(!is.null(obj[[col.ctype]])){
	    obj[[sprintf("%s.old",col.ctype)]] <- obj[[sprintf("%s",col.ctype)]]
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="BC" ] <- "BRCA"
	    ##obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="LUNG" ] <- "NSCLC"
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="LUNG" ] <- "LC"
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="Melanoma" ] <- "MELA"
	    obj[[col.ctype]][ obj[[sprintf("%s.old",col.ctype)]]=="BM" & obj[[col.dataset]]=="AML.PeterVanGalen2019" ] <- "AML"
    }

    #### change dataset here
    if(!is.null(obj[[col.dataset]])){
	    obj[[sprintf("%s.old",col.dataset)]] <- obj[[sprintf("%s",col.dataset)]]
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CRC.ZiyiLi10X" ] <- "CRC.LeiZhang2020_10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.YaoHeSS2" ] <- "HCC.QimingZhang2019_SS2"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.YaoHe10X" ] <- "HCC.QimingZhang2019_10X"
	    ######## ???
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.YaoHe10X" ] <- "HCC.QimingZhang2019_10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.LichunMa2019" ] <- "LIHC.LichunMa2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.LichunMa2019" ] <- "LIHC.LichunMa2019"
	    #########
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.Diether2018" ] <- "NSCLC.DietherLambrechts2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.QianqianSong2019" ] <- "NSCLC.QianqianSong2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.RapolasZilionis2019" ] <- "NSCLC.RapolasZilionis2019"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.LivnatJerby-Arnon2018" ] <- "MELA.LivnatJerby-Arnon2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.MosheSade-Feldman2018" ] <- "MELA.MosheSade-Feldman2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="Melanoma.HanjieLi2018" ] <- "MELA.HanjieLi2019"
	    ##obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="STAD.BoxiKang2019" ] <- "STAD.BoxiKang2020"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="STAD.BoxiKang2019" ] <- "STAD.BoxiKang2021"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="HCC.zhangLabSS2" ] <- "HCC.ChunhongZheng2017"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="LUNG.zhangLabSS2" ] <- "NSCLC.XinyiGuo2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CRC.zhangLabSS2" ] <- "CRC.LeiZhang2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Elham2018.10X" ] <- "BRCA.ElhamAzizi2018_10X"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Elham2018.Indrop" ] <- "BRCA.ElhamAzizi2018_InDrop"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.Peter2018" ] <- "BRCA.PeterSavas2018"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="BC.zhangLab5P" ] <- "BRCA.zhangLab5P"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="CHOL.zhangLabSS2" ] <- "CHOL.thisStudy"
	    obj[[col.dataset]][ obj[[sprintf("%s.old",col.dataset)]]=="NPC.XiliangWang2019" ] <- "NPC.YangLiu2021"
	    obj[[col.dataset]] <- gsub("zhangLab5P","thisStudy",obj[[col.dataset]])
    }
    
    return(obj)
}


correctCellInfo <- function(cellInfo.tb)
{
	#cellInfo.tb[,sampleID:=""]
	#cellInfo.tb[,treatment:=""]
	cellInfo.tb$sampleID <- ""
	cellInfo.tb$treatment <- ""
	cellInfo.tb[dataset=="AML.PeterVanGalen2019",sampleID:=batchV]
	cellInfo.tb[dataset=="AML.PeterVanGalen2019" & grepl("-D0",batchV),treatment:="baseline"]
	cellInfo.tb[dataset=="AML.PeterVanGalen2019" & !grepl("-D0",batchV),treatment:="post.treatment"]
	cellInfo.tb[dataset=="AML.PeterVanGalen2019" & grepl("^BM",batchV),treatment:="normal"]

	cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset),sampleID:=libraryID]
	cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset) & grepl("\\.pre\\b",sampleID),
				treatment:="baseline"]
	cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset) & grepl("\\.post\\b",sampleID),
				treatment:="post.treatment"]

	cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset),sampleID:=libraryID]
	cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset) & grepl("\\.pre\\b",sampleID),
				treatment:="baseline"]
	cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset) & grepl("\\.post\\b",sampleID),
				treatment:="post.treatment"]

	cellInfo.tb[grepl("^NSCLC.RapolasZilionis2019",dataset),sampleID:=sprintf("%s.%s",patient,loc)]
	cellInfo.tb[grepl("^NSCLC.RapolasZilionis2019",dataset) & patient=="p2",treatment:="post.treatment"]
	cellInfo.tb[grepl("^NSCLC.RapolasZilionis2019",dataset) & patient!="p2",treatment:="baseline"]

	cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset),sampleID:=patient]
	cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset) & 
				patient %in% c("C25","C39","C56","C60","C66",
							   "H21","H23","H28","H30","H38","H34","H65"),
				treatment:="baseline"]
	cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset) &
				patient %in% c("C26","C29","C35","C42","C46","H18","H37"),
				treatment:="post.treatment"]
	########### ???
	#f.dataset <- grepl("^LIHC.LichunMa2019",cellInfo.tb$dataset) & cellInfo.tb$cancerType=="CHOL"
	#cellInfo.tb$dataset[f.dataset] <- "CHOL.LichunMa2019"
	#f.dataset <- grepl("^LIHC.LichunMa2019",cellInfo.tb$dataset) & cellInfo.tb$cancerType=="HCC"
	#cellInfo.tb$dataset[f.dataset] <- "HCC.LichunMa2019"
	###########

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Melanoma.HanjieLi2018/Melanoma.HanjieLi2018.cellInfo.txt")
	tmp.id.mapping <- structure(tmp.tb$sampleID,names=tmp.tb$cellID)
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset),sampleID:=tmp.id.mapping[cellID] ]
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset) & 
				patient %in% c("p1","p11","p12","p13","p15","p16","p17",
							   "p18","p19","p21","p23","p24","p25","p26","p3"),
			treatment:="baseline"]
	cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset) & 
				patient %in% c("p10","p2","p20","p27","p28","p4","p5","p6","p8","p9"),
			treatment:="post.treatment"]

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Melanoma.LivnatJerby-Arnon2018/Melanoma.LivnatJerby-Arnon.cellInfo.txt")
	tmp.tb[treatment.group=="treatment.naive",treatment.group:="baseline"]
	tmp.treat.mapping <- structure(tmp.tb$treatment.group,names=tmp.tb$cellID)
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),sampleID:=patient]
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),treatment:=tmp.treat.mapping[cellID]]
	cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset) & patient %in% c("Mel129pa","Mel129pb"),patient:="Mel129"]

	tmp.tb <- fread("/lustre1/zeminz_pkuhpc/zhenglt/work/panC/data/Melanoma.MosheSade-Feldman2018/sampleInfo.txt")
	tmp.id.mapping <- structure(tmp.tb$oID,names=tmp.tb$title)
	tmp.tb[,treatment.group:=""]
	tmp.tb[Timepoint=="Pre",treatment.group:="baseline"]
	tmp.tb[Timepoint=="Post",treatment.group:="post.treatment"]
	tmp.treat.mapping <- structure(tmp.tb$treatment.group,names=tmp.tb$title)
	cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),sampleID:=tmp.id.mapping[cellID] ]
	cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),treatment:=tmp.treat.mapping[cellID]]

	######## treatment naive only samples
	f.cell <- cellInfo.tb[,(grepl("^(HNSCC.SidharthVPuram2017|CHOL.thisStudy|CRC.LeiZhang2018|HCC.ChunhongZheng2017|NSCLC.XinyiGuo2018)",dataset))]
	##cellInfo.tb[f.cell ,libraryID:=sprintf("%s.%s",patient,loc)]
	cellInfo.tb[f.cell, sampleID:=sprintf("%s%s",patient,loc)]
	cellInfo.tb[f.cell, treatment:="baseline"]

	f.cell <- cellInfo.tb[, grepl("^BRCA.ElhamAzizi2018",dataset) & patient %in% c("BC9","BC10","BC11")]
	cellInfo.tb[f.cell,libraryID:=sprintf("%s%s",patient,loc)]
							
	f.cell <- cellInfo.tb[,grepl("^(BRCA.ElhamAzizi2018|BRCA.PeterSavas2018|CRC.LeiZhang2020.10X|HCC.QimingZhang2019.10X|HCC.QimingZhang2019|NPC.YangLiu2021|NSCLC.DietherLambrechts2018|NSCLC.QianqianSong2019|PACA.JunyaPeng2019|RC.MatthewDYoung2018|STAD.BoxiKang2020)",dataset) | 
			( grepl("thisStudy",dataset) & dataset!="CHOL.thisStudy" ) ]
	cellInfo.tb[f.cell, sampleID:=libraryID]
	cellInfo.tb[f.cell, treatment:="baseline"]

	#### correct cancerType
	cellInfo.tb[dataset=="HCC.QimingZhang2019.10X" & patient=="D20171215",cancerType:="CHOL"]
	####### ???
	##cellInfo.tb[dataset=="HCC.QimingZhang2019.10X" & patient=="D20171215",dataset:="CHOL.QimingZhang2019.10X"]
	#######

	#### add stype
	cellInfo.tb[,stype:=sapply(strsplit(as.character(meta.cluster),"\\."),"[",1)]
	#### add patient.uid
	cellInfo.tb[,patient.uid:=sprintf("%s.%s",dataset,patient)]
	#### samples for frequency analysis
	cellInfo.tb[,usedForFreq:="Y"]
	cellInfo.tb[dataset=="CRC.LeiZhang2020.10X" & patient %in% c("P0408","P0613","P1025","P1026"),usedForFreq:="Y"]
	cellInfo.tb[dataset=="CRC.LeiZhang2020.10X" & patient %in% c("P0410","P0104"),usedForFreq:="SortingMyeloid"]
	cellInfo.tb[dataset=="HCC.QimingZhang2019.SS2" & patient=="D20171109", usedForFreq:="Dup.CrossPlatformBenchmark"]
	cellInfo.tb[grepl("^(CHOL.thisStudy|HCC.ChunhongZheng2017|NSCLC.XinyiGuo2018|CRC.LeiZhang2018)",dataset) &
				stype=="CD4",usedForFreq:="SortingTreg"]
	#cellInfo.tb[dataset=="SCC.KathrynEYost2019" & patient=="su010", patient:="" ]
	cellInfo.tb[usedForFreq!="Y",.N,by=c("dataset","patient","stype","usedForFreq","cancerType")][order(dataset,stype),]
## CRC.LeiZhang2020.10X, samples below contain all immune cells (CD45+): "CRC.P0408", "CRC.P0613", "CRC.P1025", "CRC.P1026"; other samples were sorted by other enrichment criteria

## D20171109 present in both HCC.QimingZhang2019.10X and HCC.QimingZhang2019.SS2. It was used to platform benchmark in the original study. Remove D20171109 in HCC.QimingZhang2019.SS2.
## cellInfo.tb[patient=="D20171109",table(dataset)]

## su010 present in both BCC.KathrynEYost2019 and SCC.KathrynEYost2019
## cellInfo.tb[patient=="su010",table(dataset)]

	{
		cellInfo.tb[grepl("AML.PeterVanGalen2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^BCC.KathrynEYost2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^SCC.KathrynEYost2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^NSCLC.RapolasZilionis2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment),]
		cellInfo.tb[grepl("^LIHC.LichunMa2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]
		cellInfo.tb[grepl("^MELA.HanjieLi2019",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]
		cellInfo.tb[grepl("^MELA.LivnatJerby-Arnon2018",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]
		cellInfo.tb[grepl("^MELA.MosheSade-Feldman2018",dataset),.N,
					by=c("patient","sampleID","treatment","cancerType","dataset","loc")
					][order(treatment,cancerType,patient),]

	}

	return(cellInfo.tb)
}


plotNightingaleRose <- function(dat.plot.NightingaleRose,empty_bar=2,
					    y.pretty=c(0,20000,40000,60000),
					    y.lim.min=-20000,
					    exp.tick=0.5,
					    my.title=NULL,
					    y.colum="numOfCell",
					    x.colum="cancerType",
					    group.colum="pub",
					    pallette.group=structure(scales::viridis_pal(option = "viridis")(5)[c(2,3)],
								     names=c("published","thisStudy")))
{
    dat.plot.NightingaleRose$y.value <- dat.plot.NightingaleRose[[y.colum]]
    dat.plot.NightingaleRose$x.value <- dat.plot.NightingaleRose[[x.colum]]
    #empty_bar <- 1
    to_add <- matrix(NA, empty_bar, ncol(dat.plot.NightingaleRose))
    colnames(to_add) <- colnames(dat.plot.NightingaleRose)
    dat.plot.NightingaleRose <- rbind(dat.plot.NightingaleRose, to_add)
    dat.plot.NightingaleRose[,id:=as.integer(x.value)]
    dat.plot.NightingaleRose[is.na(x.value),
			     id:=seq(max(dat.plot.NightingaleRose$id,na.rm=T)+1,
				     max(dat.plot.NightingaleRose$id,na.rm=T)+empty_bar)]
    dat.plot.NightingaleRose[,id.factor:=as.factor(id)]

    ##label_data <- dat.plot.NightingaleRose[,.(numOfCell=sum(.SD$numOfCell)),by=c("cancerType","id")]
    label_data <- dat.plot.NightingaleRose[,.(y.value=sum(.SD$y.value)),by=c(x.colum,"id")]
    label_data$x.value <- label_data[[x.colum]]
    number_of_bar <- nrow(label_data)
    #label_data[,id:=number_of_bar-as.integer(cancerType)+1]
    #label_data[,id:=as.integer(cancerType)]
    angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar
    #angle <-  90 - 360 * (label_data$id+0.5) /number_of_bar
    label_data$hjust<-ifelse( angle < -90, 1, 0)
    label_data$angle<-ifelse(angle < -90, angle+180, angle)
    if(is.null(y.pretty)){
	y.pretty <- label_data[,pretty(y.value)]
    }
    if(is.null(y.lim.min)){
	y.lim.min <- -y.pretty[2]
    }
    nticks <- length(y.pretty)

    p <- ggbarplot(dat.plot.NightingaleRose,x="id.factor",y="y.value",group=group.colum,
			   fill=group.colum,
			   color=NA,width=0.8,
			   title=my.title,
			   #ylab="number of cells",xlab="cancer types",
			   legend="right") + 
	#facet_wrap(~tech.cate,scales="free_x") +
	#coord_flip() +
	scale_fill_manual(values=pallette.group)+
	### axis
	geom_segment(data=data.table(x=rep(max(dat.plot.NightingaleRose$id)+0.5,1),
				     y.start=0,
				     y.end=y.pretty[nticks]),
		     aes(x = x, y = y.start, xend = x, yend = y.end),
		     colour = "grey", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
	### ticks
	geom_segment(data=data.table(x.start=rep(max(dat.plot.NightingaleRose$id)+0.5,nticks),
				     x.end=rep(max(dat.plot.NightingaleRose$id)+0.5,
					       nticks)-seq(0.3,0.2,-0.1/(nticks-1))*exp.tick,
				     y.start=y.pretty,
				     y.end=y.pretty),
		     aes(x = x.start, y = y.start, xend = x.end, yend = y.end),
		     colour = "grey", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
	### axis.y
	annotate("text", x = rep(max(dat.plot.NightingaleRose$id)+0.5,nticks),
		 y = y.pretty,
		 label = y.pretty,
		 color="grey", size=5 , angle=0, hjust=1.2,vjust=0.5) +
	ylim(y.lim.min,y.pretty[nticks]) +
	coord_polar(start=0,direction=1) +
	### label
	geom_text(data=label_data, aes(x=id, y=y.value+y.pretty[2]/10,
				       label=x.value,
				       angle=angle,
				       hjust=hjust),
		  color="black", alpha=0.6, size=4.2,
		  #angle= label_data$angle,
		  inherit.aes = FALSE ) +
	theme(
	      plot.title = element_text(hjust = 0.5),
	      axis.ticks=element_blank(),axis.line=element_blank(),
	      panel.grid = element_blank(),plot.margin = unit(rep(0,4), "cm"),
	      axis.text = element_blank(), axis.title = element_blank())
    return(p)
}

do.tissueDist <- function(cellInfo.tb,out.prefix,pdf.width=3,pdf.height=5,verbose=0)
{
    library("Startrac")
    dir.create(dirname(out.prefix),F,T)

    cellInfo.tb[,meta.cluster:=as.character(meta.cluster)]
    loc.avai.vec <- unique(cellInfo.tb[["loc"]])
    loc.avai.vec <- intersect(c("P","N","T"),loc.avai.vec)
    count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
    freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)

    {
	count.dist.melt.ext.tb <- test.dist.table(count.dist)
	p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
	OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
	OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
	rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
    }

    startrac.dist <- unclass(calTissueDist(cellInfo.tb,colname.cluster="meta.cluster"))
    startrac.dist <- startrac.dist[,loc.avai.vec]
    
    cuts <- c(0, 0.8, 1.2,Inf)
    startrac.dist.bin.values <- factor(c("-", "+/-", "+"),levels=c("-", "+/-", "+"))
    startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
							       ncol=ncol(startrac.dist))
    colnames(startrac.dist.bin) <- colnames(startrac.dist)
    rownames(startrac.dist.bin) <- rownames(startrac.dist)

#	sscClust:::plot.matrix.simple(freq.dist.bin,
#								  col.ht=rev(structure(colorRampPalette(brewer.pal(9,name="Blues"))(10),
#												   names=0:9 )),
#								  par.legend=list(labels=rev(sprintf("%s%%~%s%%",10*(0:9),c(10*(1:9),100) )),
#												  at=0:9),
#								  out.prefix=sprintf("%s.freq.dist",out.prefix),
#								  show.number=F,clust.row=T,exp.name=expression(italic(Freq)),
#								  #palatte=(brewer.pal(n = 7,name = "Blues")),
#								  #palatte=viridis::viridis(7),
#								  pdf.width = 4.5, pdf.height = pdf.height)

#	sscClust:::plot.matrix.simple(startrac.dist.bin,
#								  col.ht=rev(structure(viridis::viridis(3),
#													   names=levels(startrac.dist.bin.values))),
#								  out.prefix=sprintf("%s.startrac.dist.bin",out.prefix),
#								  show.number=F,clust.row=T,exp.name=expression(italic(R)[o/e]),
#								  pdf.width = pdf.width, pdf.height = pdf.height)

    sscVis::plotMatrix.simple(startrac.dist,
							      out.prefix=sprintf("%s.startrac.dist",out.prefix),
							      show.number=F,
							      clust.row=T,
							      #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
							      exp.name=expression(italic(R)[o/e]),
							      z.hi=2,
							      #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
							      palatte=viridis::viridis(7),
							      pdf.width = 4, pdf.height = pdf.height)

    sscVis::plotMatrix.simple(OR.dist.mtx,
							      out.prefix=sprintf("%s.OR.dist",out.prefix),
							      show.number=F,
							      #clust.row=T,
							      #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
							      waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
							      exp.name=expression(italic(OR)),
							      z.hi=4,
							      #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
							      palatte=viridis::viridis(7),
							      pdf.width = 4, pdf.height = pdf.height)
    if(verbose==1){
	return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
		    "p.dist.tb"=p.dist.tb,
		    "OR.dist.tb"=OR.dist.tb,
		    "OR.dist.mtx"=OR.dist.mtx))
    }else{
	return(OR.dist.mtx)
    }

}

test.dist.table <- function(count.dist,min.rowSum=0)
{
    count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
    sum.col <- colSums(count.dist)
    sum.row <- rowSums(count.dist)
    count.dist.tb <- as.data.frame(count.dist)
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid","cid","count")
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
						   this.row <- count.dist.melt.tb$rid[i]
						   this.col <- count.dist.melt.tb$cid[i]
						   this.c <- count.dist.melt.tb$count[i]
						   other.col.c <- sum.col[this.col]-this.c
						   this.m <- matrix(c(this.c,
								      sum.row[this.row]-this.c,
								      other.col.c,
								      sum(sum.col)-sum.row[this.row]-other.col.c),
								    ncol=2)
						   res.test <- fisher.test(this.m)
						   data.frame(rid=this.row,
							      cid=this.col,
							      p.value=res.test$p.value,
							      OR=res.test$estimate)
					       }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
								    by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
    #count.dist.melt.ext.tb[adj.p.value < 0.05,]
    #count.dist.melt.ext.tb[p.value < 0.05 & OR > 0,]
    
    return(count.dist.melt.ext.tb)

}

fetchMetaClusterID2CusterFullName <- function(col.use="cluster.name.full")
{
    require("data.table")
	ret.tb <- fread("../data/metaInfo/name.conversion.txt")
	ret.vec <- structure(ret.tb[[col.use]],names=ret.tb$meta.cluster)
	return(ret.vec)
}


calProliferationScore <- function(obj,assay.name,gene.prol,out.prefix=NULL,method="mean")
{
    f.gene <- rowData(obj)[,"display.name"] %in% gene.prol
    exp.sub <- as.matrix(assay(obj[f.gene,],assay.name))
    f.zero <- apply(exp.sub,1,function(x){ all(x==0) })
    if(sum(f.zero) > 0){
	    cat(sprintf("Number of gene with value zero in all cells: %d\n",sum(f.zero)))
    }
    exp.sub <- exp.sub[!f.zero,]

    if(method=="mean")
    {
        require("sscClust")
        score.prol <- colMeans(exp.sub)
        dat.score <- classify.outlier(score.prol,out.prefix=out.prefix)
        out.tb <- as.data.table(dat.score$score.cls.tb)
        colnames(out.tb)[1] <- "cellID"
    }else if(method=="AUCell")
    {
        require("AUCell")
        #####
        #f.gene <- rowData(obj)[,"display.name"] %in% gene.prol
        #exp.sub <- as.matrix(assay(obj,assay.name))
        exp.sub <- (assay(obj,assay.name))
        rownames(exp.sub) <- rowData(obj)[,"display.name"]
        f.zero <- apply(exp.sub,1,function(x){ all(x==0) })
        if(sum(f.zero) > 0){
            cat(sprintf("Number of gene with value zero in all cells: %d\n",sum(f.zero)))
        }
        exp.sub <- exp.sub[!f.zero,]

        #####
        pdf(sprintf("%s.buildRankings.1.pdf",out.prefix),width=7,height=7)
        cells_rankings <- AUCell_buildRankings(exp.sub)
        dev.off()

        geneSets <- list("prol"=intersect(rownames(exp.sub),gene.prol))
        #### geneSets <- GSEABase::GeneSet(genes, setName="geneSet1") # alternative
        cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)

        pdf(sprintf("%s.exploreThresholds.1.pdf",out.prefix),width=7,height=7)
        #par(mfrow=c(3,3))
        set.seed(123)
        cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
        dev.off()

        if("L_k2" %in% rownames(cells_assignment$prol$aucThr$thresholds)){
            th.prol <- cells_assignment$prol$aucThr$thresholds["L_k2","threshold"]
        }else{
            th.prol <- 0.09
        }

        ##geneSetName <- rownames(cells_AUC)[grep("prol", rownames(cells_AUC))]
        pdf(sprintf("%s.exploreThresholds.2.pdf",out.prefix),width=7,height=7)
        AUCell_plotHist(cells_AUC["prol",], aucThr=th.prol)
        abline(v=th.prol)
        dev.off()
        
        out.tb <- data.table(cellID=colnames(cells_AUC),
                     proliferationScore.bin=as.integer(getAUC(cells_AUC)["prol",]>th.prol),
                     proliferationScore=getAUC(cells_AUC)["prol",])
        out.tb$classification <- out.tb$proliferationScore.bin
        dat.score <- cells_AUC

        my.dat.score <- classify.outlier(getAUC(cells_AUC)["prol",],out.prefix=out.prefix)
        out.tb[,myCls:=my.dat.score$score.cls.tb[cellID,"classification"]]

    }
    return(list("out.tb"=out.tb,"detail"=dat.score))
}

#### visualization of clonotype
vis.clonotype <- function(out.prefix,
			  clone.x,sce.list,cinfo.clone.tb,
			  gene.sig.list,
			  sig.pretty=seq(-3,3,1),
			  gene.use=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4",
				     "IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
			  mcls.plot=c("CD8.c15.ISG.IFIT1",
				      "CD8.c12.Tex.CXCL13"),
			  mcls.sig=NULL,
              stype="CD8",
			  pdf.width=16.5,pdf.height=6,
			  prob.mat=NULL,sortByProb=F)
{
    require("circlize")
    #####x <- "UCEC-P20190312_C001978:59"
    ##clone.perCell.tb <- cinfo.clone.tb[clone.id==clone.x,]
    if(is.null(mcls.plot)){
        clone.perCell.tb <- cinfo.clone.tb[clone.id==clone.x,]
        mcls.plot <- as.character(sort(unique(clone.perCell.tb$meta.cluster)))
    }else{
        clone.perCell.tb <- cinfo.clone.tb[clone.id==clone.x & meta.cluster %in% mcls.plot,]
    }
    ##data.id <- clone.perCell.tb$dataset.old[1]
    data.id <- clone.perCell.tb$dataset[1]
    cat(sprintf("Processing clonotype %s (%s)\n",clone.x,data.id))

    if(is.null(sce.list[[data.id]])){
        if(data.id=="NSCLC.XinyiGuo2018") { t.id <- "LC.XinyiGuo2018" 
        }else if(data.id=="STAD.BoxiKang2020"){ t.id <- "STAD.BoxiKang2019"
        }else{ t.id <- data.id }
        sce.list[[data.id]] <- readRDS(sprintf("../data/expression/%s/byDataset/%s.sce.rds",stype,t.id))
    }

    #### mean   
    ##f.gene.zz <- match(unique(gene.sig.tb$geneID),rowData(sce.list[[data.id]])$display.name)
    f.gene.zz <- match(unique(c(unname(do.call(c,gene.sig.list)),gene.use)),
                       rowData(sce.list[[data.id]])$display.name)
    f.gene.zz <- f.gene.zz[!is.na(f.gene.zz)]    

    sce.z <- ssc.assay.zscore(obj=sce.list[[data.id]][f.gene.zz,],
						      assay.name="norm_exprs")
    sce.plot.avg <- NULL
    if(length(gene.sig.list) > 0){
        sce.plot.avg <- ssc.average.gene(sce.z,assay.name="norm_exprs.z",gene=gene.sig.list)
        sce.plot.avg <- sce.plot.avg[,clone.perCell.tb$Cell_Name]
        print("all(clone.perCell.tb$Cell_Name==colnames(sce.plot.avg))")
        print(all(clone.perCell.tb$Cell_Name==colnames(sce.plot.avg)))
        colData(sce.plot.avg) <- DataFrame(clone.perCell.tb[,c("Cell_Name","meta.cluster","clone.id"),
                           with=F])
        sce.plot.avg$meta.cluster <- as.character(sce.plot.avg$meta.cluster)
        colnames(sce.plot.avg) <- sce.plot.avg$Cell_Name
        sce.plot.avg$cluster.name <- fetchMetaClusterID2CusterFullName()[sce.plot.avg$meta.cluster]
    }

    #sig.pretty <- seq(-1,1,0.5)
    #####
    if(is.null(mcls.sig)){
        mcls.sig <- mcls.plot
    }
    sig.plot.vec <- fetchMetaClusterID2CusterFullName()[mcls.sig]
    #print(str(sig.plot.vec))
    #print(sce.plot.avg)
    
    sce.plot <- sce.z
    sce.plot <- sce.plot[,clone.perCell.tb$Cell_Name]
    print("all(clone.perCell.tb$Cell_Name==colnames(sce.plot))")
    print(all(clone.perCell.tb$Cell_Name==colnames(sce.plot)))
    colData(sce.plot) <- cbind(colData(sce.plot),clone.perCell.tb[,c("meta.cluster","clone.id"),with=F])
    sce.plot$meta.cluster <- as.character(sce.plot$meta.cluster)
    sce.plot$cluster.name <- fetchMetaClusterID2CusterFullName()[sce.plot$meta.cluster]

    if(!is.null(sce.plot.avg)){
        print("all(colnames(sce.plot)==colnames(sce.plot.avg))")
        print(all(colnames(sce.plot)==colnames(sce.plot.avg)))
        for(a.sig in names(sig.plot.vec)){
            #print(sig.plot.vec[a.sig])
            a.dat <- assay(sce.plot.avg,"exprs")[sig.plot.vec[a.sig],]
            #a.dat <- assay(sce.plot.avg,"exprs")[a.sig,]
            a.dat[a.dat < sig.pretty[1] ] <- sig.pretty[1]
            a.dat[a.dat > sig.pretty[length(sig.pretty)] ] <- sig.pretty[length(sig.pretty)]
            colData(sce.plot)[[sprintf("sig.%s",a.sig)]] <- a.dat
            
        }
    }

    for(a.sig in names(sig.plot.vec)){
        ### add cluster assignment probability
        if(!is.null(prob.mat)){
            colData(sce.plot)[[sprintf("prob.%s",a.sig)]] <- prob.mat[colnames(sce.plot),a.sig]
        }
    }
    
    {
        dat.for.sort <- data.table(cellID=colnames(sce.plot), meta.cluster=sce.plot$meta.cluster)

        if(!is.null(sce.plot.avg)){
            dat.for.sort <- cbind(dat.for.sort,as.data.table(t(assay(sce.plot.avg[sig.plot.vec,]))))
        }
        #print("ssssssss")
        if(!is.null(prob.mat)){
            .tmp.colData <- colData(sce.plot)
            #sce.plot.debug <<- sce.plot
            #dat.for.sort.debug <<- dat.for.sort
            dat.for.sort <- as.data.table(cbind(dat.for.sort,
                                                as.data.table(.tmp.colData[,grepl("prob.",colnames(.tmp.colData))])))
        }
        #print("ttt")
        #print(head(dat.for.sort))

        if(!is.null(sce.plot.avg)){
            dat.tmp <- dat.for.sort[!meta.cluster %in% names(sig.plot.vec),]
            for(mcls.xx in mcls.sig)
            {
                dat.tmp <- rbind(dat.tmp,setorderv(dat.for.sort[meta.cluster==mcls.xx,],mcls.xx))
            }
            colnames(dat.for.sort)[3:(2+length(sig.plot.vec))] <- names(sig.plot.vec)
        }

        if(!is.null(prob.mat) & sortByProb){
            dat.tmp <- dat.for.sort[!meta.cluster %in% names(sig.plot.vec),]
            for(mcls.xx in mcls.sig)
            {
                dat.tmp <- rbind(dat.tmp,setorderv(dat.for.sort[meta.cluster==mcls.xx,],sprintf("prob.%s",mcls.xx)))
            }
        }
        dat.for.sort <- dat.tmp
    }

    #gene.use <- intersect(c("HAVCR2","CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "CTLA4","PDCD1","IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
#			  rowData(sce.plot)$display.name)

    if(is.null(gene.use)){
        gene.use <- unique(unname(do.call(c,gene.sig.list)))
    }
    print(str(gene.use))

    gene.use <- intersect(gene.use,
			  rowData(sce.plot)$display.name)
    ###print(str(gene.use))
    f.gene.use <- match(gene.use,rowData(sce.plot)$display.name)
    sce.plot.comb <- sce.plot[f.gene.use,]

    annotation_legend_param <- NULL
    if(!is.null(sce.plot.avg)){
        sce.plot.comb <- sce.plot[f.gene.use,
                      colnames(ssc.assay.hclust(sce.plot.avg[(sig.plot.vec),],
                                "exprs",order.col=T,
                                clustering.distance="cosine",
                                clustering.method="ward.D2"))]
    
        annotation_legend_param <- llply(mcls.sig,function(x){ list(at=sig.pretty) })
        names(annotation_legend_param) <- sprintf("sig.%s",mcls.sig)
    }

    annotation_legend_param_prob <- llply(mcls.sig,function(x){
                                              list(at=seq(0,1,0.2),color_bar = "continuous",
                                                   legend_direction = "horizontal",
                                                   legend_width = unit(4, "cm"), legend_height = unit(2,"cm")) })
    names(annotation_legend_param_prob) <- sprintf("prob.%s",mcls.sig)
    colSet.plot <- g.colSet
    if(!is.null(prob.mat)){
        annotation_legend_param <- c(annotation_legend_param,annotation_legend_param_prob)
        colSet.prob <- llply(mcls.sig,function(x){
                                 colorRamp2(seq(0,1,length=9),sscVis:::getColorPaletteFromNameContinuous("Greys"),space="LAB") })
        names(colSet.prob) <- sprintf("prob.%s",mcls.sig)
        colSet.plot <- c(colSet.plot,colSet.prob)
    }
    #print(sprintf("%s.example.%s.comb",out.prefix,gsub(":","_",clone.x,perl=T))) 
    #print((sce.plot.comb[,dat.for.sort$cellID]))
    #sce.plot.comb.debug <<- sce.plot.comb
    col.use <- if(is.null(prob.mat)) c("cluster.name",sprintf("sig.%s",names(sig.plot.vec)))
                else c("cluster.name",sprintf("sig.%s",names(sig.plot.vec)),sprintf("prob.%s",names(sig.plot.vec)))
    col.use <- intersect(col.use,colnames(colData(sce.plot.comb)))
    ssc.plot.heatmap(sce.plot.comb[,dat.for.sort$cellID],assay.name="norm_exprs.z",
		     out.prefix=sprintf("%s.example.%s.comb",out.prefix,gsub(":","_",clone.x,perl=T)),
		     columns=col.use,
		     columns.order=c("cluster.name"),
		     clustering.distance="cosine",
		     clustering.method="ward.D2",
		     colSet=colSet.plot,do.scale=F,
		     do.clustering.row=F,
		     do.clustering.col=F,
		     par.heatmap=list(cex.row=1.2,cex.column=0),
		     #palette.name="RdBu",
		     palette.name="cividis",
		     palette.ann.numeric="RdBu",
		     Y.level.ann.numeric=sig.pretty,
		     annotation_legend_param=annotation_legend_param,
		     z.lo = -3.0, z.hi = 3.0, z.step = 1,
		     ###pdf.width = 16.5, pdf.height = 6,
		     pdf.width = pdf.width, pdf.height = pdf.height,
		     mytitle = clone.x)

}

ana.clonotypeAcrossMcls.moreThanTwo <- function(object,
			    in.dat,
			    out.prefix,
			    aid="PathTrm",
			    lim.LLR=10,
                stype="CD8",
			    par.barplot=list("pdf.width.byPatientF"=7,"pdf.height.byPatientF"=7.5,
					     "pdf.width.byPatientT"=7,"pdf.height.byPatientT"=8),
			    par.ggbarplot=NULL,
                show.N=F,
                l.colSet=g.colSet,
			    mcls.moi=c("CD8.c02.Tm.IL7R","CD8.c10.Trm.ZNF683","CD8.c12.Tex.CXCL13"))
{

    dir.create(dirname(out.prefix),F,T)

    dat.block <- object@clonotype.dist.cluster[,mcls.moi]
    ### clone.info.flt.tb
    dat.block.clone.index <- Startrac:::mrow.entropy(dat.block)
    dat.block.clone.index[is.na(dat.block.clone.index)] <- 0
    clone.info.index.tb <- data.table(cloneID=names(dat.block.clone.index),
				      pIndex.tran=dat.block.clone.index,
                      N.mcls=rowSums(dat.block))
    f.crossAll <- rowSums(dat.block > 0)==ncol(dat.block)
    clone.info.index.tb <- clone.info.index.tb[pIndex.tran > 0 & f.crossAll,]
    clone.info.patient.tb <- in.dat[,.N,by=c("cancerType","dataset","dataset.old","patient","cloneID")]
    clone.info.flt.tb <- merge(clone.info.index.tb,clone.info.patient.tb,by="cloneID")

    clone.LLR.tb <- ldply(seq_len(nrow(clone.info.flt.tb)),function(i){
			calCloneLLR(in.dat[majorCluster %in% mcls.moi,],
			    prob.mat=prob.svm.mat.list[[ clone.info.flt.tb$dataset.old[i] ]],
			    cloneID.x=clone.info.flt.tb$cloneID[i],verbose=F)
			})
    clone.info.flt.tb <- merge(clone.info.flt.tb,clone.LLR.tb,by="cloneID")
    clone.info.flt.tb <- clone.info.flt.tb[order(-LLR,pIndex.tran),]
    clone.info.flt.tb <- clone.info.flt.tb[LLR > 1 & G.best==G.obs,]

    if(nrow(clone.info.flt.tb) < 0){
        ret.list <- list("dat.block.flt"=NULL,
                 "dat.block.flt.freq"=NULL,
                 "clone.info.flt.tb"=NULL)
        warning(sprintf("no clonotype found !"))
        return(ret.list)
    }

    ## dat.block
    dat.block.flt <- dat.block[clone.info.flt.tb$cloneID,,drop=F]
    dat.block.flt.freq <- sweep(dat.block.flt,1,rowSums(dat.block.flt),"/")


    do.plot.bar <- function(dat.block.flt.freq,clone.info.flt.tb,
			    out.prefix,
			    pdf.width=4,
			    pdf.height=9,
			    splitByPatient=F)
    {

        dat.plot.tb <- cbind(data.table(cloneID=rownames(dat.block.flt.freq)),
                     dat.block.flt.freq)
        dat.plot.tb <- merge(clone.info.flt.tb,dat.plot.tb)
        dat.plot.tb <- dat.plot.tb[order(LLR,-pIndex.tran),]
        dat.plot.melt <- melt(dat.plot.tb[,c("cloneID",mcls.moi,"patient","cancerType","N.mcls"),with=F],
                      variable.name = "meta.cluster", value.name = "freq",
                      id.vars=c("cloneID","patient","cancerType","N.mcls"))
        dat.plot.melt[,cloneID.fac:=factor(cloneID,levels=dat.plot.tb$cloneID)]
        #dat.plot.melt[,meta.cluster:=factor(meta.cluster,levels=sort(unique(as.character(meta.cluster)),
        #							     decreasing=T))]
        dat.plot.melt[,meta.cluster:=factor(meta.cluster,levels=rev(mcls.moi))]
        dat.plot.tb[,cloneID.fac:=factor(cloneID,levels=cloneID)]
        if(!is.null(lim.LLR)){
            dat.plot.tb[LLR>lim.LLR,LLR:=lim.LLR]
        }

        ### freq
        p1 <- do.call(ggbarplot,c(list(data=dat.plot.melt,x="cloneID.fac",y="freq",
                           color=NA, fill="meta.cluster"),
                      par.ggbarplot))
        p1 <- p1 + labs(x="",y="Frequency")
        if(show.N){
            dat.plot.N <- clone.info.flt.tb[cloneID %in% dat.plot.melt$cloneID,]
            dat.plot.N[,cloneID.fac:=factor(cloneID,levels=levels(dat.plot.melt$cloneID.fac))]
            p1 <- p1 + geom_text(data=dat.plot.N,
                                 aes(x=cloneID.fac,label=N.mcls),y=1,
                                 color="white",size=3.34,hjust=1.1,
                                 inherit.aes=F)
        }
        p1 <- p1 + scale_fill_manual(values=l.colSet$meta.cluster) +
            theme(axis.text.y=element_blank(),legend.position="right",
                  axis.text.x=element_text(size=12),
                  strip.placement = "outside",
                  strip.background=element_blank(),
                  strip.text.y.left=element_text(angle=0)) +
            coord_flip()
        if(splitByPatient){
            p1 <- p1 + facet_grid(patient~.,scales = "free", space = "free",switch="y")
        }

        ### pIndex.tran
        p2 <- do.call(ggbarplot,c(list(data=dat.plot.tb,x="cloneID.fac",y="pIndex.tran",
                           color=NA, fill="steelblue"),
                      par.ggbarplot))
        p2 <- p2 + 
            labs(x="",y="Entropy")+
            geom_hline(yintercept=0.1,linetype="dashed") +
            theme(
                  axis.text.y=element_blank(),
                  axis.text.x=element_text(size=12),
                  legend.position="right",
                  #strip.placement = "outside",
                  strip.background=element_blank(),
                  strip.text.y=element_blank()
                  #strip.text.y.left=element_text(angle=0)
                  ) +
            coord_flip()
        if(splitByPatient){
            #p2 <- p2 + theme(axis.text.y=element_blank())
            p2 <- p2 + facet_grid(patient~.,scales = "free", space = "free",switch="y")
        }

        ### LLR
        p3 <- do.call(ggbarplot,c(list(data=dat.plot.tb,x="cloneID.fac",y="LLR",
                           color=NA, fill="#82B446"),
                      par.ggbarplot))
        p3 <- p3 +
            labs(x="")+
            geom_hline(yintercept=1,linetype="dashed") +
            theme(
                  axis.text.y=element_blank(),
                  axis.text.x=element_text(size=12),
                  legend.position="right",
                  #strip.placement = "outside",
                  strip.background=element_blank(),
                  strip.text.y=element_blank()
                  #strip.text.y.left=element_text(angle=0)
                  ) +
            coord_flip()
        if(splitByPatient){
            ##p3 <- p3 + theme(axis.text.y=element_blank())
            p3 <- p3 + facet_grid(patient~.,scales = "free", space = "free",switch="y")
        }

        p00 <- cowplot::plot_grid(p1+theme(legend.position="none"),
                    p2,
                    p3,
                    align="h",nrow=1,
                    rel_widths=if(splitByPatient) c(1,0.4,0.4) else c(1,0.4,0.4)
                    )
        pout <- cowplot::plot_grid(p00,get_legend(p1),ncol=1,rel_heights=c(0.7,0.3))
        #ggsave(sprintf("%s.test.pdf",out.prefix),plot=pout,width=pdf.width,height=pdf.height)

        p00 <- cowplot::plot_grid(p1+theme(legend.position="none"),
                    p3,
                    align="h",nrow=1,
                    rel_widths=if(splitByPatient) c(4.5,1.78) else c(4.5,1.78)
                    )
        pout <- cowplot::plot_grid(p00,get_legend(p1),ncol=1,rel_heights=c(0.7,0.3))
        ###ggsave(sprintf("%s.pdf",out.prefix),plot=pout,width=6.28,height=pdf.height)
        ggsave(sprintf("%s.pdf",out.prefix),plot=pout,width=pdf.width,height=pdf.height)
        cat(sprintf("----- width: %4.2f, height: %4.2f\n",pdf.width,pdf.height))
    }

#    do.plot.bar(dat.block.flt.freq,
#		clone.info.flt.tb,
#		out.prefix=sprintf("%s.bar.freq.byPatientF",out.prefix),
#		pdf.width=par.barplot$pdf.width.byPatientF,pdf.height=par.barplot$pdf.height.byPatientF,splitByPatient=F)
    do.plot.bar(dat.block.flt.freq,
		clone.info.flt.tb,
		out.prefix=sprintf("%s.bar.freq.byPatientT",out.prefix),
		pdf.width=par.barplot$pdf.width.byPatientT,pdf.height=par.barplot$pdf.height.byPatientT,splitByPatient=T)

    #### specific clonotype
    if(F)
    {
        clone.example.vec <- head(clone.info.flt.tb$cloneID,n=10)

        gene.sig.slim.list <- list("CD8.c02.Tm.IL7R"=c("IL7R","ZFP36L2","CXCR4","ZFP36","ANXA1","GPR183"),
                       "CD8.c10.Trm.ZNF683"=c("ZNF683","HOPX","CAPG"),
                       "CD8.c05.Tem.CXCR5"=c("GZMK", "CD74", "CST7"),
                       "CD8.c06.Tem.GZMK"=c("GZMK", "CD74", "CST7"),
                       "CD8.c11.Tex.PDCD1"=c("GZMK", "CD74", "CST7"),
                       "CD8.c15.ISG.IFIT1"=c("IFIT1","RSAD2","IFIT3","IFI44L","MX1"),
                       "CD8.c12.Tex.CXCL13"=c("CXCL13", "TNFRSF9", "LAYN", "ENTPD1", "HAVCR2", "CTLA4"),
                       "CD8.c14.Tex.TCF7"=c("CD200","GNG4","CXCL13","TNFRSF4","TCF7"),
                       "CD8.c17.Tm.NME1"=c("BIRC5","RRM2","CDC20","NME1","CENPW"),
                       "CD4.c12.Tem.GZMK"=c("GZMK", "CCL4", "GZMA", "CST7", "NKG7", "EOMES", "IFNG", "TNFRSF9"),
                       "CD4.c16.Tfh.CXCR5"=c("CXCL13", "PDCD1", "GNG4", "NMB", "NR3C1", "CD200", "IGFL2","IL21",
                                             "CXCR5","BCL6"),
                       "CD4.c17.TfhTh1.CXCL13"=c("CXCL13", "PDCD1",
                                                 "GNG4", "NR3C1", "CD200", "IL21",
                                                 "LAG3","HAVCR2","CCL4","KRT86","IFNG","GZMB","BHLHE40","CTLA4","ENTPD1"),
                       "CD4.c20.Treg.TNFRSF9"=c("FOXP3","TNFRSF18","TIGIT","IL2RA","LAYN","TNFRSF9","CTLA4","BATF","CCR8")
                       )
        gene.sig.slim.list <- gene.sig.slim.list[mcls.moi]
        names(gene.sig.slim.list) <- fetchMetaClusterID2CusterFullName()[names(gene.sig.slim.list)]
        print(str(gene.sig.slim.list))

        l_ply(seq_along(clone.example.vec),function(i){
            dataset.old.id <- clone.info.flt.tb[cloneID==clone.example.vec[i],][["dataset.old"]][1]
            vis.clonotype(out.prefix=sprintf("%s.%s", out.prefix,"exampleClone"),
                  clone.x=clone.example.vec[i],
                  sce.list=sce.list,
                  cinfo.clone.tb=in.dat[stype==stype,],
                  gene.sig.list=gene.sig.slim.list, ### manually picked
                  ##gene.sig.list=gene.sig.list, ### top 50
                  sig.pretty=seq(-1,1,0.5),
                  gene.use=unique(do.call(c,unname(gene.sig.slim.list))),
                  prob.mat=prob.svm.mat.list[[dataset.old.id]],
                  mcls.plot=mcls.moi)
                     })
    }
    
    ret.list <- list("dat.block.flt"=dat.block.flt,
		"dat.block.flt.freq"=dat.block.flt.freq,
		"clone.info.flt.tb"=clone.info.flt.tb)
    saveRDS(ret.list,sprintf("%s.ret.list.rds",out.prefix))
    return(ret.list)
	
}


makeFig.ExampleGeneBarplot <- function(out.prefix,gene.to.plot,
				       gene.long.tb,
				       gene.long.collapsed.tb,
				       gene.desc.top,
				       mcls.plot,mod.sort=3,
				       th.dprime=0.15,colSet.cancerType=NULL)
{
    prepare.data.for.plot <- function(dat.long,gene.desc.top,mcls,a.gene,mod.sort=3)
    {
        mapping.dataset.tb <- unique(data.table(dataset=gene.long.tb$dataset,dataset.old=gene.long.tb$dataset.old))
        mapping.dataset.tb[dataset.old=="CHOL.YaoHe10X",dataset:="CHOL.QimingZhang2019.10X"]
        mapping.dataset.tb[dataset.old=="CHOL.LichunMa2019",dataset:="CHOL.LichunMa2019"]
        mapping.dataset.vec <- structure(mapping.dataset.tb$dataset,names=mapping.dataset.tb$dataset.old)

        dat.plot <- dat.long[meta.cluster==mcls & geneID==a.gene,]
        dat.meta <- gene.desc.top[meta.cluster==mcls & geneSymbol==a.gene,]

        if(mod.sort==1 || mod.sort==2){
            dat.plot[,dataset:=gsub("\\.CD[48].+$","",aid)]
            dat.plot[,dataset:=mapping.dataset.vec[dataset]]
            dat.plot <- dat.plot[,c("aid","geneID","P.Value","adj.P.Val","sig",
                                    "dprime","vardprime","dataset",
                                    "meta.cluster","cancerType"),with=F]
            if(mod.sort==1){
                tmp.order <- dat.plot[order(-meta.cluster,dprime),]
            }else if(mod.sort==2){
                tmp.order <- dat.plot[order(-meta.cluster,-dprime),]
                dat.plot[,cancerType:=factor(cancerType,levels=rev(unique(tmp.order$cancerType)))]
                tmp.order <- dat.plot[order(-meta.cluster,cancerType,dprime),]
            }
            dat.plot[,dataset:=factor(dataset,levels=c(unique(tmp.order$dataset),"combined"))]
            dat.plot[,dataset.cate:="perStudies"]
            dat.meta.plot <- data.table(aid="combined",geneID=a.gene,
                                    P.Value=dat.meta$comb.p,
                                    adj.P.Val=dat.meta$comb.padj,
                                    sig=dat.meta$sig,
                                    dprime=dat.meta$comb.ES,
                                    vardprime=dat.meta$comb.ES.sd^2,
                                    dataset="combined",
                                    meta.cluster=mcls,
                                    cancerType="panC",
                                    dataset.cate="combined")
        }else if(mod.sort==3){
            dat.plot <- dat.plot[,c("geneID","cancerType","dprime","vardprime",
                                    "P.Value","adj.P.Val","sig",
                                    "meta.cluster"),with=F]
            tmp.order <- dat.plot[order(-meta.cluster,sig,dprime),]
            dat.plot[,cancerType:=factor(cancerType,levels=c(unique(tmp.order$cancerType),"panC"))]
            dat.meta.plot <- data.table(geneID=a.gene,
                                        cancerType="panC",
                                        dprime=dat.meta$comb.ES,
                                        vardprime=dat.meta$comb.ES.sd^2,
                                        P.Value=dat.meta$comb.p,
                                        adj.P.Val=dat.meta$comb.padj,
                                        sig=dat.meta$sig,
                                        meta.cluster=mcls)
        }
        dat.plot <- rbind(dat.plot,dat.meta.plot)
        dat.plot[,lower:=dprime-sqrt(vardprime)]
        dat.plot[,upper:=dprime+sqrt(vardprime)]
        dat.plot[,sig:=as.logical(sig)]
        #### for display purpose
        dat.plot[sig==F,adj.P.Val:=1]
        dat.plot[sig==T & adj.P.Val > 0.01,adj.P.Val:=0.01]
        return(dat.plot)
    }

    ### gene by gene
    dat.fig.list <- llply(gene.to.plot,function(a.gene){
	    if(mod.sort==3){
		    dat.plot <- prepare.data.for.plot(gene.long.collapsed.tb,gene.desc.top,mcls.plot,a.gene,mod.sort=mod.sort)
	    }else{
		    dat.plot <- prepare.data.for.plot(gene.long.tb,gene.desc.top,mcls.plot,a.gene,mod.sort=mod.sort)
	    }
	    ###dat.plot.debug <<- dat.plot
	    vline.x <- NULL
	    if(mod.sort==3){
		    vline.x <- dat.plot[sig==T,min(as.integer(cancerType))]
	    }

	    
	    ###print(vline.x)
	    p <- ggplot(dat.plot, aes_string(if(mod.sort==3) "cancerType" else "dataset", "dprime")) + 
		    geom_col(aes(fill=adj.P.Val),color=NA,position = "dodge2",width=0.8) +
		    scale_fill_distiller(palette = "Blues",breaks=c(0,0.002,0.004,0.006,0.008,0.01),
							     limits=c(0,0.01),na.value="lightgray",
							     guide = guide_colorbar(label.theme = element_text(angle = 45, hjust = 1),
													    direction = "horizontal")) +
		    ####geom_col(aes(fill=sig),color=NA,position = "dodge2") +
		    ####scale_fill_manual(values=c("TRUE"="#F8766D","FALSE"="#00BFC4")) +
		    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
		    geom_hline(yintercept=th.dprime,linetype="dashed",alpha=0.8,color="lightgray") +
		    theme_classic() +
		    #####labs(x="dataset",y="Effect Size",title=a.gene) +
		    ##labs(x=if(mod.sort==3) "cancerType" else "dataset",y="Effect Size",title=a.gene) +
		    labs(x="",y="Effect Size",title=a.gene) +
		    theme(plot.title = element_text(hjust = 0.5))
		    ##theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5),
		    #theme(axis.text.x=element_custom(fill=colSet$cancerType[as.character(levels(dat.plot$cancerType))]),
	    if(mod.sort==3){
		    p <- p +theme(axis.text.x=element_text(angle = 90,size=10, hjust = 1,vjust=0.5),
					      axis.text.y=element_text(size=10))
	    }else{
		    mapping.d2c.tb <- unique(dat.plot[,c("dataset","cancerType"),with=F])
		    mapping.d2c.tb[,dataset:=as.character(dataset)]
		    mapping.d2c.tb[,cancerType:=as.character(cancerType)]
		    mapping.d2c <- structure(mapping.d2c.tb$cancerType,names=mapping.d2c.tb$dataset)
		    par.fill <- colSet.cancerType[mapping.d2c[as.character(intersect(levels(dat.plot$dataset),unique(dat.plot$dataset)))]]
		    #print(par.fill)
		    #print(str(par.fill))
		    p <- p + theme(axis.text.x=element_custom(size=10, fill=par.fill))
	    }
	    if(!is.null(vline.x) && is.finite(vline.x)){
		    #p <- p + geom_vline(xintercept=vline.x-0.5,linetype="dashed",color="lightgray",alpha=0.8)
		    p <- p + geom_vline(xintercept=vline.x-0.5,linetype="dashed",color="red",alpha=0.8)
	    }
	    #ggsave(sprintf("%s.cmp.gene.example.%s.pdf",out.prefix,a.gene),
	    #	   width=if(mod.sort==3) 6.5 else 7,height=if(mod.sort==3) 2.5 else 4.0)
	    return(p)
			      },.parallel=T)
    names(dat.fig.list) <- gene.to.plot
    saveRDS(dat.fig.list,file=sprintf("%s.cmp.gene.example.fig.rds",out.prefix))

}

run.clusterProfiler <- function(out.prefix, gene.oi, gene.bg,
								my.title="",db.name="GO",es=NULL,semData=NULL,
                                verbose=T,
								pdf.width=c(9,7,7,18),gset.hl=NULL,...)
{
	require("clusterProfiler")
	require("org.Hs.eg.db")
	require("GOSemSim")

	ret.enrich <- NULL
	if(db.name=="GO")
	{
	    ret.enrich.raw <- clusterProfiler::enrichGO(gene.oi,universe=gene.bg,...)
	    ret.enrich <- clusterProfiler::simplify(clusterProfiler::gofilter(ret.enrich.raw,c(4,5)),
						    cutoff=0.6, by="pvalue", select_fun=min, measure="Jiang",
						    semData=semData)
	}else if(db.name=="wikipathways")
	{
	    wpgmtfile <- "../data/external/wikipathways-20200210-gmt-Homo_sapiens.gmt"
	    wp2gene <- read.gmt(wpgmtfile)
	    wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
	    wpid.id = bitr(wp2gene$gene, toType="SYMBOL", fromType="ENTREZID", OrgDb="org.Hs.eg.db")
	    wpid.id.vec <- structure(wpid.id$SYMBOL,names=wpid.id$ENTREZID)
	    wp2gene$symbol <- wpid.id.vec[wp2gene$gene]
	    wpid2gene <- wp2gene %>% dplyr::filter(!is.na(symbol)) %>% dplyr::select(wpid, symbol) #TERM2GENE
	    wpid2name <- wp2gene %>% dplyr::filter(!is.na(symbol)) %>% dplyr::select(wpid, name) #TERM2NAME

	    ret.enrich <- enricher(gene.oi, universe=gene.bg,
						       qvalueCutoff=0.1,
						       TERM2GENE=wpid2gene,TERM2NAME=wpid2name)
	}else if(db.name=="KEGG")
	{
	    eg.oi = bitr(gene.oi, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    eg.bg = bitr(gene.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    ret.enrich <- clusterProfiler::enrichKEGG(eg.oi$ENTREZID, organism = 'hsa',
							     universe=eg.bg$ENTREZID,
							     qvalueCutoff=0.1)
	    ret.enrich <- setReadable(ret.enrich,OrgDb=org.Hs.eg.db,keyType="ENTREZID")

	}else if(db.name=="MKEGG")
	{
	    eg.oi = bitr(gene.oi, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    eg.bg = bitr(gene.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	    ret.enrich <- enrichMKEGG(eg.oi$ENTREZID, organism = 'hsa',
							      universe=eg.bg$ENTREZID,
							      qvalueCutoff=0.1)
	    ret.enrich <- setReadable(ret.enrich,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
	}
	##write.table(ret.enrich@result,file=sprintf("%s.enrich.txt",out.prefix),row.names=F,sep="\t",quote=F)
    if(verbose==T){
        write.table(cbind(data.table(aid=my.title),as.data.frame(ret.enrich)),
                    file=sprintf("%s.enrich.txt",out.prefix),row.names=F,sep="\t",quote=F)
    }

	if(nrow(as.data.frame(ret.enrich))==0){
	    return(list("plot"=NULL,
			"enrichResult"=ret.enrich))
	}

	p1 <- dotplot(ret.enrich, showCategory=30) +
		ggtitle(sprintf("%s",my.title)) +
		scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
	ggsave(sprintf("%s.dotplot.pdf",out.prefix),width=pdf.width[1],height=7,useDingbats=F)
    #ret.enrich.debug <<- ret.enrich

    if(verbose==T){
	    p2 <- emapplot(ret.enrich, 
                   #### version difference
                   #pie_scale=1.5,
                   node_scale=1.5,
                   layout="kk") +
		scale_colour_distiller(palette = "Reds",breaks=c(0,0.01,0.05,0.1),limits=c(0,0.1),na.value="lightgray")
	    ggsave(sprintf("%s.emaplot.pdf",out.prefix),width=pdf.width[2],height=7,useDingbats=F)
    }else{
        p2 <- NULL
    }


	if(!is.null(es))
	{
	    if(!is.null(gset.hl)){
		    ret.enrich.hl <- clusterProfiler.dplyr::filter(ret.enrich,Description %in% gset.hl)
		    rownames(ret.enrich.hl@result) <- ret.enrich.hl@result$ID
	    }else{
		    ret.enrich.hl <- ret.enrich
	    }
        if(verbose==T){
            p3 <- cnetplot(ret.enrich.hl, foldChange=es) +
                scale_colour_gradient(name = "ES", low="blue",high="red",
                                      na.value = "#E5C494",
                                      limits=c(0,2),breaks=c(0.15,0.5,1,2),
                                      guide = guide_colorbar(title = "ES"))
                ggsave(sprintf("%s.cnet.pdf",out.prefix),width=pdf.width[3],height=7,useDingbats=F)
        }else{
            p3 <- NULL
        }

	    options(stringsAsFactors=T)
        if(verbose==T){
            p4 <- clusterProfiler::heatplot(ret.enrich, foldChange=es) +
                    theme(axis.text.x=element_text(size=3)) +
                    geom_tile(aes_(fill = ~foldChange), color = NA) +
                    scale_fill_continuous(low="blue", high="red", name = "ES",
                                          breaks=c(0.15,0.5,1,2),limits=c(0,2),
                                          na.value="lightgray",
                                          guide = guide_colorbar(title = "ES"))
                ggsave(sprintf("%s.heatmap.pdf",out.prefix),width=pdf.width[4],height=6)
        }else{
            p4 <- NULL
        }
	    options(stringsAsFactors=F)
	}
	return(list("plot"=list("p1"=p1,"p2"=p2,"p3"=p3,"p4"=p4),
				"enrichResult"=ret.enrich))
}

do.plot.freq.heatmap <- function(dat.plot.a,colSet,mapping.vec,group.var="meta.cluster",
				 dat.plot.b=NULL,out.prefix,k.column=6,k.row=1,
				 fltBySD=F,TH.cor=0.35,ann.bar.r=1,k.corMat=NULL,
				 must.include=NULL,
				 pdf.width.cor=8,pdf.height.cor=7.7,
				 pdf.width=8,pdf.height=8,...)
{
    library("ComplexHeatmap")
    library("dendextend")
    library("sscClust")
    library("RColorBrewer")

    ht.a.tb <- dcast(dat.plot.a,group.var~donor.var,value.var="freq",fill=0)
    ht.a <- as.matrix(ht.a.tb[,-1])
    rownames(ht.a) <- ht.a.tb[[1]]
    print(head(ht.a[,1:3]))
    f.rm.a <- colSums(ht.a) < 1
    cat(sprintf("samples with meta.clusters' frequency sum < 1 (data A):\n"))
    print(colnames(ht.a)[f.rm.a])
    ht.a <- ht.a[,!f.rm.a]
    dat.in.a <- ht.a
    dat.in <- dat.in.a

    if(!is.null(dat.plot.b))
    {
        ht.b.tb <- dcast(dat.plot.b,group.var~donor.var,value.var="freq",fill=0)
        ht.b <- as.matrix(ht.b.tb[,-1])
        rownames(ht.b) <- ht.b.tb[[1]]
        print(head(ht.b[,1:3]))
        f.rm.b <- colSums(ht.b) < 1
        cat(sprintf("samples with meta.clusters' frequency sum < 1 (data B):\n"))
        print(colnames(ht.b)[f.rm.b])
        ht.b <- ht.b[,!f.rm.b]
        f.samples <- intersect(colnames(dat.in.a),colnames(ht.b))
        dat.in.a <- dat.in.a[,f.samples,drop=F]
        dat.in.b <- ht.b[,f.samples,drop=F]
        dat.in <- rbind(dat.in.a,dat.in.b)
        mcls.color.b <- colSet[[group.var]][rownames(ht.b)]
    }


    if(fltBySD)
    {
        rsd <- rowSds(dat.in)
        th.sd <- quantile(rsd,0.5)
        rn.flt <- rownames(dat.in)[ rsd >= th.sd ]
        dat.in <- dat.in[rn.flt,,drop=F]
        dat.in.a <- dat.in.a[rownames(dat.in.a) %in% rn.flt,,drop=F]
        if(!is.null(dat.plot.b)){
            dat.in.b <- dat.in.b[rownames(dat.in.b) %in% rn.flt,,drop=F]
        }
    }

    mcls.color <- colSet[[group.var]][rownames(dat.in)]
    mcls.color.a <- colSet[[group.var]][rownames(dat.in.a)]
    if(!is.null(dat.plot.b)){
	    mcls.color.b <- colSet[[group.var]][rownames(dat.in.b)]
    }
    
    mcls.color.debug <<- mcls.color
    dat.in.debug <<- dat.in

    #### try PCA
    {
        sce.f <- ssc.build(dat.in)
        rowData(sce.f)$rsd <- rowSds(dat.in)
        rowData(sce.f)$gene.all <- T
        rowData(sce.f)
        sce.f$cancerType <- mapping.vec[colnames(sce.f)]
        sce.f <- ssc.reduceDim(sce.f,method="pca",method.vgene="gene.all",pca.npc=5,
                               tSNE.perplexity=30,
                               seed=12345)
        p <- ssc.plot.tsne(sce.f,columns="cancerType",colSet=colSet,reduced.name="pca")
        ggsave(file=sprintf("%s.pca.cancerType.png",out.prefix),width=6,height=4)
        #p <- ssc.plot.tsne(sce.f,gene=c("CD8.c12.Tex.CXCL13","CD8.c10.Trm.ZNF683","CD8.c15.ISG.IFIT1"),
        #				   par.geneOnTSNE = list(scales = "free",pt.alpha=1),
        #				   colSet=colSet,reduced.name="pca")
        #ggsave(file=sprintf("%s.pca.gene.png",out.prefix),width=10,height=3)
        #p <- ssc.plot.tsne(sce.f,gene=c("CD8.c12.Tex.CXCL13","CD8.c10.Trm.ZNF683","CD8.c15.ISG.IFIT1"),
        #				   par.geneOnTSNE = list(scales = "free",pt.alpha=1),
        #				   colSet=colSet,reduced.name="pca.tsne")
        #ggsave(file=sprintf("%s.pca.gene.tsne.png",out.prefix),width=10,height=3)
    }

    ######## correlation and select variables
    {
        cor.dat <- cor(t(dat.in),method="spearman")
        ###cor.dat <- cor(t(dat.in),method="pearson")
        #cor.hclust <- run.cutreeDynamic(dat.in,method.distance="cosine",deepSplit=2, minClusterSize=2)
        #cor.hclust <- run.cutreeDynamic(dat.in,method.distance="spearman",deepSplit=2, minClusterSize=2)
        cor.hclust <- run.cutree(dat.in,method.distance="spearman",k=9)
        cor.hclust$branch <- dendextend::set(cor.hclust$branch,"branches_lwd", 1.5)
        sscVis:::plotMatrix.simple(cor.dat,show.dendrogram=T,
                       out.prefix=sprintf("%s.cor",out.prefix),
                       show.number=F,
                       z.lo = -0.4, z.hi=0.4,
                       clust.column=cor.hclust$branch,
                       clust.row=cor.hclust$branch,
                       exp.name=expression(Corr),
                       row_dend_width = unit(2.8, "cm"),
                       par.heatmap=list(
                                column_dend_height = unit(2.8, "cm")
                                ),
                       palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                       pdf.width=9,pdf.height=9)

        #### filter by TH.cor
        mcls.cor.max <- apply(cor.dat,1,function(x){ max(abs(x[x!=1])) })
        f.cor <- mcls.cor.max >= TH.cor
        if(!is.null(must.include)){
            f.cor <- f.cor | (rownames(cor.dat) %in% must.include)
        }
        cat(sprintf("mcls.cor.max:"))
        print(sort(mcls.cor.max))
        cor.dat.flt <- cor.dat[f.cor,f.cor]
        if(is.null(k.corMat)){
            cor.hclust.flt <- run.cutreeDynamic(dat.in[f.cor,],method.distance="cosine",deepSplit=2, minClusterSize=2)
        }else{
            #cor.hclust.flt <- run.cutree(dat.in[f.cor,],method.distance="cosine",k=k.corMat)
            cor.hclust.flt <- run.cutree(dat.in[f.cor,],method.distance="spearman",k=k.corMat)
        }
        cor.hclust.flt$branch <- dendextend::set(cor.hclust.flt$branch,"branches_lwd", 1.5)
        sscVis::plotMatrix.simple(cor.dat.flt,show.dendrogram=T,
                      out.prefix=sprintf("%s.cor.flt",out.prefix),
                      show.number=F,
                      clust.column=cor.hclust.flt$branch,
                      clust.row=cor.hclust.flt$branch,
                      exp.name=expression(Corr),
                      #z.lo = -0.5, z.hi=0.5,
                      #par.legend=list(at=seq(-0.5,0.5,0.25)),
                      z.lo = -0.4, z.hi=0.4,
                      par.heatmap=list(row_dend_width = unit(1.5, "cm"),
                               row_names_gp=gpar(fontsize=10),
                               column_names_gp=gpar(fontsize=10),
                               column_dend_height = unit(1.5, "cm")),
                      palatte=rev(brewer.pal(n = 7, name = "RdBu")),
                      pdf.width=pdf.width.cor,pdf.height=pdf.height.cor)

    }

    if(!is.null(dat.plot.b)){
	axis_param.column <- default_axis_param("column")
	axis_param.column$gp <- gpar(fontsize=12)
	ha.col <- HeatmapAnnotation(cancerType=mapping.vec[colnames(dat.in)],
				    CD4=anno_barplot(t(dat.in.b),gp = gpar(fill = mcls.color.b,col=NA),
						     axis_param=axis_param.column,border=F,
						     bar_width=1, height = unit(3.0*ann.bar.r, "cm")),
				    CD8=anno_barplot(t(dat.in.a), gp = gpar(fill = mcls.color.a,col=NA),
						     axis_param=axis_param.column,border=F,
						     bar_width=1, height = unit(3.0*ann.bar.r, "cm")),
				    col = colSet["cancerType"],
				    #annotation_legend_param = list(),
				    show_legend = T,
				    ###simple_anno_size = unit(1, "cm"))
				    simple_anno_size = unit(0.5, "cm"))
    }else{
	ha.col <- HeatmapAnnotation(cancerType=mapping.vec[colnames(dat.in)],
				    CD8=anno_barplot(t(dat.in.a), gp = gpar(fill = mcls.color.a,col=NA),
						     bar_width=1, height = unit(2.5, "cm")),
				    col = colSet["cancerType"],
				    #annotation_legend_param = list(),
				    show_legend = T,
				    simple_anno_size = unit(1, "cm"))
    }

    makeHTPlot <- function(dat.in,out.prefix)
    {

        ###obj.hclust.col <- run.cutree(t(dat.in),k=k.column)
        ##obj.hclust.col <- run.cutree(t(dat.in),k=k.column,method.distance="spearman")
        ###obj.hclust.col <- run.cutreeDynamic(t(dat.in),deepSplit=1, minClusterSize=2)
        obj.hclust.col <- run.cutree(t(dat.in),k=k.column,method.distance="")
        obj.hclust.row <- run.cutree(dat.in,k=k.row)
        #obj.hclust.col <- run.cutree(t(dat.in),k=k.column,method.distance="cosine")
        #obj.hclust.row <- run.cutree(dat.in,k=k.row,method.distance="cosine")

        obj.hclust.col$branch <- dendextend::set(obj.hclust.col$branch,"branches_lwd", 2)
        obj.hclust.row$branch <- dendextend::set(obj.hclust.row$branch,"branches_lwd", 2.5)
        ##sscClust:::plot.matrix.simple(t(scale(t(dat.in))), show.dendrogram=T,
        dat.in.bin <- floor(floor(dat.in * 100)/5)
        dat.in.bin[ dat.in.bin > 9 ] <- 9
        dat.in.bin[1:3,1:4]
        #print(str(dat.in))
        dat.in.z <- t(scale(t(dat.in)))
        dat.in.z[dat.in.z > 3] <- 3
        dat.in.z[dat.in.z < -3] <- -3

        obj.hclust.row.debug <<- obj.hclust.row
        ha.row = rowAnnotation(overall=anno_boxplot(dat.in,outline = FALSE,width = unit(4, "cm"),
                                gp = gpar(fill = mcls.color[rownames(dat.in)])))

        sscVis::plotMatrix.simple(
                      #dat.in.z,
                      #palatte=rev(brewer.pal(n = 7,name = "RdBu")),
                      #par.legend=list(ncol = 1,at=seq(-3,3,1.5)),
                      dat.in.bin,
                      col.ht=rev(structure(viridis::magma(10),names=0:9 )),
                      par.legend=list(ncol = 1,labels=rev(sprintf("%s%%~%s%%",5*(0:9),c(5*(1:9),100) ))),
                      ####col.ht=rev(structure(colorRampPalette(brewer.pal(9,name="Blues"))(10),names=0:9 )),
                      show.dendrogram=T,
                      out.prefix=sprintf("%s",out.prefix),
                      show.number=F,
                      clust.row=obj.hclust.row$branch,
                      clust.column=obj.hclust.col$branch,
                      exp.name=expression(Freq),
                      top_annotation = ha.col,
                      column.split=if(k.column>=2) k.column else NULL,
                      par.heatmap=list(
                               #column_split=if(k.column>=2) k.column else NULL,
                               row_names_gp=gpar(fontsize=12),
                               right_annotation = ha.row,
                               column_gap = unit(0.8, "mm"),border=FALSE,
                               row_dend_width = unit(1.5, "cm"),
                               column_dend_height = unit(1.5, "cm")
                               ),
                      #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                      #palatte=viridis::magma(7),
                      #palatte=rev(brewer.pal(n = 7,name = "RdBu")),
                      pdf.width = pdf.width, pdf.height = pdf.height,...)
        return(list("obj.hclust.col"=obj.hclust.col,"obj.hclust.row"=obj.hclust.row,"cor.dat.flt"=cor.dat.flt))
    }

    ret.list.allRow <- makeHTPlot(dat.in,out.prefix)
    ret.list.fltRow <- makeHTPlot(dat.in[f.cor,,drop=F],sprintf("%s.flt",out.prefix))

    return(list("all"=ret.list.allRow,"flt"=ret.list.fltRow))
}

sigGeneVennPlot <- function(v.list,background.list,out.prefix,col.venn=NULL,fill.venn=NULL,venn.alpha=c(1,0.7,0.7))
{

    for(i in seq_along(v.list)){
	v.list[[i]] <- intersect(v.list[[i]],background.list)
    }
    n.sample <- length(unique(background.list))
    cat(sprintf("n.sample: %d\n",n.sample))

    library("VennDiagram")
    library("gridBase")
    pdf(sprintf("%s.venn.pdf",out.prefix),width=7,height=7)
    if(length(v.list)==2){
        opar <- par(mar=c(7,7,7,7),cex.lab=1.5,cex.main=1.5,xpd=T)
    }else{
        opar <- par(mar=c(2,2,2,2),cex.lab=1.5,cex.main=1.5,xpd=T)
    }
    plot.new()
    #title(main="venn",sub=sprintf("p=%s",p.value),cex=1.2)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    ##venn.plot <- venn.diagram(v.list,filename = NULL, cat.cex=1.5,cex=1.5,
    venn.plot <- venn.diagram(v.list,filename = NULL, cat.cex=1.5,cex=1.5,
                              hyper.test=T,total.population=n.sample,lower.tail=F,sub.pos=c(0.5,0.25),sub.cex=1.5,
							  col=if(is.null(col.venn)) "black" else col.venn,fill=col.venn,
							  alpha=venn.alpha[seq_len(length(v.list))],
                              margin=0.2,cat.dist=0.15,na = "remove")
    grid.draw(venn.plot)
    par(opar)
    dev.off()
}

sigGeneVennTable <- function(gene.tb,cmp,only.sig=T)
{
	##col.out <- c("comb.ES","comb.ES.sd","comb.Z","comb.p","comb.padj","sig","sig.cate")
	col.out <- c("comb.ES","comb.ES.sd","comb.p","comb.padj","sig","sig.cate")
	mapping.TF.tb <- unique(gene.tb[,c("geneID","geneSet.TF"),with=F])
	mapping.TF.vec <- structure(mapping.TF.tb$geneSet.TF,names=mapping.TF.tb$geneID)
	out.tb <- dcast(gene.tb[meta.cluster %in% cmp,c("geneID","meta.cluster",col.out),with=F],
					geneID~meta.cluster,value.var=col.out)
	out.tb[,geneSet.TF:=mapping.TF.vec[geneID]]
	f.gene <- apply(out.tb[,grepl("^sig_",colnames(out.tb),perl=T),with=F],1,
					function(x){ any(x==T,na.rm=T) })
	if(only.sig==T){
		out.tb <- out.tb[f.gene,]
	}
	return(out.tb)
}

run.nichenet <- function(gene.oi,gene.bg,out.prefix,
						 comb.width=12,comb.height=7,comb.rel.height=c(10,2),comb.rel.width=c(0.1,0.9),
						 lr.width=5,lr.height=7,n.top=20,pearson.max=0.2,
						 ligands_all=NULL,targets_all=NULL,es.df=NULL,
						 do.eval=T,eval.k.fold=5,eval.n.round=10)
{
	library("nichenetr")
	library("tidyverse")
    library("nichenetr")
    library("ggpubr")
    library("ggplot2")
    library("cowplot")

	weighted_networks <- readRDS("../data/external/Nichnet/weighted_networks.rds")
	ligand_target_matrix <- readRDS("../data/external/Nichnet/ligand_target_matrix.rds")
	lr_network <- readRDS("../data/external/Nichnet/lr_network.rds")
	str(ligand_target_matrix)
	lr_network_expressed <- lr_network %>% filter(to %in% gene.bg)
	potential_ligands  <-  lr_network_expressed %>% pull(from) %>% unique()
	ligand_activities <- predict_ligand_activities(geneset = gene.oi,
						       background_expressed_genes = gene.bg,
						       ligand_target_matrix = ligand_target_matrix,
						       potential_ligands = potential_ligands)
	ligand_activities %>% arrange(-pearson)
	best_upstream_ligands <- ligand_activities %>% top_n(n.top, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
	###
	p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) +
		geom_histogram(color="black", fill="darkorange")  +
		# geom_density(alpha=.1, fill="orange") +
		geom_vline(aes(xintercept=min(ligand_activities %>% top_n(n.top, pearson) %>% pull(pearson))),
				   color="red", linetype="dashed", size=1) +
		labs(x="ligand activity (PCC)", y = "# ligands") +
		theme_classic()
	ggsave(sprintf("%s.dist.pearson.pdf",out.prefix),width=5,height=4)
	###
	active_ligand_target_links_df <- best_upstream_ligands %>%
	    lapply(get_weighted_ligand_target_links,geneset = gene.oi,
		   ligand_target_matrix=ligand_target_matrix, n=250) %>%
	    bind_rows()
	active_ligand_target_links_df <- active_ligand_target_links_df %>% filter(!is.na(weight))
	active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df=active_ligand_target_links_df,
									  ligand_target_matrix = ligand_target_matrix,
									  cutoff = 0.25)
	order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
	order_targets = active_ligand_target_links_df$target %>% unique()
	order_targets = intersect(order_targets,rownames(active_ligand_target_links))
	active_ligand_target_links.debug <<- active_ligand_target_links
	order_targets.debug <<- order_targets
	order_ligands.debug <<- order_ligands
	vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
	colnames(vis_ligand_target) <- make.names(colnames(vis_ligand_target))
	gene.unexplained <- setdiff(gene.oi,colnames(active_ligand_target_links))
	vis_ligand_target.debug <<- vis_ligand_target
	vis_ligand_target[vis_ligand_target>0.01] <- 0.01

	p_ligand_target_network = vis_ligand_target %>%
				    make_heatmap_ggplot("Prioritized ligands","genes in receiver cells",
							color = "purple",legend_position = "top",
							x_axis_position = "top",
							legend_title = "Regulatory potential") +
					scale_fill_gradient2(low = "whitesmoke",  high = "purple",
							     limits=c(0,0.01), breaks = c(0,0.005,0.01)) +
					theme(axis.text.x = element_text(face = "italic"),legend.key.width=unit(2,"cm"))
	#ggsave(file=sprintf("%s.ligand_target_network.tmp.pdf",out.prefix),width=8,height=5)

	####
	ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>%
				    magrittr::set_rownames(ligand_activities$test_ligand)

	vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
	vis_ligand_pearson.debug <<- vis_ligand_pearson
	vis_ligand_pearson[vis_ligand_pearson > pearson.max ] <- pearson.max
	p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity",
							color = "darkorange",legend_position = "top",
							x_axis_position = "top",
							legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") +
				scale_fill_gradient2(low = "whitesmoke",  high = "darkorange", limits=c(0,pearson.max))
	p_ligand_pearson <- p_ligand_pearson + theme(legend.key.width=unit(2,"cm"))
	#if(!is.null(keywidth)){
	#	p_ligand_pearson <- p_ligand_pearson + theme(legend.key.width=keywidth)
	#}
	#ggsave(file=sprintf("%s.ligand_pearson.pdf",out.prefix),width=5,height=5)

	{
	####
	figures_without_legend = plot_grid(p_ligand_pearson +
					   theme(legend.position = "none", axis.ticks = element_blank()) +
					   theme(axis.title.x = element_text()),
				       p_ligand_target_network + theme(legend.position = "none",
								       axis.text.y=element_blank(),
								       axis.title.y = element_blank(),
								       axis.ticks = element_blank()),
					   #### donnot let align contain "v" here
					   align = "h", nrow = 1,
					   rel_widths = comb.rel.width,
					   rel_heights = c(nrow(vis_ligand_pearson),
							   nrow(vis_ligand_target) + 3))

	legends = plot_grid(as_ggplot(get_legend(p_ligand_pearson)),
			    as_ggplot(get_legend(p_ligand_target_network)),
			    nrow = 2, align = "h")
	pp <- plot_grid(figures_without_legend, legends, rel_heights = comb.rel.height, nrow = 2, align = "hv")
	ggsave(file=sprintf("%s.comb.pdf",out.prefix),width=comb.width,height=comb.height)
	}

	#### ligand-receptor network
	{
	    # get the ligand-receptor network of the top-ranked ligands
	    lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% gene.bg) %>% distinct(from,to)
	    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

	    # get the weights of the ligand-receptor interactions as used in the NicheNet model
	    lr_network_top_df = weighted_networks$lr_sig %>%
							    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

	    # convert to a matrix
	    lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
	    lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

	    # perform hierarchical clustering to order the ligands and receptors
	    dist_receptors = dist(lr_network_top_matrix, method = "binary")
	    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
	    order_receptors = hclust_receptors$labels[hclust_receptors$order]

	    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
	    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
	    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

	    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
	    p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>%
									    make_heatmap_ggplot("Prioritized ligands","Receptors expressed",
														    color = "mediumvioletred",
														    x_axis_position = "top",
														    legend_title = "Prior interaction potential")
	    ggsave(file=sprintf("%s.lr.matrix.pdf",out.prefix),width=lr.width,height=lr.height)

	}

	if(!is.null(ligands_all)){
		if(is.null(targets_all)){
			#targets_all <- active_ligand_target_links_df %>% filter(ligand %in% ligands_all) %>% pull(target)
			targets_all <- as.data.frame(t(active_ligand_target_links)) %>% rownames_to_column("ligand") %>%
				gather(target,weight,-ligand) %>% 
				filter(ligand %in% ligands_all & weight>0) %>% pull(target) %>% unique()
		}
		ligand_tf_matrix <- readRDS("../data/external/Nichnet/ligand_tf_matrix.rds")
		sig_network <- readRDS("../data/external/Nichnet/signaling_network.rds")
		gr_network <- readRDS("../data/external/Nichnet/gr_network.rds")
		active_signaling_network <- get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix,
															 ligands_all = ligands_all,
															 targets_all = targets_all,
															 weighted_networks = weighted_networks)

		# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
		active_signaling_network_min_max = active_signaling_network
		active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
		active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

		graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max,
														  ligands_all = ligands_all,
														  targets_all = targets_all,
														  sig_color = "indianred", gr_color = "steelblue")

		###gg <- DiagrammeR::render_graph(graph_min_max, layout = "tree")
		###DiagrammeR::export_graph(graph_min_max,file_name=sprintf("%s.test.png",out.prefix),file_type="png",width=7,height=7)
		data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network,
														   lr_network = lr_network,
														   sig_network = sig_network,
														   gr_network = gr_network)
		head(data_source_network) 
		####
		{
		    write_output = TRUE # change to TRUE for writing output

		    # weighted networks ('import network' in Cytoscape)
		    if(write_output){
			    bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"),
					      active_signaling_network$gr %>% mutate(layer = "regulatory")) %>%
				    write_tsv(sprintf("%s.weighted_signaling_network.txt",out.prefix))
		    }

		    # networks with information of supporting data sources ('import network' in Cytoscape)
		    if(write_output){
			    data_source_network %>% write_tsv(sprintf("%s.data_source_network.txt",out.prefix))
		    }

		    # Node annotation table ('import table' in Cytoscape)
		    specific_annotation_tbl = bind_rows(
											    tibble(gene = ligands_all, annotation = "ligand"),
											    tibble(gene = targets_all, annotation = "target"),
											    tibble(gene = c(data_source_network$from, data_source_network$to) %>%
													    unique() %>% setdiff(c(targets_all,ligands_all)) %>%
													    intersect(lr_network$to %>% unique()), annotation = "receptor"),
											    tibble(gene = c(data_source_network$from, data_source_network$to) %>%
												       unique() %>% setdiff(c(targets_all,ligands_all)) %>%
												       intersect(gr_network$from %>% unique()) %>%
												       setdiff(c(data_source_network$from, data_source_network$to) %>%
														       unique() %>% intersect(lr_network$to %>% unique())),
												       annotation = "transcriptional regulator"))
		    non_specific_annotation_tbl = tibble(gene = c(data_source_network$from, data_source_network$to) %>%
											     unique() %>% setdiff(specific_annotation_tbl$gene),
										     annotation = "signaling mediator")

		    node.annotation.tbl <- bind_rows(specific_annotation_tbl,non_specific_annotation_tbl)
		    if(!is.null(es.df)){
			    node.annotation.tbl <- left_join(node.annotation.tbl,es.df,by="gene")
		    }
		    if(write_output){
			     write_tsv(node.annotation.tbl,sprintf("%s.annotation_table.txt",out.prefix))
		    }

		}

	}
	
	### assess how well top-ranked ligands can predict a gene set of interest
	eval.out.list <- NULL
	if(do.eval)
	{
		# change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
		#eval.k.fold = 5
		#eval.n.round = 2
		eval_gene_predictions_top20_list = seq(eval.n.round) %>% lapply(assess_rf_class_probabilities,
															 folds = eval.k.fold, 
															 geneset = gene.oi,
															 background_expressed_genes = gene.bg,
															 ligands_oi = best_upstream_ligands,
															 ligand_target_matrix = ligand_target_matrix)
		target_prediction_performances_cv = eval_gene_predictions_top20_list %>%
													lapply(classification_evaluation_continuous_pred_wrapper) %>%
													bind_rows() %>%
													mutate(round=seq(1:nrow(.)))

		target_prediction_performances_discrete_cv = eval_gene_predictions_top20_list %>%
													lapply(calculate_fraction_top_predicted, quantile_cutoff = 0.95) %>%
													bind_rows() %>% ungroup() %>%
													mutate(round=rep(1:length(eval_gene_predictions_top20_list), each = 2))
		target_prediction_performances_discrete_fisher = eval_gene_predictions_top20_list %>%
															lapply(calculate_fraction_top_predicted_fisher,
																   quantile_cutoff = 0.95)

		top_predicted_genes = seq(length(eval_gene_predictions_top20_list)) %>%
								lapply(get_top_predicted_genes,eval_gene_predictions_top20_list) %>%
								reduce(full_join, by = c("gene","true_target"))
		top_predicted_genes %>% filter(true_target)
		f.gene.consensus <- apply(top_predicted_genes[,grep("^predicted_top_target_round",
															colnames(top_predicted_genes),value=T)],1,
								  function(x){ mean(x,na.rm=T) >= 0.5 })
		gene.target.missing <- top_predicted_genes[f.gene.consensus,] %>% pull(gene) %>% setdiff(gene.oi,.)
		eval.out.list <- list(target_prediction_performances_cv=target_prediction_performances_cv,
							  target_prediction_performances_discrete_cv=target_prediction_performances_discrete_cv,
							  top_predicted_genes=top_predicted_genes,
							  gene.target.missing=gene.target.missing,
							  avg.auroc=target_prediction_performances_cv$auroc %>% mean(),
							  avg.aupr=target_prediction_performances_cv$aupr %>% mean(),
							  avg.pearson=target_prediction_performances_cv$pearson %>% mean(),
							  avg.true_target.fraction_positive_predicted= target_prediction_performances_discrete_cv %>%
								  filter(true_target) %>% .$fraction_positive_predicted %>% mean(),
							  avg.false_target.fraction_positive_predicted= target_prediction_performances_discrete_cv %>%
								  filter(!true_target) %>% .$fraction_positive_predicted %>% mean(),
							  avg.target_prediction_performances_discrete_fisher=
								  target_prediction_performances_discrete_fisher %>% unlist() %>% mean()
							  )

	}
	
	return(list("p_ligand_target_network"=p_ligand_target_network,"p_ligand_pearson"=p_ligand_pearson,
				"p_ligand_receptor_network"=p_ligand_receptor_network,
				"gene.oi"=gene.oi,
				"gene.bg"=gene.bg,
				"eval.list"=eval.out.list,
				"gene.unexplained"=gene.unexplained))

}

#### also call run.nichenet()
run.venn.nicheNet.g2 <- function(mcls,sname,out.prefix,gene.bg,...)
{
    require("tictoc")
    sig.list.all <- list(gene.core.tb[meta.cluster==mcls[1] & sig==T,][["geneSymbol"]],
             gene.core.tb[meta.cluster==mcls[2] & sig==T,][["geneSymbol"]])
    names(sig.list.all) <- sname
    sigGeneVennPlot(sig.list.all,
            background.list=unique(gene.core.tb$geneSymbol),
            col.venn=colSet$meta.cluster[mcls],
            fill.venn=colSet$meta.cluster[mcls],
            out.prefix=sprintf("%s.venn.%s.%s.g2",out.prefix,sname[1],sname[2]),
            venn.alpha=c(0.7,0.7))

    sig.list.TF <- list(gene.core.tb[meta.cluster==mcls[1] & sig==T & geneSet.TF==T,][["geneSymbol"]],
            gene.core.tb[meta.cluster==mcls[2] & sig==T & geneSet.TF==T,][["geneSymbol"]])
    names(sig.list.TF) <- sname
    sigGeneVennPlot(sig.list.TF,
            background.list=unique(gene.core.tb[geneSet.TF==T,][["geneSymbol"]]),
            col.venn=colSet$meta.cluster[mcls],
            fill.venn=colSet$meta.cluster[mcls],
            out.prefix=sprintf("%s.venn.%s.%s.g2.TF",out.prefix,sname[1],sname[2]),venn.alpha=c(0.7,0.7))

    ###return(NULL)

    sig.TF.cmp2.tb <-  sigGeneVennTable(gene.core.tb,
                    cmp=mcls,
                    only.sig=T)
    conn <- gzfile(sprintf("%s.cmp.%s.%s.g2.txt.gz",out.prefix,sname[1],sname[2]),"w")
    write.table(sig.TF.cmp2.tb,file=conn,row.names=F,sep="\t",quote=F)
    close(conn)

    gene.cmp.g2.tb <- sig.TF.cmp2.tb
    #### comp2
    sig.name <- grep("^sig_",colnames(gene.cmp.g2.tb),value=T)
    f.gene <- rowSums(gene.cmp.g2.tb[,sig.name,with=F])==length(sig.name)
    gene.oi <- gene.cmp.g2.tb[f.gene,][["geneID"]]
    es.tb <- dcast(gene.core.tb[meta.cluster %in% mcls,
                 c("geneID","meta.cluster","comb.ES")],geneID~meta.cluster,value.var="comb.ES")
    es.avg.df <- tibble(gene=es.tb$geneID,avg.comb.ES=rowMeans(es.tb[,-("geneID")]))
    print(str(gene.bg))
    print(str(gene.oi))
    print(head(es.avg.df))
    #dat.debug <<- list("gene.bg"=gene.bg,"gene.oi"=gene.oi,"es.avg.df"=es.avg.df)

    tic(sprintf("run.nichenet (%s.%s.g2)",sname[1],sname[2]))
    res.nichenet.g2 <- run.nichenet(gene.oi,
                    gene.bg,out.prefix=sprintf("%s.nichenet.cmp.%s.%s.g2",
                                   out.prefix,sname[1],sname[2]),
                    es.df=es.avg.df,...)
    toc()
    return(res.nichenet.g2)
    
}

gen.gsea.script <- function(gene.desc.top,sh.dir,out.prefix,db.file,
							group.var="meta.cluster",rank.var="comb.ES",
							gsea.scoring.scheme="classic",
							gsea.max=500)
{
	for(mcls in unique(sort(gene.desc.top[[group.var]])) ){
		filename.geneSymbol <- sprintf("%s.geneSymbol.%s.rnk",out.prefix,mcls)
		dir.create(dirname(filename.geneSymbol),F,T)
		f.mcls <- gene.desc.top[[group.var]]==mcls
		write.table(gene.desc.top[f.mcls,c("geneSymbol",rank.var),with=F],
					file = filename.geneSymbol,
					quote = F,sep = "\t",row.names = F,col.names = F)

		odir <- sprintf("%s.%s.%s.GSEA",out.prefix,mcls,gsea.scoring.scheme)
		dir.create(odir,showWarnings = F,recursive = T)
		dir.create(sh.dir,showWarnings = F,recursive = T)
		sink(file = sprintf("%s/external.GSEA.%s.job",sh.dir,mcls))
		cat("#!/bin/bash\n")
		cat("#SBATCH -p q.all\n")
		cat("#SBATCH -N 1\n")
		cat("#SBATCH --ntasks-per-node=1\n")
		cat("#SBATCH --mem=2G\n")
		cat(sprintf("#SBATCH -o gsea.%s.%%j.out\n",mcls))
		cat(sprintf("#SBATCH -e gsea.%s.%%j.err\n",mcls))
		cat("#SBATCH --no-requeue\n")
		cmd.str <- sprintf("%s GSEAPreranked -gmx %s -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -create_svgs true -include_only_symbols true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max %s -set_min 5 -zip_report false -out %s \n",
                           system2("which","gsea-cli.sh", stdout = T),
						   db.file,
						   filename.geneSymbol,
						   gsea.scoring.scheme,
						   sprintf("%s",mcls),
						   gsea.max,odir)
		cat(cmd.str)
		#ret <- system(cmd.str)
		sink()
	}
}

