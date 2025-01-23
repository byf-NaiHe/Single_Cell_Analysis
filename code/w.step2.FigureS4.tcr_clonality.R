suppressPackageStartupMessages({
	library("reshape2")
	library("plyr")
	library("dplyr")
	library("data.table")
	library("ggpubr")
	library("ggsci")
	library("webr")
	library("doParallel")
})
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores=8)
options(stringsAsFactors=F)

oDir = "./OUT_FigS4"
dir.create(oDir, F, T)

nam.conv = fread(sprintf("%s/../../data/metaInfo/name.conversion.txt",oDir),sep="\t",stringsAsFactors=F,header=T)
nam.conv = as.data.frame(nam.conv)
rownames(nam.conv) = nam.conv$meta.cluster

tcr = readRDS(sprintf("%s/../../data/tcr/byCell/tcr.zhangLab.comb.flt.rds",oDir))

## calculation
tcr$cloneSize.group = tcr$cloneSize
tcr$cloneSize.group = ifelse(tcr$cloneSize.group>=3, ">=3", tcr$cloneSize.group)
tcr$cloneSize.group = factor(as.character(tcr$cloneSize.group), levels=c(">=3","2","1"))

stat.byClone = unique(tcr[,c("cloneID","meta.cluster","cloneSize.group")])[,c("meta.cluster","cloneSize.group")] %>% table %>% as.data.frame
stat.byCell = tcr[,c("meta.cluster","cloneSize.group")] %>% table %>% as.data.frame

## plot
FUNC_PLOT = function(df, ylab, outpre){
	### name conversion
	df$meta.cluster = nam.conv[as.character(df$meta.cluster),"cluster.name.full"]
	###
	# rank
	df = df %>% group_by(meta.cluster) %>% mutate(Sum=sum(Freq),Percent=Freq/sum(Freq)) %>% as.data.frame
    meta2num = df[,c("meta.cluster","Sum")] %>% unique %>% arrange(desc(Sum))
	df$meta.cluster = factor(as.character(df$meta.cluster), levels=as.character(meta2num$meta.cluster))
	# plot
	height = 5
	width = 0.2 + 0.2*length(unique(df$meta.cluster))
	p = ggplot(df, aes(x=meta.cluster, y=Freq, fill=cloneSize.group)) +
        geom_bar(stat="identity", width=0.8) +
        theme_classic() +
        theme(axis.text.x=element_text(size=8, vjust=1.0, hjust=1.1, angle=60),
			  plot.margin=margin(t=1, r=2, b=1, l=2, unit="cm"),
              legend.position="bottom",
              legend.text=element_text(size=8)) +
        scale_fill_manual(values=get_palette("locuszoom", length(unique(df$cloneSize.group)))) +
		ylab(ylab) +
		guides(fill= guide_legend(reverse=T))
	ggsave(p, file=sprintf("%s.barplot.pdf", outpre), width=width, height=height)
	
	### overall pie chart
	for (stype in c("CD4","CD8")){
		df.sub = df[grepl(stype, as.character(df$meta.cluster)),c("cloneSize.group","Freq")]
		df.sub = aggregate(Freq~cloneSize.group, df.sub, sum)
		p = PieDonut(df.sub, aes(cloneSize.group, count=Freq), r0=0, showPieName=F, explode=1) +
				scale_fill_manual(values=get_palette("locuszoom", length(unique(df.sub$cloneSize.group))))
		ggsave(p, file=sprintf("%s.pie_%s.pdf", outpre, stype), width=5, height=5)
	}
}


FUNC_PLOT(stat.byClone, "Num of Clones", sprintf("%s/TCR.clonality.byClone",oDir))
FUNC_PLOT(stat.byCell,  "Num of Cells",  sprintf("%s/TCR.clonality.byCell",oDir))







