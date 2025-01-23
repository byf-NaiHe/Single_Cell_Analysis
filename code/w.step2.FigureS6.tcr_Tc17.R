suppressPackageStartupMessages({
	library("reshape2")
	library("plyr")
	library("dplyr")
	library("data.table")
	library("ggpubr")
	library("ggsci")
	library("doParallel")
})
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores=8)
options(stringsAsFactors=F)

oDir = "./OUT_FigS6"
dir.create(oDir, F, T)

tcr = readRDS(sprintf("%s/../../data/tcr/byCell/tcr.zhangLab.comb.flt.rds",oDir))

# 1. preprocess
# patch: get the most frequent A1 and B1, because some have same freq of these two types.
tcr = ddply(tcr, .(cloneID), function(df){
        df$cloneID = gsub(":.*$","",as.character(df$cloneID))
        df$cloneSize = nrow(df)
        df$cloneSize.CD4 = nrow(df[grepl("^CD4",df$meta.cluster),])
        df$cloneSize.CD8 = nrow(df[grepl("^CD8",df$meta.cluster),])
        # patch: get the most frequent CDR3, because some have same freq of A1A2/B1G2
        A1 = table( df$Identifier.A1 )
		B1 = table( df$Identifier.B1 )
        A1 = names( A1[which.max(A1)] )
        B1 = names( B1[which.max(B1)] )
		##
		A1V = gsub("^.*(TRAV[0-9\\-]+).*$", "\\1", A1)
		A1J = gsub("^.*(TRAJ[0-9\\-]+).*$", "\\1", A1)
		B1V = gsub("^.*(TRBV[0-9\\-]+).*$", "\\1", B1)
		B1J = gsub("^.*(TRBJ[0-9\\-]+).*$", "\\1", B1)
		##
		df$A1 = A1
		df$A1V = A1V
		df$A1J = A1J
		df$A1VJ = sprintf("%s+%s",A1V,A1J)
		df$B1 = B1
		df$B1V = B1V
		df$B1J = B1J
		df$B1VJ = sprintf("%s+%s",B1V,B1J)
        #
        return(df)
}, .parallel=T)

tcr$meta.cluster = as.character(tcr$meta.cluster)
tcr.Tc17 = tcr[tcr$meta.cluster=="CD8.c16.MAIT.SLC4A10",]


# 2. plot
FUNC_PLOT = function(df, prefix, n.top, n.min){
	# to clone level
	df = df[,c("cloneID", "A1V", "A1J", "A1VJ", "B1V", "B1J", "B1VJ")] %>% unique

	### part1: TRA-VJpair
	dat = table(df$A1VJ) %>% as.data.frame
	dat$Percent = 100 * dat$Freq / sum(dat$Freq)
	dat = dat[order(dat$Freq,decreasing=T),]
	dat = dat[dat$Freq>=n.min,] %>% head(n.top)
	colnames(dat)[1] = "Gene"
	dat$Gene = factor(dat$Gene, levels=dat$Gene)
	# plot
	width = 1 + 0.15*length(dat$Gene)
	height = 3
	p = ggplot(dat, aes(x=Gene, y=Percent)) +
		geom_bar(stat="identity", width=0.8, fill="darkgrey") +
		theme_classic() +
		theme(axis.text.x=element_text(size=7, vjust=1.0, hjust=1.1, angle=60), legend.position="none") +
		labs(x="", y="Percentage (%)", title="TCR alpha genes pairs' usage")
	ggsave(p, file=sprintf("%s.TRA_VJ.pdf", prefix), width=width, height=height)
	
	### part2: TRB-VJpair
	dat = table(df$B1VJ) %>% as.data.frame
	dat$Percent = 100 * dat$Freq / sum(dat$Freq)
	dat = dat[order(dat$Freq,decreasing=T),]
	dat = dat[dat$Freq>=n.min,] %>% head(n.top)
	colnames(dat)[1] = "Gene"
	dat$Gene = factor(dat$Gene, levels=dat$Gene)
	# plot
	width = 1 + 0.15*length(dat$Gene)
	height = 3
	p = ggplot(dat, aes(x=Gene, y=Percent)) +
		geom_bar(stat="identity", width=0.8, fill="darkgrey") +
		theme_classic() +
		theme(axis.text.x=element_text(size=7, vjust=1.0, hjust=1.1, angle=60), legend.position="none") +
		labs(x="", y="Percentage (%)", title="TCR beta genes pairs' usage")
	ggsave(p, file=sprintf("%s.TRB_VJ.pdf", prefix), width=width, height=height)

}


n.top = 20
n.min = 1
FUNC_PLOT(tcr.Tc17, sprintf("%s/TCR.Tc17.geneUsage", oDir), n.top, n.min)




