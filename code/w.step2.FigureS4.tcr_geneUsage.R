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

oDir = "./OUT_FigS4"
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



# 2. plot
FUNC_PLOT = function(df, prefix, n.top, n.min){
	# to clone level
	df = df[,c("cloneID", "A1V", "A1J", "A1VJ", "B1V", "B1J", "B1VJ")] %>% unique

	### part1: AV, AJ
	freq.AV = table(df$A1V) %>% as.data.frame
	freq.AV$Percent = 100 * freq.AV$Freq / sum(freq.AV$Freq)
	freq.AV$type = "V"
	freq.AV = freq.AV[order(freq.AV$Freq,decreasing=T),]
	freq.AV = freq.AV[freq.AV$Freq>=n.min,] %>% head(n.top)
	#
	freq.AJ = table(df$A1J) %>% as.data.frame
	freq.AJ$Percent = 100 * freq.AJ$Freq / sum(freq.AJ$Freq)
	freq.AJ$type = "J"
	freq.AJ = freq.AJ[order(freq.AJ$Freq,decreasing=T),]
	freq.AJ = freq.AJ[freq.AJ$Freq>=n.min,] %>% head(n.top)
	#
	dat = rbind(freq.AV,freq.AJ)
	colnames(dat)[1] = "Gene"
	dat$Gene = factor(dat$Gene, levels=dat$Gene)
	# plot
	width = 1 + 0.15*length(dat$Gene)
	height = 3
	p = ggplot(dat, aes(x=Gene, y=Percent, fill=type)) +
        geom_bar(stat="identity", width=0.8) +
        theme_classic() +
        theme(axis.text.x=element_text(size=7, vjust=1.0, hjust=1.1, angle=60), legend.position="none") +
        scale_fill_manual(values=get_palette("npg", length(unique(dat$type)))) +
		labs(x="", y="Percentage (%)", title="TCR alpha genes' usage")
	ggsave(p, file=sprintf("%s.TRA_V_J.pdf", prefix), width=width, height=height)
	
	### part2: BV, BJ
	freq.BV = table(df$B1V) %>% as.data.frame
	freq.BV$Percent = 100 * freq.BV$Freq / sum(freq.BV$Freq)
	freq.BV$type = "V"
	freq.BV = freq.BV[order(freq.BV$Freq,decreasing=T),]
	freq.BV = freq.BV[freq.BV$Freq>=n.min,] %>% head(n.top)
	#
	freq.BJ = table(df$B1J) %>% as.data.frame
	freq.BJ$Percent = 100 * freq.BJ$Freq / sum(freq.BJ$Freq)
	freq.BJ$type = "J"
	freq.BJ = freq.BJ[order(freq.BJ$Freq,decreasing=T),]
	freq.BJ = freq.BJ[freq.BJ$Freq>=n.min,] %>% head(n.top)
	#
	dat = rbind(freq.BV,freq.BJ)
	colnames(dat)[1] = "Gene"
	dat$Gene = factor(dat$Gene, levels=dat$Gene)
	# plot
	width = 1 + 0.15*length(dat$Gene)
	height = 3
	p = ggplot(dat, aes(x=Gene, y=Percent, fill=type)) +
        geom_bar(stat="identity", width=0.8) +
        theme_classic() +
        theme(axis.text.x=element_text(size=7, vjust=1.0, hjust=1.1, angle=60), legend.position="none") +
        scale_fill_manual(values=get_palette("npg", length(unique(dat$type)))) +
		labs(x="", y="Percentage (%)", title="TCR beta genes' usage")
	ggsave(p, file=sprintf("%s.TRB_V_J.pdf", prefix), width=width, height=height)
}


n.top = 20
n.min = 1
FUNC_PLOT(tcr, sprintf("%s/TCR.geneUsage", oDir), n.top, n.min)






