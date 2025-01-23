suppressPackageStartupMessages({
	library("plyr")
	library("dplyr")
	library("ggpubr")
	library("reshape2")
	library("data.table")
	library("doParallel")
	library("ggrastr")
	library("meta")
	library("esc")
})
options(stringsAsFactors = FALSE)
pdf.options(reset=TRUE, onefile=FALSE)

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores=20)

oDir="./OUT_FigS32"
dir.create(sprintf("%s/",oDir), F, T)

metaInfo = read.table(sprintf("%s/../../data/metaInfo/TCGA.sig_scores.txt.gz",oDir), header=T, stringsAsFactors=F, check.names=F, sep="\t")
rownames(metaInfo) = metaInfo$sample
clusters = colnames(metaInfo)[22:ncol(metaInfo)]
metaInfo$cancerType = metaInfo$cancer.type.abbreviation
metaInfo$cancerType = ifelse(metaInfo$cancerType=="COAD"|metaInfo$cancerType=="READ", "CRC", metaInfo$cancerType)  # patch: merge COAD and READ as CRC

########################
# 1. TMB
dat = metaInfo
dat$TMB = dat$Non.silent.per.Mb
dat = dat[!is.na(dat$TMB),]
tmb.cutoff = 10
dat$TMB_group = ifelse(dat$TMB>tmb.cutoff, "high", "low")
dat$TMB_group = factor(dat$TMB_group, levels=c("low","high"))

## label typing biased cancers
tmp = dat[,c("cancerType","TMB_group")] %>% table %>% as.data.frame
bad_cancers = tmp[tmp$Freq < 10,"cancerType"] %>% unique %>% as.character
dat$biased = ifelse(dat$cancerType %in% bad_cancers, T, F )


## run
#test.clus = clusters
test.clus = "CD4.c17.TfhTh1.CXCL13"
l_ply ( test.clus, function(clu){
	this.dat = dat[,c("sample","cancerType", "TMB", "TMB_group", "biased", clu)]
	this.dat$Score = this.dat[,clu]
	
	title = sprintf("%s", clu)
	xlab = sprintf("TMB cutoff: %d", tmb.cutoff)
	
	### panCan - boxplot
	width = 4.2
	height = 4
	p = ggboxplot(this.dat, x="TMB_group", y="Score", outlier.shape=NA) +
		geom_point_rast(aes(color=cancerType), alpha=0.4, size=0.5, shape=16, raster.dpi=300, position=position_jitter(width=0.2, height=0)) + 
		theme_bw() +
    	        stat_compare_means(paired=F, size=4, method="t.test") + 
		scale_colour_manual(values=get_palette("aaas", length(unique(dat$cancerType)))) +
		ggtitle(title) + xlab(xlab)
	ggsave(p, file=sprintf("%s/TCGA.TMB.boxplot._%s_.pdf", oDir, clu), width=width, height=height, limitsize=F, useDingbats=F)

	### sepCan - meta
	this.dat.flt = this.dat[! this.dat$biased,]
	# t.test for each
	test = data.frame(cancerType=NULL, n1=NULL, n2=NULL, t=NULL, es=NULL, se=NULL, p=NULL)
	for (i in unique(this.dat.flt$cancerType)){
	    test.dat.1 = this.dat.flt[as.character(this.dat.flt$cancerType) == i & as.character(this.dat.flt$TMB_group) == "high", "Score"]
	    test.dat.2 = this.dat.flt[as.character(this.dat.flt$cancerType) == i & as.character(this.dat.flt$TMB_group) == "low", "Score"]
	    n1 = length(test.dat.1)
	    n2 = length(test.dat.2)
	    t = t.test(test.dat.1, test.dat.2) 
	    t.stat = t$statistic %>% unname
	    t.p = t$p.value
	    #
	    es.obj = esc_t(t=t.stat, grp1n=n1, grp2n=n2, es.type="g")
	    es = es.obj$es
	    se = es.obj$se
	    tmp = data.frame(cancerType=i, n1=n1, n2=n2, t=t.stat, es=es, se=se, p=t.p)
	    test = rbind(test, tmp, stringsAsFactors=F)
	}
	test$adj.p = p.adjust(test$p, method="BH")
	test$sig = ifelse(test$adj.p <0.05, T, F)
	test = arrange(test, cancerType)
	# combine using Hedges’g indexe
	meta = meta::metagen( TE=test$es, seTE=test$se, studlab=test$cancerType, data=test, 
                      comb.fixed=F, comb.random=T, prediction=F, sm="SMD")
	#
	width = 14
	height = 12
	pdf(sprintf("%s/TCGA.TMB.forrest._%s_.pdf", oDir, clu), width=width, height=height)
	meta::forest(meta, test.overall.random=T, digits.pval=4, colgap.forest.left="5cm", zero.pval=T, leftcols=c("studlab", "n1","n2","adj.p"))
	dev.off()

}, .parallel=T) 


########################
# 2. Specific Mutation
tcga.mut = fread(sprintf("%s/../../data/external/TCGA/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz",oDir), header=T, stringsAsFactors=F, check.names=F, sep="\t")
tcga.mut = as.data.frame(tcga.mut)
tcga.mut = tcga.mut[!is.na(tcga.mut$sample),]
rownames(tcga.mut) = tcga.mut$sample

#test.genes = c("TP53", "PIK3CA", "PTEN", "MUC16", "MUC4", "BRAF", "EGFR", "LRP1B", "FAT1", "KMT2D",
#  "NOTCH1", "KMT2C", "KRAS", "PBRM1", "PDE4DIP", "TP53", "PIK3CA", "PTEN", "MUC16",
#  "MUC4", "BRAF", "EGFR", "LRP1B", "FAT1", "KMT2D", "NOTCH1", "KMT2C", "KRAS", "PBRM1", "PDE4DIP")
test.genes = c("FAT1")
#test.clusters = clusters
test.clus = "CD4.c20.Treg.TNFRSF9"

common.patients = intersect(metaInfo$sample, colnames(tcga.mut)[2:ncol(tcga.mut)])
tcga.mut = t(tcga.mut[test.genes, common.patients])

dat = metaInfo[common.patients, c("sample","cancerType", test.clus)]
dat = cbind(dat, tcga.mut, stringsAsFactors=F)


# run
l_ply (  test.genes, function(gene){
  for (clu in test.clus){
	this.dat = dat[,c("sample","cancerType", gene, clu)]
	this.dat$Score = this.dat[,clu]
	this.dat$Group = this.dat[,gene]
	this.dat$Group = ifelse(this.dat$Group==0, "WT", "MT")
	this.dat$Group = factor(this.dat$Group, levels=c("WT","MT"))
	info = sprintf("%s-%s", clu, gene)
	
	## panCan - boxplot
	width = 4.2
	height = 4
	p = ggboxplot(this.dat, x="Group", y="Score", outlier.shape=NA) +
		geom_point_rast(aes(color=cancerType), alpha=0.4, size=0.5, shape=16, raster.dpi=300, position=position_jitter(width=0.2, height=0)) +
		theme_bw() +
		stat_compare_means(paired=F, size=4, method="t.test") + 
		scale_colour_manual(values=get_palette("aaas", length(unique(dat$cancerType)))) +
		ggtitle(info) + xlab("")
	ggsave(p, file=sprintf("%s/TCGA.mutGene.boxplot._%s_.pdf", oDir, info), width=width, height=height, limitsize=F, useDingbats=F)
  }
}, .parallel=T) 









