suppressPackageStartupMessages({
	library("SingleCellExperiment")
	library("GSVA")
	library("parallel")
	library("plyr")
	library("dplyr")
	library("data.table")
	library("ggpubr")
	library("ggsci")
})
options(stringsAsFactors=F)

oDir="./OUT_FigS38"
dir.create(oDir, F, T)

sce = readRDS(sprintf("%s/../../data/external/Hugo_ICB.bulkRNA.GSE78220.sce.rds",oDir))
sce = sce[,sce$Timepoint=="Pre" & ! is.na(sce$Response) & (sce$Response == "R" | sce$Response == "NR")]

geneSets = fread(sprintf("%s/../../data/metaInfo/immuneTyping.RF_gene_importance.txt.gz",oDir), header=T, stringsAsFactors=F)
geneSets = geneSets %>% arrange(immuneType, desc(importance), gene) %>% group_by(immuneType) %>% slice(1:50)  #top 50 for each type
inpList = list()
for (i in unique(geneSets$immuneType)){
	tmp = geneSets[geneSets$immuneType == i,]
	inpList[[i]] = tmp$gene
}

###### 1. calculate score
mat = as.matrix(assay(sce,"exprs"))
gsva.res = gsva(mat, inpList, method="gsva", kcdf="Gaussian", mx.diff=T, parallel.sz=8)
gsva.res = t(gsva.res)
info = cbind( colData(sce[,rownames(gsva.res)])[,c("Response"),drop=F], DataFrame(gsva.res) )
info = as.data.frame(info)
	
###### 2.barplot
pdat = info[,c("Response", "C1C2", "C3.8")]
pdat$id = rownames(info)
pdat$score = pdat$C1C2 - pdat$C3.8
##
cutoff = 0.06
pdat$flag = ifelse(pdat$score>0, "pos", "neg")
pdat.sub.good = pdat[abs(pdat$score)>=cutoff,]
test = fisher.test(table(pdat.sub.good[,c("Response","flag")]), workspace=1e4)
pdat.sub.bad = pdat[abs(pdat$score)<cutoff,]
##
width =  7.5
height = 3
p = ggbarplot(pdat, x="id", y="score", fill="Response", color="white", 
         palette="jco", sort.val="desc", sort.by.groups=FALSE, 
         x.text.angle=60, ylab="Score.C1C2 - Score.C3-8", legend.title="Response") +
	theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
		  plot.title=element_text(size=8, face="bold")) +
	ggtitle(sprintf("Fisher, p=%.3f", test$p.value)) +
	geom_bar(dat=pdat.sub.bad, fill="grey", stat="identity") +
	geom_hline(yintercept=cutoff, linetype="dashed", color="grey", size=0.6) + 
	geom_hline(yintercept=-cutoff, linetype="dashed", color="grey", size=0.6) 
ggsave(p, file=sprintf("%s/Hugo_ICB.immueType.barplot.pdf", oDir), width=width, height=height, limitsize=F)
	


