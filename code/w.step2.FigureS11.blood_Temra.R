suppressPackageStartupMessages({
library("reshape2")
library("stringr")
library("plyr")
library("dplyr")
library("ggpubr")
library("ggsci")
library("ggrastr")
library("Matrix")
library("data.table")
library("sscClust")
library("future")
library("Seurat")
library("sctransform")
library("R.utils")
library("doParallel")
library("tictoc")
library("cowplot")
library("AUCell")
})
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores=8)
plan(strategy="multiprocess", workers=8)
options(future.globals.maxSize=4*1024^3)
options(stringsAsFactors=F)
options(bitmapType='cairo')

oDir = "./OUT_FigS11/"
dir.create(oDir, F, T)

### 1. prepare
for (stype in  c("CD4","CD8")){

# read sce objs and merge them all
healthy.files = list.files(sprintf("%s/../../data/external/healthy_donor_blood_Tcells/%s/",oDir, stype), full.names=T)
patient.files = list.files(sprintf("%s/../../data/expression/%s/byDataset/",oDir, stype), pattern="sce.rds", full.names=T)

keep_cols = c("cellID","cellID.uniq","patient","patient.uid","batchV","libraryID","dataset","TCR","stype","cancerType","loc")

tmp.sce = readRDS(healthy.files[1])
gene_over = rownames(tmp.sce)
mat = assay(tmp.sce, "norm_exprs")
meta.data = colData(tmp.sce)[,keep_cols]

for(i in c(healthy.files[2:length(healthy.files)], patient.files)){
    this.sce = readRDS(i)
    this.sce = this.sce[,this.sce$loc=="P"]
    if( ncol(this.sce)<30 ){ next }
    this.sce = this.sce[!duplicated(rowData(this.sce)$display.name),]
    rownames(this.sce) = rowData(this.sce)$display.name
    colnames(this.sce) = this.sce$cellID.uniq
    this.sce$patient.uid = paste0(this.sce$dataset, ".", this.sce$patient)
    this.sce$stype = stype
    #
    gene_over = intersect(rownames(mat), rownames(this.sce))
    mat = mat[gene_over,]
    this.sce = this.sce[gene_over,]
    #
    mat = cbind(mat, assay(this.sce,"norm_exprs"))
    meta.data = rbind( meta.data, colData(this.sce)[,keep_cols])
}

# construct
sce = ssc.build(mat, assay="norm_exprs")
colData(sce) = meta.data
saveRDS(sce, sprintf("%s/health_and_deseased_blood.%s.sce.rds",oDir,stype))

}

### 2.AUCell
for (stype in  c("CD4","CD8")){

sce = readRDS(sprintf("%s/health_and_deseased_blood.%s.sce.rds",oDir,stype))
sig = fread(sprintf("%s/../../data/metaInfo/signature_genes.txt.gz", oDir), header=F, stringsAsFactors=F, sep="\t") %>% as.data.frame()
colnames(sig) = c("meta.cluster", "geneSymbol")
sig = sig[grepl("Temra", as.character(sig$meta.cluster)),]
sig = sig[grepl(stype, as.character(sig$meta.cluster)), "geneSymbol"]
geneSets = list(Temra=sig)

# run
cells_rankings = AUCell_buildRankings(assay(sce, "norm_exprs"))
cells_AUC = AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment = AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)

# info
all(colnames(cells_AUC@assays@data$AUC) == colnames(sce))
sce$AUC_Temra = cells_AUC@assays@data$AUC[1,]
sce$AUC_Temra.pred = ifelse( colnames(sce) %in% cells_assignment$Temra$assignment, "Temra", "Other")
info = as.data.frame(colData(sce))

# freq
if (stype=="CD4"){
    ## these datasets have enriched Treg 
    rm_datasets = c("CHOL.zhangLabSS2","HCC.zhangLabSS2","LUNG.zhangLabSS2","CRC.zhangLabSS2")
    info = info[! info$dataset %in% rm_datasets,]
}
freq = ddply(info, .(patient.uid,loc), function(.df){
	n_total = nrow(.df)
	n_temra = nrow(.df[.df$AUC_Temra.pred=="Temra",])
	freq = n_temra/n_total
	#
	df2 = unique(.df[,c("patient","patient.uid","dataset","cancerType","loc"),drop=F])
	df2$freq = freq
	df2$N.all = n_total
    #
    return(df2)
}, .parallel=T)

write.table(freq, file=sprintf("%s/health_and_deseased_blood.%s.freq.txt",oDir,stype), sep="\t", quote=F, col.names=T, row.names=F)

}


### 3. plot
for (stype in  c("CD4","CD8")){
min.N = 30

freq = read.table(sprintf("%s/health_and_deseased_blood.%s.freq.txt",oDir,stype), header=T, stringsAsFactors=F)
freq = freq[freq$N.all>min.N, ]
freq$type = ifelse(freq$cancerType=="HD","healthy","diseased")
freq$cancerType = ifelse(freq$cancerType=="HD", "0.Healthy", freq$cancerType)

#
keep_cols = c("patient","cancerType","dataset","type","freq")
pdat = freq[,keep_cols] 
#
p = ggboxplot(pdat, x="type", y="freq", color="type",outlier.shape=NA) +
    geom_jitter(aes(color=type),alpha=0.7, shape=16, width=0.2, size=3) +
    stat_compare_means(comparisons=list(c("healthy","diseased")), size=3, label.y.npc=0.5 ) +
    scale_colour_manual(values=get_palette("npg",2)) +
    xlab("") +
    ggtitle(sprintf("%s.Temra",stype)) +
    theme_classic2()
ggsave(sprintf("%s/health_and_deseased_blood.%s.freq.boxplot.pdf",oDir,stype), p, width=3.5, height=4)

}



