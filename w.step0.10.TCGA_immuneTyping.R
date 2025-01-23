suppressPackageStartupMessages({
        library("SingleCellExperiment")
        library("plyr")
        library("dplyr")
        library("ggpubr")
        library("data.table")
        library("scran")
        library("ranger")
})

oDir="./"


################
# 1. prepare inputs
### tcga data
sce.tcga = readRDS(sprintf("%s/../data/external/TCGA/TCGA_Toil_Recomp.rds",oDir))
sce.tcga = sce.tcga[,sce.tcga$sampleType=="Primary.Solid.Tumor" | sce.tcga$sampleType=="Primary.Blood.Derived.Cancer-Peripheral.Blood"]

### ref matrix
mtx.ref = fread(sprintf("%s/../data/expression/bulk/panC.bulkRNA.TPM.txt.gz",oDir), header=T, sep="\t", stringsAsFactors=F) %>% as.data.frame
rownames(mtx.ref) = mtx.ref[,2]
mtx.ref = mtx.ref[,grepl("\\.T$",colnames(mtx.ref))]
colnames(mtx.ref) = gsub("\\.T$", "", colnames(mtx.ref))
mtx.ref = as.matrix(mtx.ref)

## common genes
gene.common = intersect(rownames(sce.tcga), rownames(mtx.ref))   ## 18214 genes
mtx.ref = mtx.ref[gene.common,]
sce.tcga = sce.tcga[gene.common,]

mtx.tcga = assay(sce.tcga, "Toil_TPM") %>% as.matrix

### ref immuneType
immune.type = read.table(sprintf("%s/../data/metaInfo/panC.freq.all.sample.T.donorCluster.txt",oDir), header=T, stringsAsFactors=F, sep="\t", comment.char="")
immune.type$cancerType = as.character( gsub("^([^\\.]+).*", "\\1", immune.type$donorID) )
rownames(immune.type) = immune.type$donorID
immune.type = immune.type[grep("(thisStudy|CRC.LeiZhang2018|NSCLC.XinyiGuo2018|HCC.ChunhongZheng2017|STAD.BoxiKang2020)", immune.type$donorID), ]
immune.type$patient = immune.type$donorID
immune.type$patient = as.character( gsub("^[^\\.]+\\.[^\\.]+\\.(.*)", "\\1", immune.type$patient) )
rownames(immune.type) = immune.type$patient

immune.type$immuneType = paste("C", immune.type$displayCluster, sep="")
immune.type$immuneType = ifelse(immune.type$immuneType=="C1"|immune.type$immuneType=="C2", "C1C2", immune.type$immuneType)
immune.type$immuneType = ifelse(immune.type$immuneType!="C1C2", "C3-8", immune.type$immuneType)
immune.type = immune.type[,c("patient","immuneType")]

## common patients
patient.common = intersect(immune.type$patient, colnames(mtx.ref))
immune.type = immune.type[patient.common,]
mtx.ref = mtx.ref[,patient.common]

table(immune.type$immuneType)
#C1C2 C3-8
#  17   24


##############
# 2. run RF
rf.model = ranger::ranger(x=t(mtx.ref), y=factor(immune.type$immuneType), probability=T, seed=99, num.trees=1000, importance="permutation", local.importance=T, num.threads=8, classification=T) 
rf.pred = predict(rf.model, t(mtx.tcga), num.threads=8)
#
pred = c()
types = colnames(rf.pred$predictions)
for (i in 1:nrow(rf.pred$predictions)){
	tmp = rf.pred$predictions[i,]
	ind = which(tmp==max(tmp))
	pred = c(pred, types[ind])
}
names(pred) = colnames(mtx.tcga)
rf.pred$predictions.type = pred


##############
# 3. save
immuneType = rf.pred$predictions.type
prob = as.data.frame(rf.pred$predictions)
res = cbind(prob, immuneType, stringsAsFactors=F)
#
sce.tcga = sce.tcga[,rownames(res)]
info = cbind( colData(sce.tcga[,rownames(res)]), DataFrame(res) )
write.table(as.data.frame(info), file=sprintf("%s/../data/metaInfo/TCGA.immuneTypes.txt",oDir), sep="\t", quote=F, col.names=T, row.names=F)


## get importance of genes for each immuneType
## for each category
ref.pred = c()
ref.types = colnames(rf.model$predictions)
for (i in 1:nrow(rf.model$predictions)){
	tmp = rf.model$predictions[i,]
	ind = which(tmp==max(tmp))
	ref.pred = c(ref.pred, ref.types[ind])
}
#
importance.each = as.data.frame(matrix(NA, nrow=0, ncol=3))
for (i in unique(ref.pred)){
	tmp = rf.model$variable.importance.local[ ref.pred==i, ]
	tmp = colMeans(tmp)
	tmp = tmp[ order(tmp,decreasing=T) ]
	this = data.frame(immuneType=i, gene=names(tmp), importance=unname(tmp))
	importance.each = rbind(importance.each, this, stringsAsFactors=F)
}
importance.each = importance.each[importance.each$importance>0.0001,]
write.table(importance.each, file=sprintf("%s/../data/metaInfo/immuneTyping.RF_gene_importance.txt",oDir), sep="\t", quote=F, col.names=T, row.names=F)


