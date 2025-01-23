suppressPackageStartupMessages({
	library("SingleCellExperiment")
	library("GSVA")
	library("parallel")
	library("plyr")
	library("dplyr")
	library("data.table")
})

oDir="./"
dir.create(sprintf("%s/",oDir), F, T)

### bulk data
sce.tcga = readRDS(sprintf("%s/../data/external/TCGA/TCGA_Toil_Recomp.rds",oDir))
sce.tcga = sce.tcga[,sce.tcga$sampleType=="Primary.Solid.Tumor" | sce.tcga$sampleType=="Primary.Blood.Derived.Cancer-Peripheral.Blood"]
sce.tcga$cancer.type.abbreviation = as.character(sce.tcga$cancer.type.abbreviation)

# signatures
sigs = fread(sprintf("%s/../data/metaInfo/signature_genes.txt.gz", oDir), header=F, stringsAsFactors=F, sep="\t")
colnames(sigs) = c("cluster", "gene")

inpList = list()
for (i in unique(sigs$cluster)){
	tmp = sigs[sigs$cluster == i,]
	inpList[[i]] = tmp$gene
}

### run
gsva.res = matrix(nrow=0, ncol=length(unique(sigs$cluster))) 
for (cancer in unique(sce.tcga$cancer.type.abbreviation)){
	sce.sub = sce.tcga[,sce.tcga$cancer.type.abbreviation==cancer]
	mat = as.matrix(assay(sce.sub,"Toil_TPM"))
	#
	tmp.res = gsva(mat, inpList, method="gsva", kcdf="Gaussian", mx.diff=T, parallel.sz=8)
	tmp.res = t(tmp.res)
	#
	gsva.res = rbind(gsva.res, tmp.res)
}

### save
info = cbind( colData(sce.tcga[,rownames(gsva.res)]), DataFrame(gsva.res) )
write.table(as.data.frame(info), file=sprintf("%s/../data/metaInfo/TCGA.sig_scores.txt",oDir), sep="\t", quote=F, col.names=T, row.names=F)




