# !Note: using R4.0 and Seurat 4.0 here

################################################################

suppressPackageStartupMessages({
library("reshape2")
library("stringr")
library("plyr")
library("dplyr")
library("ggpubr")
library("ggsci")
#library("ggrastr")
library("cowplot")
library("patchwork")
library("data.table")
library("R.utils")
library("Matrix")
library("SingleCellExperiment")
library("Seurat")
library("sctransform")
library("harmony")
})
#
stype = "CD8"
oDir = "./OUT_FigS17/"
dir.create(oDir, F, T)


# 0. read data
sce_all = readRDS(sprintf("%s/../../data/expression/CD8/integration/int.CD8.S35.sce.merged.rds",oDir))
sce_ref = readRDS(sprintf("%s/../../data/expression/CD8/integration/int.CD8.S35.HVG.continumOnly.v1.sce.Path.rds",oDir))
#all(rownames(sce_all) == rownames(sce_ref))
gene_used = read.table(sprintf("%s/../../data/metaInfo/int.CD8_Tex.genes.txt",oDir), header=F, stringsAsFactors=F)$V1
colSet = readRDS(sprintf("%s/../../data/metaInfo/panC.colSet.list.rds",oDir))


## 1. tidy inputs
### ref
mat_ref = assay(sce_ref, "exprs")
meta_ref = as.data.frame(colData(sce_ref))
meta_ref = meta_ref[,!grepl("^RNA_",colnames(meta_ref))]
seu_ref = CreateSeuratObject(mat_ref, project="ref", meta.data=meta_ref)
seu_ref = SetAssayData(seu_ref, "scale.data", mat_ref)

seu_ref = RunPCA(seu_ref, features=gene_used, npcs=15, verbose=F)
seu_ref = RunUMAP(seu_ref, reduction="pca",dims=1:15, return.model=T)
seu_ref = RunHarmony(seu_ref, c("dataset"))
seu_ref = RunUMAP(seu_ref, reduction="harmony", reduction.name="harmony.umap", dims=1:15, return.model=T)


### query
sce_qry = sce_all[,as.character(sce_all$meta.cluster)=="CD8.c15.ISG.IFIT1"]
mat_qry = assay(sce_qry, "exprs")
meta_qry = as.data.frame(colData(sce_qry))
meta_qry = meta_qry[,!grepl("^RNA_",colnames(meta_qry))]
seu_qry = CreateSeuratObject(mat_qry, project="ref", meta.data=meta_qry)
seu_qry = SetAssayData(seu_qry, "scale.data", mat_qry)


## 2. projection
anchors = FindTransferAnchors(reference=seu_ref, query=seu_qry, reduction="pcaproject", reference.reduction="pca", dims=1:15, features=gene_used) 

seu_qry = TransferData(anchorset=anchors, reference=seu_ref, query=seu_qry, refdata="meta.cluster", dims=1:15)
seu_qry = IntegrateEmbeddings(anchorset=anchors, reference=seu_ref, query=seu_qry, new.reduction.name="ref.pca", dims.to.integrate=1:15)
seu_qry = ProjectUMAP(query=seu_qry, query.reduction="ref.pca", query.dims=1:15, reference=seu_ref, reference.reduction="pca", reference.dims=1:15, reduction.model="harmony.umap")

save(seu_qry, meta_qry, anchors, colSet, file=sprintf("%s/%s_ISG.mapping.Rdata", oDir, stype))
#load(sprintf("%s/%s_ISG.mapping.Rdata", oDir, stype))


## 3.plot
### tidy
nam.conv = fread(sprintf("%s/../../data/metaInfo/name.conversion.txt",oDir),sep="\t",stringsAsFactors=F,header=T)
nam.conv = as.data.frame(nam.conv)
rownames(nam.conv) = nam.conv$meta.cluster
seu_qry$predicted.cluster.name = nam.conv[as.character(seu_qry$predicted.id),"cluster.name"]
#
seu_qry$predicted.cluster.name = ifelse(seu_qry$predicted.id.score > 0.5, as.character(seu_qry$predicted.cluster.name), "belowThres")
seu_qry_flt = seu_qry[,seu_qry$predicted.cluster.name!="belowThres"]


### stat
pdat = as.data.frame(table(seu_qry$predicted.cluster.name))
colnames(pdat)[1] = "predicted.cluster.name"
pdat$Percent = pdat$Freq / sum(pdat$Freq) * 100
pdat = arrange(pdat, desc(pdat$Freq))

orders = c(as.character(pdat$predicted.cluster.name)[as.character(pdat$predicted.cluster.name)!="belowThres"], "belowThres")
pdat$predicted.cluster.name = factor(as.character(pdat$predicted.cluster.name), levels=orders)

pdat$group = ifelse(as.character(pdat$predicted.cluster.name)=="belowThres", "bad", "good")
colors = structure(c("#D3D3D3","#459F7E"), names=c("bad", "good"))

p = ggplot(pdat, aes(x=predicted.cluster.name, y=Percent, fill=group)) +
                        geom_bar(stat="identity", position="stack", width=0.8) +
                        theme_classic() +
                        theme(axis.text.x=element_text(size=8, vjust=1.0, hjust=1.1, angle=60),
                                        plot.margin=margin(t=1, r=2, b=1, l=2, unit="cm"),
                                        legend.position="none",
                                        legend.text=element_text(size=8)) +
						scale_fill_manual(values=colors) +
						xlab("Predicted labels for ISG+ cells")
ggsave(p, file=sprintf("%s/%s_ISG.mapping.stat.pdf", oDir, stype), height=3.5, width=5)


### umap with ref
tmp1 = cbind(data.frame(Embeddings(seu_qry_flt, "ref.umap")), seu_qry_flt$predicted.cluster.name)
tmp2 = cbind(data.frame(reducedDim(sce_ref, "recal.harmony.umap")), "ref cells")
colnames(tmp1) = c("ref_UMAP1","ref_UMAP2","group")
colnames(tmp2) = c("ref_UMAP1","ref_UMAP2","group")
dat = rbind(tmp1, tmp2)
dat$group = factor(as.character(dat$group), levels=c("Tn", "IL7R+ Tm", "GZMK+ early Tem", "GZMK+ Tem", "GZMK+ Tex",
									   "ZNF683+CXCR6+ Trm", "terminal Tex","ref cells"))

colors = c(colSet$cluster.name.short.CD8, "lightgrey")
names(colors)[length(colors)] = "ref cells"

p = ggplot(dat[as.character(dat$group)!="ref cells",], aes(x=ref_UMAP1, y=ref_UMAP2, color=group)) + 
		geom_point(dat=dat[as.character(dat$group)=="ref cells",],size=0.7, shape=16) +
		geom_point(size=0.7, shape=16) +
		theme_classic2() +
		scale_color_manual(values=colors[unique(as.character(dat$group))]) +
		guides(color=guide_legend(override.aes=list(size=3))) +
		theme(axis.text=element_text(size=10), axis.title=element_text(size=12),plot.title=element_text(size=14)) +
		coord_fixed(ratio=1, xlim=c(-10, 10),ylim = c(-10, 10),expand=F) +
		ggtitle("Predicted labels for ISG+ cells")
print(p)
ggsave(p, file=sprintf("%s/%s_ISG.mapping.umap.pdf",oDir, stype), width=5.5, height=8, useDingbats=F)




