suppressPackageStartupMessages({
library("reshape2")
library("stringr")
library("plyr")
library("dplyr")
library("ggpubr")
library("ggsci")
library("ggrastr")
library("data.table")
library("sscClust")
library("Seurat")

})
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores=8)
options(stringsAsFactors=F)

oDir = "./OUT_FigS22"
dir.create(oDir, F, T)
colSet = readRDS(sprintf("%s/../../data/metaInfo/panC.colSet.list.rds",oDir))


### 1. find cells in Tex cluster with MAIT TCR
tcr = readRDS(sprintf("%s/../../data/tcr/byCell/tcr.zhangLab.comb.flt.rds",oDir))

# get CD4/CD8 size, CDR3 seq
tcr = ddply(tcr, .(cloneID), function(df){
	df$cloneID = gsub(":.*$","",as.character(df$cloneID))
	df$cloneSize = nrow(df)
	df$cloneSize.CD8 = nrow(df[grepl("^CD8",df$meta.cluster),])
	# patch: get the most frequent CDR3, because some have same freq of A1A2/B1G2
	CDR3 = table(sprintf("%s+%s", df$CDR3.A1, df$CDR3.B1))
	A1 = table( df$Identifier.A1 )
	B1 = table( df$Identifier.B1 )
	##
	CDR3 = names(CDR3[which.max(CDR3)])
	A1 = names( A1[which.max(A1)] )
	B1 = names( B1[which.max(B1)] )
        ##
	A1V = gsub("^.*(TRAV[0-9\\-]+).*$", "\\1", A1)
	A1J = gsub("^.*(TRAJ[0-9\\-]+).*$", "\\1", A1)
	B1V = gsub("^.*(TRBV[0-9\\-]+).*$", "\\1", B1)
	B1J = gsub("^.*(TRBJ[0-9\\-]+).*$", "\\1", B1)
        ##
	df$CDR3 = CDR3
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


# get cells
tar.cells = tcr[grepl("TRAV1-2",tcr$A1)  &
	  (grepl("TRAJ33", tcr$A1) |
	   grepl("TRAJ20", tcr$A1) |
	   grepl("TRAJ12", tcr$A1)) &
	   as.character(tcr$meta.cluster)=="CD8.c12.Tex.CXCL13" &
	   as.character(tcr$loc)=="T"
	   , ]

tar.cells = ddply(tar.cells, .(cloneID), function(df){
    df$Size2 = nrow(df)
    return(df)
}, .parallel=T)

tar.cells = tar.cells[tar.cells$Size2>1 & as.character(tar.cells$meta.cluster)=="CD8.c12.Tex.CXCL13",]
tar.cells = tar.cells[, c("Cell_Name","cloneID")]

write.table(tar.cells, file=sprintf("%s/MAIT_Tex.cellID.list",oDir), sep="\t", quote=F, col.names=T, row.names=F)


### 2. construct sc merged sce obj
files = list.files(sprintf("%s/../../data/expression/CD8/byDataset/",oDir), full.names=T)
files = files[grepl("zhangLab5P|zhangLabSS2|STAD.BoxiKang2019",files)]
dataset.id = gsub("^.*/", "", files)
dataset.id = gsub(".sce.rds$", "", dataset.id)
dataset.id = gsub(".mod$", "", dataset.id)

keep.info = c("cellID","patient","cancerType","loc","batchV","dataset","cellID.uniq","meta.cluster") 
used.genes = c("CXCL13","PDCD1","LAYN","HAVCR2","ENTPD1", "CTLA4",  # 6
               "SLC4A10", "KLRB1", "ZBTB16", "NCR3", "RORC","RORA", # 6
               "CCR6", "TMIGD2",  "IL23R", "IL17A", "IL17F", "IL26",  "IL18RAP", "IL17RE")  # 8

seu.list = list()
for (i in 1:length(files)){
    sce = readRDS(files[i])
    rownames(sce) = rowData(sce)$seu.id
    meta = as.data.frame(colData(sce))[,keep.info]
    seu = CreateSeuratObject(counts=assay(sce,"norm_exprs")[used.genes,], meta.data=meta, min.cells=0, min.features=0, project=dataset.id[i])
    seu.list = c(seu.list, seu)
}

seu.merge = Seurat:::merge.Seurat(x=seu.list[[1]], y=seu.list[2:length(seu.list)], merge.data=T)
seu.merge =  ScaleData(seu.merge, features=rownames(seu), vars.to.regress=c("batchV"))

sce.merge = ssc.build(GetAssayData(seu.merge,slot="data"), assay="norm_exprs")
assay(sce.merge, "norm_exprs.scaled") = GetAssayData(seu.merge, slot="scale.data")
colData(sce.merge) = DataFrame( seu.merge@meta.data[,4:ncol(seu.merge@meta.data)] )
rowData(sce.merge)$display.name = rownames(sce.merge)

saveRDS(sce.merge, file=sprintf("%s/MAIT_Tex.expr.sce.rds",oDir))



### 3. plot
source(sprintf("%s/../source_FigS22_func.R",oDir))

sce.t = sce.merge[used.genes, tar.cells$Cell_Name]
rowData(sce.t)['geneGroup'] = c(rep("Tex signature",6), rep("Tc17 signature",14))
sce.t$cloneID = tar.cells$cloneID
colSet[['cloneID']] = structure( get_palette("Set2",length(unique(sce.t$cloneID))), names=unique(sce.t$cloneID) )
colSet[['geneGroup']] = structure( c("dimgray","lightgrey"), names=unique(rowData(sce.t)$geneGroup) )

heatmap_sm.2(sce.t, assay.name="norm_exprs.scaled", colSet=colSet,
                 out.prefix=sprintf("%s/MAIT_Tex.expr.heatmap", oDir),
                 columns=c("cloneID"),  rows=c("geneGroup"),
                 columns.order=c("cloneID"),
                 show_column_names=F, show_row_names=T,
                 row_names_side="right",
                 row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
                 border=T, ann.bar.height=0.5, palette.name="RdBu",
                 pdf.width=14, pdf.height=9, do.scale=F,
                 z.lo=-0.5, z.hi=0.5, z.step=0.1, 
                 do.clustering.row=T, do.clustering.col=F,
                 dend.row=T, dend.col=F,
                 row_dend_width = unit(2, "cm"),
                 column_dend_height = unit(2, "cm")
            )










