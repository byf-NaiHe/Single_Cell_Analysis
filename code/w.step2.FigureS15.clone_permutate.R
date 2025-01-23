suppressPackageStartupMessages({
    library("reshape2")
    library("plyr")
    library("dplyr")
    library("data.table")
    library("ggpubr")
    library("ggsci")
    library("doParallel")
    library("ComplexHeatmap")
    library("circlize")
    library("MASS")
})
RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores=8)

oDir = "./OUT_FigS15/"
dir.create(oDir, F, T)
options(stringsAsFactors=F)

colSet = readRDS(sprintf("%s/../../data/metaInfo/panC.colSet.list.rds",oDir))
tcr = readRDS(sprintf("%s/../../data/tcr/byCell/tcr.zhangLab.comb.flt.rds",oDir))

### 0. tidy and create data
tcr = tcr[stype=="CD8",]
tcr$meta.cluster = tcr$meta.cluster %>% as.character %>% as.factor
tcr$cluster.name = tcr$cluster.name %>% as.character %>% as.factor

## get terminal Tex clones
tex.clones = ddply(tcr, .(cloneID), function(df){
    thisID = df$cloneID[1]
    thisNum = nrow(df)
    thisTb = df$meta.cluster %>% table
    if (thisNum>=3 && thisTb['CD8.c12.Tex.CXCL13']!=0){
        return( data.frame(cloneID=thisID))
    }
}, .parallel=T)
tex.clones = tex.clones$cloneID

## create dump data
dump.dat = tcr[,c("Cell_Name","cloneID","meta.cluster")] %>% as.data.frame
colnames(dump.dat)[3] = "Raw"
dump.dat$Raw = as.character(dump.dat$Raw)
#
set.seed(1)
n=1000
for (i in 1:n){
    this.id = sprintf("N%.6d",i)
    this.vec = sample(dump.dat$Raw, length(dump.dat$Raw), F)
    dump.dat = cbind(dump.dat, this.vec, stringsAsFactors=F)
    colnames(dump.dat)[ncol(dump.dat)] = this.id
}


### 1. calculate
dump.dat.slim = dump.dat[dump.dat$cloneID %in% tex.clones,]

processONE = function(.vec){
    P1.clus = c("CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK","CD8.c11.Tex.PDCD1")
    #P1.clus = c("CD8.c11.Tex.PDCD1")
    P2.clus = c("CD8.c10.Trm.ZNF683")
    end.clus = c("CD8.c12.Tex.CXCL13")

    .vec = factor(.vec, levels=unique(c(P1.clus, P2.clus, end.clus, .vec)))
    #
    tb = .vec %>% table %>% unclass
    this.in.P1 = tb[P1.clus] %>% sum
    this.in.P2 = tb[P2.clus] %>% sum
    this.each = tb[c(P1.clus,P2.clus,end.clus)]
    #
    value = ifelse( (this.in.P1+this.in.P2)>0, this.in.P1/(this.in.P1+this.in.P2), NA)
    res = structure(c(this.in.P1, this.in.P2, value), names=c("num.in.P1","num.in.P2","value") )
    res = c(this.each,res)
    return(res)
}

test.res = ddply(dump.dat.slim, .(cloneID), function(df){
    this.ID = df$cloneID[1]
    this.Size = nrow(df)
    raw.res = processONE(df$Raw)
    ## create distribution
    dump.values = c()
    for (i in 4:ncol(df)){
        this.value = processONE(df[,i])[['value']]
        if (!is.na(this.value)){
            dump.values = c(dump.values,this.value)
        }
    }
    fit = fitdistr(dump.values,"normal")
    ## two side P value (normal distribution)
    p = ifelse( raw.res[['value']] > fit$estimate[['mean']],
                2*pnorm(raw.res[['value']], mean=fit$estimate[['mean']], sd=fit$estimate[['sd']],  lower.tail=F, log.p=F),
                2*pnorm(raw.res[['value']], mean=fit$estimate[['mean']], sd=fit$estimate[['sd']],  lower.tail=T, log.p=F))

    ## return
    ret.df = as.data.frame(t(raw.res))
    ret.df$cloneID = this.ID
    ret.df$cloneSize = this.Size
    ret.df$perm.mean = fit$estimate[['mean']]
    ret.df$perm.sd = fit$estimate[['sd']]
    ret.df$p.value = p
    ret.df = ret.df %>% dplyr::select("cloneID", "cloneSize",  everything())
    return(ret.df)
},.parallel=T)

test.res$adj.p.value = p.adjust(test.res$p.value, method="BH")

saveRDS(test.res, file=sprintf("%s/clone.perm_test.rds",oDir))


### 2.plot (slim style 3)
test.res = test.res[!is.na(test.res$p.value) & test.res$adj.p.value<0.05,] 
test.res$category = ""
test.res$category = ifelse( test.res$value > test.res$perm.mean, "P1", test.res$category)
test.res$category = ifelse( test.res$value < test.res$perm.mean, "P2", test.res$category)


test.res$rank = ifelse(test.res$category=="P1",
                       test.res$num.in.P1/(test.res$num.in.P1+test.res$num.in.P2),
                       test.res$num.in.P2/(test.res$num.in.P1+test.res$num.in.P2))
test.res = arrange(test.res, category, desc(rank))

test.res$log.adj.p = -log10(test.res$adj.p.value)
test.res = arrange(test.res, category, desc(log.adj.p))
test.res$log.adj.p[test.res$log.adj.p>10] = 10

# mat
mat = test.res[,c("num.in.P1","num.in.P2")] %>% as.matrix
colnames(mat) = c("CD8.c05/06/11(GZMK+ T cells)","CD8.C10(ZNF683+CXCR6+ Trm)")
rownames(mat) = test.res$cloneID
mat = sweep(mat, 1, test.res$num.in.P1+test.res$num.in.P2, "/")

# row annotation
tcr = readRDS(sprintf("%s/../../data/tcr/byCell/tcr.zhangLab.comb.flt.rds",oDir))
tcr = tcr[,c("cloneID","cancerType")] %>% unique %>% as.data.frame %>% tibble::column_to_rownames("cloneID")

Category.cols = structure(c("#3cb371","#ffd700"), names=c("P1","P2") )
row_annot = rowAnnotation(
                            N.log10.adj.P = anno_lines( test.res$log.adj.p, smooth=F, height=unit(2, "cm"), axis_param=list(direction="reverse")),
                            Category = test.res$category,
                            cancerType = tcr[test.res$cloneID,"cancerType"],
                            col=list(Category=Category.cols,
                                     cancerType=colSet[['cancerType']])
                            )

# draw
set.seed(1)
ht = ComplexHeatmap::Heatmap(
                             mat,
                             name="Ratio",
                             #col = mat_color_fun,
                             col = RColorBrewer::brewer.pal(n=7, name="Reds"),

                             cluster_rows=F,
                             cluster_columns=F,
                             show_row_dend=F,
                             show_column_dend=F,
                             
                             show_row_names=F,
                             show_column_names=T,
                             left_annotation=row_annot,
                             
                             row_split=test.res$category,
                             row_gap = unit(2, "mm"),
                             border = T,

                             use_raster = T
        )

pdf(sprintf("%s/clone.perm_test.heatmap.pdf",oDir), width=3, height=8)
print(ht)
dev.off()




