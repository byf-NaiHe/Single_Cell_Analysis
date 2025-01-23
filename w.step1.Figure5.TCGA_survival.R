#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library("plyr")
	library("dplyr")
	library("ggpubr")
	library("ggsci")
	library("reshape2")
    library("data.table")
	library("survival")
	library("survminer")
})
pdf.options(reset=TRUE, onefile=FALSE)

oDir="./OUT_Fig5"
dir.create(oDir, F, T)

dat = fread(sprintf("%s/../../data/metaInfo/TCGA.immuneTypes.txt.gz", oDir), header=T, stringsAsFactors=F, check.names=F, sep="\t") %>% as.data.frame()

###
dat = dat[!is.na(dat$OS.time),]
dat = dat[!is.na(dat$age_at_initial_pathologic_diagnosis),]
dat$age = dat$age_at_initial_pathologic_diagnosis
dat$cancerType = dat$cancer.type.abbreviation
dat$OS.time = round(dat$OS.time/30, 2)
dat = dat[!is.na(dat$ajcc_pathologic_tumor_stage),]
dat$stage = dat$ajcc_pathologic_tumor_stage
dat = dat[grepl("Stage I",dat$stage),]
dat$stage = gsub("[ABCD]$","",dat$stage)

### patch: 
### 1) merge COAD and READ as CRC; 2) remove typing biased cancers; 3) remove TGCT becaused of low death number
dat$cancerType = ifelse(dat$cancerType=="COAD"|dat$cancerType=="READ", "CRC", dat$cancerType)
tmp = dat[,c("cancerType","immuneType")] %>% table %>% as.data.frame
bad_cancers = tmp[tmp$Freq<5,"cancerType"] %>% unique %>% as.character
dat = dat[! dat$cancerType %in% c(bad_cancers, "TGCT"),]
###

dat$immuneType = factor(dat$immuneType)
cols = c("sample", "gender", "age", "cancerType", "OS", "OS.time", "immuneType", "stage")
dat = dat[,cols]


############ Cox regression
do.cox = function(df){
    cox=c()
	if(length(levels(df$gender)>1)){
		cox = coxph(Surv(OS.time,OS) ~ immuneType + gender + age + stage, data=df)
	}else{
		cox = coxph(Surv(OS.time,OS) ~ immuneType + age + stage, data=df)
	}
    cox.summary = summary(cox)
    nvar = length(unique(df$immuneType)) - 1
    #
    HR = round(cox.summary$conf.int[1:nvar,1], 2)
    HR.lower = round(cox.summary$conf.int[1:nvar,3], 2)
    HR.upper = round(cox.summary$conf.int[1:nvar,4], 2)
    HR.range = sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
    coef = cox.summary$coefficients[1:nvar,1]
    coef.se = cox.summary$coefficients[1:nvar,3]
    Pval = round(cox.summary$coefficients[1:nvar,5], 4)
    immuneType = gsub("immuneType","",rownames(cox.summary$conf.int)[1:nvar])
    return(data.frame(immuneType=immuneType, HR=HR, HR.range=HR.range, coef=coef, coef.se=coef.se, Pval=Pval))
}   

cdat = ddply(dat, c("cancerType"), do.cox)
cdat = cdat %>% group_by(immuneType) %>% mutate(adj.Pval=p.adjust(Pval, method="BH"))

### meta
d_ply(cdat, .(immuneType), function(.df){
    this_set = unique(.df$immuneType)
    #
    meta = meta::metagen(TE=.df$coef, seTE=.df$coef.se, studlab=.df$cancerType,
                         comb.fixed=F, comb.random=T, prediction=F, sm="HR")
    #
    pdf(sprintf("%s/TCGA.immuneType_%s.forrest.pdf",oDir,this_set), width=14, height=10)
    meta::forest(meta,  test.overall.random=T, digits.pval=4,
                 colgap.forest.left="5cm", zero.pval=T)
    dev.off()
}, .parallel=F)


############## KM curve
# plt.cancers = c("panCan", unique(dat$cancerType))
plt.cancers = c("KIRC","KIRP","LIHC","LUAD")
for (cancer in plt.cancers){
	pdat = dat
	if(cancer!="panCan"){
		pdat = dat[dat$cancerType==cancer,]
	}
	fit = survfit( Surv(OS.time,OS) ~ immuneType, data=pdat)
	#
	width = 5.5
	height = 6 + 0.1*length(unique(pdat$immuneType))
    pdf(sprintf("%s/TCGA.FigureS38.immuneType.KM.%s.pdf", oDir, cancer))
	p = ggsurvplot(fit, pdat, size=0.4, vlegend.labs=unique(pdat$immuneType),
				surv.median.line="none", pval=T, conf.int=F,
				risk.table=T, risk.table.y.text.col=T,
				palette=get_palette("jco", length(unique(pdat$immuneType))),
				legend.title="", ggtheme=theme_bw()) + 
			xlab("Months") + ggtitle(sprintf("data: %s", cancer))
    print(p)
	#ggsave(filename=sprintf("%s/TCGA.FigureS38.immuneType.KM.%s.pdf", oDir, cancer), p$plot, width=width, height=height, limitsize=F)
    dev.off()
}


