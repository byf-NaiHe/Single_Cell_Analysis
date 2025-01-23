rm(list=ls())
# 设置镜像：
options()$repos
options()$BioC_mirror
#options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos
options()$BioC_mirror

# 方法一：
options()$repos
install.packages('WGCNA')
install.packages(c("FactoMineR","factoextra"))
install.packages(c("ggplot2","pheatmap","ggpubr"))
library("FactoMineR")
library("factoextra")