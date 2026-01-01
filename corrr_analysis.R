#install.packages("tidyverse")
#install.packages("corrr")


#引用包
library(limma)
library(corrr)
library(tidyverse)

cutoff=0.4                  #相关系数过滤条件
expFile="normalize.txt"     #表达数据文件
geneFile="hubGene.csv"      #基因列表文件
setwd("C:\\180geoImmune(FGR)\\13.corrr")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件,提取核心基因的表达量
geneRT=read.csv(geneFile, header=T, sep=",", check.names=F)
data=data[as.vector(geneRT[,1]),,drop=F]

#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
rt=data[,group=="Treat",drop=F]
rt=as.data.frame(t(rt))

#绘制相关性的图形
cor_df <- correlate(rt, method="spearman") %>% rearrange()
pdf(file="cor.pdf", width=7, height=6)
rplot(cor_df, colors = c( "skyblue1", "white", "indianred2"))
dev.off()

#绘制相关性的网络图
pdf(file="network.pdf", width=7, height=6)
network_plot(cor_df, min_cor = cutoff, colors = c( "skyblue1", "white", "indianred2"))
dev.off()

