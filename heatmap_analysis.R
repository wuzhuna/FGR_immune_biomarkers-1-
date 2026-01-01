#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)

expFile="normalize.txt"      #表达数据文件
geneFile="hubGene.csv"       #基因列表文件
setwd("C:\\180geoImmune(FGR)\\11.heatmap")     #设置工作目录

#读取表达数据文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取核心基因的表达量
geneRT=read.csv(geneFile, header=T, sep=",", check.names=F)
data=data[as.vector(geneRT[,1]),]

#提取样品的分组信息(对照组和实验组)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
names(Type)=colnames(data)
Type=as.data.frame(Type)

#绘制核心基因的热图
pdf(file="heatmap.pdf", width=6.5, height=4)
pheatmap(data, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()
