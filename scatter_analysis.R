#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#引用包
library(limma)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggExtra)

corFile=0.6                 #相关系数过滤条件
pFilter=0.001               #相关性检验pvalue的过滤条件
expFile="normalize.txt"     #表达数据文件
geneFile="hubGene.csv"      #核心基因的列表文件
setwd("C:\\180geoImmune(FGR)\\14.scatter")      #设置工作目录

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

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#相关性分析
outTab=data.frame()
#对Gene1进行循环
for(i in 1:(nrow(data)-1)){
	#获取Gene1的表达量
	x=as.numeric(data[i,])
	gene1=row.names(data)[i]
	#对Gene2进行循环，进行相关性检验
	for(j in (i+1):nrow(data)){
		gene2=row.names(data)[j]
	    #提取Gene2的表达量
	    y=as.numeric(data[j,])
		corT=cor.test(x, y, method = 'spearman')
		cor=corT$estimate
		pvalue=corT$p.value
		#对满足条件的基因对进行可视化
		if((abs(cor) > corFile) & (pvalue<pFilter)){
		    outTab=rbind(outTab, cbind(Gene1=gene1, Gene2=gene2, cor, pvalue))
			#可视化
			df1=as.data.frame(cbind(x,y))
			p1=ggplot(df1, aes(x, y)) + 
				xlab(paste0(gene1, " expression"))+ ylab(paste0(gene2, " expression"))+
				geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
				stat_cor(method = 'spearman', aes(x =x, y =y))
			p2=ggMarginal(p1, type="histogram", xparams=list(fill="#008B45FF"), yparams=list(fill="#631879FF"))
			pdf(file=paste0("cor.", gene1, "_", gene2, ".pdf"), width=5, height=4.6)
			print(p2)
			dev.off()
		}
	}
}

#输出相关性分析的表格
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)
