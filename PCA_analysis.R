#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)

expFile="normalize.txt"       #表达数据文件
geneFile="interGene.txt"      #基因列表文件
setwd("C:\\180geoImmune(FGR)\\19.PCA")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取交集特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

#获取样品的分组信息(对照组和实验组)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

#PCA分析
data=as.matrix(data)
data.class=rownames(data)
data.pca=prcomp(data)
pcaPredict=predict(data.pca)

#定义数据框
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], group=group)

#定义椭圆的函数
PCA.mean=aggregate(PCA[,1:2], list(group=PCA$group), mean)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$group))){
	df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
                  veganCovEllipse(cov.wt(cbind(PC1,PC2),
                  wt=rep(1/length(PC1),length(PC1)))$cov,
                  center=c(mean(PC1),mean(PC2))))),group=g))
}

#绘制PCA的图形
pdf(file="PCA.pdf", width=6, height=4.6)
ggplot(data=PCA, aes(PC1, PC2)) + geom_point(aes(color=group)) +
	scale_colour_manual(name="Type",  values=c("blue", "red"))+
    geom_path(data=df_ell, aes(x=PC1, y=PC2,colour=group), size=1.2, linetype=2)+
    annotate("text", x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$group)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

