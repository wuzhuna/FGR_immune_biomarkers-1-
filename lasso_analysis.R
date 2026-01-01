#install.packages("glmnet")


#引用包
set.seed(12345)
library(limma)
library(glmnet)

expFile="normalize.txt"     #表达数据文件
geneFile="hubGene.csv"      #核心基因的列表文件
setwd("C:\\180geoImmune(FGR)\\15.lasso")      #设置工作目录

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
rt=t(data)

#构建模型
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
#绘制Lasso回归的图形
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit)
dev.off()
#绘制交叉验证的图形
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf", width=6, height=5.5)
plot(cvfit)
dev.off()

#输出疾病的特征基因
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
