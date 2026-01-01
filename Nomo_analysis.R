#install.packages("rms")
#install.packages("rmda")


#引用包
library(rms)
library(rmda)

inputFile="normalize.txt"     #表达数据文件
geneFile="interGene.txt"      #基因列表文件
setwd("C:\\180geoImmune(FGR)\\22.Nomo")      #设置工作目录

#读取表达数据文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
#row.names(data)=gsub("-", "_", row.names(data))

#读取基因列表文件, 提取交集特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#获取样品分组信息(对照组和实验组)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)

#数据打包
ddist=datadist(rt)
options(datadist="ddist")

#构建模型，绘制列线图
lrmModel=lrm(Type~ ., data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.001,0.1,0.3,0.5,0.7,0.9,0.99),
	lp=F, funlabel="Risk of Disease")
#输出列线图
pdf(file="Nomo.pdf", width=8, height=6)
plot(nomo)
dev.off()

#绘制校准曲线
cali=calibrate(lrmModel, method="boot", B=1000)
pdf(file="Calibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()

