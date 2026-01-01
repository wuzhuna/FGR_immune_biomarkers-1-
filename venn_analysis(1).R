#install.packages("VennDiagram")


library(VennDiagram)      #引用包
setwd("C:\\180geoInImmune(FGR)\\17.venn")     #设置工作目录
geneList=list()

#读取LASSO回归的特征基因
rt=read.table("LASSO.gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对基因取unique
geneList[["LASSO"]]=uniqGene

#读取SVM机器学习的特征基因
rt=read.table("SVM-RFE.gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对基因取unique
geneList[["SVM-RFE"]]=uniqGene

#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集特征基因的名称
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGene.txt", sep="\t", quote=F, col.names=F, row.names=F)

