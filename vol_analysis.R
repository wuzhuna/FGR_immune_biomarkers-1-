#install.packages("ggplot2")
#install.packages("ggrepel")


#引用包
library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=0.585           #logFC过滤条件
P.Val.Filter=0.05       #矫正后的p值过滤条件
diffFile="all.txt"          #所有基因差异分析的结果文件
geneFile="hubGene.csv"      #基因列表文件
setwd("C:\\180geoImmune(FGR)\\12.vol")       #设置工作目录

#读取差异分析的结果文件
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
row.names(rt)=rt[,1]
#定义显著性
Sig=ifelse((rt$P.Value<P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#绘制火山图
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(P.Value)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "grey","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))

#在图形中标注基因的名称
geneRT=read.csv(geneFile, header=T, sep=",", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(rt))
showData=rt[sameGene,]
p1=p+geom_label_repel(data=showData,
                    box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                    size=2.5, aes(label=id)) + theme_bw()
#输出图形
pdf(file="vol.pdf", width=5.25, height=4.5)
print(p1)
dev.off()


