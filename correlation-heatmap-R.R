data3<-read.csv("ESCC-gene-list-exp.csv",sep=",",header=T,row.names=1)
data3[1:6,1:5]
M<-cor(data3)
library(corrplot)
col <- colorRampPalette(c("darkblue", "white", "darkred"))(201) # 生成201种颜色的渐变  

pdf("TCGA-相关性热图.pdf")
corrplot(M,method="color",type = "full",col=col)
dev.off()



gse<-read.csv("GSE53625.csv",header=T,sep=",")
dim(gse)
gse[1:6,1:6]
genelist2<-read.csv("genelist.csv",header=T,sep=",")
head(genelist2)
dim(genelist2)
genelist22<-as.matrix(genelist2[,1:2])
colnames(genelist22)<-c("gene","ID_REF")
datay<-merge(genelist22,gse,by="ID_REF")
dim(datay)
datay[1:5,1:5]
data2y<-datay[,-1]
data2y[1:5,1:5]
data22y<-as.data.frame(t(data2y))
data22y[1:5,1:5]
write.table(data22y,"GSE53625.list2.csv",sep=",",quote=F,row.names=T,col.names=F)
data3y<-read.csv("GSE53625.list2.csv",sep=",",header=T,row.names=1)
data3yy<-data3y[seq(1,nrow(data3y),2),]
library(corrplot)
My<-cor(data3yy)
col <- colorRampPalette(c("darkblue", "white", "darkred"))(201) 
pdf("gse53625-相关性.pdf")
corrplot(My,method="color",type = "full",col=col)
dev.off()


