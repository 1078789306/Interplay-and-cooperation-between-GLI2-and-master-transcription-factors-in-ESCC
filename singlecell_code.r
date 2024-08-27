library(Seurat)

library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(readxl)

library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(Seurat) 
library(msigdbr) 
library(GSVA) 
library(tidyverse) 
library(clusterProfiler) 
library(patchwork) 
library(limma) 
library(readxl)

#####基因表达热土
escc_P6<-readRDS("D:\\ESCC\\GSE188900\\P6患者单独合并\\0data\\ESCC_P6_integrate合并_注释后.rds")
escc_P6@meta.data$sample[escc_P6$orig.ident=="GSM5691647"|escc_P6$orig.ident=="GSM5691648"]<-"Tumor"
escc_P6@meta.data$sample[escc_P6$orig.ident=="GSM5691649"]<-"Normal"
normal<-subset(escc_P6,subset = sample=="Normal")
GSE188900<-readRDS("D:\\ESCC\\GSE188900\\tumor\\escc_scissor_tumor0.005.rds")
GSE188900<-merge(GSE188900,normal)
GSE188900_metadata<-GSE188900@meta.data
GSE188900_metadata <- GSE188900_metadata %>% 
  mutate(sample = ifelse(is.na(sample), "tumor", sample))
GSE188900<-AddMetaData(GSE188900,GSE188900_metadata)

Idents(GSE188900)<-"SCT"
sq_epi<-subset(GSE188900,subset = c(main_cell_type=="Squamous epithelium"|main_cell_type=="Glandular epithelium"))

#计算平均表达量
gene_cell_exp <- AverageExpression(sq_epi,
                                   features = TF_26$V1,
                                   group.by = 'sample',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('tumor'="#9ECABE",
                                                  'Normal'="#F6F5B4")))#颜色设置
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("#98F5FF","white","#FF4040"))(5),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)






GSE188900<-readRDS("D:\\ESCC\\GSE188900\\tumor\\escc_scissor_tumor0.005.rds")

####全部的correlation
sq_epi<-subset(GSE188900,subset = c(main_cell_type=="Epithelial cells"))


sq_epi_data<-as.data.frame(sq_epi@assays$RNA@data)
sq_epi_data$gene<-rownames(sq_epi_data)
TF_6<-read.table("D:\\ESCC_LIU\\List1-26-TFs.txt")

colnames(TF_6)<-"gene"
sq_epi_data<-merge(sq_epi_data,TF_6,by="gene")
rownames(sq_epi_data)<-sq_epi_data[,1]
sq_epi_data<-sq_epi_data[,-1]

corr<-round(cor(t(sq_epi_data)),3)
library(ggcorrplot)
library(corrplot)
p.mat<-round(cor_pmat(t(sq_epi_data)),3)
col2 <- colorRampPalette( c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF","#FDDBC7", "#F4A582", "#D6604D", "#B2182B","#67001F"))(100)
# 绘制相关性热图




corr<-corr[c(22,14,3,1:2,4:13,15:21,23:26),]

corr<-corr[,c(22,14,3,1:2,4:13,15:21,23:26)]
corrplot(corr, method = "square", col = col2)

remove(sq_epi)
gc()

TF_6<-read.table("D:\\ESCC_LIU\\Gene list 2.txt")

colnames(TF_6)<-"gene"
sq_epi_data<-merge(sq_epi_data,TF_6,by="gene")
rownames(sq_epi_data)<-sq_epi_data[,1]
sq_epi_data<-sq_epi_data[,-1]

corr<-round(cor(t(sq_epi_data)),3)



corr<-corr[c(6,1,3,2,4,5,7),]

corr<-corr[,c(6,1,3,2,4,5,7)]
corrplot(corr, method = "square", col = col2)


library(ggstatsplot)

exprSet<-as.data.frame(t(as.data.frame(sq_epi@assays$RNA@data)))
ggscatterstats(data = exprSet,
               y = TP63,
               x = RUNX1,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram
               title = "Relationship between TP63 and RUNX1")




#####只有恶性细胞的correlation
GSE188900<-readRDS("D:\\ESCC\\GSE188900\\tumor\\escc_scissor_tumor0.005.rds")





save(GSE188900,file="F:\\GSE188900.RData")
table(GSE188900$main_cell_type)
DimPlot(GSE188900, reduction = 'tsne', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
Idents(GSE188900)<-GSE188900$scissor
VlnPlot(GSE188900,features = TF_6$gene)
table(GSE188900$main_cell_type)
sq_epi<-subset(GSE188900,subset = c(main_cell_type=="Squamous epithelium"|main_cell_type=="Glandular epithelium"))


coptkat<-read.table("D:\\ESCC_LIU\\test_copykat_prediction.txt",sep="\t",header = T)
metadata<-sq_epi@meta.data
colnames(coptkat)<-c("cell","copykat")
metadata$cell<-rownames(metadata)
metadata<-merge(metadata,coptkat,by="cell")
rownames(metadata)<-metadata$cell
sq_epi<-AddMetaData(sq_epi,metadata)

Idents(sq_epi)<-sq_epi$copykat
DimPlot(sq_epi,reduction = "tsne")

table(sq_epi$copykat)
sq_epi_ex<-subset(sq_epi,subset=copykat=="aneuploid")
sq_epi_data<-as.data.frame(sq_epi_ex@assays$RNA@data)
sq_epi_data$gene<-rownames(sq_epi_data)

TF_6<-read.table("D:\\ESCC_LIU\\List1-26-Tumor2.txt")

colnames(TF_6)<-"gene"
sq_epi_data<-merge(sq_epi_data,TF_6,by="gene")
rownames(sq_epi_data)<-sq_epi_data[,1]
sq_epi_data<-sq_epi_data[,-1]

corr<-round(cor(t(sq_epi_data)),3)
library(ggcorrplot)
library(corrplot)
p.mat<-round(cor_pmat(t(sq_epi_data)),3)
col2 <- colorRampPalette( c("#053061","#2166AC","#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF","#FDDBC7", "#F4A582", "#D6604D", "#B2182B","#67001F"))(100)
# 绘制相关性热图




corr<-corr[c(20,11,5,15,7,18,1:4,6,8:10,12:14,16:19,21:24),]

corr<-corr[,c(20,11,5,15,7,18,1:4,6,8:10,12:14,16:19,21:24)]
corrplot(corr, method = "square", col = col2)

remove(sq_epi)
gc()
sq_epi_data<-as.data.frame(sq_epi_ex@assays$RNA@data)
sq_epi_data$gene<-rownames(sq_epi_data)

TF_6<-read.table("D:\\ESCC_LIU\\Gene list 2.txt")

colnames(TF_6)<-"gene"
sq_epi_data<-merge(sq_epi_data,TF_6,by="gene")
rownames(sq_epi_data)<-sq_epi_data[,1]
sq_epi_data<-sq_epi_data[,-1]

corr<-round(cor(t(sq_epi_data)),3)



corr<-corr[c(6,1,3,2,4,5,7),]

corr<-corr[,c(6,1,3,2,4,5,7)]
corrplot(corr, method = "square", col = col2)



library(ggstatsplot)

exprSet<-as.data.frame(t(as.data.frame(sq_epi_ex@assays$RNA@data)))
ggscatterstats(data = exprSet,
               y = TP63,
               x = GLI2,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram", # #类型可以换成density,boxplot,violin,densigram
               title = "Relationship between TP63 and GLI2")






library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
saveRDS(moncel2,"D:\\ESCC\\GSE188900\\tumor\\sq_epi\\拟时分析.rds")
moncle<-sq_epi
table(moncle$copykat)
moncle<-subset(moncle,subset=c(copykat=="diploid"|copykat=="aneuploid"))
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
celltype<-moncle$copykat
celltype<-as.data.frame(celltype)
celltype$cellname<-rownames(celltype)

celltype$sample<-substr(celltype$cellname,1,10)
celltype2<-unite(celltype,"cellnamenew",c("celltype","sample"), sep="-", remove = F)
moncel2<-moncle
moncel2$cellname<-celltype2$cellnamenew


expr_matrix <- as(as.matrix(moncel2@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 

p_data <- moncel2@meta.data
#p_data$celltype <- pbmc@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
f_data <- data.frame(gene_short_name = row.names(moncel2),row.names = row.names(moncel2))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
#构建CDS对象(一定要确保identical(rownames(fd),rownames(expr_matrix))为TURE
gc()

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
rownames(fd)
rownames(expr_matrix)
expr_matrix<- expr_matrix[rownames(fd), ]



identical(rownames(fd),rownames(expr_matrix))

remove(moncle)

#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

remove(expr_matrix)
library(data.table)
gc()
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
gc()
dim(cds)
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。
#这一步输入的expressed_genes来自于步骤4。
#️️后续分析使用的是该方法
#也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(pbmc) 
####diff这一流程里面的fullModelFormulaStr后输入的应该是你自己的细胞分类或者注释结果

###########选择轨迹基因
##使用seurat选择的高变基因⚠️
escc_sq<-FindVariableFeatures(escc_sq)
express_genes <- VariableFeatures(escc_sq)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

###########选择轨迹基因
##使用seurat选择的高变基因⚠️
escc_sq<-FindVariableFeatures(escc_sq)
express_genes <- VariableFeatures(escc_sq)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)


##使用clusters差异表达基因
Idents(moncel2)<-moncel2$copykat
DefaultAssay(moncel2)
deg.cluster <- FindAllMarkers(moncel2)
express_genes <- subset(deg.cluster,p_val<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)


diff <- differentialGeneTest(cds[disp.genes,],fullModelFormulaStr="~copykat",cores=1) 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
write.table(deg,file="train.monocle.DEG.xlsx",col.names=T,row.names=F,sep="\t",quote=F)
DimPlot(moncel, reduction = "umap",group.by = "integrated_snn_res.0.1")
## 轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

gc()
cds <- orderCells(cds)
#️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
cds <- orderCells(cds, root_state = 6) #把State5设成拟时间轴的起始点

plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 

plot_cell_trajectory(cds,color_by="main_cell_type", size=1,show_backbone=TRUE)

plot_cell_trajectory(cds, color_by = "copykat",size=1,show_backbone=TRUE)
cds$orig.ident
cds$SCT_snn_res.0.1
plot_cell_trajectory(cds, color_by = "SCT_snn_res.0.1",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)


exprData <- cds@assayData$exprs
exprData <- LogNormalize(exprData) #对数据normaliza一下，会比较显著，记得引用seurat包
cds$GLI2 <- exprData["GLI2",] #输入想看的基因，cds是monocle的对象
plot_cell_trajectory(cds, color_by = "TP63",cell_size=1,show_backbone=TRUE)+scale_colour_gradient(low="grey",high = "red")

Idents(moncel2)<-moncel2$copykat
DefaultAssay(moncel2)<-"SCT"
library(ggpubr)
a<-VlnPlot(moncel2,features = c("TP63"))+
NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif") +
  ylim(0,2.5)

b<-VlnPlot(moncel2,features = c("GLI2"))+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif") +
  ylim(0,1.5)

c<-VlnPlot(moncel2,features = c("RUNX1"))+
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif") +
  ylim(0,2)
cowplot::plot_grid(a,b,c, ncol = 3)
moncel2_count<-as.data.frame(t(as.data.frame(moncel2@assays$SCT@counts)))
moncel2_metadata<-moncel2@meta.data
TP63<-as.data.frame(moncel2_count$TP63)
TP63$cell<-rownames(moncel2_count)
colnames(TP63)[1]<-"TP63_expression"
moncel2_metadata$cell<-rownames(moncel2_metadata)
metadata<-merge(TP63,moncel2_metadata,by="cell")

library(ggsignif)
table(metadata$copykat)

comparisons<-list(c("aneuploid","diploid"))

p<-ggplot(metadata, aes(x = reorder(copykat,TP63_expression),TP63_expression)) # x分组变量，y表达变量

p+geom_violin() #画出violin plot
p+geom_violin(aes(fill = copykat))+#geom_boxplot(width = 0.1)+
  geom_signif(comparisons = comparisons,map_signif_level = T, textsize = 6, test = t.test, step_increase = 0.2) +
  guides(fill = "none") + xlab(NULL) + theme_classic()+ggtitle("MAR_TP63")+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.85))#按组别填充颜色

