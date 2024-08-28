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

gene_cell_exp <- AverageExpression(sq_epi,
                                   features = TF_26$V1,
                                   group.by = 'sample',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)


library(ComplexHeatmap)

df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('tumor'="#9ECABE",
                                                  'Normal'="#F6F5B4")))

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
               marginal.type = "densigram", 
               title = "Relationship between TP63 and RUNX1")





GSE188900<-readRDS("D:\\ESCC\\GSE188900\\tumor\\escc_scissor_tumor0.005.rds")





save(GSE188900,file="F:\\GSE188900.RData")
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
               marginal.type = "densigram", 
               title = "Relationship between TP63 and GLI2")






library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
moncle<-sq_epi
table(moncle$copykat)
moncle<-subset(moncle,subset=c(copykat=="diploid"|copykat=="aneuploid"))

celltype<-moncle$copykat
celltype<-as.data.frame(celltype)
celltype$cellname<-rownames(celltype)

celltype$sample<-substr(celltype$cellname,1,10)
celltype2<-unite(celltype,"cellnamenew",c("celltype","sample"), sep="-", remove = F)
moncel2<-moncle
moncel2$cellname<-celltype2$cellnamenew


expr_matrix <- as(as.matrix(moncel2@assays$RNA@counts), 'sparseMatrix')


p_data <- moncel2@meta.data

f_data <- data.frame(gene_short_name = row.names(moncel2),row.names = row.names(moncel2))

gc()

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
rownames(fd)
rownames(expr_matrix)
expr_matrix<- expr_matrix[rownames(fd), ]



identical(rownames(fd),rownames(expr_matrix))

remove(moncle)


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
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。



disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)




diff <- differentialGeneTest(cds[disp.genes,],fullModelFormulaStr="~copykat",cores=1) 

head(diff)


deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)


write.table(deg,file="train.monocle.DEG.xlsx",col.names=T,row.names=F,sep="\t",quote=F)
DimPlot(moncel, reduction = "umap",group.by = "integrated_snn_res.0.1")

ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  

plot_ordering_genes(cds)
dev.off()

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

gc()
cds <- orderCells(cds)

cds <- orderCells(cds, root_state = 6) #把State5设成拟时间轴的起始点

plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 

plot_cell_trajectory(cds,color_by="main_cell_type", size=1,show_backbone=TRUE)

plot_cell_trajectory(cds, color_by = "copykat",size=1,show_backbone=TRUE)
cds$orig.ident
cds$SCT_snn_res.0.1
plot_cell_trajectory(cds, color_by = "SCT_snn_res.0.1",size=1,show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)


exprData <- cds@assayData$exprs
exprData <- LogNormalize(exprData) 
cds$GLI2 <- exprData["GLI2",] 
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

