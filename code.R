library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(RColorBrewer)
library("Matrix")

bc <- readRDS("Breast-cancer-metastasis-SCTransform-annotation.rds")
Epithelial_bc <- readRDS("Breast-cancer-epithelial-SCTransform.rds")

################fig 1 a####################
plot <- DimPlot(bc, reduction = "tsne", pt.size = 0.1,cols= c("B cells"="#FFA500","Plasma B cells"="#F9F900",
                                                                        "CD4 T cells"="#FF9797","CD8 T cells"="#00F5FF",
                                                                        "NK cells"="#5CADAD","Myeloid cells"="#FFBFFF",
                                                                        "CAFs cells"="#66B3FF","Perivascular like cells"="#d0d0d0",
                                                                        "Endothelial cells"="#AB82FF","Epithelial cells"="#53FF53"),raster=FALSE)
ggsave(plot,filename="fig1a-tsne-cell_type-annotation.png",width=7,height=5)
###
plot <- FeaturePlot(bc, reduction = "umap", features = c("CD79A",#B cells
                                                   "MZB1","SDC1","JCHAIN",###Plasma B cells
                                                   "CD3D","CD4",##CD 4 T cells
                                                   "CD8A",##CD 8 T cells
                                                   "GNLY",##NK cells
                                                   "LYZ",##Myeloid cells
                                                   "PDGFRA",##CAF cells
                                                   "RGS5",##PVL cells
                                                   "PLVAP",##Endothelial cells
                                                   "EPCAM","KRT19"##Epithelial cells
), cols = c('#d8dcd6','red3'), ncol = 4, pt.size = 0.5, raster=FALSE)
ggsave(plot,filename="fig1a-dimplot-marker-genes-cell_type-annotation.png",width=23,height=20)
###
plot <- VlnPlot(bc, features = c("CD79A",#B cells
                           "MZB1","SDC1","JCHAIN",###Plasma B cells
                           "CD3D", "CD4",##CD 4 T cells
                           "CD8A",##CD 8 T cells
                           "GNLY",##NK cells
                           "LYZ",##Myeloid cells
                           "PDGFRA",##CAF cells
                           "RGS5",##PVL cells
                           "PLVAP",##Endothelial cells
                           "EPCAM","KRT19"##Epithelial cells
),  log = TRUE, ncol = 4,pt.size=0, raster=FALSE)
ggsave(plot,filename="fig1a-vlnplot-marker-genes-cell_type-annotation.png",width=25,height=20)

################fig 1 b####################
plot <- DoHeatmap(bc,disp.min = -1, features = c("CD79A",#B cells
                                          "MZB1","SDC1","JCHAIN",###Plasma B cells
                                          "CD3D", "CD4",##CD 4 T cells
                                          "CD8A",##CD 8 T cells
                                          "GNLY",##NK cells
                                          "LYZ",##Myeloid cells
                                          "PDGFRA",##CAF cells
                                          "RGS5",##PVL cells
                                          "PLVAP",##Endothelial cells
                                          "EPCAM","KRT19"##Epithelial cells
                                          ))+scale_fill_gradientn(colors = c("#000080", "#00CD00", "#F9F900"))
ggsave(plot,filename="fig1b-heatmap-cell_type-annotation.png",width=20,height=10)
#NoLegend()
################fig 1 C percent####################
library(tidyverse)
allcolour=c("#FFA500","#F9F900","#FF9797","#00F5FF","#5CADAD","#FFBFFF","#66B3FF","#d0d0d0","#AB82FF","#53FF53")

table(Idents(bc))
prop.table(table(Idents(bc)))
table(bc@meta.data$condition)
View(bc@meta.data)
as.data.frame(prop.table(table(Idents(bc), bc@meta.data$orig.ident), margin = 2))-> pdf -> td
as.data.frame(prop.table(table(Idents(bc), bc@meta.data$condition), margin = 2))-> pdf -> td
plt<- ggplot(td,aes(x=td[,2],y=td[,3],fill=td[,1]))+
  scale_x_discrete()+
  coord_flip()+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="replicate",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=allcolour)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))
ggsave(plt,filename="fig1c-condition-celltype-persent.png",width=8,height=3)
###
as.data.frame(table(Idents(bc), bc@meta.data$orig.ident), margin = 2)-> pdf -> td
as.data.frame(table(Idents(bc), bc@meta.data$condition), margin = 2)-> pdf -> td
plt<- ggplot(td,aes(x=td[,2],y=td[,3],fill=td[,1]))+
  scale_x_discrete()+
  coord_flip()+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="replicate",y="Cell number")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=allcolour)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))
ggsave(plt,filename="fig1c-condition-celltype-number.png",width=8,height=4)
###
as.data.frame(prop.table(table(Idents(bc), bc@meta.data$condition), margin = 1))-> pdf -> td
as.data.frame(table(Idents(bc), bc@meta.data$condition), margin = 2)-> pdf -> td

allcolour=c("#AB82FF","#FFDD55")
td$Var1 <- factor(td$Var1, levels=c("Epithelial cells",
                                    "Endothelial cells",
                                    "Perivascular like cells",
                                    "CAFs cells",
                                    "Myeloid cells",
                                    "NK cells",
                                    "CD8 T cells",
                                    "CD4 T cells",
                                    "Plasma B cells",
                                    "B cells"))

plt<- ggplot(td,aes(x=td[,1],y=td[,3],fill=td[,2]))+
  scale_x_discrete()+
  coord_flip()+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="replicate",y="Cell number")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=allcolour,limits=c("PT","LNM"))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))
ggsave(plt,filename="fig1c-condition-condition-percent.png",width=8,height=4)
###########################fig 1d epithelial######################################
Epithelial_bc@meta.data$condition <- factor(Epithelial_bc@meta.data$condition,levels = c("PT","LNM"))

allcolour=c("#AB82FF","#FFDD55")
plot <- DimPlot(Epithelial_bc, group.by = "condition",reduction = "tsne",cols=allcolour,  label = F, pt.size = 0.5,raster=FALSE)
ggsave(plot,filename="fig1d-tsne-epithlial-condition.png",width=5.5,height=5)
#######
Epithelial_bc@meta.data$id <- c(rep("id1",258),rep("id2",2768),rep("id3",52),rep("id4",5400),rep("id5",4330),rep("id6",502),rep("id7",1435),rep("id8",1494),rep("id1",60),rep("id2",3),rep("id3",1),rep("id4",0),rep("id5",86),rep("id6",113),rep("id7",63),rep("id8",1062))
allcolour=c("#00F5FF","#FFBFFF","#66B3FF","#FF8888","#53FF53","#AB82FF","#F9F900","#FFA500")
plot <- DimPlot(Epithelial_bc, group.by = "id",reduction = "tsne",cols=allcolour , label = F, pt.size = 0.5,raster=FALSE)
ggsave(plot,filename="fig1d-tsne-epithlial-id.png",width=5.5,height=5)
###########################fig 1e######################################
allcolour=c("#AB82FF","#FFDD55")
plot <- DimPlot(Epithelial_bc, reduction = "tsne",  label = F, pt.size = 0.5,raster=FALSE)
ggsave(plot,filename="fig1e-tsne-epithlial-cluster.png",width=5.5,height=5)
###########################fig 1g######################################
plot <- DoHeatmap(Epithelial_bc,disp.min = -1, features = c("SCGB1D2","TFF1","IGFBP5",#0
                                                   "MUCL1","APOD","S100A8",#1
                                                   "STC2","EFHD1","BAMBI",#2
                                                   "SFRP1","KRT15","SAA1",#3
                                                   "IRX2","SSFA2","UGT2B4",#4
                                                   "SRGN","IL32","CXCR4",#5
                                                   "MRPS30-DT","MRPS30","CNTNAP2",#6
                                                   "PNMT","AAMDC","KRT81",#7
                                                   "SCGB3A1","TFPI2","CYP4X1",#8
                                                   "MRPS30-DT","MX1","MRPS30","IFI27",#9
                                                   "UGT2B4","CLEC3A","SSFA2",#10
                                                   "HIST3H2A","SLC7A2","HIST1H1C",#11
                                                   "PIP","CYP2A6","COX6C",#12
                                                   "MALAT1","NEAT1","ZNF587",#13
                                                   "ACTA2","TAGLN","KRT14",#14
                                                   "MKI67","NUSAP1","ZWINT",#15
                                                   "AC005150.1","LRP1B","SCGB2A1","FGF10",#16
                                                   "GOLGA8A","GSDMB","USP54",#17
                                                   "C2orf40","ACTA2","CCL2"#18
    ))+scale_fill_gradientn(colors = c("#000080", "#00CD00", "#F9F900"))
ggsave(plot,filename="fig1g-heatmap-cluster.png",width=20,height=10)
###############################fig 1 h infercnv#######################################
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(AnnoProbe) 
library(infercnv)
bc <- readRDS("E:/A-analysis/BRCA/Nature communication-scRNAseq-lymph-metastasis/Breast-cancer-metastasis-SCTransform-annotation.rds")

scRNA_harmony = bc[, Idents(bc) %in% c( "B cells","Epithelial cells" )]
mono_sce = scRNA_harmony[,Idents(scRNA_harmony)== "Epithelial cells"]
scRNA_harmony$cell_type = as.character(Idents(scRNA_harmony))
mono_sce$cell_type = mono_sce$orig.ident
mono_sce$cell_type = as.character(mono_sce$cell_type)
scRNA_harmony$cell_type[match(colnames(mono_sce),colnames(scRNA_harmony))] =  mono_sce$cell_type
scRNA_harmony@meta.data$cell_type <- factor(scRNA_harmony@meta.data$cell_type, levels=c("B cells","CA1","LN1",
                                                                                        "CA2","LN2",
                                                                                        "CA3","LN3",
                                                                                        "CA4","LN4",
                                                                                        "CA5","LN5",
                                                                                        "CA6","LN6",
                                                                                        "CA7","LN7",
                                                                                        "CA8","LN8"))

dat <- GetAssayData(scRNA_harmony,assay = "RNA",slot = "counts")
dat <- as.data.frame(dat)

geneInfor=annoGene(rownames(dat),"SYMBOL",'human')  
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
dat=dat[match(geneInfor[,1], rownames(dat)),]  
rownames(geneInfor) <- geneInfor$SYMBOL   
geneInfor <- geneInfor[,-1]   

geneInfor$numbers <- as.numeric(gsub("[^0-9]", "", geneInfor$chr))
geneInfor <- geneInfor[order(geneInfor$numbers), ]
geneInfor <- geneInfor[,-which(names(geneInfor) == "numbers")]

#制作mate信息
meta <- subset(scRNA_harmony@meta.data,select = c("cell_type"))  

identical(colnames(dat),rownames(meta))  
identical(rownames(dat),rownames(geneInfor))

options(scipen = 100)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_names=c("B cells")) 
getwd()
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, #cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="B-cells-epithelial/", 
                             cluster_by_groups=T,
                             write_expr_matrix = T,
                             analysis_mode="subclusters",
                             denoise=T,
                             HMM=F,
                             num_threads = 30) 
saveRDS(infercnv_obj, file = "B-epithlial-infercnv_obj.rds")
infercnv_obj <- readRDS(file = "E:/A-analysis/BRCA/Nature communication-scRNAseq-lymph-metastasis/fig 1/B-epithlial-infercnv_obj.rds")

#######cnv scores##########################
cnv_table <- read.table("E:/A-analysis/BRCA/Nature communication-scRNAseq-lymph-metastasis/fig 1/B-cells-epithelial/infercnv.observations.txt", header=T,check.names = F)
data <- cnv_table
expr=data %>% as.matrix()
expr.scale =scale(t(expr))
tmp1=sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2=apply(expr.scale, 2, max) - apply(expr.scale,2,min)
expr_1=t(2*sweep(tmp1, 2, tmp2, "/")-1)
cnv_score=as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
cnv_score=rownames_to_column(cnv_score, var='cell')
###fig 1 i############
Anno <- subset(scRNA_harmony@meta.data,select = c("condition"))
Anno$cell=rownames(Anno)
colnames(Anno)[1]="cluster"
test=merge(cnv_score,Anno,by="cell", all=F)
test$cluster <- factor(test$cluster, levels=c("PT","LNM"))
write.table(test,file="cnv-scores-condition.xls",sep="\t",row.names=T,quote=F)

plot <- ggplot2::ggplot(test,aes(x=cluster,y=cnv_score))+
  geom_violin(aes(fill=cluster),cex=1.2)+ 
  scale_fill_manual(values = c('#00F5FF','#FB5554','#42F203','#579ABB','#B978AE',"#f0f0f4","#868B31","#FFBFFF","#66B3FF","#FF8888""#00F5FF","#FFBFFF","#66B3FF","#FF8888","#53FF53","#AB82FF"","#53FF53","#AB82FF","#F9F900","#00F5FF","#FFBFFF","#66B3FF","#FF8888","#53FF53","#AB82FF","#F9F900"))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none')
ggsave(plot,filename="fig1i-cnv-scores-condition.png",width=3,height=4)
###fig 1 j############
Anno <- subset(scRNA_harmony@meta.data,select = c("cell_type"))
Anno$cell=rownames(Anno)
colnames(Anno)[1]="cluster"
test=merge(cnv_score,Anno,by="cell", all=T)
write.table(test,file="cnv-scores-cell_type.xls",sep="\t",row.names=T,quote=F)

plot <- ggplot2::ggplot(test,aes(x=cluster,y=cnv_score))+
  geom_violin(aes(fill=cluster),cex=1.2)+ 
  scale_fill_manual(values = c("#E7298A",'#FB5554','#42F203',"#4DAF4A",'#00F5FF',"#FFD92F","#FCCDE5","#BEAED4",
                               "#D9D9D9","#666666","#FFFF99","#E6AB02","#386CB0","#7570B3"))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none')
ggsave(plot,filename="fig1i-cnv-scores-orig-ident.png",width=8,height=4)



