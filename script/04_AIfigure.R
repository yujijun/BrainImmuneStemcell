############ load all package ##### 
setwd("/mnt/data/yjj/BrainImmune")
library(SingleR)
library(celldex)
library(AnnotationHub)
library(BiocParallel)
library(Seurat)
library(M3Drop)
library(ggplot2)
library(stringr)
library(dplyr)
library(Matrix)
library(gdata)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(ggthemes)
source("./reference/reference_script/scRNA_color.R")
extrafont::loadfonts()
############ Hyperparameter #### 
base_path <- "/mnt/data/yjj/BrainImmune/"
setwd(base_path)
output_path <- "./output/figure/"
cluster_color.v2 <- c("#0a62ae","#139277","#91cbbb","#7d91bc",
                      "#de7f07","#ba2d1b","#fdde86","#6791a7",
                      "#7d8aaf","#48b1d3","#e43983","#ef9572",
                      "#dd3e23","#8c741a","#99a2c0","#0a60ad",
                      "#1169b4","#4d4d4a","#6f68a7","f09f7c")

anno.final.color <- c(
  "#ba2d1b","#1169b4","#91cbbb","#139277","#6B8E23",
  "#48b1d3","#7d91bc","#4d4d4a","#de7f06","#BC8F8F",
  "#fdde86","#e43983","#8C549C")

anno.stemcell <- c(
  "#ba2d1b",rep("#DCDCDC",11)
)
mytheme <- theme(axis.title = element_text(face = "bold",size = 12,family = "Arial"),
                 axis.text = element_text(face = "bold",size = 10,family = "Arial"),
                 legend.text = element_text(face = "bold",size = 10,family = "Arial"))
prefix = "r4_"
############ Input  ####
load("/mnt/data/yjj/BrainImmune/output/data/r1_03_control.seuobj.singlet.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r3_control.seuobj.singlet.lastmeta.RData")
control.seuobj.singlet@meta.data <- control.seuobj.singlet.lastmeta

############ AI figure ######
#### Dimplot for sum.v2 ####
library(RColorBrewer)
plot1.2 <- DimPlot(control.seuobj.singlet,label = T,repel = T,reduction = "umap",
                  label.size = 5,cols = anno.final.color,pt.size = 0.3) + mytheme
ggsave(filename = paste0(output_path,prefix,"/01.2_Dimplot_sum.ident.v2.pdf"),plot = plot1.2,
       device = "pdf",height = 20,width = 23,units = "cm")

#### Dimplot for stemcell ####
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
plot1.3 <- DimPlot(control.seuobj.singlet,repel = F,reduction = "umap",
                   label.size = 5,cols = anno.stemcell,pt.size = 0.3) + 
  mytheme + NoLegend() 
ggsave(filename = paste0(output_path,prefix,"/01.3_Dimplot_Stemcell.pdf"),plot = plot1.3,
       device = "pdf",height = 10,width = 10,units = "cm")

#### Dotplot #### 
Allcellfeatures <- c("Kit","Cd34","Cd79a","Cd79b","Cd3g","Cd8b1","Cd4",
                     "Il7r","Nkg7","Gzma","Gata3","S100a8","S100a9","Camp","Ngp",
                     "Cpa3","Gata2",
                     "Ctss","Ms4a6c","Ly6c2","Lgals3",
                     "Cd209a","H2-Aa",
                     "Apoe",
                     "C1qa","Ms4a7",
                     "Hba-a1","Car2")
plot2 <- DotPlot(control.seuobj.singlet,features = Allcellfeatures,
                 cols = c("lightgrey", "#FF9900")) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) + 
  mytheme
ggsave(filename = paste0(output_path,prefix,"/02.1_DotPlot_allcellfeature.pdf"),plot = plot2,
       device = "pdf",height = 10,width = 22,units = "cm")

#### Stem cell #### 
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
StemSubset.3.final <- subset(control.seuobj.singlet,idents = "Stem cells")
StemSubset.3.final <- subset(StemSubset.3.final,
                             UMAP_1 > 0)
load("./output_20210922/StemSubset.final.RData")
Idents(StemSubset.3.final) <- StemSubset.3.final$singlerCluster.v4
StemCell_mark <- c("Lin28a","Kit","Ly6a","Flt3","Cd34",
                   "Slamf1","Cd48","Fcgr4","Csf1r")
#StemCell.subset <- subset(StemSubset.3.final,Lin28a == 0)
plot3 <- FeaturePlot(StemSubset.3.final,features = StemCell_mark,
                     reduction = "umap") + mytheme
ggsave(filename = paste0(output_path,prefix,"/03_FeaturePlot_stemcell_v2.pdf"),plot = plot3,
       device = "pdf",height = 15,width = 20,units = "cm")

#### B cell ####
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
Bcells <- subset(control.seuobj.singlet,idents = "B cells")
# Neutrophils.inBcells <- colnames(Bcells)[which(Bcells$singlerCluster == "Neutrophils (GN)")]
# control.seuobj.singlet$sum.v2[Neutrophils.inBcells] <- "Neutrophils"
Idents(Bcells) <- Bcells$singlerCluster
plot4 <- DimPlot(Bcells,label = T,label.box = T,
                 reduction = "umap",cols = col_cluster,
                 repel = T) + mytheme
ggsave(filename = paste0(output_path,prefix,"/04_DimPlot_Bcell.pdf"),plot = plot4,
       device = "pdf",height = 15,width = 20,units = "cm")

#### Bcell feature ####
Bcell_features <- c("Cd79a","Ighm","Ebf1",
                    "H2-DMb2","Ms4a1","Cd74","Ighd",
                    "Cd72","Vpreb3",
                    "Dntt","Igll1",
                    "Cd34")
plot5 <- FeaturePlot(Bcells,features = Bcell_features,reduction = "umap") + mytheme
ggsave(filename = paste0(output_path,prefix,"/05_FeaturePlot_Bcell.pdf"),plot = plot5,
       device = "pdf",height = 20,width = 30,units = "cm")

# Neutrophils.inBcells.v1 <- colnames(Bcells)[which(Bcells$singlerCluster == "B cells (proB.CLP)")]
# control.seuobj.singlet$sum.v2[Neutrophils.inBcells.v1] <- "Neutrophils"
# Bcell_subset <- subset(Bcells,tSNE_1 > -2)

#### cell ratio ####
orig <- c(rep("scRNAseq",4),rep("FACS",5))
celltype <- c("Myeloid cells","B cell","T cells","Erythrocytes",
              "Myeloid cells","B cell","T cells","Erythrocytes","Megakaryocyte")
celltype <- factor(celltype,levels = c("Myeloid cells","B cell","T cells","Erythrocytes","Megakaryocyte"))
cellfill <- c("#7d91bc","#1169b4","#91cbbb","#e43983",
              "#7d91bc","#1169b4","#91cbbb","#e43983","#F39B7FB2")
cellratio <- c(56.78,24.81,8.95,4.58,51.17,29.95,6.29,1.73,0.82)
df.cellratio <- data.frame(Method = orig,
                           Celltype = celltype,
                           CellRatio = cellratio,
                           cellfill = cellfill)
plot6 <- df.cellratio %>% ggbarplot(x = "Method",
                           y = "CellRatio",
                           fill = "Celltype",
                           width = 0.4,
                           label = T,
                           lab.pos = "in") + 
  theme(legend.position = "right") + 
  scale_fill_manual(values = cellfill) + 
  mytheme + 
  theme_base() 
ggsave(filename = paste0(output_path,prefix,"/06_ggbarplot_ratio.pdf"),plot = plot6,
       device = "pdf",height = 15,width = 15,units = "cm")

save(control,file = "AI_20210926/control.RData")

#### DoHeatmap #### 
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
control.seuobj.singlet.markers <- FindAllMarkers(control.seuobj.singlet,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
control.seuobj.singlet.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
plot7 <- DoHeatmap(control.seuobj.singlet, features = top10$gene) + 
  NoLegend() + mytheme
ggsave(filename = paste0(output_path,prefix,"/07_DoHeatmap.pdf"),plot = plot7,
       device = "pdf",height = 35,width = 25,units = "cm")

#### DC #### 
DC <- subset(control.seuobj.singlet,idents = "DC")
FeaturePlot(DC,features = c("Siglech","Itgam","Xcr1"))
pDC <- subset(DC,UMAP_1 > 4.5)
control.seuobj.singlet$sum.v2 <- control.seuobj.singlet$sum.v2
control.seuobj.singlet$sum.v2 <- as.character(control.seuobj.singlet$sum.v2)
control.seuobj.singlet$sum.v2[colnames(pDC)] <- "pDC"
