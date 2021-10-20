##------ Update Time: Wed Oct 20 11:14:46 2021 ------##
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
output_fpath <- "./output/figure/"
anno.final.color <- c(
  "#ba2d1b","#1169b4","#91cbbb","#139277","#6B8E23",
  "#48b1d3","#7d91bc","#4d4d4a","#de7f06","#BC8F8F",
  "#fdde86","#e43983","#8C549C")

anno.stemcell <- c(
  "#ba2d1b",rep("#DCDCDC",12)
)
mytheme <- theme(axis.title = element_text(face = "bold",size = 12,family = "Arial"),
                 axis.text = element_text(face = "bold",size = 10,family = "Arial"),
                 legend.text = element_text(face = "bold",size = 10,family = "Arial"))
prefix = "r4_"
height = 178
width = 178
############ Input  ####
load("/mnt/data/yjj/BrainImmune/output/data/r1_03_control.seuobj.singlet.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r4_control.seuobj.singlet.lastmeta.RData")
control.seuobj.singlet@meta.data <- control.seuobj.singlet.lastmeta
load("/mnt/data/yjj/BrainImmune/output/data/r4_Bcell.RData")
############ AI figure ######
#### 01_Dimplot for sum.v2 ####
library(RColorBrewer)
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
plot1.2 <- DimPlot(control.seuobj.singlet,label = T,repel = T,reduction = "umap",
                  label.size = 5,cols = anno.final.color,pt.size = 0.3) + mytheme
ggsave(filename = paste0(output_fpath,prefix,"1_Dimplot_sum.v2.pdf"),plot = plot1.2,
       device = "pdf",height = height,width = width+20,units = "mm")

#### 02_Dimplot for only stemcell ####
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
plot1.3 <- DimPlot(control.seuobj.singlet,repel = F,reduction = "umap",
                   label.size = 5,cols = anno.stemcell,pt.size = 0.3) +
  mytheme + NoLegend()
ggsave(filename = paste0(output_fpath,prefix,"2_Dimplot_Stemcell.pdf"),plot = plot1.3,
       device = "pdf",height = height/3,width = width/3,units = "mm")

#### 03_Dotplot for all clusters ####
Allcellfeatures <- c("Kit","Cd34","Ifitm1","Cd79a","Cd79b","Cd3g","Cd8b1","Cd4",
                     "Il7r","Nkg7","Klrd1","Siglech",
                     "Gata3","S100a8","S100a9","Camp","Ngp",
                     "Cpa3","Gata2",
                     "Ctss","Ms4a6c","Ly6c2","Lgals3",
                     "H2-Aa",
                     "Apoe",
                     "C1qa","Ms4a7",
                     "Hba-a1","Car2",
                     "Gngt2","Ear2")
plot2 <- DotPlot(control.seuobj.singlet,features = Allcellfeatures,
                 cols = c("lightgrey", "#FF9900")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  mytheme
ggsave(filename = paste0(output_fpath,prefix,"3_DotPlot_allcellfeature.pdf"),plot = plot2,
       device = "pdf",height = height,width = width*3/2,units = "mm")

#### 04_Featurepot of Stem cells ####
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
StemSubset.3.final <- subset(control.seuobj.singlet,idents = "Precursor cells")
#StemSubset.3.final <- subset(StemSubset.3.final,
#                             UMAP_1 > 0)
#load("./output_20210922/StemSubset.final.RData")
#Idents(StemSubset.3.final) <- StemSubset.3.final$singlerCluster.v4
StemCell_mark <- c("Lin28a","Kit","Ly6a","Flt3","Cd34",
                   "Slamf1","Cd48","Fcgr4","Csf1r")
#StemCell.subset <- subset(StemSubset.3.final,Lin28a == 0)
plot3 <- FeaturePlot(StemSubset.3.final,features = StemCell_mark,
                     reduction = "umap") + mytheme
ggsave(filename = paste0(output_fpath,prefix,"4_FeaturePlot_stemcell_v2.pdf"),plot = plot3,
       device = "pdf",height = height,width = width*3/2,units = "mm")

#### 05_Dimplot of B cells ####
Idents(Bcells) <- Bcells$singlerCluster
plot4 <- DimPlot(Bcells,label = T,label.box = T,
                 reduction = "umap",cols = col_cluster,
                 repel = T) + mytheme
ggsave(filename = paste0(output_fpath,prefix,"4_DimPlot_Bcell.pdf"),plot = plot4,
       device = "pdf",height = 15,width = 20,units = "cm")

#### 06_Featureplot of B cells ####
Bcell_features <- c("Cd79a","Ighm","Ebf1",
                    "Dntt","Igll1",
                    "Vpreb3","Cd72",
                    "Cd74","Ms4a1","Ighd","H2-DMb2","Igha")
plot5 <- FeaturePlot(Bcells,features = Bcell_features,reduction = "umap") + mytheme
ggsave(filename = paste0(output_fpath,prefix,"5_FeaturePlot_Bcell.pdf"),plot = plot5,
       device = "pdf",height = height,width = width*3/2,units = "mm")

#### 07_cell proportion between scRNAseq and FACS ####
orig <- c(rep("scRNAseq",4),rep("FACS",5))
celltype <- c("Myeloid cells","B cell","T cells","Erythrocytes",
              "Myeloid cells","B cell","T cells","Erythrocytes","Megakaryocyte")
celltype <- factor(celltype,levels = c("Myeloid cells","B cell","T cells","Erythrocytes","Megakaryocyte"))
cellratio <- c(55.85,24.77,9.14,4.59,51.17,29.95,6.29,1.73,0.82)
df.cellratio <- data.frame(Method = orig,
                           Celltype = celltype,
                           CellRatio = cellratio)
cellfill <- c("#53A85F","#1169b4","#91cbbb","#e43983","#D6E7A3",
              "#53A85F","#1169b4","#91cbbb","#e43983","#D6E7A3")
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
ggsave(filename = paste0(output_fpath,prefix,"6_ggbarplot_ratio.pdf"),plot = plot6,
       device = "pdf",height = height,width = width*2/3,units = "mm")


#### 08_Top20 DoHeatmap of all cells #### 
load("/mnt/data/yjj/BrainImmune/output/data/r3_final.markers.sum.v2.1.RData")
final.markers.sum.v2.1 %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
plot7 <- DoHeatmap(control.seuobj.singlet, features = top20$gene,
                   angle = 90,size = 3,group.colors = anno.final.color) + 
  NoLegend() + mytheme
ggsave(filename = paste0(output_fpath,prefix,"7.2_DoHeatmap.pdf"),plot = plot7,
       device = "pdf",height = height,width = width,units = "mm")



#### 09_Balloonplot between seurat and RaceID3 ####
library(gplots)
library(tidyverse)
load("./software/RaceID3_StemID2-master/RaceID3_1019/cluster.RData")
tab.1 <- table(cluster$cpart,cluster$anno)
pdf(paste0(output_fpath,prefix,"balloonplot.pdf"),
    width = 12,height = 10)
balloonplot(tab.1)
dev.off()

#### 10_ pie plot for cell proportion
tab <- table(control.seuobj.singlet$sum.v2)
names = names(tab)
piepercent = paste0(round(100*tab/sum(tab),2), "%")
names(tab) <- paste0(names,"_",piepercent)
pdf(paste0(output_fpath,prefix,"10_pie.pdf"),
    width = 10,height = 10)
pie(tab,labels = names(tab),col = anno.final.color)
dev.off()

#### 11_All RaceID3_StemID2 result ####
#All of the cold in RaceID3_StemID2_sample_main_sampling_FaceID.R  which matched result was in `RaceID3_1019` folder.
# also filter one outlier cell version in filteroncells.R, which match result was in `RaceID3_1020` folder
