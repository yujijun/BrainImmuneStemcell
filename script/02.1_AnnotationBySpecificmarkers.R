#### Description #### 
# Annotation of CD45+ immune cells in Meninges by specific markers. 
# Update time:Sat Oct 16 17:20:06 2021

#### Library and Hyperparameter ####
library(DoubletFinder)
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
library(readxl)
extrafont::loadfonts()
# Self-defined functions
convertHumanGeneList <- function(x){
  require(biomaRt)
  load(paste0(base_path,"reference_geneset/human_mouse.ensemblid.RData"))
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
Celltype_marker <- function(top20=top20,Cell_type, cluster){
  # return Find specific markers from cellmarker reference dataset
  # top20, is top20 DEgenes which filtered from FindAllmarkers result
  # Cell_type, is a character of cell type name which could be matched in mouse.immune.cellmarker reference dataset.
  # cluster, a cluster object in top20
  tmp_marker <- 
    mouse.immune.cellmarker[which(mouse.immune.cellmarker$cellName == Cell_type),]
  tmp_marker <- str_split(tmp_marker$cellMarker,pattern = ",")
  tmp_marker <- unlist(tmp_marker)
  tmp_marker <- str_remove_all(tmp_marker,pattern = " ")
  intersect(top20$gene[which(top20$cluster == cluster)],tmp_marker)
}
Celltype_marker_vall <- function(marker_matrix, Cell_type, cluster){
  # return, Find specific markers from cellmarker reference dataset
  # marker_matrix, is Whole DEgenes matrix which came from FindAllmarkers result
  # Cell_type, is a character of cell type name which could be matched in mouse.immune.cellmarker reference dataset.
  # cluster, a cluster object in top20
  tmp_marker <- 
    mouse.immune.cellmarker[which(mouse.immune.cellmarker$cellName == Cell_type),]
  tmp_marker <- str_split(tmp_marker$cellMarker,pattern = ",")
  tmp_marker <- unlist(tmp_marker)
  tmp_marker <- str_remove_all(tmp_marker,pattern = " ")
  intersect(marker_matrix$gene[which(marker_matrix$cluster == cluster)],tmp_marker)
}
#Hyperparameter
base_path = "/mnt/data/yjj/BrainImmune"
setwd(base_path)
output_dpath = "./output/data/"
output_fpath = "./output/figure/"
prefix = "r2_"
pcadim = 20 
clusterResolution = 0.9
outputseurat_name = "01_ControlPreprocess.RData"
mytheme = theme(axis.title = element_text(face = "bold",size = 12,family = "Arial"),
                axis.text = element_text(face = "bold",size = 10,family = "Arial"),
                legend.text = element_text(face = "bold",size = 10,family = "Arial"))
source("./reference/reference_script/scRNA_color.R")
#Load referece genesets 
mouse.immune.cellmarker <- read_xlsx("./reference/reference_geneset/single.cellmarker.mouse.xlsx")
load("./reference/reference_geneset/03_部分免疫细胞19_较为准确.rda")
load("./reference/reference_geneset/02_部分免疫细胞_比较明显的区分_from常识.rda")
#### Input dataset #### 
load("./output/data/r1_03_control.seuobj.singlet.RData")
Idents(control.seuobj.singlet) <-control.seuobj.singlet$RNA_snn_res.0.9
#### Visualization about cluster situation #### 
ClusterPlot <- DimPlot(control.seuobj.singlet,
        cols = my36colors,label = T)
ggsave(filename = paste0(output_fpath,prefix,"04_ClusterPlot.pdf"),
       plot = ClusterPlot,
       device = "pdf",height = 20,width = 25,units = "cm")

#### Find and visualization of top20 markers of cluster ####
con.seu.markers <- FindAllMarkers(control.seuobj.singlet,
                                  only.pos = T)
top20 <- con.seu.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

Top20Heatmap <- DoHeatmap(control.seuobj.singlet,
                         features = top20$gene,
                         angle = 90,
                         size = 3) + mytheme
ggsave(filename = paste0(output_fpath,prefix,"05_Top20Heatmap.pdf"),
       plot = Top20Heatmap,
       device = "pdf",height = 50,width = 40,units = "cm")

#### Visualization of Specific markers####
specificmarkers <- c("Kit","Cd34","Cd79a","Cd79b","Cd3g","Cd8b1","Cd4",
                     "Il7r","Nkg7","Gzma","Gata3","S100a8","S100a9","Camp","Ngp",
                     "Cpa3","Gata2",
                     "Ctss","Ms4a6c","Ly6c2","Lgals3",
                     "Cd209a","H2-Aa",
                     "Apoe",
                     "C1qa","Ms4a7",
                     "Hba-a1","Car2")
BiomarkersHeatmap <- DoHeatmap(control.seuobj.singlet,
                         features = specificmarkers,
                         angle = 90,
                         size = 3) + mytheme
ggsave(filename = paste0(output_fpath,prefix,"06_BiomarkersHeatmap.pdf"),
       plot = BiomarkersHeatmap,
       device = "pdf",height = 20,width = 30,units = "cm")



