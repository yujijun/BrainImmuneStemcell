##------Update Time: Wed Oct 20 14:49:26 2021 ------##
#### load all package ##### 
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
#### Hyperparameter #### 
base_path <- "/mnt/data/yjj/BrainImmune/"
setwd(base_path)
output_dpath <- "./output/data/"
mytheme <- theme(axis.title = element_text(face = "bold",size = 12,family = "Arial"),
                 axis.text = element_text(face = "bold",size = 10,family = "Arial"),
                 legend.text = element_text(face = "bold",size = 10,family = "Arial"))
prefix = "r4_"
height = 178
width = 178
#### Input control seuobj ####
load("/mnt/data/yjj/BrainImmune/output/data/r1_03_control.seuobj.singlet.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r4_control.seuobj.singlet.lastmeta.RData")
control.seuobj.singlet@meta.data <- control.seuobj.singlet.lastmeta
#### B cell selected ####
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
Bcells <- subset(control.seuobj.singlet,idents="B cells")

#### B cell reclusters ####
Bcells <- NormalizeData(Bcells, normalization.method = "LogNormalize", scale.factor = 10000)
Bcells <- FindVariableFeatures(Bcells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Bcells)
Bcells <- ScaleData(Bcells, features = all.genes)
Bcells <- RunPCA(Bcells, features = VariableFeatures(object = Bcells))
Bcells <- FindNeighbors(Bcells, dims = 1:20)
res.used <- seq(0.1,1,by=0.1) #Set different resolutions
for(i in res.used){
  Bcells <- FindClusters(object = Bcells, verbose = F, resolution = res.used)
}
library(clustree)
clus.tree.out <- clustree(Bcells) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
plot(clus.tree.out)
Bcells <- RunUMAP(Bcells, dims = 1:20)
Idents(Bcells) <- Bcells$RNA_snn_res.0.9
save(Bcells,file = paste0(output_dpath,prefix,"Bcell.RData"))

#### B cell annotation by singleR ####
load("./output/data/r4_Bcell.RData")
load("./reference/reference_geneset/immgen_ensemblID_20210505.RData")
DefaultAssay(Bcells) <- "RNA"
sce.test <- as.SingleCellExperiment(Bcells)
print("ready to calculate singleR!!!")
cluster_annotation.clusters <- SingleR(test = sce.test, ref = immgen,
                                       labels = immgen$label.fine, 
                                       assay.type.test=1,
                                       BPPARAM=MulticoreParam(8),
                                       clusters = Bcells$RNA_snn_res.0.9)

save(cluster_annotation.clusters,file = paste0(output_path,
                                               prefix,
                                               "Bcells_annotation.clusters.RData"))
final_cluster <- cluster_annotation.clusters$pruned.labels
current.ident <- levels(Bcells$RNA_snn_res.0.9)
Bcells@meta.data$singlerCluster = plyr::mapvalues(x = Bcells$RNA_snn_res.0.9,
                                                  from = current.ident,
                                                  to = final_cluster)

cluster_annotation.singlecells <- SingleR(test = sce.test, ref = immgen,
                                       labels = immgen$label.fine, 
                                       assay.type.test=1,
                                       BPPARAM=MulticoreParam(8))
Bcells@meta.data$singlecells <- cluster_annotation.singlecells$pruned.labels

save(cluster_annotation.singlecells, file = paste0(output_path,prefix,"Bcells_annotation.singlecells.RData"))
save(Bcells,file = paste0(output_dpath,prefix,"Bcells.last.RData"))
#### visualization of B cells ####
plotB <- DimPlot(Bcells,label = T,label.size = 4,cols = col_cluster) + NoLegend() + mytheme
ggsave(filename = paste0(output_path,prefix,"1_Bcell_Dimplot.pdf"),plot = plotB,
       device = "pdf",height = height,width = width,units = "mm")

Bcell_features <- c("Cd79a","Ighm","Ebf1",
                    "Dntt","Igll1",
                    "Vpreb3","Cd72",
                    "Cd74","Ms4a1","Ighd","H2-DMb2","Igha"
                    )
plotB.feature <- FeaturePlot(Bcells,features = Bcell_features) + mytheme
ggsave(filename = paste0(output_path,prefix,"2_Bcell_feature.pdf"),plot = plotB.feature,
       device = "pdf",height = height,width = width*3/2,units = "mm")
