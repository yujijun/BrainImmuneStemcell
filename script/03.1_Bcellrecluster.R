# B cell reclusters ####
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

############ B cell selected ####
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
Bcells <- subset(control.seuobj.singlet,idents="B cells")

# ############  B cell reclusters ####
# Bcells <- NormalizeData(Bcells, normalization.method = "LogNormalize", scale.factor = 10000)
# Bcells <- FindVariableFeatures(Bcells, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(Bcells)
# Bcells <- ScaleData(Bcells, features = all.genes)
# Bcells <- RunPCA(Bcells, features = VariableFeatures(object = Bcells))
# Bcells <- FindNeighbors(Bcells, dims = 1:20)
# res.used <- seq(0.1,1,by=0.1) #Set different resolutions
# for(i in res.used){
#   Bcells <- FindClusters(object = Bcells, verbose = F, resolution = res.used)
# }
# library(clustree)
# clus.tree.out <- clustree(Bcells) +
#   theme(legend.position = "bottom") + 
#   scale_color_brewer(palette = "Set1") +
#   scale_edge_color_continuous(low = "grey80", high = "red")
# plot(clus.tree.out)
# Bcells <- RunUMAP(Bcells, dims = 1:20)
# Idents(Bcells) <- Bcells$RNA_snn_res.0.9
# DimPlot(Bcells)
# 
# #### Specific markers visualization ##### 
# Bcell_features <- c("Cd79a","Ighm","Ebf1",
#                     "H2-DMb2","Ms4a1","Cd74","Ighd","Igha",
#                     "Cd72","Vpreb3",
#                     "Dntt","Igll1",
#                     "Cd34")
# FeaturePlot(Bcells,features = Bcell_features)
# 
# Bcells.markers <- FindAllMarkers(Bcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Bcells.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(Bcells, features = top10$gene) + NoLegend()

#### B cell annotation by singleR ####
load("./reference/reference_geneset/immgen_ensemblID_20210505.RData")
DefaultAssay(Bcells) <- "RNA"
sce.test <- as.SingleCellExperiment(Bcells)
print("ready to calculate singleR!!!")
cluster_annotation.clusters <- SingleR(test = sce.test, ref = immgen,
                                       labels = immgen$label.fine, 
                                       assay.type.test=1,
                                       BPPARAM=MulticoreParam(8),
                                       clusters = Bcells$RNA_snn_res.0.9)

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

save(Bcells, file = paste0(output_path,prefix,"Bcells.RData"))

