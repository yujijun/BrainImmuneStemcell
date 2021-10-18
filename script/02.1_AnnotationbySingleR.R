# SingleR cluster ####
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

### input #### 
load("./output/data/r2_control2125.RData")
load("./reference/reference_geneset/immgen_ensemblID_20210505.RData")
DefaultAssay(control2125) <- "RNA"
sce.test <- as.SingleCellExperiment(control2125)

#### singler annotation ####
cluster_annotation.singlecell <- SingleR(test = sce.test, ref = immgen,
                                         labels = immgen$label.fine, 
                                         assay.type.test=1,
                                         BPPARAM=MulticoreParam(8))
control2125$singler <- cluster_annotation.singlecell$pruned.labels
save(cluster_annotation.singlecell,file	= "./output/data/r2_cluster_annotation.singlecell.RData")
cluster_annotation.clusters <- SingleR(test = sce.test, ref = immgen,
                                       labels = immgen$label.fine, 
                                       assay.type.test=1,
                                       BPPARAM=MulticoreParam(8),
                                       clusters = control2125$RNA_snn_res.1)

final_cluster <- cluster_annotation.clusters$pruned.labels
current.ident <- levels(control2125$RNA_snn_res.1)
control2125@meta.data$singlerCluster = plyr::mapvalues(x = control2125$RNA_snn_res.1,
                                                   from = current.ident, 
                                                   to = final_cluster)
save(control2125,file="./output/data/r2_control2125.singleleR.RData")

Idents(control2125) <- control2125$singler
DimplotSinglecell <- DimPlot(control2125,label = T,label.size = 2) + NoLegend()
ggsave(filename = paste0(output_fpath,prefix,"14_control2125.singleRcell.pdf"),
       plot = DimplotSinglecell,
       device = "pdf",height = 25,width = 25,units = "cm")
Idents(control2125) <- control2125$singlerCluster
DimplotSingleRcluster <- DimPlot(control2125,label = T)
ggsave(filename = paste0(output_fpath,prefix,"15_control2125.singleRcluster.pdf"),
       plot = DimplotSingleRcluster,
       device = "pdf",height = 20,width = 20,units = "cm")

final_id <- c("Precursor cells",
              "Precursor cells",
              "B cells",
              "Mast cells")
current_id <- levels(control2125$RNA_snn_res.0.9)
control2125$sum.v2 <- plyr::mapvalues(
  x = control2125$RNA_snn_res.0.9,
  from = current_id,
  to = final_id
)
Idents(control2125) <- control2125$sum.v2
Dimplot.control2125.final <- DimPlot(control2125,label = T)
ggsave(filename = paste0(output_fpath,prefix,"16_control2125.final.pdf"),
       plot = Dimplot.control2125.final,
       device = "pdf",height = 20,width = 25,units = "cm")

control2125.final.meta <- control2125@meta.data
save(control2125.final.meta,file = paste0(output_dpath,
                                          prefix,
                                          "control2125.final.meta.RData"))




#### final all annotation 