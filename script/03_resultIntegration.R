############# Library and Hyperparameter ####
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
base_path = "/mnt/data/yjj/BrainImmune"
setwd(base_path)
output_dpath = "./output/data/"
output_fpath = "./output/figure/"
prefix = "r3_"
mytheme = theme(axis.title = element_text(face = "bold",size = 12,family = "Arial"),
                axis.text = element_text(face = "bold",size = 10,family = "Arial"),
                legend.text = element_text(face = "bold",size = 10,family = "Arial"))

############# input ####
load("/mnt/data/yjj/BrainImmune/output/data/r2_control.seuobj.singlet.initmeta.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r2_control2125.final.meta.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r2_control9.meta.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r1_03_control.seuobj.singlet.RData")

############# integration and visualization ####
control.seuobj.singlet$sum.v2 <- control.seuobj.singlet.initmeta$sum.v1
control.seuobj.singlet$sum.v2 <- as.character(control.seuobj.singlet$sum.v2)
control2125.final.meta$sum.v1.1 <- control2125.final.meta$sum
control2125.final.meta$sum.v1.1 <- plyr::mapvalues(
  x = control2125.final.meta$sum.v1.1,
  from = "B plasma cells",
  to = "B cells"
)
control.seuobj.singlet$sum.v2[rownames(control2125.final.meta)] <- as.character(control2125.final.meta$sum.v1.1)
control.seuobj.singlet$sum.v2[rownames(control9.meta)] <- as.character(control9.meta$sum)
control.seuobj.singlet$sum.v2 <- factor(control.seuobj.singlet$sum.v2,
                                        levels = c(
                                          'Precursor cells',
                                          'B cells',
                                          'T cells',
                                          'NK cells',
                                          'pDCs',
                                          "ILC2",
                                          "Neutrophils",
                                          "Mast cells",
                                          "Monocytes",
                                          "cDCs",
                                          "Macrophages",
                                          "Erythrocytes",
                                          "Fibroblast"
                                        ))
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2

Dimplot.control.singlet.sum.v2 <- DimPlot(control.seuobj.singlet,label = T) + NoLegend()
ggsave(filename = paste0(output_fpath,prefix,"1_Dimplot.control.singlet.sum.v2.pdf"),
       plot = Dimplot.control.singlet.sum.v2,
       device = "pdf",height = 20,width = 20,units = "cm")

control.seuobj.singlet.lastmeta <- control.seuobj.singlet@meta.data
save(control.seuobj.singlet.lastmeta,file = 
       paste0(output_dpath,prefix,"control.seuobj.singlet.lastmeta.RData"))


############## calculation final top20 marker genes  #############
Idents(control.seuobj.singlet) <- control.seuobj.singlet$sum.v2
final.markers.sum.v2.1 <- FindAllMarkers(control.seuobj.singlet,
                                   only.pos = T)
save(final.markers.sum.v2.1,file = paste0(output_dpath,prefix,
                                    "final.markers.sum.v2.1.RData"))

top20 <- final.markers.sum.v2.1 %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

Top20Heatmap <- DoHeatmap(control.seuobj.singlet,
                          features = top20$gene,
                          angle = 90,
                          size = 3) + mytheme
ggsave(filename = paste0(output_fpath,prefix,"3_Top20Heatmap.control.singlet.sum.v2.1.pdf"),
       plot = Top20Heatmap,
       device = "pdf",height = 60,width = 30,units = "cm")

Celltype_marker(top20,"Precursor cells","Precursor cells")
