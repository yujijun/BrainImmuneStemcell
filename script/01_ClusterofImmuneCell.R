#### Description #### 
# This is a QC and cluster about CD45+ immune cells in Meninges. 
# Update Time: Sat Oct 16 10:33:01 2021

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
input_path = "./input/"
output_dpath = "./output/data/"
output_fpath = "./output/figure/"
prefix = "r1_"
CCS_regression = FALSE
mt_regression = FALSE
pcadim = 20 
clusterResolution = 1
nFeature_min = 200 
nFeature_max = 7000
outputseurat_name = "01_ControlPreprocess.RData"
mytheme = theme(axis.title = element_text(face = "bold",size = 12,family = "Arial"),
                 axis.text = element_text(face = "bold",size = 10,family = "Arial"),
                 legend.text = element_text(face = "bold",size = 10,family = "Arial"))

#### Create seurat object ####
control <- Read10X(data.dir = paste0(input_path,'Mectrl/raw_feature_bc_matrix'))
Rp_genes <- grep(rownames(control),pattern = "^Rp")
control.seu <- CreateSeuratObject(counts = control,
                                  min.cells = 3,
                                  min.features = 100,
                                  project = "control")
control.seu[["percent.mt"]] <- PercentageFeatureSet(control.seu, pattern = "^mt-")
control.seu[["percent.Rp"]] <- PercentageFeatureSet(control.seu, pattern = "^Rp")
plotQC <- VlnPlot(control.seu,features = c("percent.mt","percent.Rp",
                                 "nFeature_RNA","nCount_RNA"),ncol = 4) +
  mytheme
ggsave(filename = paste0(output_fpath,prefix,"01_plotQC.pdf"),plot = plotQC,
       device = "pdf",height = 20,width = 20,units = "cm")
control <- subset(control.seu,subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt < 5)
#### Find Doublet #####
#Normalize
control <- NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)
#FindVariableFeatures
control <- FindVariableFeatures(control, selection.method = "vst", nfeatures = 2000)
#Scale
control <- ScaleData(control, features = rownames(control))
#PCA
control <- RunPCA(control, features = VariableFeatures(control),npcs = 50)
#cluster
control <- FindNeighbors(control, dims = 1:pcadim)
control <- FindClusters(control, resolution = clusterResolution)
control <- RunUMAP(control, dims = 1:20)
control <- RunTSNE(control, dims = 1:20)
# calculating parameters
sweep.res.list <- paramSweep_v3(control, PCs = 1:pcadim, sct = FALSE)
save(sweep.res.list,file = paste0(output_dpath,prefix,"sweep.res.list.RData"))
for(i in 1:length(sweep.res.list)){
  if(length(sweep.res.list[[i]]$pANN[is.nan(sweep.res.list[[i]]$pANN)]) != 0){
    if(i != 1){
      sweep.res.list[[i]] <- sweep.res.list[[i - 1]]
    }else{
      sweep.res.list[[i]] <- sweep.res.list[[i + 1]]
    }
  }
}
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
# defined desired numbers of doublet cells
nExp_poi <- round(0.1*length(colnames(control)))
# Calculating doublet cells
control <- doubletFinder_v3(control, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(control@meta.data)[ncol(control@meta.data)]="DoubletFinder"
DimplotDoublet <- DimPlot(control,reduction = "umap",pt.size = 0.5,group.by = "DoubletFinder") + 
  mytheme
ggsave(filename = paste0(output_fpath,prefix,"01_plotDoublet.pdf"),plot = DimplotDoublet,
       device = "pdf",height = 20,width = 20,units = "cm")
#### Save proprocessing result ####
control.seuobj <- control
save(control.seuobj,file = paste0(output_dpath,prefix,"01_control.seuobj.beforedoublet.RData"))
control.doublet.meta <- control@meta.data
save(control.doublet.meta,file = paste0(output_dpath,prefix,"01_control.doublet.meta.RData"))
