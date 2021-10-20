# BrainImmuneStemcell
Update Time: Sat Oct 16 09:59:41 2021

## Desciription: 
This is a repository about Single cell analysis of CD45+ immune cells in meninges

## Work Space:

R version 4.0.3 (2020-10-10)

Platform: x86_64-pc-linux-gnu (64-bit)

Running under: Ubuntu 20.04.1 LTS

other attached packages:  
 [1] DoubletFinder_2.0.3        
 [2] ggpubr_0.4.0               
 [3] RColorBrewer_1.1-2         
 [4] uwot_0.1.10                
 [5] scran_1.18.7               
 [6] scater_1.18.6 
 [7] SingleCellExperiment_1.12.0
 [8] gdata_2.18.0               
 [9] Matrix_1.3-4               
[10] dplyr_1.0.7                
[11] stringr_1.4.0              
[12] ggplot2_3.3.5              
[13] M3Drop_1.16.0              
[14] numDeriv_2016.8-1.1        
[15] SeuratObject_4.0.2         
[16] Seurat_4.0.4               
[17] BiocParallel_1.24.1        
[18] AnnotationHub_2.22.1       
[19] BiocFileCache_1.14.0       
[20] dbplyr_2.1.1               
[21] celldex_1.0.0              
[22] SingleR_1.4.1              
[23] SummarizedExperiment_1.20.0
[24] Biobase_2.50.0             
[25] GenomicRanges_1.42.0       
[26] GenomeInfoDb_1.26.7        
[27] IRanges_2.24.1             
[28] S4Vectors_0.28.1           
[29] BiocGenerics_0.36.1        
[30] MatrixGenerics_1.2.1       
[31] matrixStats_0.61.0        

## Documentation 

### Input of files 

Input files shoule be in `input/Mectrl/raw_feature_bc_matrix` folder 
### Script files

01_ClusterofImmuneCell.R  #perform clustering 

02.1_AnnotationBySpecificmarkers.R  #Annotation by specific markers methods

02.2_AnnotationbySingleR.R  #Annotated by SingleR

03_resultIntegration.R # integration results from 01,02.1,02.2. 

```
load("/mnt/data/yjj/BrainImmune/output/data/r2_control.seuobj.singlet.initmeta.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r2_control2125.final.meta.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r2_control9.meta.RData") #Tcell and NK cells
load("/mnt/data/yjj/BrainImmune/output/data/r1_03_control.seuobj.singlet.RData")
```

03.1_Bcellrecluster.R  # for B cell recluster. Attention: recluster before B recluster.

04_AIfigure.R # visualization 

### Output files
├── data
│   ├── r1_01_control.seuobj.withdoublet.RData
│   ├── r1_02_control.doublet.meta.RData
│   ├── r1_03_control.seuobj.singlet.RData # final Count dataset
│   ├── r1_control.seuobj.singlet.lastmeta.RData
│   ├── r1_sweep.res.list.RData
│   ├── r2_cluster_annotation.singlecell.RData
│   ├── r2_con2125.seu.markers.RData
│   ├── r2_con9.seu.markers.RData
│   ├── r2_con.seu.singlet.0.9.marker.RData
│   ├── r2_control2125.final.meta.RData
│   ├── r2_control2125.RData
│   ├── r2_control2125.singleleR.RData
│   ├── r2_control9.meta.RData
│   ├── r2_control.seuobj.singlet.initmeta.RData
│   ├── r2_control.seuobj.singlet.lastmeta.RData
│   ├── r3_control.seuobj.singlet.lastmeta.RData
│   ├── r3_final.markers.sum.v2.1.RData # final markers
│   ├── r3_final.markers.sum.v2.RData
│   ├── r4_Bcell.last.RData
│   ├── r4_Bcells_annotation.singlecells.RData
│   ├── r4_Bcells_annotation.clusters.RData
│   └── r4_control.seuobj.singlet.lastmeta.RData #final metadata version 
└── figure
    ├── r1_01_plotDoublet.pdf
    ├── r1_01_plotQC.pdf
    ├── r2_01_clustertree.pdf
    ├── r2_04_ClusterPlot.pdf
    ├── r2_05_Top20Heatmap.con.singlet.0.9.pdf
    ├── r2_06_BiomarkersHeatmap.con.singlet.0.9.pdf
    ├── r2_07_Featureplot.con.singlet.0.9.pdf
    ├── r2_08_Mast.pdf
    ├── r2_09_Monocyte.pdf
    ├── r2_10_cluster2125.pdf
    ├── r2_11_Cluster2125.pdf
    ├── r2_12_Top20Heatmap.con2125.singlet.0.9.pdf
    ├── r2_13_control2125.sum.pdf
    ├── r2_14_control2125.singleRcell.pdf
    ├── r2_15_cluster9.pdf
    ├── r2_15_control2125.singleRcluster.pdf
    ├── r2_16_Cluster9.pdf
    ├── r2_16_control2125.final.pdf
    ├── r2_17_Top20Heatmap.con9.singlet.0.9.pdf
    ├── r2_1_Dimplot.control.singlet.sum.v2.pdf
    ├── r3_1_Dimplot.control.singlet.sum.v2.pdf
    ├── r3_2_Top20Heatmap.control.singlet.sum.v2.pdf
    ├── r3_3_Top20Heatmap.control.singlet.sum.v2.1.pdf
    ├── r4_10_pie.pdf
    ├── r4_1_Bcell_Dimplot.pdf
    ├── r4_1_Dimplot_sum.v2.pdf
    ├── r4_2_Bcell_feature.pdf
    ├── r4_2_Dimplot_Stemcell.pdf
    ├── r4_3_DotPlot_allcellfeature.pdf
    ├── r4_4_FeaturePlot_stemcell_v2.pdf
    ├── r4_5_FeaturePlot_Bcell.pdf
    ├── r4_6_ggbarplot_ratio.pdf
    ├── r4_7.1_DoHeatmap.pdf
    ├── r4_7.2_DoHeatmap.pdf
    └── r4_7_DoHeatmap.pdf

All r4_* was the final version 