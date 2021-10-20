#procedure description for all scripts and files
## running log 
./script folders with all script of analysis 

-----------
**01_ClusterofImmuneCell.R**  was to perform seurat object,find doublet cells, cluster.

**Input:** ./input/Mectrl/raw_feature_bc_matrix/*

**Output:**

result was saved into `output/data` folder; 

r1_01_control.seuobj.withdoublet.RData # control seurat object after QC but with doublet

r1_02_control.doublet.meta.RData # doublet metadata from control object[r1_01_control.seuobj.withdoublet.RData]

r1_03_control.seuobj.singlet.RData #control seurat object after QC and filtering doublet cells.

r1_sweep.res.list.RData #This is a process backup file of DoubletFinder.

Attention: 注意在filter Sinelet 之后，重新进行降维和聚类；
---------------

---------
script:
** 02.1_AnnotationBySpecificmarkers.R **
Input:
r1_03_control.seuobj.singlet.RData 
Output:


--------


# > table(control.choosen$sum.v2[colnames(sc@fdata)])
# 
# Precursor cells         B cells         T cells        NK cells 
#              49              94              97              99 
#            pDCs            ILC2     Neutrophils      Mast cells 
#              75              98              70               8 
#       Monocytes            cDCs     Macrophages    Erythrocytes 
#             100              86              54              96 
#      Fibroblast 
#              33 



