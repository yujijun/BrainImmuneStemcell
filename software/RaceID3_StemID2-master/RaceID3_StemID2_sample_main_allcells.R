#### load all library and dataset ####
work_path <- "/mnt/data/yjj/BrainImmune/software/RaceID3_StemID2-master"
setwd(work_path)
source("RaceID3_StemID2_class.R")
library(tidyverse)
library(Seurat)
load("/mnt/data/yjj/BrainImmune/output/data/r1_03_control.seuobj.singlet.RData")
load("/mnt/data/yjj/BrainImmune/output/data/r4_control.seuobj.singlet.lastmeta.RData")
control.seuobj.singlet@meta.data <- control.seuobj.singlet.lastmeta

#### hyperparameter #### 
output_path = "RaceID3_1019_allcells/"
outminc.use=10
outlg.use=10
pdishuf.use = 10000
prefix = "outminc10.outlg10_"

# # ## input data
# choosen_sample <- control.seuobj.singlet@meta.data %>%
#     rownames_to_column() %>%
#     filter(sum.v2 != "pDCs" &
#                sum.v2 != "Mast cells" &
#                sum.v2 != "Fibroblast",
#            sum.v2 != "Precursor cells") %>%
#     group_by(sum.v2) %>%
#     sample_n(size = 100) %>%
#     distinct() %>%
#     ungroup()
# 
# choosen_sample.v1 <- control.seuobj.singlet@meta.data %>%
#     rownames_to_column() %>%
#     filter(sum.v2 == "pDCs" | sum.v2 == "Mast cells"|
#                sum.v2  == "Fibroblast" |sum.v2 == "Precursor cells")
# choosen_sample.all <- rbind(choosen_sample,choosen_sample.v1)
# choosen_cell <- choosen_sample.all$rowname
# control.choosen <- subset(control.seuobj.singlet,cells = choosen_cell)

#### create SCseq object and clustering ####
x <- control.seuobj.singlet@assays$RNA@counts
# prdata: data.frame with transcript counts for all genes (rows) in all cells (columns); with rownames == gene ids; remove ERCC spike-ins 
#prdata <- x[grep("ERCC",rownames(x),invert=TRUE),-1]
prdata <- as.data.frame(x)
## RaceID2
# initialize SCseq object with transcript counts
sc <- SCseq(prdata)

# filtering of expression data
sc <- filterdata(sc, mintotal=3000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=FALSE, sfn=FALSE, hkn=FALSE, dsn=1, rseed=17000, CGenes=NULL, FGenes=NULL, ccor=.4)

# regress out the batch effect
# optional:
#vars <- data.frame(row.names=names(sc@fdata),batch=sub("(_|)\\d.+","",names(sc@fdata)))
#sc@fdata <- varRegression(sc@fdata,vars)

# correct for cell cycle, proliferation, and expression of degradation markers by PCA
# optional:
require(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
g   <- rownames(sc@fdata)
k   <- getBM(attributes = c("external_gene_name", "go_id","name_1006"),filters="external_gene_name",values=g,mart=mart)
gCC <- name2id( k$external_gene_name[k$name_1006 == "cell cycle"],rownames(sc@fdata)) 
gCP <- name2id( k$external_gene_name[k$name_1006 == "cell proliferation"],rownames(sc@fdata))
vset <- list(gCC,gCP)
x <- CCcorrect(sc@fdata,vset=vset,CGenes=NULL,ccor=.4,nComp=NULL,pvalue=.05,quant=.01,mode="pca")
# number of principal components that have been removed
x$n
# loadings of the first principal component that has been removed
y <- x$pca$rotation[,x$n[1]]
# genes from vset are either enriched in the head or the tail of this list
tail(y[order(y,decreasing=TRUE)],10)
# reassign the corrected expression matrix to sc@fdata
sc@fdata <- x$xcor

# k-medoids clustering
sc <- clustexp(sc,clustnr=30,bootnr=50,metric="pearson",
               do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",
               SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids",FSelect=TRUE)
# compute t-SNE map
sc <- comptsne(sc,rseed=15555,sammonmap=FALSE,initial_cmd=TRUE,fast=TRUE,perplexity=30)
# detect outliers and redefine clusters
sc <- findoutliers(sc, outminc=outminc.use,outlg=outlg.use,probthr=1e-3,thr=2**-(1:40),outdistquant=.95)
# reassign clusters based on random forest
sc <- rfcorrect(sc,rfseed=12345,final=TRUE,nbfactor=5)
sc.name = paste0(output_path,prefix,"sc.RData")
save(sc,file = sc.name)


#### save all the figure #####
pdf(file = paste0(output_path,prefix,"_plotsaturation.v1.pdf"),
    width = 10,
    height = 10)
plotsaturation(sc,disp=TRUE)
# plot change of the within-cluster dispersion as a function of the cluster number: only if sat == TRUE
dev.off()

pdf(file = paste0(output_path,prefix,"_plotsaturation.v2.pdf"),
    width = 10,
    height = 10)
plotsaturation(sc)
dev.off()

# silhouette of k-medoids clusters
pdf(file = paste0(output_path,prefix,"_plotsilhouette.pdf"),
    width = 10,
    height = 20)
plotsilhouette(sc)
dev.off()

# Jaccard's similarity of k-medoids clusters
pdf(file = paste0(output_path,prefix,"_plotjaccard.pdf"),
    width = 10,
    height = 10)
plotjaccard(sc)
dev.off()

# barchart of outlier probabilities
pdf(file = paste0(output_path,prefix,"_plotoutlierprobs.pdf"),
    width = 10,
    height = 10)
plotoutlierprobs(sc)
dev.off()

# regression of background model
pdf(file = paste0(output_path,prefix,"_plotbackground.pdf"),
    width = 10,
    height = 10)
plotbackground(sc)
dev.off()

# dependence of outlier number on probability threshold (probthr)
pdf(file = paste0(output_path,prefix,"_plotsensitivity.pdf"),
    width = 10,
    height = 10)
plotsensitivity(sc)
dev.off()

# heatmap of k-medoids cluster
pdf(file = paste0(output_path,prefix,"_clustheatmap.kmedoids.pdf"),
    width = 10,
    height = 10)
clustheatmap(sc,final=FALSE,hmethod="single")
dev.off()

# heatmap of final cluster
pdf(file = paste0(output_path,prefix,"_clustheatmap.final.pdf"),
    width = 10,
    height = 10)
clustheatmap(sc,final=TRUE,hmethod="single")
dev.off()

# highlight k-medoids clusters in t-SNE map
pdf(file = paste0(output_path,prefix,"_plottsne.kmedoids.pdf"),
    width = 10,
    height = 10)
plottsne(sc,final=FALSE)
dev.off()

# highlight final clusters in t-SNE map
pdf(file = paste0(output_path,prefix,"_plottsne.final.pdf"),
    width = 10,
    height = 10)
plottsne(sc,final=TRUE)
dev.off()

# highlight cell labels in t-SNE map
scLabel <- control.seuobj.singlet@meta.data[colnames(sc@fdata),]$sum.v2
pdf(file = paste0(output_path,prefix,"_plotlabelstsne.pdf"),
    width = 10,
    height = 10)
plotlabelstsne(sc,labels=scLabel)
dev.off()

# highlight groups of cells by symbols in t-SNE map
pdf(file = paste0(output_path,prefix,"_plotsymbolstsne.pdf"),
    width = 10,
    height = 10)
plotsymbolstsne(sc,types=scLabel)
dev.off()

# highlight transcirpt counts of a set of genes in t-SNE map, e. g. all Apoa genes
pdf(file = paste0(output_path,prefix,"_plotexptsne.pdf"),
    width = 10,
    height = 10)
g <- c("Kit")
plotexptsne(sc,g,n="Kit genes",logsc=TRUE)
dev.off()
## StemID2

# initialization
ltr <- Ltree(sc)
# computation of the entropy
ltr <- compentropy(ltr)
# computation of the projections for all cells
ltr <- projcells(ltr,cthr=2,nmode=FALSE)
# computation of the projections for all cells after randomization
ltr <- projback(ltr,pdishuf=pdishuf.use,nmode=FALSE,fast=FALSE,rseed=17000)
# assembly of the lineage tree
ltr <- lineagetree(ltr,pthr=0.05,nmode=FALSE,fast=FALSE)
# determination of significant differentiation trajectories
ltr <- comppvalue(ltr,pethr=0.05,nmode=FALSE,fast=FALSE)

save(ltr,file = paste0(output_path,prefix,"ltr.RData"))

## diagnostic plots
# histogram of ratio between cell-to-cell distances in the embedded and the input space
pdf(file = paste0(output_path,prefix,"_plotdistanceratio.pdf"),
    width = 10,
    height = 10)
plotdistanceratio(ltr)
dev.off()

# t-SNE map of the clusters with more than cthr cells including a minimum spanning tree for the cluster medoids
pdf(file = paste0(output_path,prefix,"_plotmap.pdf"),
    width = 10,
    height = 10)
plotmap(ltr)
dev.off()

# visualization of the projections in t-SNE space overlayed with a minimum spanning tree connecting the cluster medoids
pdf(file = paste0(output_path,prefix,"_plotmapprojections.pdf"),
    width = 10,
    height = 10)
plotmapprojections(ltr)
dev.off()

# lineage tree showing the projections of all cells in t-SNE space
pdf(file = paste0(output_path,prefix,"_plottreewithcells.pdf"),
    width = 10,
    height = 10)
plottree(ltr,showCells=TRUE,nmode=FALSE,scthr=0)
dev.off()

# lineage tree without showing the projections of all cells
pdf(file = paste0(output_path,prefix,"_plottree.withoutcells.pdf"),
    width = 10,
    height = 10)
plottree(ltr,showCells=FALSE,nmode=FALSE,scthr=0)
dev.off()

# heatmap of the enrichment p-values for all inter-cluster links
pdf(file = paste0(output_path,prefix,"_plotlinkpv.pdf"),
    width = 10,
    height = 10)
plotlinkpv(ltr)
dev.off()

# heatmap of the link score for all inter-cluster links
pdf(file = paste0(output_path,prefix,"_plotlinkscore.pdf"),
    width = 10,
    height = 10)
plotlinkscore(ltr)
dev.off()

# heatmap showing the fold enrichment (or depletion) for significantly enriched or depleted links
pdf(file = paste0(output_path,prefix,"_projenrichment.pdf"),
    width = 10,
    height = 10)
projenrichment(ltr)
dev.off()

## extract projections onto all links for all cells in a given cluster i
x <- getproj(ltr,i=1)
pdf(file = paste0(output_path,prefix,"_pheatmap_project.cluster1.v1.pdf"),
    width = 10,
    height = 10)
# heatmap of all projections for cluster i
pheatmap(x$pr)
dev.off()


pdf(file = paste0(output_path,prefix,"_pheatmap_project.cluster1.v2.pdf"),
    width = 10,
    height = 10)
# heatmap of z-score for all projections for cluster i
pheatmap(x$prz)
dev.off()

## extracting all cells on two branches sharing the same cluster and computing differentially expressed genes between these two branches
# x <- branchcells(ltr,list("1.3","1.2"))
# # z-scores for differentially expressed genes
# head(x$diffgenes$z)
# # plotting the cells on the two branches as additional clusters in the t-SNE map
# plottsne(x$scl)

## computing the StemID score
x <- compscore(ltr,nn=1)
#plotting the StemID score
pdf(file = paste0(output_path,prefix,"_plotscore.pdf"),
    width = 10,
    height = 15)
plotscore(ltr,1)
dev.off()
print("well done !!! all is good.")

# # retrieve cells from branch in pseudo-temporal order as inferred by the projection coordinates
# n <- cellsfromtree(ltr,c(6,2))
# # filter out lowly expressed genes
# fs  <- filterset(ltr@sc@ndata,n=n$f,minexpr=2,minnumber=1)
# # compute self organizing map (SOM) of co-expressed genes
# s1d <- getsom(fs,nb=1000,k=5,locreg=TRUE,alpha=.5)
# ps  <- procsom(s1d,corthr=.85,minsom=3)
# # coloring scheme for clusters ( vector with colours)
# fcol <- ltr@sc@fcol
# y    <- ltr@sc@cpart[n$f]
# # plot average z-score for all modules derived from the SOM
# plotheatmap(ps$nodes.z,xpart=y,xcol=fcol,ypart=unique(ps$nodes),xgrid=FALSE,ygrid=TRUE,xlab=T)
# # plot z-score profile of each gene
# plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE,ylab = T)
# # plot normalized expression profile of each gene
# plotheatmap(ps$all.e,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
# # plot binarized expression profile of each gene (z-score < -1, -1 < z-score < 1, z-score > 1)
# plotheatmap(ps$all.b,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
# 
# # extract all genes from node 1 of the SOM
# g <- names(ps$nodes)[ps$nodes == 1]
# # plot average expression profile of these genes along the trajectory
# plotexpression(fs,y,g,n$f,k=25,col=fcol,name="Node 1",cluster=FALSE,locreg=TRUE,alpha=.5,types=NULL)
# # plot expression of a single gene
# plotexpression(fs,y,"Cd79b",n$f,k=25,col=fcol,cluster=FALSE,locreg=TRUE,alpha=.5,types=NULL)
# # plot average expression profile of these genes along the trajectory, highlighting batch origin
# plotexpression(fs,y,g,n$f,k=25,col=fcol,name="Node 1",cluster=FALSE,locreg=TRUE,alpha=.5,types=sub("\\_\\d+","",n$f))
# 







