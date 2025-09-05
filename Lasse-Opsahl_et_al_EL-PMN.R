#This script recreates Figures 2B-N, S2A-D, S5D-F
#Raw data files for the novel dataset generated in this manuscript are available through the NIH Gene Expression Omnibus (GEO), accession number GSE292712.
#Control lungs (GSM8864044, GSM8864046, GSM8864048, GSM8864052,	GSM8864056)
#iKRAS 16wk ON lungs (GSM8864043, GSM8864045, GSM8864047)
#iKRAS 16wk ON + 1wk OFF lungs (GSM9123780, GSM9123781)

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#version.string R version 4.4.0 (2024-04-24) --- Puppy Cup 
#Seurat Version 5.3.0

#Load packages
library(CoGAPS)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(magrittr)
library(devtools)
library(janitor)
library(scCustomize)
library(readxl)
library(fgsea)
library(msigdbr)
library(slurmR)
library(cowplot)
library(tidyverse)
library(pheatmap)
library(scGSVA)
library(enrichplot)
library(viridis)

#### Pre-processing ----------------------------------------------------------------------------------------------------####

#Read in metadata file
Samples <- read_excel("IntegrateEL-PMNiKRASscRNAseqSamples.xlsx")
Samples$Description <- NULL

#keep only Lung samples
SamplesLung <- subset(Samples, Tissue %in% c("Lung"))

samples <- as.character(as.list(SamplesLung$Sample))
samples_list <- vector("list", length = length(samples))

#Make seurat objects and add metadata
for (x in 1:nrow(SamplesLung)){
  print(SamplesLung[x,"SampleID"])
  sc <- Read_CellBender_h5_Mat(as.character(SamplesLung[x,"file"]))
  samples_list[[x]] <- CreateSeuratObject(sc, min.cells = 3, min.features = 100)
  #Metadata
  metaDataGroups <- colnames(SamplesLung)
  metaData <- as.character(as.list(SamplesLung[x,]))
  for (y in 1:length(metaData)){
    group <- metaDataGroups[[y]]
    z <- metaData[[y]]
    samples_list[[x]][[group]] <- z
    samples_list[[x]]@meta.data[["file"]] <- NULL
  }
}
names(samples_list) <- samples

#merge seurat objects
Lung <- Merge_Seurat_List(samples_list, add.cell.ids = samples)

#Normalize
Lung <- NormalizeData(object = Lung, normalization.method = "LogNormalize", scale.factor = 10000)

#Check percent mitochondrial genes to check for doublets and poor cell quality
Idents(Lung) <- "Sample"
Lung[["percent.mt"]] <- PercentageFeatureSet(object = Lung, pattern = "^mt-")

#Keep cells with between 800 and 100,000 RNA reads and less than 15% mitochondrial genes
Lung <- subset(x = Lung, subset = nCount_RNA > 800 & nCount_RNA < 100000 & percent.mt < 15)

#Identify variable genes
Lung <- FindVariableFeatures(Lung, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(Lung)
Lung <- ScaleData(Lung, verbose = T, features = all.genes)
Lung <- RunPCA(Lung, npcs = 30, verbose = T)

#Integrate
Lung <- IntegrateLayers(object = Lung, method = RPCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- Lung@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#17

#Find Neighbors and cluster cells
Lung <- FindNeighbors(object = Lung, dims = 1:varFinal) #change dims to 1:x (number of PCs)
Lung <- FindClusters(object = Lung, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
Lung <- RunUMAP(Lung, reduction = "integrated.rpca", dims = 1:varFinal, verbose = F) #change dims to 1:x (number of PCs)

#Join seurat layers
Lung <- JoinLayers(Lung)

#Identify clusters based on marker expression
DotPlot(Lung, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdgfrb","Rgs5","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                         "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Try4","Ctrb1","Kras",
                         "Trp53","EGFP","Sec14l3","Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Cd200r3","Fcer1a","Mcam","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Label Manual Clusters:
Idents(Lung) <- "seurat_clusters"
Lung <- RenameIdents(Lung, 
                     "0" = "B Cell", 
                     "1" = "CD4+ T Cell", 
                     "2" = "Macrophage", 
                     "3" = "B Cell", 
                     "4" = "Neutrophil", 
                     "5" = "Endothelial", 
                     "6" = "CD8+ T Cell", 
                     "7" = "Neutrophil", 
                     "8" = "Alveolar Epithelial",
                     "9" = "NK Cell", 
                     "10" = "CD8+ T Cell", 
                     "11" = "Treg", 
                     "12" = "B Cell", 
                     "13" = "Interstitial Macrophage", 
                     "14" = "Monocyte", 
                     "15" = "CD4+ T Cell", 
                     "16" = "Alveolar Macrophage", 
                     "17" = "Endothelial",
                     "18" = "RBC", 
                     "19" = "Endothelial", 
                     "20" = "Interstitial Macrophage", 
                     "21" = "CD8+ T Cell",
                     "22" = "B Cell", 
                     "23" = "Fibroblast", 
                     "24" = "RBC", 
                     "25" = "Macrophage", 
                     "26" = "mDC", 
                     "27" = "Endothelial", 
                     "28" = "Fibroblast",
                     "29" = "Endothelial",
                     "30" = "NK Cell",
                     "31" = "Pericyte",
                     "32" = "B Cell",
                     "33" = "Endothelial",
                     "34" = "Fibroblast", 
                     "35" = "Neutrophil", 
                     "36" = "B Cell",
                     "37" = "Endothelial",
                     "38" = "Pericyte",
                     "39" = "Basophil", 
                     "40" = "pDC",
                     "41" = "Proliferating T Cell", 
                     "42" = "cDC1", 
                     "43" = "Proliferating Interstitial Macrophage",
                     "44" = "Plasma Cell",
                     "45" = "Fibroblast",
                     "46" = "Endothelial",
                     "47" = "Alveolar Epithelial",
                     "48" = "Neutrophil",
                     "49" = "Monocyte",
                     "50" = "Mesothelial",
                     "51" = "Bronchial Epithelial",
                     "52" = "Neutrophil",
                     "53" = "Proliferating B Cell",
                     "54" = "Proliferating T Cell", 
                     "55" = "Alveolar Epithelial", 
                     "56" = "Endothelial", 
                     "57" = "Active DC", 
                     "58" = "Neutrophil",
                     "59" = "Proliferating cDC1",
                     "60" = "Endothelial", 
                     "61" = "Alveolar Epithelial", 
                     "62" = "Treg", 
                     "63" = "Alveolar Macrophage", 
                     "64" = "B Cell", 
                     "65" = "B Cell", 
                     "66" = "B Cell", 
                     "67" = "Alveolar Macrophage", 
                     "68" = "B Cell",
                     "69" = "Cardiomyocyte", 
                     "70" = "Alveolar Macrophage")

Lung[["manual_clusters"]] <- Lung@active.ident

#Simple Clusters
Idents(Lung) <- "manual_clusters"
Lung <- RenameIdents(Lung,"Alveolar Epithelial" = "Epithelial",
                     "Bronchial Epithelial" = "Epithelial",
                     "Endothelial" = "Endothelial",
                     "Mesothelial" = "Mesothelial",
                     "Fibroblast" = "Fibroblast",
                     "Macrophage" = "Macrophage",
                     "Alveolar Macrophage" = "Macrophage",
                     "Interstitial Macrophage" = "Macrophage",
                     "Active DC" = "DC",
                     "cDC1" = "DC",
                     "pDC" = "DC",
                     "mDC" = "DC",
                     "Monocyte" = "Monocyte",
                     "Basophil" = "Granulocyte",
                     "Neutrophil" = "Granulocyte",
                     "CD8+ T Cell" = "T Cell",
                     "Treg" = "T Cell",
                     "CD4+ T Cell" = "T Cell",
                     "NK Cell" = "NK Cell",
                     "B Cell" = "B Cell",
                     "Proliferating cDC1" = "Proliferating",
                     "Proliferating B Cell" = "Proliferating",
                     "Proliferating T Cell" = "Proliferating",
                     "Proliferating Interstitial Macrophage" = "Proliferating",
                     "Pericyte" = "Pericyte",
                     "RBC" = "RBC",
                     "Cardiomyocyte" = "Cardiomyocyte")
Lung[["simple_clusters"]] <- Lung@active.ident

#Re-order clusters
Idents(Lung) <- "simple_clusters"

new_order <- c("Fibroblast",
               "Pericyte",
               "Endothelial",
               "Macrophage",
               "Monocyte",
               "Granulocyte",
               "DC",
               "T Cell",
               "NK Cell",
               "B Cell",
               "Plasma Cell",
               "Epithelial",
               "Mesothelial",
               "Proliferating",
               "Cardiomyocyte",
               "RBC")

Lung@active.ident <- factor(Lung@active.ident, levels = new_order)
Lung[["simple_clusters"]] <- Lung@active.ident

#Re-order groups
Lung$Status <- factor(Lung@meta.data[["Status"]], levels = c("Healthy","EL-PMN","OFF"))

#### Pre-processing Fibroblasts ----------------------------------------------------------------------------------------------------####

#subset out the fibroblasts from the larger object
LungFB <- subset(Lung, idents = "Fibroblast")

#split seurat layers by sample
LungFB[["RNA"]] <- split(LungFB[["RNA"]], f = LungFB$Sample)

#Identify variable genes
LungFB <- FindVariableFeatures(LungFB, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungFB)
LungFB <- ScaleData(LungFB, verbose = T, features = all.genes)
LungFB <- RunPCA(LungFB, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungFB@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#17

#Find Neighbors and cluster cells
LungFB <- FindNeighbors(object = LungFB, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungFB <- FindClusters(object = LungFB, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungFB <- RunUMAP(LungFB, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-fibroblasts from the object
LungFB <- subset(LungFB, idents = c("5","8","17","18","20","21"), invert = T) #remove contamination

#Repeat pipeline after removing contaminating cells (Part 2)
#Identify variable genes
LungFB <- FindVariableFeatures(LungFB, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungFB)
LungFB <- ScaleData(LungFB, verbose = T, features = all.genes)
LungFB <- RunPCA(LungFB, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungFB@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#17

#Find Neighbors and cluster cells
LungFB <- FindNeighbors(object = LungFB, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungFB <- FindClusters(object = LungFB, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungFB <- RunUMAP(LungFB, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-fibroblasts from the object
LungFB <- subset(LungFB, idents = c("8","13","16","20"), invert = T) #remove contamination

#Repeat pipeline after removing contaminating cells (Part 3)
#Identify variable genes
LungFB <- FindVariableFeatures(LungFB, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungFB)
LungFB <- ScaleData(LungFB, verbose = T, features = all.genes)
LungFB <- RunPCA(LungFB, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungFB@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#19

#Find Neighbors and cluster cells
LungFB <- FindNeighbors(object = LungFB, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungFB <- FindClusters(object = LungFB, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungFB <- RunUMAP(LungFB, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#rejoin seurat layers
LungFB <- JoinLayers(LungFB)

#Run CoGAPS
LungFB_list <- SplitObject(LungFB, split.by = 'Sample')
LungFB_list <- LungFB_list[unlist(lapply(LungFB_list, function(x){dim(x)[[2]] > 5}))]
LungFB_list <- lapply(LungFB_list, function(x){
  message(unique(x$Sample))
  x <- NormalizeData(x)
  x <- ScaleData(x)
  x <- FindVariableFeatures(x)
})
features <- lapply(LungFB_list, VariableFeatures) %>% reduce(union)
mtx <- GetAssayData(LungFB, assay = 'RNA', layer = 'data')[features,]

params <- CogapsParams(nPatterns=5, nIterations=50000, seed=42, sparseOptimization=TRUE)
params <- setDistributedParams(params, nSets=5)
res <- CoGAPS(as.matrix(mtx), params)

saveRDS(res, "LungFB_cogaps_result.Rds")

patterns_in_order <-t(res@sampleFactors[colnames(mtx),])

LungFB[["CoGAPS"]] <- CreateAssayObject(counts = patterns_in_order)

save(LungFB, file = 'LungFB_CoGAPS_May2025.RData')

#Name fibroblast seurat clusters based on their highest expressed CoGAPS Pattern

DefaultAssay(LungFB) <- "CoGAPS"
pattern_names = rownames(LungFB@assays$CoGAPS)

VlnPlot(LungFB, features = pattern_names, group.by = "seurat_clusters", pt.size = 0.1, ncol = 5)&stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)

Idents(LungFB) <- "seurat_clusters"
LungFB <- RenameIdents(LungFB,
                       "0" = "FB2",  
                       "1" = "FB4", 
                       "2" = "FB3", 
                       "3" = "FB5", 
                       "4" = "FB5", 
                       "5" = "FB5", 
                       "6" = "FB2", 
                       "7" = "FB5",
                       "8" = "FB4",
                       "9" = "FB5",
                       "10" = "FB1",  
                       "11" = "FB5", 
                       "12" = "FB1", 
                       "13" = "FB1", 
                       "14" = "FB2",
                       "15" = "FB4",
                       "16" = "FB2",
                       "17" = "FB5",
                       "18" = "FB5",
                       "19" = "FB2",
                       "20" = "FB2",  
                       "21" = "FB2")

LungFB[["CoGAPS_clusters"]] <- LungFB@active.ident

#Re-order clusters
Idents(LungFB) <- "CoGAPS_clusters"

new_order <- c("FB1",
               "FB2",
               "FB3",
               "FB4",
               "FB5")

LungFB@active.ident <- factor(LungFB@active.ident, levels = new_order)
LungFB[["CoGAPS_clusters"]] <- LungFB@active.ident


#### Pre-processing Macrophages and Monocytes ----------------------------------------------------------------------------------------------------####

#subset out the fibroblasts from the larger object
LungMac <- subset(Lung, idents = c("Macrophage","Monocyte","Interstitial Macrophage",
                                   "Alveolar Macrophage","Proliferating Interstitial Macrophage"))

#split seurat layers by sample
LungMac[["RNA"]] <- split(LungMac[["RNA"]], f = LungMac$Sample)

#Identify variable genes
LungMac <- FindVariableFeatures(LungMac, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungMac)
LungMac <- ScaleData(LungMac, verbose = T, features = all.genes)
LungMac <- RunPCA(LungMac, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungMac@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#14

#Find Neighbors and cluster cells
LungMac <- FindNeighbors(object = LungMac, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungMac <- FindClusters(object = LungMac, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungMac <- RunUMAP(LungMac, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungMac, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                            "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                            "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-macrophages&monocytes from the object
LungMac <- subset(LungMac, idents = c("20","21","22","23","24","27","29","30","32"), invert = T) #remove contamination

##Repeat pipeline after removing contaminating cells (Part 2)

#Identify variable genes
LungMac <- FindVariableFeatures(LungMac, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungMac)
LungMac <- ScaleData(LungMac, verbose = T, features = all.genes)
LungMac <- RunPCA(LungMac, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungMac@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#14

#Find Neighbors and cluster cells
LungMac <- FindNeighbors(object = LungMac, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungMac <- FindClusters(object = LungMac, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungMac <- RunUMAP(LungMac, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungMac, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","C1qa","Apoe", 
                            "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                            "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-macrophages&monocytes from the object
LungMac <- subset(LungMac, idents = c("12","17","22","27","28"), invert = T) #remove contamination

##Repeat pipeline after removing contaminating cells (Part 3)

#Identify variable genes
LungMac <- FindVariableFeatures(LungMac, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungMac)
LungMac <- ScaleData(LungMac, verbose = T, features = all.genes)
LungMac <- RunPCA(LungMac, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungMac@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#14

#Find Neighbors and cluster cells
LungMac <- FindNeighbors(object = LungMac, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungMac <- FindClusters(object = LungMac, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungMac <- RunUMAP(LungMac, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungMac, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","Arg1","C1qa","Apoe", 
                            "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                            "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-macrophages&monocytes from the object
LungMac <- subset(LungMac, idents = c("17","21"), invert = T) #remove contamination

##Repeat pipeline after removing contaminating cells (Part 4)

#Identify variable genes
LungMac <- FindVariableFeatures(LungMac, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(LungMac)
LungMac <- ScaleData(LungMac, verbose = T, features = all.genes)
LungMac <- RunPCA(LungMac, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- LungMac@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#14

#Find Neighbors and cluster cells
LungMac <- FindNeighbors(object = LungMac, dims = 1:varFinal) #change dims to 1:x (number of PCs)
LungMac <- FindClusters(object = LungMac, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
LungMac <- RunUMAP(LungMac, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(LungMac, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","C1qa","Apoe", 
                            "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                            "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#rejoin layers
LungMac <- JoinLayers(LungMac)

Idents(LungMac) <- "seurat_clusters"
LungMac <- RenameIdents(LungMac,
                        "0" = "Macrophage",  
                        "1" = "Macrophage", 
                        "2" = "Monocyte", 
                        "3" = "Alveolar Macrophage", 
                        "4" = "Classical Monocyte", 
                        "5" = "Macrophage", 
                        "6" = "Macrophage", 
                        "7" = "Interstitial Macrophage",
                        "8" = "Interstitial Macrophage",
                        "9" = "Interstitial Macrophage",
                        "10" = "Classical Monocyte",  
                        "11" = "Interstitial Macrophage", 
                        "12" = "Activated Macrophage", 
                        "13" = "Alveolar Macrophage", 
                        "14" = "Activated Macrophage",
                        "15" = "Monocyte",
                        "16" = "Alveolar Macrophage",
                        "17" = "Monocyte",
                        "18" = "Alveolar Macrophage",
                        "19" = "Alveolar Macrophage",
                        "20" = "Monocyte",  
                        "21" = "Activated Macrophage", 
                        "22" = "Alveolar Macrophage", 
                        "23" = "Interstitial Macrophage", 
                        "24" = "Activated Macrophage", 
                        "25" = "Interstitial Macrophage", 
                        "26" = "Classical Monocyte")

LungMac[["Myeloid_clusters"]] <- LungMac@active.ident

#Re-order clusters
Idents(LungMac) <- "Myeloid_clusters"

new_order <- c("Monocyte",
               "Classical Monocyte",
               "Alveolar Macrophage",
               "Interstitial Macrophage",
               "Macrophage",
               "Activated Macrophage")

LungMac@active.ident <- factor(LungMac@active.ident, levels = new_order)
LungMac[["Myeloid_clusters"]] <- LungMac@active.ident

#### Figures ----------------------------------------------------------------------------------------------------####

####Figure 2B ####
Idents(Lung) <- "simple_clusters"
DimPlot(Lung, cols = c("#D7A8FF","#9990FF","#6A5D99","#FF9360",
                       "#FFB2D8","#73D279","#FFE18F","#52E7ED",
                       "#FFB300","#94435F","#009193","#D266A1",
                       "#7FA4CA","#528AFF","#DF4945","lightgray"), raster = F, split.by = "Status") #16 colors

#### Supplementary Figure S2A ####
#Clustering DotPlot
Idents(Lung) <- "simple_clusters"
DotPlot(Lung,split.by = "Status", features=c("Pdgfra", "Clec3b", "Col1a1","Col1a2", #fibroblast
                                             "Pdgfrb","Rgs5", #Pericyte
                                             "Pecam1", "Cdh5","Flt1", #Endothelial
                                             "Ptprc","H2-Eb1","Itgam","Cd68","Adgre1","Mrc1","Siglecf","Itgax", #macrophage
                                             "Ly6c2","Cd14",  #monocyte
                                             "Fcgr3","Ly6g", #granulocyte
                                             "Flt3", "Batf3", #DC
                                             "Cd3e", "Cd4","Cd8a", #T Cell
                                             "Nkg7","Klrb1c","Ccr5", #NK Cell
                                             "Cd79a", "Cd19", "Ms4a1", #B Cell
                                             "Igha","Ighm", #Plasma cell
                                             "Lamp3","Cxcl15","Aqp5", "Sftpd","Krt8", "Krt18", #Epithelial
                                             "Msln", "Wt1", #Mesothelial
                                             "Ccna2", "Mki67", #Proliferating
                                             "Myl7", "Tnni3","Tnnt2", #Cardiomyocyte
                                             "Hbb-bt" #RBC
), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#### Supplementary Figure S2B ####
Idents(Lung) <- "Sample"
groups <-levels(factor(unique(Idents(Lung))))

Idents(Lung) <- "simple_clusters"
pt <- table(Idents(Lung), Lung$orig.ident)
pt <- as.data.frame(pt)

i=1
for (i in 1:length(groups)) {
  print(groups[i])
  Idents(Lung) <- "Sample"
  x <- subset(Lung, idents = groups[i])
  Idents(x) <- "simple_clusters"
  pt1 <- as.data.frame(table(Idents(x), x$orig.ident))
  if (nrow(pt) != length(pt1)){
    z <- c(pt$Var1,pt1$Var1)
    a <- z[!(duplicated(z)|duplicated(z, fromLast=TRUE))]
    row <- which(!(duplicated(z)|duplicated(z, fromLast=TRUE)), arr.ind = TRUE)
    pt1<-pt1 %>% add_row(Var1 = a, Var2 = "SeuratProject", Freq = 0, .before = row)
  }
  list <- pt1$Freq
  pt$y <- list
  names(pt)[names(pt) == 'y'] <- groups[i]
  i=i+1
}

#save table
pt <- pt %>%adorn_totals("row")
write.csv(pt, "Lung_CellAbundance_AllCells_perSample_5.18.25.csv", row.names = T)

pt_per <-adorn_percentages(pt, "col")
write.csv(pt_per, "Lung_CellPercentages_AllCells_perSample_5.18.25.csv", row.names = T)

#Relative abundance bar graph
pt_per <- pt_per %>% filter(row_number() <= n()-1)
pt_per <- pt_per[-c(2,3)]

pt_plot <- data.frame()
for (i in 1:length(groups)) {
  x <- pt_per[c("Var1",groups[i])]
  names(x)[names(x) == groups[i]] <- "Frequency"
  names(x)[names(x) == 'Var1'] <- "Cell_Types"
  x$Group <- groups[i]
  pt_plot <- rbind(pt_plot, x)
}

ggplot(pt_plot, aes(fill = Cell_Types, x = Group, y = Frequency))+
  scale_x_discrete(limits = c("Healthy_16wk_ON_lung_1","Healthy_16wk_ON_lung_2","Healthy_16wk_ON_lung_3",
                              "Healthy_20wk_ON_lung_1","Healthy_20wk_ON_lung_2",
                              "iKras_16wk_ON_lung_1","iKras_16wk_ON_lung_2","iKras_16wk_ON_lung_3",
                              "iKras_16wk_ON_1wk_OFF_lung_1","iKras_16wk_ON_1wk_OFF_lung_2")) +
  geom_bar(position = "stack", stat = "identity")+
  ggtitle("Lung Relative Cell Type Abundance")+
  theme(plot.title = element_text(hjust = 0.5))+
  RotatedAxis()+
  theme(axis.text.x=element_text(size=11), axis.title.x = element_blank())+
  scale_fill_manual(values = c("Fibroblast" = "#D7A8FF",
                               "Pericyte" = "#9990FF",
                               "Endothelial" = "#6A5D99",
                               "Macrophage" = "#FF9360",
                               "Monocyte" = "#FFB2D8",
                               "Granulocyte" = "#73D279",
                               "DC" = "#FFE18F",
                               "T Cell" = "#52E7ED",
                               "NK Cell" = "#FFB300",
                               "B Cell" = "#94435F",
                               "Plasma Cell" = "#009193",
                               "Epithelial" = "#D266A1",
                               "Mesothelial" = "#7FA4CA",
                               "Proliferating" = "#528AFF",
                               "Cardiomyocyte" = "#DF4945",
                               "RBC" = "lightgray"))

#### Figure 2C ####
VlnPlot(LungFB, features = "cogaps_Pattern-1", group.by = "Status", pt.size = 0.1, c("lightgray","#D7A8FF", "grey40"))+NoLegend()+
  theme(axis.title.x = element_blank())+stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)

#### Supplementary Figure S2C ####
DefaultAssay(LungFB) <- "CoGAPS"
pattern_names = rownames(LungFB@assays$CoGAPS)

VlnPlot(LungFB, features = pattern_names, group.by = "Status", pt.size = 0.1, c("lightgray","#D7A8FF", "grey40"), ncol = 5)+NoLegend()&
  stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)

#### Figure 2D ####
h_gene_sets = msigdbr(species = "mouse", category = "H") 
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

hallmarks_ora <- getPatternGeneSet(res,
                                   gene.sets = msigdbr_list,
                                   method = "overrepresentation")

pl_pattern1 <- plotPatternGeneSet(patterngeneset = hallmarks_ora, whichpattern = 1, padj_threshold = 0.05)+
  scale_fill_continuous(low = "lightskyblue",high = "palevioletred")
pl_pattern1

#### Figure 2E ####
DimPlot(LungFB, label = T,repel = T, cols = c("lightskyblue","plum1","#7A81FF","palevioletred","navy"))+NoLegend()

#### Figure 2F ####
Idents(LungFB) <- "Status"
groups <-levels(factor(unique(Idents(LungFB))))

Idents(LungFB) <- "CoGAPS_clusters"
pt <- table(Idents(LungFB), LungFB$orig.ident)
pt <- as.data.frame(pt)

i=1
for (i in 1:length(groups)) {
  print(groups[i])
  Idents(LungFB) <- "Status"
  x <- subset(LungFB, idents = groups[i])
  Idents(x) <- "CoGAPS_clusters"
  pt1 <- as.data.frame(table(Idents(x), x$orig.ident))
  if (nrow(pt) != length(pt1)){
    z <- c(pt$Var1,pt1$Var1)
    a <- z[!(duplicated(z)|duplicated(z, fromLast=TRUE))]
    row <- which(!(duplicated(z)|duplicated(z, fromLast=TRUE)), arr.ind = TRUE)
    pt1<-pt1 %>% add_row(Var1 = a, Var2 = "SeuratProject", Freq = 0, .before = row)
  }
  list <- pt1$Freq
  pt$y <- list
  names(pt)[names(pt) == 'y'] <- groups[i]
  i=i+1
}

pt <- pt %>%adorn_totals("row")
pt_per <-adorn_percentages(pt, "col")

#Relative abundance bar graph
pt_per <- pt_per %>% filter(row_number() <= n()-1)
pt_per <- pt_per[-c(2,3)]

pt_plot <- data.frame()
for (i in 1:length(groups)) {
  x <- pt_per[c("Var1",groups[i])]
  names(x)[names(x) == groups[i]] <- "Frequency"
  names(x)[names(x) == 'Var1'] <- "Cell_Types"
  x$Group <- groups[i]
  pt_plot <- rbind(pt_plot, x)
}

ggplot(pt_plot, aes(fill = Cell_Types, x = Group, y = Frequency))+
  scale_x_discrete(limits = c("Healthy", "EL-PMN","OFF")) +
  geom_bar(position = "stack", stat = "identity")+
  ggtitle("Lung FB Relative Cell Type Abundance")+
  theme(plot.title = element_text(hjust = 0.5))+
  RotatedAxis()+
  theme(axis.text.x=element_text(size=11), axis.title.x = element_blank())+
  scale_fill_manual(values =c("FB1" = "lightskyblue",
                              "FB2" = "plum1",
                              "FB3" = "#7A81FF",
                              "FB4" = "palevioletred",
                              "FB5" = "navy"))

#### Supplementary Figure S2D ####
VlnPlot(LungFB, features = pattern_names, group.by = "CoGAPS_clusters", pt.size = 0.1,c("lightskyblue","plum1","#7A81FF","palevioletred","navy"), ncol = 5)&
  theme(axis.title.x = element_blank())&stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)

#### Figure 2G ####
c5_gene_sets = msigdbr(species = "mouse", category = "C5") 
msigdbr_list = split(x = c5_gene_sets$gene_symbol, f = c5_gene_sets$gs_name)
GOBP_FIBROBLAST_ACTIVATION<- list(msigdbr_list[["GOBP_FIBROBLAST_ACTIVATION"]])

LungFB <- AddModuleScore(object = LungFB, features = GOBP_FIBROBLAST_ACTIVATION, name = "GOBP_FIBROBLAST_ACTIVATION")
FeaturePlot(object = LungFB, features = "GOBP_FIBROBLAST_ACTIVATION1")+
  scale_color_gradientn(colours = c("lightskyblue","palevioletred"))

#### Figure 2H ####
h_gene_sets = msigdbr(species = "mouse", category = "H") 
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
IL6_JAK_STAT3<- list(msigdbr_list[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]])
TNFA_SIGNALING_VIA_NFKB<- list(msigdbr_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]])

LungFB <- AddModuleScore(object = LungFB, features = IL6_JAK_STAT3, name = "Il6_Jak_Stat3")
LungFB <- AddModuleScore(object = LungFB, features = TNFA_SIGNALING_VIA_NFKB, name = "Tnfa_via_Nfkb")

DotPlot(LungFB, c("Pdgfra","Fap","Fn1","Des","Col1a1","Acta2","Il6","Cxcl1","Il6_Jak_Stat31","Tnfa_via_Nfkb1","cogaps_Pattern-1"), dot.scale = 10)+
  scale_color_gradientn(colours = c("lightskyblue","palevioletred"))+
  RotatedAxis()+theme(axis.title = element_blank())+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))


#### Figure 2I ####
#read in gene signatures
hallmark<-buildMSIGDB(species="mouse",keytype = "SYMBOL",anntype = "HALLMARK")

#single sample/cell GSEA
res_H<-scgsva(LungFB,hallmark,assay = "RNA",method="UCell") 

#Fit linear model for each pathway
#dataframe with all pathways - t-test compares means of 2 groups and accounts for variability
#B-stat is a bayesian-derived measurement that reflects the odds that a gene is differentially 
#expressed (if B = 0 there's a 50/50 chance the gene is DE, B>0 increases chance that gene is
#DE and higher in experimental vs reference)
paths_H <- findPathway(res_H,group = "Status", ref = "EL-PMN")

#Significance testing between groups
#Performs wilcox test and Benjamini-Hochberg Procedure for adjusted p-val
#https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/wilcoxon-signed-rank-test/
#https://www.statisticshowto.com/benjamini-hochberg-procedure/
sig_H <- sigPathway(res_H,group = "Status", ref = "EL-PMN")

ridgePlot(res_H,features="IL6_JAK_STAT3_SIGNALING",group_by="Status", color = c("lightgray","#D7A8FF","grey40"))+
  stat_summary(fun = median, geom='point', size = 2, colour = "white", shape = 20)+FontSize(x.text = 11)+
  ggtitle("Lung Fibroblasts - UCell Signature Score", subtitle = "HALLMARK_IL6_JAK_STAT3_SIGNALING")+
  theme(axis.title.y = element_blank())

#### Figure 2J ####
DimPlot(LungMac, label = T,repel = T, cols = c("#F5C3CB","#FFB2D8","#A02B93","#D266A1","#FF9360","#C04F15"))+NoLegend()

#### Figure 2K ####
Idents(LungMac) <- "Status"
groups <-levels(factor(unique(Idents(LungMac))))

Idents(LungMac) <- "Myeloid_clusters"
pt <- table(Idents(LungMac), LungMac$orig.ident)
pt <- as.data.frame(pt)

i=1
for (i in 1:length(groups)) {
  print(groups[i])
  Idents(LungMac) <- "Status"
  x <- subset(LungMac, idents = groups[i])
  Idents(x) <- "Myeloid_clusters"
  pt1 <- as.data.frame(table(Idents(x), x$orig.ident))
  if (nrow(pt) != length(pt1)){
    z <- c(pt$Var1,pt1$Var1)
    a <- z[!(duplicated(z)|duplicated(z, fromLast=TRUE))]
    row <- which(!(duplicated(z)|duplicated(z, fromLast=TRUE)), arr.ind = TRUE)
    pt1<-pt1 %>% add_row(Var1 = a, Var2 = "SeuratProject", Freq = 0, .before = row)
  }
  list <- pt1$Freq
  pt$y <- list
  names(pt)[names(pt) == 'y'] <- groups[i]
  i=i+1
}

pt <- pt %>%adorn_totals("row")
pt_per <-adorn_percentages(pt, "col")

#Relative abundance bar graph
pt_per <- pt_per %>% filter(row_number() <= n()-1)
pt_per <- pt_per[-c(2,3)]

pt_plot <- data.frame()
for (i in 1:length(groups)) {
  x <- pt_per[c("Var1",groups[i])]
  names(x)[names(x) == groups[i]] <- "Frequency"
  names(x)[names(x) == 'Var1'] <- "Cell_Types"
  x$Group <- groups[i]
  pt_plot <- rbind(pt_plot, x)
}

ggplot(pt_plot, aes(fill = Cell_Types, x = Group, y = Frequency))+
  scale_x_discrete(limits = c("Healthy", "EL-PMN","OFF")) +
  geom_bar(position = "stack", stat = "identity")+
  ggtitle("Lung Macrophage & Monocyte Relative Cell Type Abundance")+
  theme(plot.title = element_text(hjust = 0.5))+
  RotatedAxis()+
  theme(axis.text.x=element_text(size=11), axis.title.x = element_blank())+
  scale_fill_manual(values =c(c("Monocyte" = "#F5C3CB",
                                "Classical Monocyte" = "#FFB2D8",
                                "Alveolar Macrophage" = "#A02B93",
                                "Interstitial Macrophage" = "#D266A1",
                                "Macrophage" = "#FF9360",
                                "Activated Macrophage" = "#C04F15")))


#### Figure 2L ####
obj <- subset(LungMac, idents = c("Macrophage","Interstitial Macrophage",
                                  "Alveolar Macrophage","Activated Macrophage"))

#Import gene sets from msigdb. Here, I'm importing all of the mouse Hallmark gene sets. You can change the species and the gene set group based on the package directions:
h_gene_sets = msigdbr(species = "mouse", category = "H") #Hallmark
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

#Run DE with 0 thresholds:
DE <- FindMarkers(obj, ident.1 ="EL-PMN", ident.2 = "OFF", group.by = "Status", logfc.threshold = 0, min.pct = 0)

#Rank your DE by fold change (people rank by adj p value or FC, I like FC):
ranks <- DE$avg_log2FC
names(ranks) <- rownames(DE)

#Run GSEA through fgsea:
fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

#Retrieve GSEA result:
result <- apply(fgseaRes,2,as.character)

#select only pathways with an adjusted p-value <= 0.05
fgseaRes_sig <- subset(fgseaRes, fgseaRes$padj <= 0.05)

waterfall_plot_sig <- function (fgseaRes, graph_title) {
  (fgseaRes %>% 
     mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
     ggplot( aes(reorder(short_name,NES), NES)) +
     geom_bar(stat= "identity", aes(fill = padj))+
     coord_flip()+
     labs(x = "Hallmark Gene Set", y = "Normalized Enrichment Score", title = graph_title)+
     theme(axis.text.y = element_text(size = 10), 
           plot.title = element_text(hjust = 1)))+scale_fill_viridis(direction = 1, option = "B")
}
waterfall_plot_sig(fgseaRes_sig, "Pathways Enriched in EL-PMN vs. OFF - Macrophages")

#### Figure 2M ####
Idents(LungMac) <- "Myeloid_clusters"
DotPlot(LungMac, features = c("Arg1"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#### Figure 2N ####
VlnPlot(LungFB, "Csf1", split.by = "Status", cols = c("lightgray","#D266A1","grey40"))+NoLegend()&
  stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)


#### Supplementary Figure S5D ####
DotPlot(Lung, features = c("Il6"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#### Supplementary Figure S5E ####
Idents(LungFB) <- "Status"
VlnPlot(LungFB,"Il6", c("lightgray","#D7A8FF", "grey40"))+
  NoLegend()+theme(axis.title.x = element_blank(),plot.title = element_text(face = "bold.italic"))

#### Supplementary Figure S5F ####
Stacked_VlnPlot(Lung, features = c("Il6ra","Il6st"),x_lab_rotate = 45, colors_use = c("#D7A8FF","#9990FF","#6A5D99","#FF9360",
                                                                                      "#FFB2D8","#73D279","#FFE18F","#52E7ED",
                                                                                      "#FFB300","#94435F","#009193","#D266A1",
                                                                                      "#7FA4CA","#528AFF","#DF4945","lightgray"))






