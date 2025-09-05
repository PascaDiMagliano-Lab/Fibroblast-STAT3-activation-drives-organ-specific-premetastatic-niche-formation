#This script recreates Figures 5F-I, S7A-E
#Raw data files for the novel dataset generated in this manuscript are available through the NIH Gene Expression Omnibus (GEO), accession number GSE292712.
#Control PMNLung (iKRAS model: GSM8864044, GSM8864046, GSM8864048, GSM8864052,	GSM8864056; KPC/KC model: GSM9123778, GSM9123779)
#PMN PMNLung (iKRAS model: GSM8864051; KPC/KC model: GSM9123787, GSM9123788, GSM9123789)
#Metastatic PMNLung (iKRAS model: GSM8864055)

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#version.string R version 4.4.0 (2024-04-24) --- Puppy Cup 
#Seurat Version 5.3.0

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
library(infercnv)
library(ggpubr)


#### Pre-processing ----------------------------------------------------------------------------------------------------####

#Read in metadata file
Samples <- read_excel("IntegratePMN_Met_scRNAseqSamples.xlsx")
Samples$Description <- NULL

#keep only PMNLung samples
Samples<- subset(Samples, Tissue %in% c("Lung"))


samples <- as.character(as.list(SamplesPMNLung$Sample))
samples_list <- vector("list", length = length(samples))

#Make seurat objects and add metadata
for (x in 1:nrow(Samples)){
  print(Samples[x,"SampleID"])
  sc <- Read_CellBender_h5_Mat(as.character(Samples[x,"file"]))
  samples_list[[x]] <- CreateSeuratObject(sc, min.cells = 3, min.features = 100)
  #Metadata
  metaDataGroups <- colnames(Samples)
  metaData <- as.character(as.list(Samples[x,]))
  for (y in 1:length(metaData)){
    group <- metaDataGroups[[y]]
    z <- metaData[[y]]
    samples_list[[x]][[group]] <- z
    samples_list[[x]]@meta.data[["file"]] <- NULL
  }
}
names(samples_list) <- samples

#merge seurat objects
PMNLung <- Merge_Seurat_List(samples_list, add.cell.ids = samples)

#Normalize
PMNLung <- NormalizeData(object = PMNLung, normalization.method = "LogNormalize", scale.factor = 10000)

#Check percent mitochondrial genes to check for doublets and poor cell quality
Idents(PMNLung) <- "Sample"
PMNLung[["percent.mt"]] <- PercentageFeatureSet(object = PMNLung, pattern = "^mt-")

#Keep cells with between 800 and 100,000 RNA reads and less than 15% mitochondrial genes
PMNLung <- subset(x = PMNLung, subset = nCount_RNA > 800 & nCount_RNA < 100000 & percent.mt < 15)

#Identify variable genes
PMNLung <- FindVariableFeatures(PMNLung, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(PMNLung)
PMNLung <- ScaleData(PMNLung, verbose = T, features = all.genes)
PMNLung <- RunPCA(PMNLung, npcs = 30, verbose = T)

#Integrate
PMNLung <- IntegrateLayers(object = PMNLung, method = RPCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- PMNLung@reductions$pca@stdev
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
PMNLung <- FindNeighbors(object = PMNLung, dims = 1:varFinal) #change dims to 1:x (number of PCs)
PMNLung <- FindClusters(object = PMNLung, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
PMNLung <- RunUMAP(PMNLung, reduction = "integrated.rpca", dims = 1:varFinal, verbose = F) #change dims to 1:x (number of PCs)

#Join seurat layers
PMNLung <- JoinLayers(PMNLung)

#Identify clusters based on marker expression
DotPlot(PMNLung, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdgfrb","Rgs5","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                         "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Try4","Ctrb1","Kras",
                         "Trp53","EGFP","Sec14l3","Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Cd200r3","Fcer1a","Mcam","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Label Manual Clusters:
Idents(PMNLung) <- "seurat_clusters"
PMNLung <- RenameIdents(PMNLung, 
                     "0" = "T Cell", 
                     "1" = "B Cell", 
                     "2" = "Neutrophil", 
                     "3" = "B Cell", 
                     "4" = "B Cell", 
                     "5" = "CD4+ T Cell", 
                     "6" = "Neutrophil", 
                     "7" = "Endothelial", 
                     "8" = "Endothelial",
                     "9" = "NK Cell", 
                     "10" = "B Cell", 
                     "11" = "Treg", 
                     "12" = "Macrophage", 
                     "13" = "Macrophage", 
                     "14" = "CD8+ T Cell", 
                     "15" = "Fibroblast", 
                     "16" = "CD8+ T Cell", 
                     "17" = "Endothelial",
                     "18" = "Alveolar Epithelial", 
                     "19" = "Endothelial", 
                     "20" = "Alveolar Epithelial", 
                     "21" = "Fibroblast",
                     "22" = "Macrophage", 
                     "23" = "RBC", 
                     "24" = "Alveolar Macrophage", 
                     "25" = "Endothelial", 
                     "26" = "Plasma Cell", 
                     "27" = "NK Cell", 
                     "28" = "Fibroblast",
                     "29" = "Alveolar Macrophage",
                     "30" = "Fibroblast",
                     "31" = "Pericyte",
                     "32" = "B Cell",
                     "33" = "Macrophage",
                     "34" = "cDC1", 
                     "35" = "Pericyte", 
                     "36" = "Mesothelial",
                     "37" = "Endothelial",
                     "38" = "Bronchial Epithelial",
                     "39" = "Pericyte", 
                     "40" = "B Cell",
                     "41" = "Proliferating T Cell", 
                     "42" = "Active DC", 
                     "43" = "Endothelial",
                     "44" = "Fibroblast",
                     "45" = "Alveolar Epithelial",
                     "46" = "Neutrophil",
                     "47" = "Endothelial",
                     "48" = "Endothelial",
                     "49" = "Monocyte",
                     "50" = "B Cell",
                     "51" = "Macrophage",
                     "52" = "Neutrophil",
                     "53" = "B Cell",
                     "54" = "Bronchial Epithelial", 
                     "55" = "Alveolar Epithelial", 
                     "56" = "Basophil", 
                     "57" = "Tumor", 
                     "58" = "Endothelial",
                     "59" = "Proliferating Myeloid",
                     "60" = "Macrophage", 
                     "61" = "B Cell", 
                     "62" = "RBC", 
                     "63" = "pDC", 
                     "64" = "Endothelial", 
                     "65" = "Endothelial", 
                     "66" = "Neutrophil", 
                     "67" = "Alveolar Macrophage", 
                     "68" = "Proliferating Alveolar Macrophage",
                     "69" = "NK Cell", 
                     "70" = "Endothelial", 
                     "71" = "Macrophage", 
                     "72" = "Endothelial", 
                     "73" = "Fibroblast", 
                     "74" = "RBC", 
                     "75" = "Fibroblast", 
                     "76" = "Neutrophil", 
                     "77" = "Alveolar Epithelial", 
                     "78" = "Cardiomyocyte",
                     "79" = "Endothelial")

PMNLung[["manual_clusters"]] <- PMNLung@active.ident

#Simple Clusters
Idents(PMNLung) <- "manual_clusters"
PMNLung <- RenameIdents(PMNLung,"Alveolar Epithelial" = "Epithelial",
                     "Bronchial Epithelial" = "Epithelial",
                     "Tumor" = "Epithelial",
                     "Endothelial" = "Endothelial",
                     "Mesothelial" = "Mesothelial",
                     "Fibroblast" = "Fibroblast",
                     "Macrophage" = "Macrophage",
                     "Alveolar Macrophage" = "Macrophage",
                     "Active DC" = "DC",
                     "cDC1" = "DC",
                     "pDC" = "DC",
                     "Monocyte" = "Monocyte",
                     "Basophil" = "Granulocyte",
                     "Neutrophil" = "Granulocyte",
                     "CD8+ T Cell" = "T Cell",
                     "Treg" = "T Cell",
                     "CD4+ T Cell" = "T Cell",
                     "T Cell" = "T Cell",
                     "NK Cell" = "NK Cell",
                     "B Cell" = "B Cell",
                     "Proliferating Myeloid" = "Proliferating",
                     "Proliferating T Cell" = "Proliferating",
                     "Proliferating Alveolar Macrophage" = "Proliferating",
                     "Pericyte" = "Pericyte",
                     "RBC" = "RBC",
                     "Cardiomyocyte" = "Cardiomyocyte")
PMNLung[["simple_clusters"]] <- PMNLung@active.ident

#Re-order clusters
Idents(PMNLung) <- "simple_clusters"

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

PMNLung@active.ident <- factor(PMNLung@active.ident, levels = new_order)
PMNLung[["simple_clusters"]] <- PMNLung@active.ident

#Re-order groups
PMNLung$Status <- factor(PMNLung@meta.data[["Status"]], levels = c("Healthy","PMN","Metastatic"))


#### Pre-processing Fibroblasts ----------------------------------------------------------------------------------------------------####

#subset out the fibroblasts from the larger object
PMNLungFB <- subset(Lung, idents = "Fibroblast")

#split seurat layers by sample
PMNLungFB[["RNA"]] <- split(PMNLungFB[["RNA"]], f = PMNLungFB$Sample)

#Identify variable genes
PMNLungFB <- FindVariableFeatures(PMNLungFB, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(PMNLungFB)
PMNLungFB <- ScaleData(PMNLungFB, verbose = T, features = all.genes)
PMNLungFB <- RunPCA(PMNLungFB, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- PMNLungFB@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#16

#Find Neighbors and cluster cells
PMNLungFB <- FindNeighbors(object = PMNLungFB, dims = 1:varFinal) #change dims to 1:x (number of PCs)
PMNLungFB <- FindClusters(object = PMNLungFB, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
PMNLungFB <- RunUMAP(PMNLungFB, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(PMNLungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-fibroblasts from the object
PMNLungFB <- subset(PMNLungFB, idents = c("15","17","19","20","21","22","24","27","28"), invert = T) #remove contamination

#Repeat pipeline after removing contaminating cells (Part 2)

#Identify variable genes
PMNLungFB <- FindVariableFeatures(PMNLungFB, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(PMNLungFB)
PMNLungFB <- ScaleData(PMNLungFB, verbose = T, features = all.genes)
PMNLungFB <- RunPCA(PMNLungFB, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- PMNLungFB@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#16

#Find Neighbors and cluster cells
PMNLungFB <- FindNeighbors(object = PMNLungFB, dims = 1:varFinal) #change dims to 1:x (number of PCs)
PMNLungFB <- FindClusters(object = PMNLungFB, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
PMNLungFB <- RunUMAP(PMNLungFB, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(PMNLungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-fibroblasts from the object
PMNLungFB <- subset(PMNLungFB, idents = c("13","21","22","26"), invert = T) #remove contamination

#Repeat pipeline after removing contaminating cells (Part 3)

#Identify variable genes
PMNLungFB <- FindVariableFeatures(PMNLungFB, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(PMNLungFB)
PMNLungFB <- ScaleData(PMNLungFB, verbose = T, features = all.genes)
PMNLungFB <- RunPCA(PMNLungFB, npcs = 30, verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- PMNLungFB@reductions$pca@stdev
  var <- st_dev^2
  xVar <- sum(var[1:x])/ sum(var)
  yVar <- abs(0.9-xVar)
  varList <-append(varList,yVar)
}
calcVar <- varList[which.min(abs(varList))]
position <- which(varList == calcVar) 
varFinal <- y[position]

varFinal 
#18

#Find Neighbors and cluster cells
PMNLungFB <- FindNeighbors(object = PMNLungFB, dims = 1:varFinal) #change dims to 1:x (number of PCs)
PMNLungFB <- FindClusters(object = PMNLungFB, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
PMNLungFB <- RunUMAP(PMNLungFB, reduction = "integrated.rpca", dims = 1:varFinal, verbose = T) #change dims to 1:x (number of PCs)

#Check for other contaminating cell types
DotPlot(PMNLungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Vim","Des","Pdgfrb","Rgs5","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

##Join seurat layers
PMNLungFB <- JoinLayers(PMNLungFB)


#### Read in top genes from EL-PMN lung fibroblast CoGAPS pattern 1 ----------------------------------------------------------------------------------------------------####
LungFB_cogaps_result <- readRDS("/LungFB_cogaps_result.Rds")

patterns<-patternMarkers(LungFB_cogaps_result)

TopPat1 <- list(patterns[["PatternMarkers"]][["Pattern_1"]])

PMNLungFB <- AddModuleScore(object = PMNLungFB, features = TopPat1, name = "TopPat")

#### Figures ----------------------------------------------------------------------------------------------------####

####Figure 5F ####
DimPlot(Lung, cols = c("#D7A8FF","#9990FF","#6A5D99","#FF9360",
                       "#FFB2D8","#73D279","#FFE18F","#52E7ED",
                       "#FFB300","#94435F","#009193","#D266A1",
                       "#7FA4CA","#528AFF","#DF4945","lightgray"), raster = F, split.by = "Status")

#### Supplementary Figure S7A ####
Idents(PMNLung) <- "simple_clusters"
DotPlot(PMNLung,features=c("Pdgfra", "Clec3b", "Col1a1","Col1a2", #fibroblast
                        "Pdgfrb","Rgs5", #Pericyte
                        "Pecam1", "Cdh5","Flt1", #Endothelial
                        "Ptprc","H2-Eb1","Itgam","Cd68","Adgre1","Mrc1", #macrophage
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

#### Supplementary Figure S7B ####
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
write.csv(pt, "Lung_CellAbundance_AllCells_perSample_9.3.25.csv", row.names = T)

pt_per <-adorn_percentages(pt, "col")
write.csv(pt_per, "Lung_CellPercentages_AllCells_perSample_9.3.25.csv", row.names = T)



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
                              "K;Hif1a_Healthy_14m_lung","KP_Healthy_4m_lung",
                              "iKras;p53_20wk_ON_lung_PDA","KPC_4m_lung_1","KPC_8m_lung","KC_17.5m_lung",
                              "iKras;p53_20wk_ON_lung_Met")) +
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

#### Supplementary Figure S7C ####
LungEpi <- subset(Lung, idents = c("Epithelial"))

#access counts data from seurat object
counts_matrix = GetAssayData(LungEpi, slot="counts")

cell_annotations <- data.frame(cell_id = colnames(counts_matrix),
                               cell_type = LungEpi$Status) # Replace 'cell_type_column' with your relevant column

write.table(cell_annotations, "LungEpi_cell_annotations.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


download.file(file.path("https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions"), 
              "mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions")


# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="LungEpi_cell_annotations.txt",
                                    delim="\t",
                                    gene_order_file="mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions",
                                    ref_group_names=c("Healthy"))

options(scipen = 100)

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T)

#### Supplementary Figure S7D ####

#take the columns from the processed expression matrix
  #count the number of of genes per cell that have a modified expression
  #divide by the number of genes
scores=apply(infercnv_obj@expr.data,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })

df <- as.data.frame(scores)
df <- t(df)

LungEpi[["score"]] <- CreateAssayObject(counts = df)

Idents(LungEpi) <- "score"
VlnPlot(LungEpi, "scores", group.by = "Status", cols = c("lightgray","#D29EBB","#D266A1"))+NoLegend()+theme(axis.title.x = element_blank())+
  stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)+
  stat_compare_means(method = "anova")& ggtitle(expression(paste(bold("inferCNV Score"))))

#### Supplementary Figure S7E ####
FeaturePlot(LungEpi, "EGFP", split.by = "Status", cols = c("lightgray","#D266A1"), pt.size =1)+ theme(legend.position = "right")

####Figure 5G ####
h_gene_sets = msigdbr(species = "mouse", category = "H") 
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
IL6_JAK_STAT3<- list(msigdbr_list[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]])
TNFA_SIGNALING_VIA_NFKB<- list(msigdbr_list[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]])

PMNLungFB <- AddModuleScore(object = PMNLungFB, features = IL6_JAK_STAT3, name = "Il6_Jak_Stat3")
PMNLungFB <- AddModuleScore(object = PMNLungFB, features = TNFA_SIGNALING_VIA_NFKB, name = "Tnfa_via_Nfkb")

DotPlot(PMNLungFB, c("Pdgfra","Fap","Fn1","Des","Col1a1","Acta2","Il6","Cxcl1","Il6_Jak_Stat31","Tnfa_via_Nfkb1","TopPat1"), split.by = "Status", dot.scale = 10)+
  scale_color_gradientn(colours = c("lightskyblue","palevioletred"))+
  RotatedAxis()+theme(axis.title = element_blank())+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))

####Figure 5H ####
Idents(PMNLungFB) <- "Status"
VlnPlot(PMNLungFB, "TopPat1", cols = c("lightgray","#D7A8FF", "#78206E"))+NoLegend()+
  stat_summary(fun = median, geom='point', size = 15, colour = "white", shape = 95)+
  theme(axis.title.x = element_blank(),plot.title = element_text(face = "bold.italic"))

####Figure 5I ####
hallmark<-buildMSIGDB(species="mouse",keytype = "SYMBOL",anntype = "HALLMARK")

#single sample/cell GSEA
res_H<-scgsva(PMNLungFB,hallmark,assay = "RNA",method="UCell") 

#Fit linear model for each pathway
#dataframe with all pathways - t-test compares means of 2 groups and accounts for variability
#B-stat is a bayesian-derived measurement that reflects the odds that a gene is differentially 
#expressed (if B = 0 there's a 50/50 chance the gene is DE, B>0 increases chance that gene is
#DE and higher in experimental vs reference)
paths_H <- findPathway(res_H,group = "Status", ref = "Healthy")

#Significance testing between groups
#Performs wilcox test and Benjamini-Hochberg Procedure for adjusted p-val
#https://www.statisticshowto.com/probability-and-statistics/statistics-definitions/wilcoxon-signed-rank-test/
#https://www.statisticshowto.com/benjamini-hochberg-procedure/
sig_H <- sigPathway(res_H,group = "Status", ref = "Healthy")

ridgePlot(res_H,features="IL6_JAK_STAT3_SIGNALING",group_by="Status", color = c("lightgray","#D7A8FF","#78206E"))+
  stat_summary(fun = median, geom='point', size = 2, colour = "white", shape = 20)+FontSize(x.text = 11)+
  ggtitle("Lung Fibroblasts - UCell Signature Score", subtitle = "HALLMARK_IL6_JAK_STAT3_SIGNALING")+
  theme(axis.title.y = element_blank())


