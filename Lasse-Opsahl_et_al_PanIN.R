#This script recreates Figures 4I-J, S5A-C
#Raw data files for the novel dataset generated in this manuscript are available through the NIH Gene Expression Omnibus (GEO), accession number GSE292712.
#Control pancreas (GSM9123775, GSM9123776, GSM9123777,GSM9162119, GSM9162120)
#iKRAS 16wk ON pancreas (GSM9123784, GSM9123785, GSM9123786)
#iKRAS 16wk ON + 1wk OFF pancreas (GSM9123782, GSM9123783)

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#version.string R version 4.4.0 (2024-04-24) --- Puppy Cup 
#Seurat Version 5.3.0

#Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(magrittr)
library(devtools)
library(scCustomize)
library(readxl)
library(msigdbr)



#### Pre-processing ----------------------------------------------------------------------------------------------------####
#Read in metadata file
Samples <- read_excel("IntegrateEL-PMNiKRASscRNAseqSamples.xlsx")
Samples$Description <- NULL

#keep only Pancreas samples
SamplesPanc <- subset(Samples, Tissue %in% c("Pancreas"))

samples <- as.character(as.list(SamplesPanc$Sample))
samples_list <- vector("list", length = length(samples))

#Make seurat objects and add metadata
for (x in 1:nrow(SamplesPanc)){
  print(SamplesPanc[x,"SampleID"])
  sc <- Read_CellBender_h5_Mat(as.character(SamplesPanc[x,"file"]))
  samples_list[[x]] <- CreateSeuratObject(sc, min.cells = 3, min.features = 100)
  #Metadata
  metaDataGroups <- colnames(SamplesPanc)
  metaData <- as.character(as.list(SamplesPanc[x,]))
  for (y in 1:length(metaData)){
    group <- metaDataGroups[[y]]
    z <- metaData[[y]]
    samples_list[[x]][[group]] <- z
    samples_list[[x]]@meta.data[["file"]] <- NULL
  }
}
names(samples_list) <- samples

#merge seurat objects
Panc <- Merge_Seurat_List(samples_list, add.cell.ids = samples)

#Normalize
Panc <- NormalizeData(object = Panc, normalization.method = "LogNormalize", scale.factor = 10000)

#Check percent mitochondrial genes to check for doublets and poor cell quality
Idents(Panc) <- "Sample"
Panc[["percent.mt"]] <- PercentageFeatureSet(object = Panc, pattern = "^mt-")

#Keep cells with between 800 and 100,000 RNA reads and less than 15% mitochondrial genes
Panc <- subset(x = Panc, subset = nCount_RNA > 800 & nCount_RNA < 100000 & percent.mt < 15)

#Identify variable genes
Panc <- FindVariableFeatures(Panc, selection.method = "vst", nfeatures = 2000)

#Scale data and run principal component analysis (PCA)  dimensionality reduction
all.genes <- rownames(Panc)
Panc <- ScaleData(Panc, verbose = T, features = all.genes)
Panc <- RunPCA(Panc, npcs = 30, verbose = T)

#Integrate
Panc <- IntegrateLayers(object = Panc, method = RPCAIntegration,
                        orig.reduction = "pca", new.reduction = "integrated.rpca",
                        verbose = T)

#Calculate 90% Variance:
y <- (10:30)
varList <- c()
for (x in 10:30){
  st_dev <- Panc@reductions$pca@stdev
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
Panc <- FindNeighbors(object = Panc, dims = 1:varFinal) #change dims to 1:x (number of PCs)
Panc <- FindClusters(object = Panc, resolution = 2) #change resolution to higher (more clusters) or lower (fewer clusters)
Panc <- RunUMAP(Panc, reduction = "integrated.rpca", dims = 1:varFinal, verbose = F) #change dims to 1:x (number of PCs)

#Join seurat layers
Panc <- JoinLayers(Panc)

#Identify clusters based on marker expression
DotPlot(Panc, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdgfrb","Rgs5","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                         "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Try4","Ctrb1","Kras",
                         "Trp53","EGFP","Sec14l3","Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Cd200r3","Fcer1a","Mcam","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Label Manual Clusters:
Idents(Panc) <- "seurat_clusters"
Panc <- RenameIdents(Panc, 
                     "0" = "B Cell", 
                     "1" = "B Cell", 
                     "2" = "Acinar", 
                     "3" = "Acinar", 
                     "4" = "CD4+ T Cell", 
                     "5" = "CD4+ T Cell", 
                     "6" = "Treg", 
                     "7" = "CD8+ T Cell", 
                     "8" = "CD8+ T Cell",
                     "9" = "B Cell", 
                     "10" = "CD8+ T Cell", 
                     "11" = "Fibroblast", 
                     "12" = "CD4+ T Cell", 
                     "13" = "Endothelial", 
                     "14" = "ADM", 
                     "15" = "Fibroblast", 
                     "16" = "Acinar", 
                     "17" = "Macrophage",
                     "18" = "CD4+ T Cell", 
                     "19" = "Treg", 
                     "20" = "Macrophage", 
                     "21" = "CD8+ T Cell",
                     "22" = "B Cell", 
                     "23" = "ADM", 
                     "24" = "ADM", 
                     "25" = "B Cell", 
                     "26" = "Neutrophil", 
                     "27" = "Ductal", 
                     "28" = "Plasma Cell",
                     "29" = "NK Cell",
                     "30" = "Fibroblast",
                     "31" = "Endothelial",
                     "32" = "Proliferating B Cell",
                     "33" = "Active DC",
                     "34" = "Proliferating T Cell", 
                     "35" = "B Cell", 
                     "36" = "pDC",
                     "37" = "B Cell",
                     "38" = "cDC1",
                     "39" = "Active DC", 
                     "40" = "RBC",
                     "41" = "Fibroblast", 
                     "42" = "Proliferating Myeloid", 
                     "43" = "mDC",
                     "44" = "Macrophage",
                     "45" = "NK Cell",
                     "46" = "Pericyte",
                     "47" = "Mesothelial",
                     "48" = "Pericyte",
                     "49" = "Treg")

Panc[["manual_clusters"]] <- Panc@active.ident

#Simple Clusters
Idents(Panc) <- "manual_clusters"
Panc <- RenameIdents(Panc,
                     "Ductal" = "Ductal",
                     "Acinar" = "Acinar",
                     "ADM" = "ADM",
                     "Endothelial" = "Endothelial",
                     "Mesothelial" = "Mesothelial",
                     "Fibroblast" = "Fibroblast",
                     "Macrophage" = "Macrophage",
                     "Active DC" = "DC",
                     "cDC1" = "DC",
                     "pDC" = "DC",
                     "mDC" = "DC",
                     "Neutrophil" = "Granulocyte",
                     "CD8+ T Cell" = "T Cell",
                     "Treg" = "T Cell",
                     "CD4+ T Cell" = "T Cell",
                     "NK Cell" = "NK Cell",
                     "B Cell" = "B Cell",
                     "Plasma Cell" = "Plasma Cell",
                     "Proliferating B Cell" = "Proliferating",
                     "Proliferating T Cell" = "Proliferating",
                     "Proliferating Myeloid" = "Proliferating",
                     "Pericyte" = "Pericyte",
                     "RBC" = "RBC")
Panc[["simple_clusters"]] <- Panc@active.ident

#Re-order clusters
Idents(Panc) <- "simple_clusters"

new_order <- c("Fibroblast",
               "Pericyte",
               "Endothelial",
               "Macrophage",
               "Granulocyte",
               "DC",
               "T Cell",
               "NK Cell",
               "B Cell",
               "Plasma Cell",
               "Ductal",
               "ADM",
               "Acinar",
               "Mesothelial",
               "Proliferating",
               "RBC")

Panc@active.ident <- factor(Panc@active.ident, levels = new_order)
Panc[["simple_clusters"]] <- Panc@active.ident

#Re-order groups
Panc$Status <- factor(Panc@meta.data[["Status"]], levels = c("Healthy","PanIN","OFF"))
DimPlot(Panc, split.by = "Status")

Idents(Panc) <- "simple_clusters"
PancFB <- subset(Panc,idents = "Fibroblast")



#### Figures ----------------------------------------------------------------------------------------------------####
####Figure 4I ####
DimPlot(Panc, cols = c("#D7A8FF","#9990FF","#6A5D99","#FF9360",
                       "#73D279","#FFE18F","#52E7ED","#FFB300",
                       "#94435F","#009193","#DF4945","#D266A1",
                       "#DAA536","#7FA4CA","#528AFF","lightgray"), raster = F, split.by = "Status")

#### Supplementary Figure S5A ####
Idents(Panc) <- "simple_clusters"
DotPlot(Panc, features=c("Pdgfra", "Clec3b", "Col1a1","Col1a2", #fibroblast
                         "Pdgfrb","Rgs5", #Pericyte
                         "Pecam1", "Cdh5","Flt1", #Endothelial
                         "Ptprc","H2-Eb1","Itgam","Cd68","Adgre1","Mrc1", #macrophage
                         "Fcgr3","Ly6g", #granulocyte
                         "Flt3", "Batf3", #DC
                         "Cd3e", "Cd4","Cd8a", #T Cell
                         "Nkg7","Klrb1c", #NK Cell
                         "Cd79a", "Cd19", "Ms4a1", #B Cell
                         "Igha","Ighm", #Plasma cell
                         "Krt8", "Krt18", #Epithelial
                         "Try4","Ctrb1", #Acinar
                         "Msln", "Wt1", #Mesothelial
                         "Ccna2", "Mki67", #Proliferating
                         "Hbb-bt" #RBC
), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

#### Supplementary Figure S5B ####
Idents(Panc) <- "Sample"
groups <-levels(factor(unique(Idents(Panc))))

Idents(Panc) <- "simple_clusters"
pt <- table(Idents(Panc), Panc$orig.ident)
pt <- as.data.frame(pt)

i=1
for (i in 1:length(groups)) {
  print(groups[i])
  Idents(Panc) <- "Sample"
  x <- subset(Panc, idents = groups[i])
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

#if out of order
pt<-pt[match(new_order, pt$Var1),]

#save table
pt <- pt %>%adorn_totals("row")
write.csv(pt, "Panc_CellAbundance_AllCells_perSample_5.20.25.csv", row.names = T)

pt_per <-adorn_percentages(pt, "col")
write.csv(pt_per, "Panc_CellPercentages_AllCells_perSample_5.20.25.csv", row.names = T)



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
  scale_x_discrete(limits = c("Healthy_16wk_ON_pancreas_1","Healthy_16wk_ON_pancreas_2","Healthy_16wk_ON_pancreas_3",
                              "Healthy_20wk_ON_pancreas_1","Healthy_20wk_ON_pancreas_2",
                              "iKras_16wk_ON_pancreas_1","iKras_16wk_ON_pancreas_2","iKras_16wk_ON_pancreas_3",
                              "iKras_16wk_ON_1wk_OFF_pancreas_1","iKras_16wk_ON_1wk_OFF_pancreas_2")) +
  geom_bar(position = "stack", stat = "identity")+
  ggtitle("Pancreas Relative Cell Type Abundance")+
  theme(plot.title = element_text(hjust = 0.5))+
  RotatedAxis()+
  theme(axis.text.x=element_text(size=11), axis.title.x = element_blank())+
  scale_fill_manual(values = c("Fibroblast" = "#D7A8FF",
                               "Pericyte" = "#9990FF",
                               "Endothelial" = "#6A5D99",
                               "Macrophage" = "#FF9360",
                               "Granulocyte" = "#73D279",
                               "DC" = "#FFE18F",
                               "T Cell" = "#52E7ED",
                               "NK Cell" = "#FFB300",
                               "B Cell" = "#94435F",
                               "Plasma Cell" = "#009193",
                               "Ductal" = "#DF4945",
                               "ADM" = "#D266A1",
                               "Acinar" = "#DAA536",
                               "Mesothelial" = "#7FA4CA",
                               "Proliferating" = "#528AFF",
                               "RBC" = "lightgray"))

#### Supplementary Figure S5C ####
DotPlot(Panc, features = c("Il6"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

####Figure 4J ####
Idents(PancFB) <- "Status"
VlnPlot(PancFB,"Il6", c("lightgray","#D7A8FF", "grey40"))+
  NoLegend()+theme(axis.title.x = element_blank(),plot.title = element_text(face = "bold.italic"))
