setwd("C:\\Work\\set6")
library("DESeq2")
library(loadings)

#seurat and accessory libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(ggrepel)

#atlas counts and info
atCounts <- read.csv("SS2_counts.csv", header=T, row.names=1)
#altered info file, less columns and different order - possible batch stuff to the end
atInfo <- read.table("SS2_pheno2.txt", header = T, sep = "\t", quote = "", row.names=1)


######################read in and parse bulk data

set6 <- read.delim("set6.count.trim.tsv")

#set gene IDs to rownames
geneID <- set6[,1]
set6 <- set6[,2:25]
rownames(set6) <- geneID

#cut out the paths
x <- strsplit(colnames(set6), "RNA.")
newNames <- character()
for(i in 1:24){
newNames[i] <- x[[i]][2]
}
x <- strsplit(newNames, ".sorted.bam")
newNames <- unlist(x)
colnames(set6) <- newNames

#fill info table using file names
cond <- c(rep("KO", 12), rep("WT", 12))
stage <- character()
for(i in 1:24){
stage[i] <- strsplit(colnames(set6)[i], "\\.")[[1]][2]
}
condStage <- factor(paste0(cond, "_", stage), levels=c("KO_G", "KO_8", "KO_16", "KO_24", "WT_G", "WT_8", "WT_16", "WT_24"))
meta <- data.frame(cond, stage, condStage)
rownames(meta) <- colnames(set6)

#sensible order
so <- c("1614.G.1_S2", "1614.G.2_S4", "1614.G.3_S6", "1614.8.2_S14", "1614.8.3_S16", "1614.8.5_S18", "1614.16.1_S2", "1614.16.2_S6", "1614.16.3_S10", "1614.24.1_S4", "1614.24.2_S8", "1614.24.3_S12", "507.G.1_S1", "507.G.2", "507.G.3_S5", "507.8.2_S13", "507.8.3_S15", "507.8.5_S17", "507.16.1_S1", "507.16.2_S5", "507.16.3_S9", "507.24.1_S3", "507.24.2_S7", "507.24.3_S11")
meta <- meta[so,]
set6 <- set6[,so]
condStage <- meta$condStage
condStageRep <- c("KO.G.1", "KO.G.2", "KO.G.3", "KO.8.2", "KO.8.3", "KO.8.5", "KO.16.1", "KO.16.2", "KO.16.3", "KO.24.1", "KO.24.2", "KO.24.3", "WT.G.1", "WT.G.2", "WT.G.3", "WT.8.2", "WT.8.3", "WT.8.5", "WT.16.1", "WT.16.2", "WT.16.3", "WT.24.1", "WT.24.2", "WT.24.3")
meta <- cbind(meta, condStageRep) 


setwd("C:\\Work\\set6\\atlas")



###########################################FILTER COUNT DATA (as before in main script but not keeping the normalisation)
#DESEQ2 for normalisation and maybe for PCA, suggested by Dario
#start with 5254 rows
dds <- DESeqDataSetFromMatrix(countData = set6, colData = meta, design = ~ cond + stage)
#keep those genes with at least 10 counts in at least 2 samples (example was 3 but this is subgroup size-1)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
#5018 now
dds <- estimateSizeFactors(dds)
#1614.G.3_S6 seems small, only 0.56 size factor to 2dp
sizeFactors(dds)
normSet6 <- varianceStabilizingTransformation(dds, blind=F)

#filter by variance here
ns6Data <- assay(normSet6)
nrow(ns6Data)
vars <- character()
varCount <- 0
varVal <- numeric()
for(i in 1:nrow(ns6Data)){
if(!is.na(var(ns6Data[i,]))){
if(var(ns6Data[i,]) > 0.05)
{
varCount <- varCount+1
vars[varCount] <- rownames(ns6Data)[i]
varVal[varCount] <- var(ns6Data[i,])
}
}
}
#5003 for PCA then?
varCount

#ordered by variance, decreasing
names(varVal) <- vars
varVal <- varVal[order(varVal, decreasing=T)]
vars <- names(varVal)

set6 <- set6[vars,]



##############################analysis - data already QC'd

#match on genes and keep only those in both
inBoth <- intersect(rownames(set6), rownames(atCounts))
set6 <- set6[inBoth,]
atCounts <- atCounts[inBoth,]

#split into datasets as mentioned in paper
idc <- atCounts[,rownames(atInfo)[which(atInfo$ShortenedLifeStage3=="Trophozoite" | atInfo$ShortenedLifeStage3=="Schizont" | atInfo$ShortenedLifeStage3=="Ring" | atInfo$ShortenedLifeStage3=="Merozoite")]]
liver <- atCounts[,rownames(atInfo)[which(atInfo$ShortenedLifeStage3=="EEF")]]
gameto <- atCounts[,rownames(atInfo)[which(atInfo$ShortenedLifeStage3=="Female" | atInfo$ShortenedLifeStage3=="Male")]]
ookoo <- atCounts[,rownames(atInfo)[which(atInfo$ShortenedLifeStage3=="Ookinete" | atInfo$ShortenedLifeStage3=="Oocyst")]]
sporo <- atCounts[,rownames(atInfo)[which(atInfo$ShortenedLifeStage3=="bbSpz" | atInfo$ShortenedLifeStage3=="sgSpz")]]

setOok <- set6[,c(4:12,16:24)] 
setGam <- set6[,c(1:3,13:15)]

ookoo <- cbind(ookoo, set6[,c(4:12,16:24)])
gameto <- cbind(gameto, set6[,c(1:3,13:15)])

#reading into Seurat objects
idc <- CreateSeuratObject(counts = idc, project = "idc")
liver <- CreateSeuratObject(counts = liver, project = "liver")
gameto <- CreateSeuratObject(counts = gameto, project = "gameto")
ookoo <- CreateSeuratObject(counts = ookoo, project = "ookoo")
sporo <- CreateSeuratObject(counts = sporo, project = "sporo")

projList <- list(idc, liver, gameto, ookoo, sporo)

# normalize and identify variable features for each dataset independently
projList <- lapply(X = projList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = projList)
#use to integrate
anchors <- FindIntegrationAnchors(object.list = projList, anchor.features = features)
# this command creates an 'integrated' data assay
intData <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(intData) <- "integrated"

intData <- ScaleData(intData, verbose = FALSE)
intData <- RunPCA(intData, verbose = FALSE)
#first 12
png("PCA elbow.png")
ElbowPlot(intData, ndims=50)
dev.off()

intData <- RunUMAP(intData, dims = 1:12, verbose = FALSE)
intData <- FindNeighbors(intData, dims = 1:12, verbose = FALSE)
intData <- FindClusters(intData, verbose = FALSE)

#plot QC metrics
png("intFeature.png", width=1000, height=1000)
FeaturePlot(intData, features=c("nFeature_RNA"))
dev.off()
png("intCount.png", width=1000, height=1000)
FeaturePlot(intData, features=c("nCount_RNA"))
dev.off()

slf3 <- atInfo$ShortenedLifeStage3
names(slf3) <- rownames(atInfo)
ShortLifeStage3 <- slf3[rownames(intData@meta.data)]

intData <- AddMetaData(intData, ShortLifeStage3, col.name="ShortLifeStage3")


class1 <- intData@meta.data$ShortLifeStage3
class1[1138:1143] <- "SET6-GAM"
class1[1537:1554] <- "SET6-OOK"

class2 <- intData@meta.data$ShortLifeStage3
class2[1138:1140] <- "SET6-KO.GAM"
class2[1141:1143] <- "SET6-WT.GAM"
class2[1537:1539] <- "SET6-KO.OOK08"
class2[1540:1542] <- "SET6-KO.OOK16"
class2[1543:1545] <- "SET6-KO.OOK24"
class2[1546:1548] <- "SET6-WT.OOK08"
class2[1549:1551] <- "SET6-WT.OOK16"
class2[1552:1554] <- "SET6-WT.OOK24"

intData <- AddMetaData(intData, class1, col.name="class1")
intData <- AddMetaData(intData, class2, col.name="class2")


# Visualization of integration
png("integAll.png", width=1000, height=1000)
DimPlot(intData, reduction = "umap", group.by = "class1", label=T)
dev.off()

png("integAll2.png", width=1000, height=1000)
DimPlot(intData, reduction = "umap", group.by = "class2", label=T)
dev.off()










#SUBSET
Idents(intData) <- "class1"
setStuff <- subset(intData, idents=c("SET6-GAM","SET6-OOK"))
setStuff2 <- subset(intData, idents=c("Ookinete","Male","Female","SET6-GAM","SET6-OOK"))
setStuff3 <- subset(intData, idents=c("Ookinete","Male","Female","SET6-GAM","SET6-OOK","Ring","Oocyst"))

Idents(intData) <- "class2"
setStuff1 <- subset(intData, idents=c("SET6-KO.GAM","SET6-WT.GAM","SET6-KO.OOK08","SET6-KO.OOK16","SET6-KO.OOK24","SET6-WT.OOK08","SET6-WT.OOK16","SET6-WT.OOK24"))

Idents(intData) <- "orig.ident"

#more integration graphs - each pairwise comparison
png("gamook.png", width=1000, height=1000)
DimPlot(setStuff, reduction = "umap", group.by = "class1")
dev.off()
png("gamook1.png", width=1000, height=1000)
DimPlot(setStuff1, reduction = "umap", group.by = "class2", pt.size=3)
dev.off()

pointSizes <- rep(1, 366)
pointSizes[c(202:207,349:366)] <- 2

png("gamookAt.png", width=1000, height=1000)
DimPlot(setStuff2, reduction = "umap", group.by = "class1", pt.size=pointSizes)
dev.off()
png("gamookrinooc.png", width=1000, height=1000)
DimPlot(setStuff3, reduction = "umap", group.by = "class1")
dev.off()




#save...
save.image(file = "AtlasWork.RData")







###########################################PCA

ns6Data <- GetAssayData(intData, "scale.data")
pcNum <-12

set6PCA <- prcomp(t(ns6Data), scale=T)

#sort(set6PCA$rotation[,1])
loadings <- set6PCA$rotation %*% diag(set6PCA$sdev)
colnames(loadings) <- colnames(set6PCA$x)
write.table(loadings[,1:pcNum], file="loadings.txt", quote=F, sep="\t")

###########################################P VALUES FOR LOADINGS

###########LEAKY FUNCTION OR REALLY TOO BIG? "Error: cannot allocate vector of size 30.9 Gb"
#pcaL <- pca_loading(set6PCA)



###########################################PCA PLOTS


pcNumH <- pcNum/2
subCols <- c("#EE6447", "#4EB3D2", "#D62F1F", "#B30000", "#7E0000", "#2B8BBD", "#0768AC", "#074081")

class3 <- intData@meta.data$class2
class3[c(1138:1140, 1537:1545)] <- "KO"
class3[c(1141:1143, 1546:1554)] <- "WT"
intData <- AddMetaData(intData, class3, col.name="class3")

class4 <- intData@meta.data$class2
class4[c(1138:1143, 1537:1554)] <- "set6"
intData <- AddMetaData(intData, class4, col.name="class4")

Idents(intData) <- "class4"

onlySet <- subset(intData, idents=c("set6"))

#for all interesting PCs
for(i in 1:pcNumH){

x1<- (i-1)*2 + 1
x2<- x1+1

Idents(intData) <- "class1"
fn <- paste0("PC", x1, ".PC", x2, "_gamook.png")
png(fn, width=1000, height=1000)
thingy <- DimPlot(object=intData, dims = c(x1, x2), reduction = "pca", label=T)
print(thingy)
dev.off()

Idents(intData) <- "class2"
fn <- paste0("PC", x1, ".PC", x2, "_subgroups.png")
png(fn, width=1000, height=1000)
thingy <- DimPlot(object=intData, dims = c(x1, x2), reduction = "pca", label=T)
print(thingy)
dev.off()

Idents(intData) <- "class3"
fn <- paste0("PC", x1, ".PC", x2, "_KOvWT.png")
png(fn, width=1000, height=1000)
thingy <- DimPlot(object=intData, dims = c(x1, x2), reduction = "pca", label=T)
print(thingy)
dev.off()

Idents(onlySet) <- "class2"
fn <- paste0("PC", x1, ".PC", x2, "_set6.png")

png(fn, width=1000, height=1000)
thingy <- DimPlot(object=onlySet, dims = c(x1, x2), reduction = "pca", pt.size=4, cols=subCols, label=T)

print(thingy)
dev.off()



}









