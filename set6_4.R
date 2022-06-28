setwd("C:\\Work\\set6\\")
library("DESeq2")
library("pcaExplorer")
library("ggplot2")
library("ggrepel")
library("graphics")
library("dendextend")
library("corrplot")
library("tidyr")
library("glmpca")
library("colorspace")
library("devtools")
library("loadings")
library("gplots")
library("data.table")




###########################################READ IN DATA, ORDER DATA, GENERATE INFO TABLE

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




###########################################FILTER AND NORMALISE COUNT DATA
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

ns6Data <- ns6Data[vars,]


###########################################PCA STUFF NOW BELOW, DONE WITH GGPLOT2 INSTEAD OF THIS WEIRD WRAPPER
#was meant to do everything but no scree part so have to do prcomp anyway, and no obvious way to change colours

#pcaplot(normSet6, intgroup = c("cond"), ntop = 5003,
#        pcX = 1, pcY = 2, title = "Set6 Cond - PC1 vs PC2",
#        ellipse = F)

#pcaplot(normSet6, intgroup = c("cond"), ntop = 5003,
#        pcX = 3, pcY = 4, title = "Set6 Cond - PC3 vs PC4",
#        ellipse = F)

#pcaplot(normSet6, intgroup = c("cond"), ntop = 5003,
#        pcX = 5, pcY = 6, title = "Set6 Cond - PC3 vs PC4",
#        ellipse = F)

#png("set6_pc12.png", width=1500, height=750)
#pcaplot(normSet6, intgroup = c("condStage"),ntop = 5003, pcX = 1, pcY = 2,
#	title = "Set6 Cond and Stage - PC1 vs PC2", ellipse = F,
#	text_labels=F) + geom_label_repel(aes(label = meta$condStage), size = 4)
#dev.off()

#png("set6_pc34.png", width=1500, height=750)
#pcaplot(normSet6, intgroup = c("condStage"),ntop = 5003, pcX = 3, pcY = 4,
#	title = "Set6 Cond and Stage - PC3 vs PC4", ellipse = F,
#	text_labels=F) + geom_label_repel(aes(label = meta$condStageRep), size = 4)
#dev.off()

#with max genes
#pcaplot(normSet6, intgroup = c("condStage"),ntop = 10000, pcX = 1, pcY = 2,
#	title = "Set6 Cond and Stage - PC1 vs PC2", ellipse = F,
#	text_labels=F) + geom_label_repel(aes(label = meta$condStageRep), size = 4)

#put in a restrictive one too i.e. back to 1000
#png("set6_pc12_1000genes.png", width=1500, height=750)
#pcaplot(normSet6, intgroup = c("condStage"), ntop = 1000, pcX = 1, pcY = 2,
#	title = "Set6 Cond and Stage - PC1 vs PC2", ellipse = F,
#	text_labels=F) + geom_label_repel(aes(label = meta$condStageRep), size = 4)
#dev.off()

#and one with the original labels for Scott to see the groups
#png("set6_pc12_exp.png", width=1500, height=750)
#pcaplot(normSet6, intgroup = c("condStage"),ntop = 5003, pcX = 1, pcY = 2,
#	title = "Set6 Cond and Stage - PC1 vs PC2", ellipse = F,
#	text_labels=F) + geom_label_repel(aes(label = meta$condStageRep), size = 4)
#dev.off()

#didn't like this, didn't seem to load properly 
#pcaExplorer(dst=normSet6)




###########################################PCA

set6PCA <- prcomp(t(ns6Data), scale=T)
png("screeplot.png", width=750, height=750)
screeplot(set6PCA, type="lines", main="Set6 PCA", npcs=20)
dev.off()
#eigenvalues
set6PCA$sdev^2
summary(set6PCA)
set6PCA$x[1:5,1:3]
plot(set6PCA$x[,1:2])
text(set6PCA$x[,1], set6PCA$x[,2], meta$condStage)
#WT.16.3 is a bit off from expected
plot(set6PCA$x[,1:2])
text(set6PCA$x[,1], set6PCA$x[,2], meta$condStageRep)

#look at loadings, maybe use PCA2GO?
#they are very small?! - multiply by sdev
#sort(set6PCA$rotation[,1])
loadings <- set6PCA$rotation %*% diag(set6PCA$sdev)
colnames(loadings) <- colnames(set6PCA$x)
write.table(loadings, file="loadings.txt", quote=F, sep="\t")

#check it, PC1 and PBANKA_1014300 are 0.9916 in the matrix
cor.test(ns6Data["PBANKA_1014300",], set6PCA$x[,1])
#works




###########################################P VALUES FOR LOADINGS

#Scott M wanted p values for the loadings
pcaL <- pca_loading(set6PCA)
#check the loadings are calculated the same and order is maintained
all.equal(pcaL$loading$R, loadings)
#output the p values - repeats columns for some reason?!
#check values vs the cor.test above
pcaL$loading$p.value["PBANKA_1014300",1]
#summary above shows smaller than a set limit, have to ask specifically for p to get the actual value
cor.test(ns6Data["PBANKA_1014300",], set6PCA$x[,1])$p.value
write.table(pcaL$loading$p.value[,1:24], file="loadingsP.txt", quote=F, sep="\t")

#generate master table
loadOut <- rownames(pcaL$loading$p.value)
for(i in 1:24){
loadOut <- cbind(loadOut, pcaL$loading$R[,i], pcaL$loading$p.value[,i])
pcNames <- c(paste0(colnames(pcaL$loading$R)[i], "_cor"), paste0(colnames(pcaL$loading$R)[i], "_p"))
pos1 <- ncol(loadOut)-1
pos2 <- ncol(loadOut)
colnames(loadOut)[pos1:pos2] <- pcNames
}
colnames(loadOut)[1] <- "gene"
write.table(loadOut, file="loadingsCorP.txt", quote=F, sep="\t", row.names=F)




###########################################PCA PLOTS

#1D plots of first 10 PCs (being generous with the elbow, in case anyone wants more)
theseCols <- c("#EE6447", "#D62F1F", "#B30000", "#7E0000", "#4EB3D2", "#2B8BBD", "#0768AC", "#074081")
pcData <- cbind(data.frame(set6PCA$x[,1:10]), condStageRep, condStage, theseCols)
pcData <- df.long <- pivot_longer(pcData, cols=1:10, names_to = "PC", values_to = "Value")
pcData <- data.frame(pcData)
pcData[,"PC"] <- factor(pcData[,"PC"], levels=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

#first ten PCs, 1D
png("10_PCs_1D.png", width=1000, height=500)
ggplot(pcData, aes(x=PC, y=Value, color=condStage)) + geom_point() + geom_text_repel(label=pcData$condStageRep) + scale_color_manual("Condition and Stage", values=theseCols) 
dev.off()

#bigger so more legible, still some labels not shown
png("10_PCs_1D_big.png", width=1500, height=1000)
ggplot(pcData, aes(x=PC, y=Value, color=condStage)) + geom_point() + geom_text_repel(label=pcData$condStageRep) + scale_color_manual("Condition and Stage", values=theseCols) 
dev.off()

#PC1 and PC2
newPCs <- cbind(data.frame(set6PCA$x[,1:10]), condStageRep, condStage, theseCols)
png("set6_pc12.png", width=750, height=750)
ggplot(newPCs, aes(x=PC1, y=PC2, color=condStage)) + geom_point() + geom_text_repel(label=newPCs$condStageRep, size=6) + scale_color_manual("Condition and Stage", values=theseCols) + theme(text = element_text(size = 15)) + geom_point(size=3)
dev.off()

#PC3 has flipped
png("set6_pc34.png", width=750, height=750)
ggplot(newPCs, aes(x=PC3, y=PC4, color=condStage)) + geom_point() + geom_text_repel(label=newPCs$condStageRep, size=6) + scale_color_manual("Condition and Stage", values=theseCols) + theme(text = element_text(size = 15)) + geom_point(size=3)
dev.off()

#do PCA for top 1000 genes by variance
set6PCA_1000 <- prcomp(t(ns6Data[1:1000,]), scale=T)
newPCs_1000 <- cbind(data.frame(set6PCA_1000$x[,1:10]), condStageRep, condStage, theseCols)
png("set6_pc12_1000.png", width=750, height=750)
ggplot(newPCs_1000, aes(x=PC1, y=PC2, color=condStage)) + geom_point() + geom_text_repel(label=newPCs$condStageRep, size=6) + scale_color_manual("Condition and Stage", values=theseCols) + theme(text = element_text(size = 15)) + geom_point(size=3)
dev.off()




###########################################CHECK GLMPCA DOESN'T LOOK ESPECIALLY DIFFERENT

#ADD IN BATCH IF BATCH EXISTS - is there a batch effect to be controlled?

#compare pca with glmpca, it didn't like being centered - "min(Y) >= 0 is not TRUE"
set.seed(202)
scaleMat <- scale(ns6Data, scale=T, center=F)
set6GLMPCA <- glmpca(scaleMat, 2, fam="poi")
print(set6GLMPCA)

#check optimizer decreased deviance
plot(set6GLMPCA$dev,type="l",xlab="iterations",ylab="Poisson deviance")

#plot factors and compare with PCA
p1 <- ggplot(cbind(as.data.frame(set6PCA$x), meta), aes(x=PC1, y=PC2, color=condStage)) + geom_point() + geom_text_repel(label=meta$condStage) + scale_color_manual(values=theseCols) + theme(legend.position="none")
p2 <- ggplot(cbind(set6GLMPCA$factors, meta), aes(x=dim1, y=dim2, color=condStage)) + geom_point() + geom_text_repel(label=meta$condStage) + scale_color_manual("Condition and Stage", values=theseCols)

png("glmp_vs_pca.png", width=1200, height=600)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()




###########################################HIERARCHICAL CLUSTERING

dist_mat <- dist(as.data.frame(scale(t(ns6Data))), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
#set the hclust labels - you need to get the order from the object

dend <- as.dendrogram(hclust_avg)
labels(dend) <- meta$condStageRep[hclust_avg$order]

theseColsD <- c(rep("#EE6447", 3), rep("#D62F1F", 3), rep("#B30000", 3), rep("#7E0000", 3), rep("#4EB3D2", 3), rep("#2B8BBD", 3), rep("#0768AC", 3), rep("#074081", 3)) 
labels_colors(dend) <- theseColsD[order.dendrogram(dend)]

png("set6_dend.png", width=1000, height=1000)
plot(dend, horiz=T)
dev.off()

#top 1000 genes
dist_mat_1000 <- dist(as.data.frame(scale(t(ns6Data[1:1000,]))), method = 'euclidean')
hclust_avg2 <- hclust(dist_mat_1000, method = 'average')
#set the hclust labels - you need to get the order from the object

dend2 <- as.dendrogram(hclust_avg2)
labels(dend2) <- meta$condStageRep[hclust_avg2$order]
labels_colors(dend2) <- theseColsD[order.dendrogram(dend2)]

png("set6_dend_1000.png", width=1000, height=1000)
plot(dend2, horiz=T)
dev.off()



###########################################LINEAR MODELLING

#fit <- lm(Y ~ 0 + X)

#Where Y is the vector of gene expressions after vst for one of the KO libraries 
#(e.g. the expressions for "KO.8.5"). X is a "reference" matrix of wild type gene expressions. 
#Each column of X could be either each WT library (so X would have columns WT.G1, WT.G.2, ..., WT.24.3)
# or (maybe better?) the average of each WT stage (so X would have columns: WT.G, WT.8, WT.16, WT.24). 
#The coefficients from the fitted model would give an indication of how much each WT library or WT 
#stage explains the gene expression of the KO library. I'm not 100% sure this is fully sensible so 
#don't overthink it and don't get too involved! The idea, later, could be to resample genes in a 
#bootstrap-like fashion to get some sense of how accurate the estimates are.

aveData <- cbind(rowMeans(ns6Data[,1:3], na.rm=TRUE), rowMeans(ns6Data[,4:6], na.rm=TRUE), rowMeans(ns6Data[,7:9], na.rm=TRUE), rowMeans(ns6Data[,10:12], na.rm=TRUE), rowMeans(ns6Data[,13:15], na.rm=TRUE), rowMeans(ns6Data[,16:18], na.rm=TRUE), rowMeans(ns6Data[,19:21], na.rm=TRUE), rowMeans(ns6Data[,22:24], na.rm=TRUE))
colnames(aveData) <- c("KO.G", "KO.8", "KO.16", "KO.24", "WT.G", "WT.8", "WT.16", "WT.24")

#most similar to G, 8, G (close to 8), 24, respectively
fitKOG <- lm(aveData[,"KO.G"] ~ 0 + aveData[,5:8])
summary(fitKOG)
fitKO8 <- lm(aveData[,"KO.8"] ~ 0 + aveData[,5:8])
summary(fitKO8)
fitKO16 <- lm(aveData[,"KO.16"] ~ 0 + aveData[,5:8])
summary(fitKO16)
fitKO24 <- lm(aveData[,"KO.24"] ~ 0 + aveData[,5:8])
summary(fitKO24)

#try top genes by variance, same story
fitKOG_1000 <- lm(aveData[1:1000,"KO.G"] ~ 0 + aveData[1:1000,5:8])
summary(fitKOG_1000)
fitKO8_1000 <- lm(aveData[1:1000,"KO.8"] ~ 0 + aveData[1:1000,5:8])
summary(fitKO8_1000)
fitKO16_1000 <- lm(aveData[1:1000,"KO.16"] ~ 0 + aveData[1:1000,5:8])
summary(fitKO16_1000)
fitKO24_1000 <- lm(aveData[1:1000,"KO.24"] ~ 0 + aveData[1:1000,5:8])
summary(fitKO24_1000)

#do the KO side of the equation 1 at a time
fitList <- list()
estMat <- matrix(nrow=4, ncol=12)
rownames(estMat) <- c("WT_G", "WT_8", "WT_16", "WT_24")
colnames(estMat) <- meta$condStageRep[1:12]
for(i in 1:12){
fitList[[i]] <- lm(ns6Data[,i] ~ 0 + aveData[,5:8])
estMat[,i] <- fitList[[i]]$coefficients
print(meta$condStageRep[i])
print(summary(fitList[[i]]))
}

png("coefficients.png", height=750, width=900)
heatmap.2(t(estMat[,seq(from=12, to=1)]), cellnote=round(t(estMat[,seq(from=12, to=1)]), 2), 
	notecol="black", dendrogram="none", Rowv=FALSE, Colv=FALSE, 
	col = hcl.colors(50, palette="Blue-Yellow3"), scale="none", key=TRUE, density.info="none",
	trace="none", cexRow=1, cexCol=1, symm=F, symkey=T, symbreaks=T)
dev.off()

#NOT ZERO-CENTERED
#heatmap(t(estMat), Rowv=NA, Colv=NA, col=hcl.colors(50, palette="Blue-Yellow3"))
#legend(x="right", legend=c("1", "0", "-1"), fill= rev(hcl.colors(3, palette="Blue-Yellow3")))

#do the same but for simple linear reg
estMat2 <- matrix(nrow=4, ncol=12)
rownames(estMat2) <- c("WT_G", "WT_8", "WT_16", "WT_24")
colnames(estMat2) <- meta$condStageRep[1:12]
for(i in 1:12){
for(j in 1:4){
j2 <- j + 4
fit <- lm(ns6Data[,i] ~ 0 + aveData[,j2])
estMat2[j,i] <- fit$coefficients
print(paste0(meta$condStageRep[i], " and ", colnames(aveData)[j2]))
print(summary(fit))
}
}

png("coefficients2.png", height=750, width=900)
heatmap.2(t(estMat2[,seq(from=12, to=1)]), cellnote=round(t(estMat2[,seq(from=12, to=1)]), 2), 
	notecol="black", dendrogram="none", Rowv=FALSE, Colv=FALSE, 
	col = hcl.colors(50, palette="Blue-Yellow3"), scale="none", key=TRUE, density.info="none",
	trace="none", cexRow=1, cexCol=1, symm=F, symkey=T, symbreaks=T)
dev.off()

############################Dario's graph suggestion
stage <- c('00h', '08h', '16h', '24h')
KO <- paste(rep(paste0('KO.', stage), each= 3), 1:3, sep= '-')
WT <- paste0('WT.', stage)

dat <- data.table(
    expand.grid(WT= WT, KO= KO)
)
dat[, WT.stage := rep(stage, 12)]
dat[, KO.stage := rep(stage, each= 12)]
dat[, coef := as.numeric(estMat)]
dat[, coef2 := as.numeric(estMat2)]
dat[, repl := rep(rep(1:3, each= 4), 4)]

png("coefficients1.1.png", height=700, width=700)
ggplot(data= dat, aes(x= WT.stage, y= coef, colour= as.character(repl), label= repl)) +
    geom_point() +
    geom_text(hjust= -1) +
    facet_wrap(~ paste('KO stage:', KO.stage), ncol= 2) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), legend.position= 'none') +
	xlab("WT Stage") + ylab("Coefficient") + ggtitle("Multiple Linear Regression") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("coefficients2.1.png", height=700, width=700)
ggplot(data= dat, aes(x= WT.stage, y= coef2, colour= as.character(repl), label= repl)) +
    geom_point() +
    geom_text(hjust= -1) +
    facet_wrap(~ paste('KO stage:', KO.stage), ncol= 2) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), legend.position= 'none') +
	xlab("WT Stage") + ylab("Coefficient") + ggtitle("Multiple Linear Regression") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()




###########################################CORRELATIONS

kocor = cor(aveData)
png("corr.png", width=1000, height=1000)
corrplot(kocor)
dev.off()
png("corrNum.png", width=1000, height=1000) 
corrplot(kocor, method = 'number') # colorful number
dev.off()

png("corrNumWT.png", width=500, height=500) 
corrplot(kocor[5:8,5:8], method = 'number') # colorful number
dev.off()


#stacked bar thing
KOs <- factor(c(rep(c("KO-G", "KO-8", "KO-16", "KO-24"), each = 4)), levels=c("KO-G", "KO-8", "KO-16", "KO-24")) 
WTs  <- factor(c(rep(c("WT-G", "WT-8", "WT-16", "WT-24"), times = 4)), levels=c("WT-G", "WT-8", "WT-16", "WT-24"))
Correlations <- round(c(kocor[1,5:8], kocor[2,5:8], kocor[3,5:8], kocor[4,5:8])*100, digits=1)
Data <- data.frame(KOs, WTs, Correlations)
#SAVE THIS? not v informative
ggplot(Data, aes(x = KOs, y = Correlations, fill = WTs, label = Correlations)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))


#colours are different now but misleading becuase of the range of the scale
heatmap(kocor[1:4,5:8], Rowv=NA, Colv=NA)

#do with ggplot2 and have key so not misleading
png("KO-WT-cor.png", width=500, height=500)
ggplot(Data, aes(WTs, KOs, fill= Correlations)) + 
  geom_tile(color = "black") +
  geom_text(aes(label = Correlations), color = "white", size = 4) +
  coord_fixed()
dev.off()

#do the same with WT-16 removed
aveData2 <- cbind(rowMeans(ns6Data[,1:3], na.rm=TRUE), rowMeans(ns6Data[,4:6], na.rm=TRUE), rowMeans(ns6Data[,7:9], na.rm=TRUE), rowMeans(ns6Data[,10:12], na.rm=TRUE), rowMeans(ns6Data[,13:15], na.rm=TRUE), rowMeans(ns6Data[,16:18], na.rm=TRUE), rowMeans(ns6Data[,19:20], na.rm=TRUE), rowMeans(ns6Data[,22:24], na.rm=TRUE))
colnames(aveData2) <- c("KO.G", "KO.8", "KO.16", "KO.24", "WT.G", "WT.8", "WT.16", "WT.24")

#
kocor2 = cor(aveData2)
Correlations2 <- round(c(kocor2[1,5:8], kocor2[2,5:8], kocor2[3,5:8], kocor2[4,5:8])*100, digits=1)
Data2 <- data.frame(KOs, WTs, Correlations2)

#
png("KO-WT-cor-outlier.png", width=500, height=500)
ggplot(Data2, aes(WTs, KOs, fill= Correlations2)) + 
  geom_tile(color = "black") +
  geom_text(aes(label = Correlations2), color = "white", size = 4) +
  coord_fixed()
dev.off()







#do corr with individual KOs
aveData3 <- cbind(ns6Data[,1:12], aveData[,5:8])
colnames(aveData3) <- c(colnames(ns6Data)[1:12], colnames(aveData)[5:8])
kocor3 = cor(aveData3)
Correlations3 <- round(c(kocor3[1,13:16], kocor3[2,13:16], kocor3[3,13:16], kocor3[4,13:16], kocor3[5,13:16], kocor3[6,13:16], kocor3[7,13:16], kocor3[8,13:16], kocor3[9,13:16], kocor3[10,13:16], kocor3[11,13:16], kocor3[12,13:16])*100, digits=1)
WTs2  <- factor(c(rep(c("WT-G", "WT-8", "WT-16", "WT-24"), times = 12)), levels=c("WT-G", "WT-8", "WT-16", "WT-24"))
KOs2 <- factor(c(rep(meta$condStageRep[1:12], each = 4)), levels=meta$condStageRep[1:12]) 
Data3 <- data.frame(KOs2, WTs2, Correlations3)

png("indKO-WT-cor.png", width=400, height=750)
ggplot(Data3, aes(WTs2, KOs2, fill= Correlations3)) + 
	geom_tile(color = "black") +
	geom_tile(data = Data3[c(1,5,9,14,17,22,25,30,34,39,43,46),], fill = NA, color = "red", size = 1) +
	geom_text(aes(label = Correlations3), color = "white", size = 4) +
	coord_fixed() +
	xlab("WT Stage") + ylab("KO Stage and Replicate") + ggtitle("Correlations") +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(legend.title=element_blank())
dev.off()

#same again but with WT-16-3 removed
aveData4 <- cbind(ns6Data[,1:12], aveData2[,5:8])
colnames(aveData4) <- c(colnames(ns6Data)[1:12], colnames(aveData)[5:8])
kocor4 = cor(aveData4)
Correlations4 <- round(c(kocor4[1,13:16], kocor4[2,13:16], kocor4[3,13:16], kocor4[4,13:16], kocor4[5,13:16], kocor4[6,13:16], kocor4[7,13:16], kocor4[8,13:16], kocor4[9,13:16], kocor4[10,13:16], kocor4[11,13:16], kocor4[12,13:16])*100, digits=1)
WTs2  <- factor(c(rep(c("WT-G", "WT-8", "WT-16", "WT-24"), times = 12)), levels=c("WT-G", "WT-8", "WT-16", "WT-24"))
KOs2 <- factor(c(rep(meta$condStageRep[1:12], each = 4)), levels=meta$condStageRep[1:12]) 
Data4 <- data.frame(KOs2, WTs2, Correlations4)

png("indKO-WT-cor-outrem.png", width=500, height=750)
ggplot(Data4, aes(WTs2, KOs2, fill= Correlations4)) + 
	geom_tile(color = "black") +
	geom_tile(data = Data3[c(1,5,9,14,17,22,25,30,34,39,43,46),], fill = NA, color = "red", size = 1) +
	geom_text(aes(label = Correlations4), color = "white", size = 4) +
	coord_fixed() +
	xlab("WT Stage") + ylab("KO Stage and Replicate") + ggtitle("Correlations") +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(legend.title=element_blank())
dev.off()


#######################and last two graphs again with Dario's code
dat[, corr3 := Data3$Correlations3]
dat[, corr4 := Data4$Correlations4]

png("indKO-WT-cor-D.png", height=700, width=700)
ggplot(data= dat, aes(x= WT.stage, y= corr3, colour= as.character(repl), label= repl)) +
    geom_point() +
    geom_text(hjust= -1) +
    facet_wrap(~ paste('KO stage:', KO.stage), ncol= 2) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), legend.position= 'none') +
	xlab("WT Stage") + ylab("rho") + ggtitle("Correlations") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


png("indKO-WT-cor-outrem-D.png", height=700, width=700)
ggplot(data= dat, aes(x= WT.stage, y= corr4, colour= as.character(repl), label= repl)) +
    geom_point() +
    geom_text(hjust= -1) +
    facet_wrap(~ paste('KO stage:', KO.stage), ncol= 2) +
    theme_light() +
    theme(strip.text= element_text(colour= 'black'), legend.position= 'none') +
	xlab("WT Stage") + ylab("rho") + ggtitle("Correlations") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()




###########################################DIFFERENTIAL EXPRESSION - should have written this bit as a function

#the data has been transformed to a log2-ish scale already so maybe the diff exp is actually just a subtraction?
#"The transformed data is on the log2 scale for large counts." - http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variance-stabilizing-transformation
#NOTE: THIS IS A 'LIMMA' STYLE FOLD CHANGE, AN APPROXIMATION THAT IS SLIGHTLY OFF FROM THE TRUE FOLD CHANGE
#i.e. it is the difference between the means of the log2 values, instead of the log2 transform of the difference between the mean of the raw values

KoWtDe <- data.frame(nrow=nrow(ns6Data), ncol=32)

#all KO vs all WT
counter <- 1
#for every KO group
for(i in 1:4){
#for every WT group
for(j in 1:4){
for(thisRow in 1:nrow(ns6Data)){

#KO
x1 <- (i-1)*3+1
x2 <- (i-1)*3+3
#WT
x3 <- (j-1)*3+13
x4 <- (j-1)*3+15

#store
KoWtDe[thisRow,counter] <- mean(ns6Data[thisRow,x1:x2]) - mean(ns6Data[thisRow,x3:x4])
KoWtDe[thisRow,counter+1] <- wilcox.test(ns6Data[thisRow,x1:x2], ns6Data[thisRow,x3:x4])$p.value
}

#column header
this1 <- strsplit(colnames(ns6Data)[x1], "_")[[1]][1]
this2 <- strsplit(colnames(ns6Data)[x3], "_")[[1]][1]
colnames(KoWtDe)[counter] <- paste0(substr(this1, 1, nchar(this1)-2), "x", substr(this2, 1, nchar(this2)-2), "_log2FC")
colnames(KoWtDe)[counter + 1] <-paste0(substr(this1, 1, nchar(this1)-2), "x", substr(this2, 1, nchar(this2)-2), "_p")
 
counter <- counter+2
}
}


#all KO ooks vs all WT ooks
ooks <- data.frame(nrow=nrow(ns6Data), ncol=2)
for(thisRow in 1:nrow(ns6Data)){

ooks[thisRow,1] <- mean(ns6Data[thisRow,4:12]) - mean(ns6Data[thisRow,16:24])
ooks[thisRow,2] <- wilcox.test(ns6Data[thisRow,4:12], ns6Data[thisRow,16:24])$p.value
}
colnames(ooks) <- c("KoOokxWtOok_log2FC", "KoOokxWtOok_p")



#all KO against all WT
all <- data.frame(nrow=nrow(ns6Data), ncol=2)
for(thisRow in 1:nrow(ns6Data)){

all[thisRow,1] <- mean(ns6Data[thisRow,1:12]) - mean(ns6Data[thisRow,13:24])
all[thisRow,2] <- wilcox.test(ns6Data[thisRow,1:12], ns6Data[thisRow,13:24])$p.value
}
colnames(all) <- c("KoxWt_log2FC", "KoxWt_p")






#each KO group again all WT
KoAllWtDe <- data.frame(nrow=nrow(ns6Data), ncol=8)
counter <- 1

#for every KO group
for(i in 1:4){
for(thisRow in 1:nrow(ns6Data)){

#KO
x1 <- (i-1)*3+1
x2 <- (i-1)*3+3

#store
KoAllWtDe[thisRow,counter] <- mean(ns6Data[thisRow,x1:x2]) - mean(ns6Data[thisRow,13:24])
KoAllWtDe[thisRow,counter+1] <- wilcox.test(ns6Data[thisRow,x1:x2], ns6Data[thisRow,13:24])$p.value
}

#column header
this1 <- strsplit(colnames(ns6Data)[x1], "_")[[1]][1]
colnames(KoAllWtDe)[counter] <- paste0(substr(this1, 1, nchar(this1)-2), "xWT", "_log2FC")
colnames(KoAllWtDe)[counter + 1] <-paste0(substr(this1, 1, nchar(this1)-2), "xWT", "_p")
 
counter <- counter+2
}



#each KO ook against all WT ooks
KoAllWtOoksDe <- data.frame(nrow=nrow(ns6Data), ncol=6)
counter <- 1

#for every ook KO group
for(i in 2:4){
for(thisRow in 1:nrow(ns6Data)){

#KO
x1 <- (i-1)*3+1
x2 <- (i-1)*3+3

#store
KoAllWtOoksDe[thisRow,counter] <- mean(ns6Data[thisRow,x1:x2]) - mean(ns6Data[thisRow,16:24])
KoAllWtOoksDe[thisRow,counter+1] <- wilcox.test(ns6Data[thisRow,x1:x2], ns6Data[thisRow,16:24])$p.value
}

#column header
this1 <- strsplit(colnames(ns6Data)[x1], "_")[[1]][1]
colnames(KoAllWtOoksDe)[counter] <- paste0(substr(this1, 1, nchar(this1)-2), "xWT", "_log2FC")
colnames(KoAllWtOoksDe)[counter + 1] <-paste0(substr(this1, 1, nchar(this1)-2), "xWT", "_p")
 
counter <- counter+2
}

ID <- rownames(ns6Data)
allDE <- cbind(ID, all, ooks, KoWtDe, KoAllWtDe, KoAllWtOoksDe)
write.table(allDE, file="diffExp.txt", quote=F, row.names=F, sep="\t")


#do again with WT outlier removed like with previous parts of the analysis?






