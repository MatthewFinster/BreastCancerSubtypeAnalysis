## Loading packages
library(cluster)
library(genefilter)
library(gplots)
library(limma)
library(qvalue)
library(survival)
library(ggplot2)
library(ggfortify)
library(org.Hs.eg.db)
library(GO.db)
library(xtable)

## --------------------------------- NOTES ----------------------------------
# Author: Matthew Finster

# Notes: There is plenty of code that is commented out throughout this R file.
#        The reason for this is that often I conduct manual checks to see that
#        the code works as expected. I like to keep these checks in the code as
#        I can conduct the checks again later and not forget the function at
#        each processing step.

# Date:  25/10/2024

## ------------------------------ PREPARATION --------------------------------

## -- Setting working directory 
setwd("C:\\Users\\anais\\Documents\\Matt\\Etudies\\13. Models for Bioinformatics\\Applied Project\\Assessment 2 - Datasets-20241010")

## -- Loading and inspecting the expression data 
expressionData <- readRDS("STA5MB_2024_BC_expression_data.RDS")
# head(expressionData)
dim(expressionData) # 22283 genes, 251 patients

## -- Loading and inspecting the clinical data
clinicalData <- readRDS("STA5MB_2024_BC_clinical_data.rds")
# clinicalData
dim(clinicalData) # 251 patients, 9 columns of clinical data

## -- Pre-processing: Removing the control probes in the expression data
filteredData <- expressionData[!startsWith(rownames(expressionData), "AFFX"), ]
# head(filteredData)
dim(filteredData) # 22215 genes, 251 patients

## -- Pre-processing: normalising the expression data
normalisedData <- normalizeBetweenArrays(filteredData, method = "scale")
dim(normalisedData)

## -- Pre-processing: scaling the expression data

# This transposes the matrix to scale each gene with 0 as the centre,
# then transposes back before assigning it to the new scaled variable
scaled.expressionData <- t(scale(t(normalisedData), center = TRUE))
# head(scaled.expressionData)

## ------------------------------- CLUSTERING ---------------------------------
# Remember: Clustering individuals to find dissimilar and similar individuals
#           and thus potentially subtypes of breast cancer

## -- Euclidean distances between each individual's gene expression data
distances <- dist(t(scaled.expressionData))
# head(distances)
# individualA <- scaled.expressionData[,151]
# individualB <- scaled.expressionData[,1]
# rem <- sqrt(sum((individualA - individualB)^2))
# distancesmatrix <- as.matrix(distances)
# distancesmatrix[1,151]

## -- Trialling various clustering techniques

# Upon initial inspection, clustering using ward.D2 linkage appears to produce 
# the most logical groupings
par(mfrow = c(2,2))
clusteringWard <- hclust(distances, method = "ward.D2")
plot(clusteringWard, cex = 0.5, main = "Ward linkage", sub = "", xlab = "", ylab = "Distances")
clusteringSingle <- hclust(distances, method = "single")
plot(clusteringSingle, cex = 0.5, main = "Single linkage", sub = "", xlab = "", ylab = "Distances")
clusteringComplete <- hclust(distances, method = "complete")
plot(clusteringComplete, cex = 0.5, main = "Complete linkage", sub = "", xlab = "", ylab = "Distances")
clusteringAverage <- hclust(distances, method = "average")
plot(clusteringAverage, cex = 0.5, main = "Average linkage", sub = "", xlab = "", ylab = "Distances")

# -- Verification via silhouette function

# This accounts for the avg distance between the individual and other individuals 
# in their cluster (a), as well as the average distance between the nearest 
# neighbours in the next closest cluster (b). Seeking distinct
# clusters so large silhouette values (after b - a) are better.
# Notes: [,3] in the silhouette function takes the silhouette width. The
# median allows us the find the average silhouette value using different
# numbers of clusters. k are the different numbers of clustered trialled

par(mfrow = c(2,2))

K <- 2:5
silW <- NULL
silS <- NULL
silC <- NULL
silA <- NULL

# Ward silhouettes
for (i in K) {
  silW <- c(silW, median(silhouette(cutree(clusteringWard, k = i),dist = distances)[,3],na.rm = T))
}
plot(K, silW, type = "l", main = "Ward linkage", xlab = "Number of clusters", ylab = "Median silhouette value")

# Single silhouettes
for (i in K) {
  silS <- c(silS, median(silhouette(cutree(clusteringSingle, k = i),dist = distances)[,3],na.rm = T))
}
plot(K, silS, type = "l", main = "Single linkage", xlab = "Number of clusters", ylab = "Median silhouette value")

# Complete silhouettes
for (i in K) {
  silC <- c(silC, median(silhouette(cutree(clusteringComplete, k = i),dist = distances)[,3],na.rm = T))
}
plot(K, silC, type = "l", main = "Complete linkage", xlab = "Number of clusters", ylab = "Median silhouette value")

# Average silhouettes
for (i in K) {
  silA <- c(silA, median(silhouette(cutree(clusteringAverage, k = i),dist = distances)[,3],na.rm = T))
}
plot(K, silA, type = "l", main = "Average linkage", xlab = "Number of clusters", ylab = "Median silhouette value")


# The silhouette analysis showed that 3 clusters produced the largest median
# silhouette value for ward.D2 linkage (approx 0.035); 2 produced the largest median 
# silhouette value for single linkage (approx 0.24); 5 produced the largest median
# silhouette value for complete linkage (approx 0.015); and 2 clusters produced 
# the largest median silhouette value for average linkage as well (approx 0.24).
# Although the silhouette analysis seems to infer that single linkage
# and average linkage produced the most well-defined clusters;
# based on the reading of the dendrograms, these methods would 
# group 250 individuals into one cluster and 1 individual (85A03) into the other.
# Not ideal but could be relevant if 85A03 has a special rare type of breast cancer.
# Given the silhouette value for 3 clusters (Ward) was higher than
# the silhouette value for 5 clusters (complete), this linkage method was thought
# to be best. However,this cannot be certain until we examine the data to see
# if the clusters have biological differences.

# -- Obtaining optimal clusters
clW <- cutree(clusteringWard, k = K[which.max(silW)])
clS <- cutree(clusteringSingle, k = K[which.max(silS)])
clC <- cutree(clusteringComplete, k = K[which.max(silC)])
clA <- cutree(clusteringAverage, k = K[which.max(silA)])

# -- PCA plots
par(bg = "white")
pcaW <- princomp(scaled.expressionData[, ])
plot(pcaW$loadings[, 1:2], col = clW)
title("Ward linkage (k = 3)")

pcaS <- princomp(scaled.expressionData[, ])
plot(pcaS$loadings[, 1:2], col = clS)
title("Single linkage (k = 2)")

pcaC <- princomp(scaled.expressionData[, ])
plot(pcaC$loadings[, 1:2], col = clC)
title("Complete linkage (k = 5)")

pcaA <- princomp(scaled.expressionData[, ])
plot(pcaA$loadings[, 1:2], col = clA)
title("Average linkage (k = 2)")


## ------ CHECKING THE BIOLOGICAL RELEVANCE OF THE CLUSTERS VIA HEATMAPS ------

# -- Finding the 5% of genes with the highest variance 
#    using the unscaled gene expression data
22215/100*5
rowVariances <- rowVars(normalisedData)
highVarianceGenes <- order(-rowVariances)[1:1111]


# -- Producing heatmap for Ward
colsW <- c("blue", "purple", "orange")
head(cbind(colnames(normalisedData), colsW))
heatmap.2(scaled.expressionData[highVarianceGenes,],labCol = clW, trace = "none", ColSideColors = colsW[clW], col = redgreen(100), Colv = as.dendrogram(clusteringWard))


# -- Producing heatmap for Single
colsS <- c("pink", "purple")
head(cbind(colnames(normalisedData), colsS))
heatmap.2(scaled.expressionData[highVarianceGenes,],labCol = clS, trace = "none", ColSideColors = colsS[clS], col = redgreen(100), Colv = as.dendrogram(clusteringSingle))

# -- Producing heatmap for Complete
colsC <- c("blue", "orange", "purple", "pink", "yellow")
head(cbind(colnames(normalisedData), colsC))
heatmap.2(scaled.expressionData[highVarianceGenes,],labCol = clC, trace = "none",ColSideColors = colsC[clC], col = redgreen(100), Colv = as.dendrogram(clusteringComplete))

# -- Producing heatmap for Average
colsA <- c("pink", "purple")
head(cbind(colnames(normalisedData), colsA))
heatmap.2(scaled.expressionData[highVarianceGenes,],labCol = clA, trace = "none",ColSideColors = colsA[clA], col = redgreen(100), Colv = as.dendrogram(clusteringAverage))


## ------ INFERENTIAL STATISTICS: DIFFERENTIALLY EXPRESSED GENES ANALYSIS ------

# -- Specifying a design matrix
design <- model.matrix(~0 + as.factor(clW))
colnames(design) <- c("A", "B", "C")

# -- Constructing DE object
DE.object <- lmFit(normalisedData, design)

# -- Making contrasts between Cluster 1 and Cluster 2,
#    then reorienting model to contrast 1 and 2
AvsB.contrast <- makeContrasts(BvsA = B - A, levels = design)
AvsB.DE.object <- contrasts.fit(DE.object, AvsB.contrast)
# Making contrasts between Cluster 2 and Cluster 3,
# then reorienting model to contrast 2 and 3
BvsC.contrast <- makeContrasts(CvsB = C-B, levels = design)
BvsC.DE.object <- contrasts.fit(DE.object, BvsC.contrast)
# Making contrasts between Cluster 1 and cluster 3,
# then reorienting model to contrast 1 and 3
AvsC.contrast <- makeContrasts(CvsA = C-A, levels = design)
AvsC.DE.object <- contrasts.fit(DE.object, AvsC.contrast)

# -- Performing empirical Bayes on each contrast
AvsB.DE.object <- eBayes(AvsB.DE.object)
BvsC.DE.object <- eBayes(BvsC.DE.object)
AvsC.DE.object <- eBayes(AvsC.DE.object)

# -- Volcano plots of each contrast
par(mfrow=c(2,2))
volcanoplot(AvsB.DE.object, main = "Cluster 1 vs Cluster 2 contrast")
abline(h = -log10(0.01), col = "red", lty = 2)
abline(v = c(-1, 1), col = "red", lty = 2)
volcanoplot(BvsC.DE.object, main = "Cluster 2 vs Cluster 3 contrast")
abline(h = -log10(0.01), col = "red", lty = 2)
abline(v = c(-1, 1), col = "red", lty = 2)
volcanoplot(AvsC.DE.object, main = "Cluster 1 vs Cluster 3 contrast")
abline(h = -log10(0.01), col = "red", lty = 2)
abline(v = c(-1, 1), col = "red", lty = 2)

# -- Obtaining DE genes with qvalue <= 0.05
AvsB.qval <- qvalue(AvsB.DE.object$p.value[,1], fdr.level = 0.05)
BvsC.qval <- qvalue(BvsC.DE.object$p.value[,1], fdr.level = 0.05)
AvsC.qval <- qvalue(AvsC.DE.object$p.value[,1], fdr.level = 0.05)

## -- Number of DE genes
sum(AvsB.qval$qvalues < 0.05)
sum(BvsC.qval$qvalues < 0.05)
sum(AvsC.qval$qvalues < 0.05)

## -- Top 100 DE genes by q-value
top_AvsB_qvals <- AvsB.DE.object[order(AvsB.qval$qvalues)[1:100], ]
top_BvsC_qvals <- BvsC.DE.object[order(BvsC.qval$qvalues), ][1:100, ]
top_AvsC_qvals <- AvsC.DE.object[order(AvsC.qval$qvalues), ][1:100, ]

## -- Finding the top 5 genes ranked by q value, and their logFC changes
top5qvalsAvsB <- sort(AvsB.qval$qvalues)[1:5]
top5qvalsBvsC <- sort(BvsC.qval$qvalues)[1:5]
top5qvalsAvsC <- sort(AvsC.qval$qvalues)[1:5]
top_AvsB_qvals$coefficients[1:5]
top_BvsC_qvals$coefficients[1:5]
top_AvsC_qvals$coefficients[1:5]

## ---------- GENE ONTOLOGY AND KEGG PATHWAYS ANALYSES ON DE GENES ------------

## -- Reading the annotation data
annotationData <- readRDS("STA5MB_2024_BC_annotations.RDS")
#dim(annotationData)
#annotationData

## -- Getting entrez IDs for top 100 DE genes by q values
# Cluster 1 vs 2
top_AvsB_probeIDs <- rownames(top_AvsB_qvals)
AvsB_ensemblIDs <- annotationData[annotationData$affy_hg_u133_plus_2 %in% top_AvsB_probeIDs, c("ensembl_gene_id")]
homo_sapiens <- org.Hs.eg.db
AvsB_entrezIDs <- select(
  homo_sapiens,
  keys = AvsB_ensemblIDs,
  columns = c("ENTREZID", "ENSEMBL"),
  keytype = "ENSEMBL"
)
AvsB_entrezIDs <- unique(na.omit(AvsB_entrezIDs$ENTREZID))

# Cluster 2 vs 3
top_BvsC_probeIDs <- rownames(top_BvsC_qvals)
BvsC_ensemblIDs <- annotationData[annotationData$affy_hg_u133_plus_2 %in% top_BvsC_probeIDs, c("ensembl_gene_id")]
BvsC_entrezIDs <- select(
  homo_sapiens,
  keys = BvsC_ensemblIDs,
  columns = c("ENTREZID", "ENSEMBL"),
  keytype = "ENSEMBL"
)
BvsC_entrezIDs <- unique(na.omit(BvsC_entrezIDs$ENTREZID))

# Cluster 1 vs 3
top_AvsC_probeIDs <- rownames(top_AvsC_qvals)
AvsC_ensemblIDs <- annotationData[annotationData$affy_hg_u133_plus_2 %in% top_AvsC_probeIDs, c("ensembl_gene_id")]
AvsC_entrezIDs <- select(
  homo_sapiens,
  keys = AvsC_ensemblIDs,
  columns = c("ENTREZID", "ENSEMBL"),
  keytype = "ENSEMBL"
)
AvsC_entrezIDs <- unique(na.omit(AvsC_entrezIDs$ENTREZID))

## -- GO for each contrast
GO_AvsB <- goana(AvsB_entrezIDs, species='Hs')
GO_BvsC <- goana(BvsC_entrezIDs, species='Hs')
GO_AvsC <- goana(AvsC_entrezIDs, species='Hs')
GO_top_AvsB <- topGO(GO_AvsB, n=20, ontology='BP')
GO_top_BvsC <- topGO(GO_BvsC, n=20, ontology='BP')
GO_top_AvsC <- topGO(GO_AvsC, n=20, ontology='BP')

## -- KEGG for each contrast
KEGG_AvsB <- kegga(AvsB_entrezIDs, species='Hs')
KEGG_BvsC <- kegga(BvsC_entrezIDs, species='Hs')
KEGG_AvsC <- kegga(AvsC_entrezIDs, species='Hs')
KEGG_top_AvsB <- topKEGG(KEGG_AvsB, n=20)
KEGG_top_BvsC <- topKEGG(KEGG_BvsC, n=20)
KEGG_top_AvsC <- topKEGG(KEGG_AvsC, n=20)

## -- Skeletons for tables in latex
go_table_AvsB <- xtable(GO_top_AvsB)
go_table_BvsC <- xtable(GO_top_BvsC)
go_table_AvsC <- xtable(GO_top_AvsC)
kegg_table_AvsB <- xtable(KEGG_top_AvsB)
kegg_table_BvsC <- xtable(KEGG_top_BvsC)
kegg_table_AvsC <- xtable(KEGG_top_AvsC)

## --------------------------- SURVIVAL ANALYSES ------------------------------

# -- Identifying significant genes common to all three contrasts
significant_genes <- intersect(intersect(rownames(AvsB.DE.object)[AvsB.qval$significant],
                                         rownames(BvsC.DE.object)[BvsC.qval$significant]),
                               rownames(AvsC.DE.object)[AvsC.qval$significant])

## -- Calculating gene scores by summing across these common significant genes for each sample
gene.scores <- colSums(normalisedData[significant_genes, ]) 

## -- Creating a survival object (using Surv_time and event from the clinical data)
Y <- Surv(time = clinicalData$Surv_time, event = clinicalData$event)

# Convergence issues with LNstatus and ERstatus when running the Cox model in
# the next step, likely due to the LN? and ER? values.
# table(clinicalData$ERstatus, clinicalData$event)
# table(clinicalData$LNstatus, clinicalData$event)
# Given that + EstrogenStatus (ERstatus) is associated with better survival
# and given that - LymphNodeMetastatisStatus (LNstatus) is assosiated with better survival,
# and the ER? values did not have any deaths and the LN? values did not have any deaths,
# these unknown ER? and LN? values were merged into the ER+ and LN- groups.
clinicalData$ERstatus[clinicalData$ERstatus == "ER?"] <- "ER+"
clinicalData$LNstatus[clinicalData$LNstatus == "LN?"] <- "LN-"

## -- Changing the reference level to cluster 2
clinicalData$clW <- relevel(as.factor(clW), ref = "2")

## -- Fitting the Cox proportional hazards model
# cox.modelSimple <- coxph(Y ~ gene_scores, data = clinicalData)
# cox.modelSimple <- coxph(Y ~ as.factor(clW), data = clinicalData)
# cox.modelRobust <- coxph(Y ~ gene_scores + age + as.factor(histgrade) + tumor_size_mm + PRstatus + ERstatus + LNstatus, data = clinicalData)
cox.modelRobust <- coxph(Y ~ gene.scores + as.factor(clW) + age + as.factor(histgrade) + tumor_size_mm + as.factor(PRstatus) + as.factor(ERstatus) + as.factor(LNstatus), data = clinicalData)
xtable(cox.modelRobust)

## -- New Cox model: combining cluster 1 and cluster 3 into a new cluster group
clinicalData$clW_grouped <- ifelse(clinicalData$clW == 1 | clinicalData$clW == 3, "1&3", as.character(clinicalData$clW))
clinicalData$clW_grouped <- relevel(as.factor(clinicalData$clW_grouped), ref = "2")
cox.modelGrouped <- coxph(Y ~ gene.scores + as.factor(clW_grouped) + age + as.factor(histgrade) + tumor_size_mm + as.factor(PRstatus) + as.factor(ERstatus) + as.factor(LNstatus), data = clinicalData)
summary(cox.modelGrouped)

## -- Kaplan-Meier survival curves
km_fit <- survfit(Y ~ clW, data = clinicalData)
km_fitGrouped <- survfit(Y ~ clW_grouped, data = clinicalData)

## -- Plotting the Kaplan-Meier curves
autoplot(km_fit)+labs(x = "\n Survival time (years) since diagnosis ", y = "Survival probabilities", 
                      title = "K-M survival curves for breast cancer patients (Clusters 1 vs 2 vs 3) \n", color = "Cluster Group", fill = "Cluster Group")
autoplot(km_fitGrouped)+labs(x = "\n Survival time (years) since diagnosis ", y = "Survival probabilities", 
                      title = "K-M survival curves for breast cancer patients (Cluster 2 vs Cluster 1&3) \n", color = "Cluster Group", fill = "Cluster Group")
