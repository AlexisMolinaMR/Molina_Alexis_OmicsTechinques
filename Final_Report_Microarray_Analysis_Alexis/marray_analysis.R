#' ---
#' title: "Microarrays analysis Marina Carrera"
#' output:
#' html_document: default
#' ---
#' ## Summary of the project
#' Analysis of microarray data of GSE112497 using learned tools during Omics course. 
#' 
#' ### Summary of the data 	
#' Analysis of gene expression in human CAFs with control (shCtrl) or shRNA against NNMT (shNNMT) expression or normal fibroblasts expressing control (Ctrl) or NNMT overexpression (NNMT) constructs. Cells were infected with lentivirus to express the indicated shRNA construct. Hypothesis is that knockdown of NNMT will affect expression of genes via regulation of histone methylation status.
#'   	
#' ### Overall design	
#' Total RNA was collected from human cancer associated fibroblasts expressing the indicated shCtrl or shNNMT constructs or from normal WI-38 fibroblasts expressing Ctrl or NNMT constructs and expression analyzed with Affymetrix microarray.
#'   	
#' 
#' ## Citations
#' Eckert MA, Coscia F, Chryplewicz A, Chang JW et al. Proteomics reveals NNMT as a master metabolic regulator of cancer-associated fibroblasts. Nature 2019 May;569(7758):723-728. PMID: 31043742
#' 
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#### Load the needed packages before using the code, those can be found in the script
## ----packages, include=FALSE---------------------------------------------
if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}

installifnot("pd.mogene.1.0.st.v1")
installifnot("mogene10sttranscriptcluster.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")

#' #### Targets summary
#' There are 12 samples which can be joined into 4 groups: ShCtrl, ShNNMT, Ctrl, NNMT
## ----Read targets--------------------------------------------------------
targets <- read.csv('./data/targets.csv', header = T,sep = ',')
targets

## ----include=FALSE-------------------------------------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir, "results")
setwd(resultsDir)

#' Reading cel files and getting into rawData to use later
## ----Read cel, include=F-------------------------------------------------
Cel <- list.celfiles(file.path(dataDir))

rawData <- read.celfiles(file.path(dataDir,Cel))


#' Getting variables for plots
## ----Predifine variables for plots---------------------------------------
sampleNames <- as.character(targets$Short_Name)
sampleColor <- as.character(targets$Colors)

pdf("RawData Plots")

#' #### Quality control by boxplot, PCA and cluster of rawData
#' Here we check the quality of our data looking for outliers, artifacts, etc
## ----Boxplot raw---------------------------------------------------------
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)



#' Data seems to be ok with some variation, but without erros or strange values, some samples seem to look more alike, we will check it with clustering and PCA
#' 
#' Now we cluster data
## ----Cluster Raw---------------------------------------------------------
clust.euclid.averageR <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.averageR, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)

#' Clusters of this plot are supported by the fact that they have samples that look alike in the boxplot.
#' 
#' Now we do PCA of the data
## ----PCA raw-------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

dev.off()
#' PCA shows as that some groups are really clustered like SHNNMT and others are pretty distinct like Ctrl or NNMT.
#' 
#' ####Data normalization and analysis
#' Here we normalize data to get rid of possible artifacts or errors and check how it changes after the procedure
#' 
## ----Normalize, include=F------------------------------------------------
eset<-rma(rawData)


pdf("Normalized data plots")

#' Redo boxplot, Clustering and PCA
## ----Boxplot normalized--------------------------------------------------
#BOXPLOT
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)


#' We can see that data has close to no differences between the observations now, so we may have lost some statistical information, but also some possible artifacts are also lost, so it seems to be worth it
#' 
#' Now we redo clustering and check for the groups
## ----Cluster normalized--------------------------------------------------
clust.euclid.averageN <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.averageN, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)

#' Now clusters are formed by different samples of the same observations, which means that the difference between the samples is not big, but between the groups is bigger. Also Sh are closer than not Sh
#' 
#' We repeat PCA
## ----PCA normalized------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(eset), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

dev.off()
#' Now samles of the same experiment are togetherm which means that the expression of genes between the types of experiments which is good.
#' 
#' ####Array quality matrices
#' Here we check for the corectness of data once more to verify the previous statements 
## ----Filtering, include=F------------------------------------------------
arrayQualityMetrics(rawData,  reporttitle="Quality Control of Data", force=TRUE)
annotation(eset) <- "org.Mm.eg.db"

eset_filtered <- nsFilter(eset, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)

#' Number of genes that are out
## ----Filter out----------------------------------------------------------
print(eset_filtered$filter.log$numLowVar)

#' 
#' Genes that are still used
## ----Filter characteristics----------------------------------------------
print(eset_filtered$eset)

#' 
#' ####Design, plots, regression
#' This is where we are going to create design matrix, contrast matrix and liner model
## ----Start analysis------------------------------------------------------
Groups <- targets$Grups
types <- factor(Groups)
matrix_design <-model.matrix(~0+types)
colnames(matrix_design)<-c("shCtrl","shNNMT","Ctrl","NNMT")
rownames(matrix_design) <- sampleNames
print(matrix_design)


#' 
## ----Contrast------------------------------------------------------------
matrix_contrasts <- makeContrasts (
   	shCtrlvsshNNMT = shCtrl-shNNMT,
  shCtrlvsCtrl = shCtrl-Ctrl,
   shCtrlvsNNMT= shCtrl-NNMT,
 shNNMTvsCtrl = shNNMT-Ctrl,
shNNMTvsNNMT  = shNNMT-NNMT,
  CtrlvsNNMT=Ctrl-NNMT ,
  levels=matrix_design)

matrix_contrasts

#' Fir linear model with our data
## ----Fit model-----------------------------------------------------------
fit1 <- lmFit(eset_filtered$eset, matrix_design)
fit.main1 <- contrasts.fit(fit1, matrix_contrasts)
fit.main1


#' 
#' Accounting for bayes to improve the model
## ----adjust bayes--------------------------------------------------------
fit.main1 <- eBayes(fit.main1)
fit.main1

#' 
#' ####Toptables of the genes
#' Once we have fitted our data into a linear model we create a topTable for each comparison to see what genes are differentially expressed between the observations
## ----topTables-----------------------------------------------------------
topTab1_shCtrlvsshNNMT <- topTable (fit.main1, number=nrow(fit.main1), coef="shCtrlvsshNNMT", adjust="fdr"); 
topTab2_shCtrlvsCtrl <- topTable (fit.main1, number=nrow(fit.main1), coef="shCtrlvsCtrl", adjust="fdr"); 
topTab3_shCtrlvsNNMT <- topTable (fit.main1, number=nrow(fit.main1), coef="shCtrlvsNNMT", adjust="fdr");
topTab4_shNNMTvsCtrl <- topTable (fit.main1, number=nrow(fit.main1), coef="shNNMTvsCtrl", adjust="fdr"); 
topTab5_shNNMTvsNNMT<- topTable (fit.main1, number=nrow(fit.main1), coef="shNNMTvsNNMT", adjust="fdr"); 
topTab6_CtrlvsNNMT <- topTable (fit.main1, number=nrow(fit.main1), coef="CtrlvsNNMT", adjust="fdr");



#' 
#' 
#' 
#' ####Volcano plot
#' Volcanoplot for the visualization of data
## ----Volcano-------------------------------------------------------------
volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(matrix_contrasts), sep="\n"))
abline(v = c(-3, 3))





pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep = "\n"))
abline(v = c(-3, 3))
dev.off()



## ------------------------------------------------------------------------
my_frame <- data.frame(exprs(eset))

HMdata <- merge(my_frame, topTab1_shCtrlvsshNNMT, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
HMdata1 <- data.matrix(HMdata, rownames.force=TRUE)

HMdata <- merge(my_frame, topTab2_shCtrlvsCtrl, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
HMdata2 <- data.matrix(HMdata, rownames.force=TRUE)

HMdata <- merge(my_frame, topTab3_shCtrlvsNNMT, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
HMdata3 <- data.matrix(HMdata, rownames.force=TRUE)

HMdata <- merge(my_frame, topTab4_shNNMTvsCtrl, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
HMdata4 <- data.matrix(HMdata, rownames.force=TRUE)

HMdata<- merge(my_frame, topTab5_shNNMTvsNNMT, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
HMdata5 <- data.matrix(HMdata, rownames.force=TRUE)

HMdata <- merge(my_frame, topTab6_CtrlvsNNMT, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
HMdata6 <- data.matrix(HMdata, rownames.force=TRUE)

#HEATMAP PLOT
my_palette <- colorRampPalette(c("yellow2", "darkblue"))(n = 299)
png("HeatMap1")
heatmap.2(HMdata1,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap1",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()
png("HeatMap2")
heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap2",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()
png("HeatMap3")
heatmap.2(HMdata3,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap3",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()
png("HeatMap4")
heatmap.2(HMdata4,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap4",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()
png("HeatMap5")
heatmap.2(HMdata5,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap5",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()
png("HeatMap6")
heatmap.2(HMdata6,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap6",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()


#Annotations



BiocManager::install("pd.clariom.s.human")



BiocManager::install("clariomshumantranscriptcluster.db")

require("clariomshumantranscriptcluster.db")

columns(clariomshumantranscriptcluster.db)




probes_tot<-rownames(unique(topTab1_shCtrlvsshNNMT))

write.csv(select(clariomshumantranscriptcluster.db,probes_tot,
                 columns = c("SYMBOL","GENENAME","PROBEID","REFSEQ","PFAM")),"ann_shCtrlvsshNNMT.csv")

probes_tot<-rownames(unique(topTab2_shCtrlvsCtrl))

write.csv(select(clariomshumantranscriptcluster.db,probes_tot,
                 columns = c("SYMBOL","GENENAME","PROBEID","REFSEQ","PFAM")),"ann_shCtrlvsCtrl.csv")


probes_tot<-rownames(unique(topTab3_shCtrlvsNNMT))

write.csv(select(clariomshumantranscriptcluster.db,probes_tot,
                 columns = c("SYMBOL","GENENAME","PROBEID","REFSEQ","PFAM",'GO')),"ann_shCtrlvsNNMT.csv")


probes_tot<-rownames(unique(topTab4_shNNMTvsCtrl))

write.csv(select(clariomshumantranscriptcluster.db,probes_tot,
                 columns = c("SYMBOL","GENENAME","PROBEID","REFSEQ","PFAM","GO")),"ann_shNNMTvsCtrl.csv")


probes_tot<-rownames(unique(topTab5_shNNMTvsNNMT))

write.csv(select(clariomshumantranscriptcluster.db,probes_tot,
                 columns = c("SYMBOL","GENENAME","PROBEID","REFSEQ","PFAM","GO")),"ann_shNNMTvsNNMT.csv")


probes_tot<-rownames(unique(topTab6_CtrlvsNNMT))

write.csv(select(clariomshumantranscriptcluster.db,probes_tot,
                 columns = c("SYMBOL","GENENAME","PROBEID","REFSEQ","PFAM","GO")),"ann_CtrlvsNNMT.csv")


#' 
#' 
#' #### Prepare data for heatmap and make the plot
#' Heatmaps are found in results, the previously observed information is confirmed
#' 
#' #### Annotations
#' Annotations are found in results
