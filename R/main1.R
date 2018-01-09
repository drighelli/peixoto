source("R/utilities.R")
source("R/filterLowCounts.R")
source("R/plotFunctions.R")
source("R/normalizationFunctions.R")
library(plotly)
library(RUVSeq)

###Importing data

countMatrix <- ReadDataFrameFromTsv(file.name.path="../data/NEW_PFCmatrix.txt")
# head(countMatrix)

designMatrix <- ReadDataFrameFromTsv(file.name.path="../design/all_samples.tsv.csv")
head(designMatrix)

filteredCountsProp <- filterLowCounts(counts.dataframe=countMatrix, 
                                      is.normalized=FALSE,
                                      design.dataframe=designMatrix,
                                      cond.col.name="gcondition",
                                      method.type="Proportion")


##Plot PCA of log unnormalized data

pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(filteredCountsProp), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="Prop-Un-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(filteredCountsProp), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="Prop-Un-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(filteredCountsProp), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="Prop-Un-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)



##Normalizations
###Upper Quartile Normalization

normPropCountsUqua <- NormalizeData(data.to.normalize=filteredCountsProp, 
                                    norm.type="uqua", 
                                    design.matrix=designMatrix)

pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)


##Negative control genes
## Using Negative Control Genes to normalize data

neg.ctrls <- ReadDataFrameFromTsv(file.name.path="../data/negative_controls.tsv", row.names.col=NULL)
neg.ctrls.ens <- as.character(neg.ctrls$To)
neg.ctrls.est <- neg.ctrls.ens[which(neg.ctrls.ens %in% rownames(filteredCountsProp))]

pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua[neg.ctrls.est,]), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua[neg.ctrls.est,]), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua[neg.ctrls.est,]), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)



###RUVg Normalization

ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=filteredCountsProp,
                                           design.matrix=designMatrix,
                                           desMatColStr="gcondition",
                                           estimated.gene.names=neg.ctrls.est,
                                           k=1,
                                           isLog=FALSE)

normExprData <- ruvedExprData@assayData$normalizedCounts
pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="RUV-uq-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="RUV-uq-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="RUV-uq-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)

###Upper Quartile + RUVg Normalization

ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=round(normPropCountsUqua),
                                           design.matrix=designMatrix,
                                           desMatColStr="gcondition",
                                           estimated.gene.names=neg.ctrls.est,
                                           k=1,
                                           isLog=FALSE)

normExprData <- ruvedExprData@assayData$normalizedCounts
pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)

###Upper Quartile + RUVs Normalization

library(RUVSeq)
groups <- makeGroups(paste0(designMatrix$genotype, designMatrix$classic))[c(1, 3),]
ruvedSExprData <- RUVs(as.matrix(round(normPropCountsUqua)), cIdx = neg.ctrls.est,
                       scIdx = groups, k = 5)

normExprData <- ruvedSExprData$normalizedCounts

pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)

pal <- RColorBrewer::brewer.pal(9, "Set1")
plotRLE(normExprData, outline=FALSE, col=pal[designMatrix$gcondition])






### edgering
source("R/edgeRFunctions.R")
source("R/VolcanoPlotFunctions.R")
desMat <- cbind(designMatrix, ruvedSExprData$W)
colnames(desMat) <- c(colnames(designMatrix), colnames(ruvedSExprData$W))

cc <- c("KOSD5 - KOHC5", "KORS2 - KOHC7", "WTSD5 - WTHC5", "WTRS2 - WTHC7")

rescList1 <- applyEdgeR(counts=normExprData, design.matrix=desMat,
                        factors.column="gcondition", 
                        weight.columns=c("W_1", "W_2", "W_3", "W_4"),
                        contrasts=cc, useIntercept=FALSE, p.threshold=1,
                        verbose=TRUE)

for(i in 1:length(rescList1))
{
    PlotVolcanoPlot(de.results=rescList1[[i]], counts.dataframe=normExprData, design.matrix=desMat,
                    show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[i], threshold=0.01)
}

