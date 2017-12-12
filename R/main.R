source("R/utilities.R")
source("R/filterLowCounts.R")
source("R/plotFunctions.R")
source("R/normalizationFunctions.R")
countMatrix <- ReadDataFrameFromTsv(file.name.path="data/PFC_counts.txt")
head(countMatrix)

designMatrix <- ReadDataFrameFromTsv(file.name.path="design/all_samples.tsv.csv")
head(designMatrix)

filteredCountsProp <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
                                    design.dataframe=designMatrix,
                                    cond.col.name="gcondition",
                                    method.type="Proportion")
head(filteredCountsProp) #15823
# 
# filteredCountsWilc <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
#                                       design.dataframe=designMatrix,
#                                       cond.col.name="gcondition",
#                                       method.type="Wilcoxon")
# head(filteredCountsWilc) #0
# 
# 
# filteredCountsCpm <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
#                                       design.dataframe=designMatrix,
#                                       cond.col.name="gcondition",
#                                       method.type="CPM") 
# head(filteredCountsCpm) ##16190

PlotPCAPlotlyFunction(counts.data.frame=log1p(filteredCountsProp), 
                      design.matrix=designMatrix, shapeColname="genotype",
                      colorColname="gcondition", save.plot=FALSE,
                      plotly.flag=TRUE, show.plot.flag=TRUE)



normPropCountsUqua <- NormalizeData(data.to.normalize=filteredCountsProp, norm.type="uqua", 
                                    design.matrix=designMatrix)


PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua), 
                      design.matrix=designMatrix, shapeColname="condition",
                      colorColname="genotype", save.plot=FALSE,
                      plotly.flag=FALSE, show.plot.flag=TRUE, prefix.plot="uqua")


neg.ctrls <- ReadDataFrameFromTsv(file.name.path="data/negative_controls.tsv", row.names.col=NULL)
head(neg.ctrls)

neg.ctrls.ens <- as.character(neg.ctrls$To)

neg.ctrls.est <- neg.ctrls.ens[which(neg.ctrls.ens %in% rownames(filteredCountsProp))]
library("RUVSeq")
ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=filteredCountsProp,
                                           design.matrix=designMatrix,
                                           desMatColStr="gcondition",
                                           estimated.gene.names=neg.ctrls.est,
                                           k=1,
                                           isLog=FALSE)

ruv.counts <- ruvedExprData@assayData$normalizedCounts

PlotPCAPlotlyFunction(counts.data.frame=log1p(ruv.counts), 
                      design.matrix=designMatrix, shapeColname="condition",
                      colorColname="genotype", save.plot=FALSE,
                      plotly.flag=FALSE, show.plot.flag=TRUE, prefix.plot="RUVg")

## uqua + ruv
ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=round(normPropCountsUqua),
                                           design.matrix=designMatrix,
                                           desMatColStr="gcondition",
                                           estimated.gene.names=neg.ctrls.est,
                                           k=1,
                                           isLog=FALSE)

ruv.counts <- ruvedExprData@assayData$normalizedCounts

PlotPCAPlotlyFunction(counts.data.frame=log1p(ruv.counts), 
                      design.matrix=designMatrix, shapeColname="condition",
                      colorColname="genotype", save.plot=FALSE,
                      plotly.flag=FALSE, show.plot.flag=TRUE, prefix.plot="uqua+RUVg",
                      xPCA="PC2", yPCA="PC3")
