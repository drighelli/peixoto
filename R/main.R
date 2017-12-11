source("R/utilities.R")
source("R/filterLowCounts.R")
source("R/plotFunctions.R")
countMatrix <- ReadDataFrameFromTsv(file.name.path="data/PFC_counts.txt")
head(countMatrix)

designMatrix <- ReadDataFrameFromTsv(file.name.path="design/all_samples.tsv.csv")
head(designMatrix)

filteredCountsProp <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
                                    design.dataframe=designMatrix,
                                    cond.col.name="gcondition",
                                    method.type="Proportion")
head(filteredCountsProp) #15823

filteredCountsWilc <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
                                      design.dataframe=designMatrix,
                                      cond.col.name="gcondition",
                                      method.type="Wilcoxon")
head(filteredCountsWilc) #0


filteredCountsCpm <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
                                      design.dataframe=designMatrix,
                                      cond.col.name="gcondition",
                                      method.type="CPM") 
head(filteredCountsCpm) ##16190

PlotPCAPlotlyFunction(counts.data.frame=filteredCountsProp, 
                      design.matrix=designMatrix, shapeColname="gcondition",
                      colorColname="gcondition", save.plot=FALSE,
                      plotly.flag=TRUE, show.plot.flag=TRUE)



