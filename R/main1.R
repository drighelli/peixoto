source("R/utilities.R")
source("R/filterLowCounts.R")
source("R/plotFunctions.R")
source("R/VolcanoPlotFunctions.R")
source("R/MAPlotFunctions.R")
source("R/normalizationFunctions.R")
library(plotly)
library(RUVSeq)

###Importing data

# firstcountMatrix <- ReadDataFrameFromTsv(file.name.path="data/PFC_counts.txt")
# dim(firstcountMatrix)
secondcountMatrix <- ReadDataFrameFromTsv(file.name.path="data/NEW_PFCmatrix.txt")
dim(secondcountMatrix)

refSeqCountMatrix <- ReadDataFrameFromTsv(file.name.path="data/refSEQ_countMatrix.txt")


designMatrix <- ReadDataFrameFromTsv(file.name.path="design/all_samples_short_names.tsv.csv")
head(designMatrix)

filteredCountsProp <- filterLowCounts(counts.dataframe=refSeqCountMatrix, 
                                      is.normalized=FALSE,
                                      design.dataframe=designMatrix,
                                      cond.col.name="gcondition",
                                      method.type="Proportion")#14531 


# rownames(filteredCountsProp)
# ind.row.entr <- which(rownames(filteredCountsProp) %in% gene.map.mm9$entrezgene)#13523
# 
# ind.gm.entr <- which(gene.map.mm9$entrezgene %in% rownames(filteredCountsProp))


# library(biomaRt)
# 
# names <- as.integer(rownames(filteredCountsProp)[order(rownames(filteredCountsProp))])
# 
# mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl', 
#                 host="may2012.archive.ensembl.org") ## using mm9
# listAttributes(mart)[grep("entrez", listAttributes(mart)[,1]),1]
# attrs <- c("ensembl_gene_id", "external_gene_id", "entrezgene")
# gene.map <- getBM(attributes=attrs, mart=mart, filters="entrezgene", 
#                   values=names)
# 
# ## gene.map
# gm.no.dp <- gene.map[-which(duplicated(gene.map$entrezgene)),]
# ind.gm <- which(gm.no.dp$entrezgene %in% names)
# gm.no.dpigm <- gm.no.dp[ind.gm,]
# gm.no.dpigm.o <- gm.no.dpigm[order(gm.no.dpigm$entrezgene),]
# 
# ## names
# ind.nm <- which(names %in% gene.map$entrezgene)
# names.gm <- names[ind.nm]
# names.gm.o <- names.gm[order(names.gm)]
# ndf.map <- cbind.data.frame(gm.no.dpigm.o,
#                     names.gm.o)
# 
# filtCountsProp.o <- filteredCountsProp[order(as.integer(rownames(filteredCountsProp))),]
# 
# ind.rn <- which(as.integer(rownames(filtCountsProp.o)) %in% ndf.map$names.gm.o)
# 
# filtCountsProp.o$map <- NA
# filtCountsProp.o$ensmap <- NA
# filtCountsProp.o$map[ind.rn] <- ndf.map$names.gm.o
# filtCountsProp.o$ensmap[ind.rn] <- ndf.map$ensembl_gene_id


# mapNames <- function(names, map, from.col, to.col)
# {
#     
# }

# gene.map.mm10 <- readRDS("data/gene_map_mm10.RDS")
# length(which(rownames(filteredCountsProp) %in% gene.map.mm10$entrezgene))#13523 
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
neg.ctrls <- ReadDataFrameFromTsv(file.name.path="data/negative_controls.tsv", row.names.col=NULL)
head(neg.ctrls)
neg.ctrls.ens <- neg.ctrls$To

neg.ctrls.est <- neg.ctrls.ens[which(neg.ctrls.ens %in% filtCountsProp.o$ensmap)]

pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua[neg.ctrls.est,]), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua[neg.ctrls.est,]), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua[neg.ctrls.est,]), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA-Norm")
plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)



###RUVg Normalization
# 
# ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=filteredCountsProp,
#                                            design.matrix=designMatrix,
#                                            desMatColStr="gcondition",
#                                            estimated.gene.names=neg.ctrls.est,
#                                            k=1,
#                                            isLog=FALSE)
# 
# normExprData <- ruvedExprData@assayData$normalizedCounts
# pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="RUV-uq-Norm")
# pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="RUV-uq-Norm")
# pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="RUV-uq-Norm")
# plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)
# 
# ###Upper Quartile + RUVg Normalization
# 
# ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=round(normPropCountsUqua),
#                                            design.matrix=designMatrix,
#                                            desMatColStr="gcondition",
#                                            estimated.gene.names=neg.ctrls.est,
#                                            k=1,
#                                            isLog=FALSE)
# 
# normExprData <- ruvedExprData@assayData$normalizedCounts
# pc1_2 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC2", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
# pc2_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC2", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
# pc1_3 <- PlotPCAPlotlyFunction(counts.data.frame=log1p(normExprData), design.matrix=designMatrix, shapeColname="classic", colorColname="gcondition", xPCA="PC1", yPCA="PC3", plotly.flag=TRUE, show.plot.flag=FALSE, prefix.plot="UQUA+RUV-Norm")
# plotly::subplot(pc1_2, pc2_3, pc1_3, nrows=2, margin = 0.1, titleX=TRUE, titleY=TRUE)

###Upper Quartile + RUVs Normalization

library(RUVSeq)
groups <- makeGroups(paste0(designMatrix$genotype, designMatrix$classic))[c(1, 3),]
ruvedSExprData <- RUVs(as.matrix(round(normPropCountsUqua)), cIdx=neg.ctrls.est,
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
source("R/MAPlotFunctions.R")
desMat <- cbind(designMatrix, ruvedSExprData$W)
colnames(desMat) <- c(colnames(designMatrix), colnames(ruvedSExprData$W))

cc <- c("KOSD5 - KOHC5", "KORS2 - KOHC7", "WTSD5 - WTHC5", "WTRS2 - WTHC7")

rescList1 <- applyEdgeR(counts=normExprData, design.matrix=desMat,
                        factors.column="gcondition", 
                        weight.columns=c("W_1", "W_2", "W_3", "W_4"),
                        contrasts=cc, useIntercept=FALSE, p.threshold=1,
                        verbose=TRUE)

    h <- PlotHistPvalPlot(rescList1[[i]], design.matrix=desMat, show.plot.flag=FALSE, 
                     plotly.flag=FALSE, save.plot=FALSE, 
                     prefix.plot=names(rescList1)[i])


for(i in 1:length(rescList1))
{
    filename <- paste0(names(rescList1)[i], "_edgeR")
    WriteDataFrameAsTsv(data.frame.to.save=rescList1[[i]], file.name.path=file.path("res", filename))
}

    
## mapping ensembl gene id using biomart
library(biomaRt)

names <- as.integer(rownames(rescList1[[1]])[order(rownames(rescList1[[1]]))])

mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl', 
                host="may2012.archive.ensembl.org") ## using mm9
listAttributes(mart)[grep("entrez", listAttributes(mart)[,1]),1]
attrs <- c("ensembl_gene_id", "external_gene_id", "entrezgene")
gene.map <- getBM(attributes=attrs, mart=mart, filters="entrezgene", 
                  values=names)

## gene.map
gm.no.dp <- gene.map[-which(duplicated(gene.map$entrezgene)),]
ind.gm <- which(gm.no.dp$entrezgene %in% names)
gm.no.dpigm <- gm.no.dp[ind.gm,]
gm.no.dpigm.o <- gm.no.dpigm[order(gm.no.dpigm$entrezgene),]

## names
ind.nm <- which(names %in% gene.map$entrezgene)
names.gm <- names[ind.nm]
names.gm.o <- names.gm[order(names.gm)]
ndf.map <- cbind.data.frame(gm.no.dpigm.o,
                    names.gm.o)

res.o <- rescList1[[1]][order(as.integer(rownames(rescList1[[1]]))),]

ind.rn <- which(as.integer(rownames(res.o)) %in% ndf.map$names.gm.o)
res.o$map <- NA
# filtCountsProp.o$ensmap <- NA
res.o$map[ind.rn] <- ndf.map$external_gene_id
# res.o$ensmap[ind.rn] <- ndf.map$ensembl_gene_idres.o
# filtCountsProp.onna <- filtCountsProp.o[-which(is.na(filtCountsProp.o$ensmap))]

# rownames(res.o) <- res.o$map

## plotting

PlotVolcanoPlot(de.results=res.o, counts.dataframe=normExprData, design.matrix=desMat,
                show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[1], threshold=0.05)
PlotMAPlotCounts(de.results=res.o, counts.dataframe=normExprData, design.matrix=desMat,
                 show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[1], threshold=0.05)


PlotVolcanoPlot(de.results=rescList1[[2]], counts.dataframe=normExprData, design.matrix=desMat,
                show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[2], threshold=0.05)
PlotMAPlotCounts(de.results=rescList1[[2]], counts.dataframe=normExprData, design.matrix=desMat,
                 show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[2], threshold=0.05)


PlotVolcanoPlot(de.results=rescList1[[3]], counts.dataframe=normExprData, design.matrix=desMat,
                show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[3], threshold=0.05)
PlotMAPlotCounts(de.results=rescList1[[3]], counts.dataframe=normExprData, design.matrix=desMat,
                 show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[3], threshold=0.05)


PlotVolcanoPlot(de.results=rescList1[[4]], counts.dataframe=normExprData, design.matrix=desMat,
                show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[4], threshold=0.05)
PlotMAPlotCounts(de.results=rescList1[[4]], counts.dataframe=normExprData, design.matrix=desMat,
                 show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, prefix.plot=names(rescList1)[4], threshold=0.05)



library(biomaRt)

mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl')#, 
                #host="may2012.archive.ensembl.org")
listAttributes(mart)[grep("entrez", listAttributes(mart)[,1]),1]
attrs <- c("ensembl_gene_id", "external_gene_id", "refseq_mrna", "entrezgene")
gene.map <- getBM(attributes=attrs, mart=mart)
head(gene.map)
gene.map[which(!is.na(gene.map$entrezgene)),]

res <- rescList1[[1]]

gene.map$refseq_mrna <- gsub(pattern="NM_", replacement="", gene.map$refseq_mrna)

length(unique(rownames(filteredCountsProp)[which( rownames(countMatrix) %in% gene.map$entrezgene)]))

length(rownames(filteredCountsProp)[which( rownames(filteredCountsProp) %in% gene.map$ensembl_gene_id )])

length(unique(gene.map$external_gene_name[which(gene.map$ensembl_gene_id %in% rownames(filteredCountsProp))]))
gene.map[which( gene.map$ensembl_gene_id %in% rownames(filteredCountsProp)),"external_gene_name"]



