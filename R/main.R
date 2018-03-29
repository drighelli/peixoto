source("R/utilities.R")
source("R/filterLowCounts.R")
source("R/plotFunctions.R")
source("R/normalizationFunctions.R")
countMatrix <- ReadDataFrameFromTsv(file.name.path="data/refSEQ_countMatrix.txt")
head(countMatrix)
dim(countMatrix)

designMatrix <- ReadDataFrameFromTsv(file.name.path="design/all_samples_short_names.tsv.csv")
head(designMatrix)


filteredCountsProp <- filterLowCounts(counts.dataframe=countMatrix, is.normalized=FALSE,
                                    design.dataframe=designMatrix,
                                    cond.col.name="gcondition",
                                    method.type="Proportion")
# head(filteredCountsProp) #14531
# dim(filteredCountsProp)


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
# 
# g=PlotPCAPlotlyFunction(counts.data.frame=log1p(filteredCountsProp), 
#                       design.matrix=designMatrix, shapeColname="genotype",
#                       colorColname="gcondition", save.plot=FALSE,
#                       plotly.flag=TRUE, show.plot.flag=TRUE)
# plotly::ggplotly(g)

normPropCountsUqua <- NormalizeData(data.to.normalize=filteredCountsProp, 
                                    norm.type="uqua", 
                                    design.matrix=designMatrix)


g=PlotPCAPlotlyFunction(counts.data.frame=log1p(normPropCountsUqua), 
                      design.matrix=designMatrix, shapeColname="condition",
                      colorColname="genotype", save.plot=FALSE,
                      plotly.flag=FALSE, show.plot.flag=TRUE, prefix.plot="uqua")

plotly::ggplotly(g)


# neg.ctrls.ens <- as.character(neg.ctrls$To)

# neg.ctrls.est <- neg.ctrls.ens[which(neg.ctrls.ens %in% rownames(filteredCountsProp))]
# library("RUVSeq")
# ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=filteredCountsProp,
#                                            design.matrix=designMatrix,
#                                            desMatColStr="gcondition",
#                                            estimated.gene.names=neg.ctrls.est,
#                                            k=1,
#                                            isLog=FALSE)
# 
# ruv.counts <- ruvedExprData@assayData$normalizedCounts
# 
# gr <- PlotPCAPlotlyFunction(counts.data.frame=log1p(ruv.counts), 
#                       design.matrix=designMatrix, shapeColname="condition",
#                       colorColname="genotype", save.plot=FALSE,
#                       plotly.flag=FALSE, show.plot.flag=TRUE, prefix.plot="RUVg")
# plotly::ggplotly(gr)
# 
# ## uqua + ruv
ruvedExprData <- RUVgNormalizationFunction(data.to.normalize=round(normPropCountsUqua),
                                           design.matrix=designMatrix,
                                           desMatColStr="gcondition",
                                           estimated.gene.names=neg.ctrls.entrez,
                                           k=1,
                                           isLog=FALSE)
# 
# ruv.counts <- ruvedExprData@assayData$normalizedCounts
# 
# gqr <- PlotPCAPlotlyFunction(counts.data.frame=log1p(ruv.counts), 
#                       design.matrix=designMatrix, shapeColname="condition",
#                       colorColname="genotype", save.plot=FALSE,
#                       plotly.flag=FALSE, show.plot.flag=TRUE, prefix.plot="uqua+RUVg",
#                       xPCA="PC2", yPCA="PC3")
# plotly::ggplotly(gqr)

## converting genes
# 
# neg.ctrls <- ReadDataFrameFromTsv(file.name.path="data/full_SD_Neg_Control_genes_BMC_genomics.csv", row.names.col=NULL)
# head(neg.ctrls)
# 

rownames <- as.integer(rownames(countMatrix))

rownames <- rownames[order(rownames)]
rownames.map <- convertGenesViaBiomart(specie="mm10", filter="entrezgene",
                    filter.values=rownames, c("external_gene_name",
                                "mgi_symbol", "entrezgene"))
# 
# countMatrix$symbols <- NA
# idx.re <- which(rownames %in% rownames.map$entrezgene)
# idx.er <- which(rownames.map$entrezgene %in% rownames)
# 
# 
# ## with mgi symb
# # countMatrix$symbols[idx.re] <- rownames.map$mgi_symbol[idx.er]
# # noNaCountMatrix <- countMatrix[-which(is.na(countMatrix$symbols)),]
# # noNaCountMatrix <- noNaCountMatrix[-which(noNaCountMatrix$symbols == ""),]
# # noNaCountMatrix$symbols[duplicated(noNaCountMatrix$symbols)]
# # dim(noNaCountMatrix) ## 20715
# 
# ## with ext gene symb
# countMatrix$symbols[idx.re] <- rownames.map$external_gene_name[idx.er]
# noNaCountMatrix <- countMatrix
# idx.na <- which(is.na(countMatrix$symbols))
# if(length(idx.na) > 0) noNaCountMatrix <- countMatrix[-idx.na,]
# idx.e <- which(noNaCountMatrix$symbols == "")
# if(length(idx.e) > 0) noNaCountMatrix <- noNaCountMatrix[-idx.e,]
# noNaCountMatrix$symbols[duplicated(noNaCountMatrix$symbols)]
# dim(noNaCountMatrix) ##20730



rownames <- as.integer(rownames(countMatrix))

rownames <- rownames[order(rownames)]
rownames.map <- convertGenesViaBiomart(specie="mm10", filter="entrezgene",
                                       filter.values=rownames, c("external_gene_name",
                                                                 "mgi_symbol", "entrezgene"))


noNaCountMatrix <- attachGeneColumnToDf(mainDf=countMatrix,
                                genesMap=rownames.map,
                                rowNamesIdentifier="entrezgene",
                                mapFromIdentifier="entrezgene",
                                mapToIdentifier="external_gene_name")
dim(noNaCountMatrix)

filteredCountsProp <- filterLowCounts(counts.dataframe=noNaCountMatrix, 
                                    is.normalized=FALSE,
                                    design.dataframe=designMatrix,
                                    cond.col.name="gcondition",
                                    method.type="Proportion")
dim(filteredCountsProp)

normPropCountsUqua <- NormalizeData(data.to.normalize=filteredCountsProp, 
                                    norm.type="uqua",
                                    design.matrix=designMatrix)


## neg controls

neg.ctrls <- ReadDataFrameFromTsv(file.name.path="data/full_SD_Neg_Control_genes_BMC_genomics.csv", row.names.col=NULL)
head(neg.ctrls)

neg.ctrls <- neg.ctrls$MGI.Symbol

neg.map <- convertGenesViaBiomart(specie="mm10", filter="mgi_symbol",
                    filter.values=neg.ctrls, c("external_gene_name",
                                "mgi_symbol", "entrezgene"))

neg.map.nna <- neg.map[-which(is.na(neg.map$entrezgene)),]

neg.map.nna <- neg.map.nna[which(neg.map.nna$entrezgene %in%  rownames(normPropCountsUqua)),]

neg.ctrls.entrez <- neg.map.nna$entrezgene



filteredCountsPropGenes <- attachGeneColumnToDf(mainDf=filteredCountsProp,
                                                genesMap=rownames.map,
                                                rowNamesIdentifier="entrezgene",
                                                mapFromIdentifier="entrezgene",
                                                mapToIdentifier="external_gene_name")

dim(filteredCountsPropGenes)
# 
# library(biomaRt)
# 
# using.names <- neg.ctrls$MGI.Symbol #2124
# which(is.na(using.names))
# 
# mart = useMart('ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl')
# listAttributes(mart)[grep("external", listAttributes(mart)[,1]),1]
# listFilters(mart)[grep("external", listFilters(mart)[,1]),1]
# attrs <- c("ensembl_gene_id", "external_gene_name", "entrezgene")
# gene.map <- getBM(attributes=attrs, mart=mart, 
#                 filters="external_gene_name", values=using.names) ## 2112x3
# 
# gene.map.u <- gene.map[-which(duplicated(gene.map$external_gene_name)), ]
# dim(gene.map.u) ## 2104
# 
# 
# entrez.names <- as.integer(rownames(filteredCountsProp))
# attrs <- c("ensembl_gene_id", "external_gene_name", "entrezgene")
# gene.map.names <- getBM(attributes=attrs, mart=mart, 
#                   filters="entrezgene", values=entrez.names) 
# 
# dim(gene.map.names)


## positive controls 
pos.ctrls <- ReadDataFrameFromTsv(file.name.path="data/SD_full_pos_control_genes_BMC_genomics.csv", row.names.col=NULL)
head(pos.ctrls)
dim(pos.ctrls) #641
length(which(is.na(pos.ctrls$MGI.Symbol))) ##48
pos.ctrls.names <- pos.ctrls$MGI.Symbol
pos.ctrls.names <- pos.ctrls.names[-which(is.na(pos.ctrls.names))] ##593

gene.map.posc <- getBM(attributes=attrs, mart=mart, 
                        filters="external_gene_name", values=pos.ctrls.names) 
dim(gene.map.posc)
gene.map.posc<-gene.map.posc[-which(duplicated(gene.map.posc$entrezgene)),]
dim(gene.map.posc) ##560


## intersections names vs pos genes
sum(entrez.names %in% gene.map.posc$entrezgene ) ## 554 (of 560)

## intersections names vs neg genes
sum(entrez.names %in% gene.map.u$entrezgene ) ## 1224 (of 2104)



