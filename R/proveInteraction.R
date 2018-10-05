source("R/includes.R")
library(plotly)
library(RUVSeq)


countMatrix <- ReadDataFrameFromTsv(file.name.path="data/refSEQ_countMatrix.txt")
# head(countMatrix)

designMatrix <- ReadDataFrameFromTsv(file.name.path="design/all_samples_short_names_noRS2HC7.tsv")
head(designMatrix)

rownames <- as.character(rownames(countMatrix))

rownames <- rownames[order(rownames)]
rownames.map <- convertGenesViaBiomart(specie="mm10", filter="entrezgene",
                                       filter.values=rownames, attrs=c("external_gene_name",
                                                                       "mgi_symbol", "entrezgene"))

noNaCountMatrix <- attachGeneColumnToDf(mainDf=countMatrix,
                                        genesMap=rownames.map,
                                        rowNamesIdentifier="entrezgene",
                                        mapFromIdentifier="entrezgene",
                                        mapToIdentifier="external_gene_name")


filteredCountsProp <- filterLowCounts(counts.dataframe=noNaCountMatrix, 
                                      is.normalized=FALSE,
                                      design.dataframe=designMatrix,
                                      cond.col.name="gcondition",
                                      method.type="Proportion")

normPropCountsUqua <- NormalizeData(data.to.normalize=filteredCountsProp, 
                                    norm.type="tmm", 
                                    design.matrix=designMatrix)

library(readxl)
rs2.ctrls <- read_excel(path="data/controls/Additional File 4 full list of BMC genomics SD&RS2.xlsx", sheet=3)
rs2.ctrls <- rs2.ctrls[order(rs2.ctrls$adj.P.Val),]
rs2.neg.ctrls <- rs2.ctrls[rs2.ctrls$adj.P.Val > 0.9, ]
rs2.neg.ctrls <- rs2.neg.ctrls$`MGI Symbol`
rs2.neg.ctrls <- rs2.neg.ctrls[-which(is.na(rs2.neg.ctrls))]
sd.ctrls <- read_excel(path="data/controls/Additional File 4 full list of BMC genomics SD&RS2.xlsx",sheet=1)
sd.ctrls <- sd.ctrls[order(sd.ctrls$adj.P.Val),]
sd.pos.ctrls <- sd.ctrls[sd.ctrls$adj.P.Val < 0.01, ]
sd.neg.ctrls <- sd.ctrls[sd.ctrls$adj.P.Val > 0.9, ]
sd.neg.ctrls <- sd.neg.ctrls$`MGI Symbol`
sd.neg.ctrls <- sd.neg.ctrls[-which(is.na(sd.neg.ctrls))]
int.neg.ctrls <- intersect(rs2.neg.ctrls, sd.neg.ctrls)
neg.map <- convertGenesViaBiomart(specie="mm10", filter="mgi_symbol",
                        filter.values=int.neg.ctrls, c("external_gene_name",
                                            "mgi_symbol", "entrezgene"))
neg.map.nna <- neg.map[-which(is.na(neg.map$entrezgene)),]
neg.ctrls.entrez <- as.character(neg.map.nna$entrezgene)
ind.ctrls <- which(rownames(normPropCountsUqua) %in% neg.ctrls.entrez)
counts.neg.ctrls <- normPropCountsUqua[ind.ctrls,]
library(RUVSeq)
neg.ctrl.list <- rownames(counts.neg.ctrls)
groups <- makeGroups(paste0(designMatrix$genotype, designMatrix$classic))[c(1, 3),]
ruvedSExprData <- RUVs(as.matrix(round(normPropCountsUqua)), cIdx=neg.ctrl.list,
                       scIdx=groups, k=5)
normExprData <- ruvedSExprData$normalizedCounts
#############################

interactionMatrixNew <- constructInteractionMatrix(design.matrix=designMatrix, 
                                                   genotype.col="genotype",
                                                   genotype.ref="WT",
                                                   condition.col="condition",
                                                   weights=ruvedSExprData$W)
# 
# mat <- designMatrix
# geno <- relevel(mat$genotype, ref="WT")
# cond <- mat$condition
# interactionMatrix <- model.matrix(~0+geno*cond+ruvedSExprData$W)
cond <- designMatrix[["condition"]]

fit <- applyEdgeRGLMFit(counts=filteredCountsProp, factors=cond, 
            design=interactionMatrixNew, is.normalized=FALSE, method="TMM",
            verbose=TRUE)

# counts <- filteredCountsProp
# dgel <- edgeR::DGEList(counts=counts, group=cond)
# dgel <- edgeR::calcNormFactors(dgel)
# edisp <- edgeR::estimateDisp(y=dgel, design=interactionMatrix)
# fit1 <- edgeR::glmFit(edisp, interactionMatrix, robust=TRUE)
# class(edisp)


lrt <- applyEdgeRLRT(fit=fit, interaction.matrix=interactionMatrixNew, 
                    verbose=TRUE)
sum(lrt$FDR<0.01)
## 5 SD5 
colnames(fit)
lrtkosd5 <- edgeR::glmLRT(fit, coef=NCOL(interactionMatrix))
top <- topTags(lrtkosd5, n=Inf)
sum(top$table$FDR<0.1)
dim(top$table)
res <- top$table
# sum(top$table$FDR<0.05)
# dim(top$table)


res.o.map <- convertGenesViaBiomart(specie="mm10", filter="entrezgene",
                        filter.values=rownames(res),
                        c("external_gene_name", "mgi_symbol", "entrezgene"))


# WriteDataFrameAsTsv(data.frame.to.save=rescList1[[1]], 
#                     file.name.path=paste0(names(rescList1)[1], "_edgeR"))

res.o <- attachGeneColumnToDf(mainDf=res,
                              genesMap=res.o.map,
                              rowNamesIdentifier="entrezgene",
                              mapFromIdentifier="entrezgene",
                              mapToIdentifier="external_gene_name")

PlotHistPvalPlot(de.results=res.o, design.matrix=designMatrix, 
                 show.plot.flag=TRUE, plotly.flag=TRUE, 
                 prefix.plot="SD5")

PlotVolcanoPlot(de.results=res.o, counts.dataframe=normExprData, 
                design.matrix=designMatrix,
                show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, 
                prefix.plot="SD5", threshold=0.05, 
                positive.ctrls.list=NULL)

## 6 RS2 
colnames(fit)
lrtkosd5 <- edgeR::glmLRT(fit, coef=6)
top <- topTags(lrtkosd5, n=Inf)
sum(top$table$FDR<0.05)
dim(top$table)
res <- top$table
# sum(top$table$FDR<0.05)
# dim(top$table)


res.o.map <- convertGenesViaBiomart(specie="mm10", filter="entrezgene",
                                    filter.values=rownames(res),
                                    c("external_gene_name", "mgi_symbol", "entrezgene"))

res.o <- attachGeneColumnToDf(mainDf=res,
                              genesMap=res.o.map,
                              rowNamesIdentifier="entrezgene",
                              mapFromIdentifier="entrezgene",
                              mapToIdentifier="external_gene_name")

PlotHistPvalPlot(de.results=res.o, design.matrix=designMatrix, 
                 show.plot.flag=TRUE, plotly.flag=TRUE, 
                 prefix.plot="SD5")

PlotVolcanoPlot(de.results=res.o, counts.dataframe=normExprData, 
                design.matrix=designMatrix,
                show.plot.flag=TRUE, plotly.flag=TRUE, save.plot=FALSE, 
                prefix.plot="SD5", threshold=0.05, 
                positive.ctrls.list=NULL)

geneProfileLucia(normalized.counts=normExprData, 
                design.matrix=design.matrix, gene.name="320184",
                show.plot=TRUE, plotly.flag=TRUE)



