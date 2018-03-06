source("R/includes.R")



designMatrix <- ReadDataFrameFromTsv("design/all_samples_short_names.tsv.csv")
exonmatrix <- ReadDataFrameFromTsv(file.name.path="data/NEW_PFC_exonMatrix.txt", 
                            row.names.col=NULL)

counts <- exonmatrix[,c(7:dim(exonmatrix)[2])]

exonDesignMatrix <- cbind(colnames(counts), designMatrix)
exonDesignMatrix$samples <- rownames(exonDesignMatrix)
exonDesignMatrix$exonsamples <- colnames(counts)
rownames(exonDesignMatrix) <- exonDesignMatrix$exonsamples

exonDesignMatrix <- exonDesignMatrix[-c(grep("RS2",exonDesignMatrix$condition)),]
exonDesignMatrix<- exonDesignMatrix[-c(grep("HC7",exonDesignMatrix$condition)),]


exonmatrix <- exonmatrix[, c(1:6, which(colnames(exonmatrix) %in% rownames(exonDesignMatrix)))]
colnames(exonmatrix)
annot <- exonmatrix[,c(1:6)]
counts <- exonmatrix[,c(7:dim(exonmatrix)[2])]

y.all <- DGEList(counts=counts, genes=annot)

y.all <- y.all[, colnames(counts)]



# colnames(counts)[designMatrix$gcondition]
# y <- sumTechReps(y.all, designMatrix$gcondition)
y <- y.all
y$samples 
y

colnames(y)

## annotation downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/

ncbi.annotation <- read.delim("data/annotation/Mus_musculus.gene_info", comment.char = "%", 
                    header=TRUE, stringsAsFactors=FALSE)

m <- match(y$genes$GeneID, ncbi$GeneID)
y$genes$Chr <- ncbi$chromosome[m]
y$genes$Symbol <- ncbi$Symbol[m]
y$genes$Strand <- NULL
head(y$genes)

library(dplyr)
library(tidyr)

long_data <- gather(exonmatrix, key = condition, value = expression, 
                S3HC5_1:WTSD5_6)

long_data %>%
    group_by(GeneID, condition) %>%
    summarize(total = sum(expression)) %>%
    group_by(GeneID) %>%
    summarize(pass = sum(total > 20) >= 5) %>%
    filter(pass) %>%
    pull(GeneID) -> keep_ids

keep <- which(exonmatrix$GeneID %in% keep_ids)
length(keep)
dim(exonmatrix)

y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y <- estimateCommonDisp(y)

head(y$pseudo.counts)
plotMDS(y)

## negative controls
library(readxl)

sd.ctrls <- read_excel(path="data/controls/Additional File 4 full list of BMC genomics SD&RS2.xlsx", sheet=1)
sd.ctrls <- sd.ctrls[order(sd.ctrls$adj.P.Val),]

sd.neg.ctrls <- sd.ctrls[sd.ctrls$adj.P.Val > 0.9, ]

sd.neg.ctrls <- sd.neg.ctrls$`MGI Symbol`
sd.neg.ctrls <- sd.neg.ctrls[-which(is.na(sd.neg.ctrls))]

tab <- table(y$genes$Symbol)
single <- names(which(tab==1))
negcon <- which(y$genes$Symbol %in% intersect(sd.neg.ctrls, single))
length(negcon)

exonDesignMatrix$gcondition <- droplevels(exonDesignMatrix$gcondition )

library(RUVSeq)
groups <- makeGroups(exonDesignMatrix$gcondition)
ruvedSExprData <- RUVs(as.matrix(round(y$pseudo.counts)), cIdx=negcon,
                       scIdx=groups, k=5, round=FALSE)

normExprData <- ruvedSExprData$normalizedCounts

plotPCA(normExprData)
plotRLE(normExprData, outline=FALSE)

design <- model.matrix(~0+exonDesignMatrix$gcondition + ruvedSExprData$W)

colnames(design) <- c(levels(exonDesignMatrix$gcondition), paste0("W", 1:5))

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)


fit <- glmFit(y, design)

fit$design


contr <- limma::makeContrasts(contrasts="KOSD5 - WTSD5", levels=design)
qlf <- glmLRT(fit, contrast=contr)
topTags(qlf, p.value=0.05, n=20)

# tt <- topTags(qlf,)
# sum( tt$table$FDR < 0.05 )

sp <- diffSpliceDGE(fit, contrast=contr, geneid="GeneID", exonid="Start")
topSpliceDGE(sp, test="Simes", FDR=0.05)

plotSpliceDGE(sp, geneid="Shank3", genecol="Symbol")
plotSpliceDGE(sp, geneid="Cntn1", genecol="Symbol")
plotSpliceDGE(sp, geneid="Rpl3", genecol="Symbol")
plotSpliceDGE(sp, geneid="Tuba1b", genecol="Symbol")

plotSpliceDGE(sp, geneid="Fabp7", genecol="Symbol")
plotSpliceDGE(sp, geneid="Hdac7", genecol="Symbol")
plotSpliceDGE(sp, geneid="Per3", genecol="Symbol")

plotSpliceDGE(sp, geneid="Kalrn", genecol="Symbol")
plotSpliceDGE(sp, geneid="Sgk1", genecol="Symbol")
plotSpliceDGE(sp, geneid="Dido1", genecol="Symbol")
