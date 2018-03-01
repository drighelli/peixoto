source("R/includes.R")



designMatrix <- ReadDataFrameFromTsv("design/all_samples_short_names.tsv.csv")
exonmatrix <- ReadDataFrameFromTsv(file.name.path="data/NEW_PFC_exonMatrix.txt", 
                            row.names.col=NULL)

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
y<-y.all
y$samples 
y

colnames(y)

## annotation downloaded from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/

ncbi.L1 <- readLines("data/annotation/Mus_musculus.gene_info", n = 1)
ncbi.colname <- unlist(strsplit(substring(ncbi.L1, 10, 234), " "))
ncbi <- read.delim("data/annotation/Mus_musculus.gene_info", skip=1, 
                    header=FALSE, stringsAsFactors=FALSE)

colnames(ncbi) <- ncbi.colname
m <- match(y$genes$GeneID, ncbi$GeneID)
y$genes$Chr <- ncbi$chromosome[m]
y$genes$Symbol <- ncbi$Symbol[m]
y$genes$Strand <- NULL
head(y$genes)
keep <- rowSums(cpm(y) > 1) >=3
summary(keep)

y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
plotMDS(y)

batch <- rep(1:5, 4)
exonDesignMatrix$gcondition<-droplevels(exonDesignMatrix$gcondition )
design <- model.matrix(~0+exonDesignMatrix$gcondition)

colnames(design) <- unique(exonDesignMatrix$gcondition)

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)


fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

fit$design


contr <- limma::makeContrasts(contrasts="KOSD5 - WTSD5", levels=design)
qlf <- glmQLFTest(fit, contrast=contr)

# tt <- topTags(qlf,)
# sum( tt$table$FDR < 0.05 )
