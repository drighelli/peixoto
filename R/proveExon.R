exonmatrix <- ReadDataFrameFromTsv(file.name.path="data/NEW_PFC_exonMatrix.txt", row.names.col=NULL)

annot <- exonmatrix[,c(1:6)]
counts <- exonmatrix[,c(7:dim(exonmatrix)[2])]

y.all <- DGEList(counts=counts, genes=annot)

y.all <- y.all[, colnames(counts)]

exonDesignMatrix <- cbind(colnames(counts), designMatrix)

