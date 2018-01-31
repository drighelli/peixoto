## mapping controls

neg.ctrls <- ReadDataFrameFromTsv(file.name.path="data/full_SD_Neg_Control_genes_BMC_genomics.csv", row.names.col=NULL)
head(neg.ctrls)

neg.ctrls.mgi <- neg.ctrls$MGI.Symbol

# neg.map <- convertGenesViaBiomart(specie="mm10", filter="mgi_symbol",
#                             filter.values=neg.ctrls, c("external_gene_name",
#                                 "mgi_symbol", "entrezgene"))


specie="mm10"
filter="mgi_symbol"
filter.values=neg.ctrls.mgi
attrs = c("external_gene_name", "mgi_symbol", "entrezgene")



switch (specie,
        "mm10" = { ds <- "mmusculus_gene_ensembl"},
        "hg38" = { ds <- "hsapiens_gene_ensembl"},
        "rnor6" = { ds <- "rnorvegicus_gene_ensembl"}
)

mart <- biomaRt::useMart("ensembl", dataset=ds)

# listAttributes(mart)[grep(external", listAttributes(mart)[,1]),1]
# listFilters(mart)[grep("external", listFilters(mart)[,1]),1]
# attrs <- c("ensembl_gene_id", "external_gene_name", "entrezgene")
gene.map.total <- biomaRt::getBM(attributes=attrs, mart=mart)
gene.map <- biomaRt::getBM(attributes=attrs, mart=mart, 
                            filters=filter, values=tolower(filter.values), uniqueRows=TRUE)

idx.on.map <- which( tolower(gene.map.total$mgi_symbol) %in% tolower(neg.ctrls.mgi))

gene.map.total


if(!is.null(filter))
{
    idx.dp <- which(duplicated(gene.map[[filter]]))
    if(length(idx.dp) > 0 ) 
    {
        gene.map <- gene.map[-idx.dp,]
    }
    gene.map <- gene.map[order(gene.map[[filter]]),]
}
idx.entr <- grep(pattern="entrezgene", x=colnames(gene.map))
if(length(idx.entr) > 0 )
{
    gene.map[,idx.entr] <- as.character(gene.map[,idx.entr])
}

length(neg.ctrls.mgi)
length(which(neg.ctrls.mgi %in% gene.map$mgi_symbol))
length(which(tolower(neg.ctrls.mgi) %in% tolower(gene.map$mgi_symbol)))

refSeqCountMatrix <- ReadDataFrameFromTsv(file.name.path="data/refSEQ_countMatrix.txt")

rownames(refSeqCountMatrix)

# designMatrix <- ReadDataFrameFromTsv(file.name.path="design/all_samples_short_names.tsv.csv")
# head(designMatrix)