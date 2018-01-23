
WriteDataFrameAsTsv <- function(data.frame.to.save, file.name.path, col.names=NA, row.names=TRUE) {
    file.name.path <- gsub(pattern = " ", replacement = "_", x = file.name.path)
    file.name <- paste0(file.name.path, ".tsv")
    write.table(x = data.frame.to.save, file = file.name, quote = FALSE, sep = "\t", col.names = col.names, row.names = row.names)
    message(file.name, " written on disk as TSV file!\n")
}

ReadDataFrameFromTsv <- function(file.name.path, row.names.col=1, header.flag=TRUE, sep="\t", quote.char="") {
    df <- read.table(file = file.name.path, sep = sep, header=header.flag, row.names = row.names.col, quote = quote.char)
    message(file.name.path, " read from disk!\n")
    return(df)
}

convertGenesViaBiomart <- function(specie=c("hg38", "mm10", "rnor6"), 
                                   attrs=NULL, filter=NULL, filter.values=NULL) 
{
    specie <- match.arg(specie)
    stopifnot((length(attrs) != 0))
    stopifnot( is.null(filter) || (!is.null(filter.values)) )
    stopifnot( (!is.null(filter)) || is.null(filter.values) )
    
    if(length(which(attrs %in% filter))==0) 
    {
        stop("please use a filter matching the attributes!")
    }
    
    switch (specie,
            "mm10" = { ds <- "mmusculus_gene_ensembl"},
            "hg38" = { ds <- "hsapiens_gene_ensembl"},
            "rnor6" = { ds <- "rnorvegicus_gene_ensembl"}
    )
    
    mart <- biomaRt::useMart("ensembl", dataset=ds)
    
    # listAttributes(mart)[grep(external", listAttributes(mart)[,1]),1]
    # listFilters(mart)[grep("external", listFilters(mart)[,1]),1]
    # attrs <- c("ensembl_gene_id", "external_gene_name", "entrezgene")
    gene.map <- biomaRt::getBM(attributes=attrs, mart=mart, 
                               filters=filter, values=filter.values)
    if(!is.null(filter)){
        idx.dp <- which(duplicated(gene.map[[filter]]))
        if(length(idx.dp) > 0 ) {
            gene.map <- gene.map[-idx.dp,]
        }
        gene.map <- gene.map[order(gene.map[[filter]]),]
    }
    
    return(gene.map)
}


attachGeneColumnToDf <- function(mainDf, genesMap, 
                                 rowNamesIdentifier=c("entrezgene", "ensembl", "symbol"),
                                 mapFromIdentifier=NULL, mapToIdentifier=NULL)
{
    match.arg(rowNamesIdentifier)
    stopifnot(!is.null(mapFromIdentifier))
    stopifnot(!is.null(mapToIdentifier))
    
    if(rowNamesIdentifier=="entrezgene") 
    {
        rownames <- as.integer(rownames(mainDf)) ## if rownames == entrezid
    } else {
        rownames <- rownames(mainDf)
    }
    
    mainDf$gene <- NA
    idx.re <- which(rownames %in% genesMap[[mapFromIdentifier]])
    idx.er <- which(genesMap[[mapFromIdentifier]] %in% rownames)
    
    mainDf$gene[idx.re] <- genesMap[[mapToIdentifier]][idx.er]
    noNaMainDf <- mainDf
    idx.na <- which(is.na(mainDf$gene))
    if(length(idx.na) > 0) noNaMainDf <- mainDf[-idx.na,]
    idx.e <- which(noNaMainDf$gene == "")
    if(length(idx.e) > 0) noNaMainDf <- noNaMainDf[-idx.e,]
    # noNaMainDf$genes[duplicated(noNaMainDf$genes)]
    return(noNaMainDf)
    
}
