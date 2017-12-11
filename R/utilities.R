
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
