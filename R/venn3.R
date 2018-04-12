Venn3de <- function(x, y, z, label1, label2, label3, 
                    title="Venn Diagram", 
                    intersection.flag=TRUE, 
                    intersection.exclusion.flag=FALSE, 
                    save.plot=FALSE,
                    plot.dir=NULL, 
                    enrich.lists.flag=FALSE, 
                    conversion.map=NULL)
{
    
    a=x
    b=y
    c=z
    plot.combined.string <- paste(label1, "_", label2, "_", label3, sep="")
    
    if(is.null(plot.dir)) {
        stop("Please provide a directory where to save VENN results!")
    }
    
    Lists <- list(a, b, c)  #put the word vectors into a list to supply lapply  
    Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
    names <- c(label1, label2, label3)                
    colnames(MAT) <- names
    rownames(MAT) <- items
    lapply(seq_along(Lists), function(i) 
    {   #fill the matrix
        MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    })
    
    outputName <- UpdateFilename(filename="VennDiagram", 
                                label1, label2, label3, 
                                extension="pdf")
    out.path.name <- UpdateFolderPath(plot.dir, "Venn3")
    prefix="venn3"
    file.path.name <- file.path(out.path.name, outputName)
    # dev.new()
    require("limma")
    limma::vennDiagram(MAT, circle.col= c("red","green","yellow"), main=title)
    
    if(save.plot)
    {
        dev.new()
        limma::vennDiagram(MAT, circle.col= c("red","green","yellow"), main=title)
        dev.print(device=pdf, file=file.path.name, width=20, height=20)
        graphics.off()
    }
    
    ab <- intersect(a, b)
    bc <- intersect(b, c)
    ac <- intersect(a, c)
    abc <- intersect(ab, bc)
    union.of.intersections <- union(union(ab, bc), union(ac, abc))
    # ret.list <- list(ab, bc, ac, abc, union.of.intersections)
    
    if(intersection.flag) 
    {
        if(is.null(plot.dir)) 
        {
            stop("Please provide a directory where to save VENN results!")
        }
        SaveInteserctionsList(gene.list=ab, conversion.map=conversion.map, 
                              root.dir=out.path.name, 
                              prefix=prefix, 
                              labels.list=c(label1, paste0("AND_",label2)), 
                              enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=bc, conversion.map=conversion.map, 
                              root.dir=out.path.name, 
                              prefix=prefix, 
                              labels.list=c(label2, paste0("AND_", label3)), 
                              enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=ac, conversion.map=conversion.map, 
                              root.dir=out.path.name, 
                              prefix=prefix, 
                              labels.list=c(label1, paste0("AND_", label3)), 
                              enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=abc, conversion.map=conversion.map, 
                              root.dir=out.path.name, 
                              prefix=prefix, 
                              labels.list=c(label1, 
                                              paste0("AND_", label2), 
                                              paste0("AND_", label3)), 
                              enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=union.of.intersections, 
                              conversion.map=conversion.map, 
                              root.dir=out.path.name, 
                              prefix=paste0(prefix, "_union_of_intersections"), 
                              labels.list=c(label1, label2, label3), 
                              enrich.lists.flag=enrich.lists.flag)
        
        ## intentionally left commented
        # union.of.lugs <- union(bc, union(ac, abc))
        # SaveInteserctionsList(gene.list=union.of.lugs, 
        #                       conversion.map=conversion.map, 
        #                       root.dir=out.path.name, 
        #                       prefix=paste0(prefix, "_union_of_LUGs"), 
        #                       labels.list=c(label1, label2, label3), 
        #                       enrich.lists.flag=enrich.lists.flag)
    }
    
    a.not.b <-  setdiff(a, b)
    b.not.a <-  setdiff(b, a)
    c.not.b <-  setdiff(c, b)
    b.not.c <-  setdiff(b, c)
    a.not.c <-  setdiff(a, c)
    c.not.a <-  setdiff(c, a)
    a.not.bc <- setdiff(a.not.b, c)
    b.not.ac <- setdiff(b.not.a, c)
    c.not.ab <- setdiff(c.not.a, b)
    # ret.list <- c(ret.list, 
    #                 a.not.b, 
    #                 b.not.a,
    #                 c.not.b,
    #                 b.not.c,
    #                 a.not.c,
    #                 c.not.a,
    #                 a.not.bc,
    #                 b.not.ac,
    #                 c.not.ab)
    ret.list <- list("XintY"=ab, 
                    "YintZ"=bc, 
                    "XintZ"=ac, 
                    "XintYintZ"=abc, 
                    "unionOfInts"=union.of.intersections,
                    "XnotY"=a.not.b,
                    "YnotX"=b.not.a,
                    "ZnotY"=c.not.b,
                    "YnotZ"=b.not.c,
                    "XnotZ"=a.not.c,
                    "ZnotX"=c.not.a,
                    "XnotYnotZ"=a.not.bc,
                    "YnotXnotZ"=b.not.ac,
                    "ZnotXnotY"=c.not.ab)
    if(intersection.exclusion.flag) 
    {
        if(is.null(plot.dir)) 
        {
            stop("Please provide a directory where to save VENN results!")
        }
        res.path <- file.path(plot.dir, "Venn3", 
                            paste0(plot.combined.string,"_gene_lists"))
        dir.create(res.path, recursive=TRUE, showWarnings=FALSE)
        SaveInteserctionsList(gene.list=a.not.b, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix, 
                            labels.list=c(label1, paste0("_NOT_",label2), 
                            enrich.lists.flag=enrich.lists.flag))
        SaveInteserctionsList(gene.list=b.not.a, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix, 
                            labels.list=c(label2, paste0("_NOT_",label1)), 
                            enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=c.not.b, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix, 
                            labels.list=c(label3, paste0("_NOT_",label2)), 
                            enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=b.not.c, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix, 
                            labels.list=c(label2, paste0("_NOT_", label3)), 
                            enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=a.not.c, conversion.map=conversion.map,
                              root.dir=out.path.name, prefix=prefix, 
                              labels.list=c(label1, paste0("_NOT_", label3)), 
                              enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=c.not.a, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix, 
                            labels.list=c(label3, paste0("_NOT_", label1)), 
                            enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=a.not.bc, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix, 
                            labels.list=c(label1, paste0("_NOT_", label2), 
                                          paste0("_NOT_", label3)), 
                            enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=b.not.ac, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label2,
                                        paste0("_NOT_", label1),
                                        paste0("_NOT_", label3)),
                            enrich.lists.flag=enrich.lists.flag)
        SaveInteserctionsList(gene.list=c.not.ab, conversion.map=conversion.map,
                            root.dir=out.path.name, prefix=prefix,
                            labels.list=c(label3, paste0("_NOT_", label1),
                                        paste0("_NOT_", label2)),
                            enrich.lists.flag=enrich.lists.flag)
        
    }
    
    # Lists <- list(a, b, c)  #put the word vectors into a list to supply lapply  
    # Lists <- lapply(Lists, function(x) as.character(unlist(x)))
    # items <- sort(unique(unlist(Lists)))   #put in alphabetical order
    # MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=3)  #make a matrix of 0's
    # names <- c(label1,label2,label3)                
    # colnames(MAT) <- names
    # rownames(MAT) <- items
    # lapply(seq_along(Lists), function(i) {   #fill the matrix
    #   MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
    # })
    # 
    # outputName=paste(label1,"_",label2,"_",label3,"_VennDiagram.pdf",sep="")
    # 
    # #b=paste(a,outputName,sep="")
    # 
    # file.path.name <- file.path(plot.dir, "Venn3")
    # dir.create(file.path.name, recursive=TRUE, showWarnings=FALSE)
    # file.path.name <- paste0(file.path.name, "/", outputName)
    # dev.new()
    # require("limma")
    # limma::vennDiagram(MAT, circle.col= c("red","green","yellow"), main=title)
    # 
    # dev.print(device=pdf, file=file.path.name, width=10, height=10)
    # graphics.off()
    return(ret.list)
  }


SaveInteserctionsList <- function(gene.list, conversion.map, root.dir, prefix, 
                                labels.list, enrich.lists.flag=FALSE) 
{
    
    for(lbl in labels.list) {
        prefix <- UpdatePrefix(prefix, lbl)
    }
    
    out.dir <- UpdateFolderPath(root.dir, prefix)
    filename <- UpdateFilename(prefix, "genes")
    
    if(!is.null(conversion.map)) {
        gene.list.df <- CreateConvertedDataframe(gene.list, conversion.map)
    } else {
        gene.list.df <- as.data.frame(gene.list)
    }
    
    WriteDataFrameAsTsv(data.frame.to.save=gene.list.df, 
                        file.name.path=file.path(out.dir, filename),
                        col.names= TRUE, row.names=FALSE)
    
    if(enrich.lists.flag) {
        enrichSuitedSingleList(de.gene.list=gene.list, 
                            functional.folder=file.path(out.dir, "functional"), 
                            filename=filename)
    }
    
}


UpdateFilename <- function(filename, ..., extension=NULL) {
    dots <- list(...)
    filename <- gsub(pattern = " ", replacement = "_", x = filename)
    if(length(dots) != 0) {
        for (str in dots) filename <- paste(filename, str, sep = "_")
    } else {
        stop("provide a string to append to ", filename)
    }
    if(!is.null(extension)) {
        filename <- paste0(filename, ".", extension)
    }
    return(filename)
}

UpdateFolderPath <- function(path, ...) {
    dots <- list(...)
    if(length(dots) != 0) {
        for (str in dots) {
            str <- gsub(pattern = " ", replacement = "_", str)
            path <- file.path(path, str)
        }
    } else {
        stop("provide a string to append to ", path)
    }
    dir.create(path, recursive=TRUE, showWarnings=FALSE)
    return(path)
}


CreateConvertedDataframe <- function(gene.list, conversion.map) {
    # de.data.frame.ord <- de.data.frame[ order(rownames(de.data.frame)),]
    gene.list.ord <- as.data.frame(gene.list[order(gene.list)])
    conversion.map.ord <- conversion.map[ order(conversion.map[,1]), ]
    ind.map <- which(conversion.map.ord[,1] %in% gene.list.ord[,1] )
    conversion.map.sub <- conversion.map.ord[ind.map,]
    if(dim(conversion.map.sub)[1]==dim(gene.list.ord)[1]) {
        gene.list.ord$symbol <- as.character(conversion.map.sub[,2])
    } else {
        stop("gene map and de.dataframe dimensions differs!")
    }
    colnames(gene.list.ord) <- colnames(conversion.map)
    
    return(gene.list.ord)
}


