

SelectGenesFromHistBreaks <- function(data.counts, breakss = 100, n.genes.per.break=1) {
    data.counts.mean <- apply(data.counts, 1, mean)
    
    data.counts.mean.filtered <- data.counts.mean[which(data.counts.mean <= summary(data.counts.mean)["3rd Qu."])]
    h <- hist(data.counts.mean.filtered, breaks = breakss, plot=FALSE)
    data.counts.mean.sorted <- sort(data.counts.mean.filtered)
    gene.names <- c()
    j<-1
    for(i in 1:length(h$counts)) {
        names <- names(data.counts.mean.sorted[ j : (j+h$counts[i]-1) ])
        if(length(names) >= n.genes.per.break) {
            n.g.names <- names[c(1:n.genes.per.break)]
        } else {
            n.g.names <- names[c(1:length(names))]
        }
        gene.names <- c(gene.names, n.g.names)
        
        j <- j+h$counts[i]
    }
    return(gene.names)
}

EstimateNegativeControlGenesForRUV <- function(de.genes, n.tail.genes=2000, counts.dataset, n.genes.per.hist.break=1, threshold=0.05) {
    de.genes.tail <- TailDeGenesByPAdj(de.genes = de.genes, n.genes = n.tail.genes, threshold = threshold)
    counts.tail <- counts.dataset[which(rownames(counts.dataset) %in% de.genes.tail),]
    estimated.genes <- SelectGenesFromHistBreaks(counts.tail, n.genes.per.break=n.genes.per.hist.break)
    return(estimated.genes)
}


RUVgNormalizationFunction <- function(data.to.normalize, 
                                      design.matrix, 
                                      desMatColStr, 
                                      estimated.gene.names, 
                                      k=1, isLog=FALSE,
                                      gene.names.col=NULL) {
    
    if( length( which(colnames(design.matrix) == desMatColStr)) == 0 ) {
        stop("Please provide a design matrix with a named ", desMatColStr, 
            " column")
    }
    # require("RUVSeq")
    if(is.null(gene.names.col))
    {
        genes <- rownames(design.matrix)
    } else {
        genes <- data.to.normalize[,"gene.names.col"]
    }
    data.to.normalize <- data.to.normalize[,
                                which(colnames(data.to.normalize) %in% colnames(design.matrix))]
    if(!isLog) {
        dataset.set <- EDASeq::newSeqExpressionSet(as.matrix(data.to.normalize),
                                    phenoData=data.frame(
                                    design.matrix[[desMatColStr]],
                                    row.names=
                                    ))
    } else {
        dataset.set <- data.to.normalize
    }
    
    
    ruved.set <- RUVSeq::RUVg(x=dataset.set, cIdx=estimated.gene.names, 
                              k=k, isLog=isLog)
    
    
    #return(as.data.frame(ruved.set@assayData$normalizedCounts))
    return(ruved.set)
}

NormalizeData <- function(data.to.normalize, norm.type = c("fqua", "uqua", "tmm", "ruvg"), design.matrix=NULL, design.matrix.factors.column=NULL, estimated.genes=NULL, is.log=FALSE) {
    ## @ norm.type can be uqua, tmm, fqua or ruvg
    
    x <- data.to.normalize
    
    if(norm.type != "fqua") {
        if(norm.type == "ruvg") {
            if(!(is.null(design.matrix) && is.null(estimated.genes))){ ## test
                normalized.data <- RUVgNormalizationFunction(
                                        data.to.normalize=data.to.normalize, 
                                        design.matrix=design.matrix, 
                                        desMatColStr=design.matrix.factors.column,
                                        estimated.gene.names=estimated.genes,
                                        isLog=is.log)
            } else {
                stop("Please select a design matrix and a list of negative",
                        " control genes for RUVg normalization")
            }
        } else {
            require("edgeR")
            x <- edgeR::DGEList(counts = x)
            if(norm.type=="uqua") {
                x <- edgeR::calcNormFactors(x, method='upperquartile')
            } 
            
            if(norm.type=="tmm") {
                x <- edgeR::calcNormFactors(x, method='TMM')
            }
            
            x <- edgeR::estimateCommonDisp(x, verbose=FALSE)
            x <- edgeR::estimateTagwiseDisp(x)
            normalized.data <- as.data.frame(x$pseudo.counts)
        }
        
    } else if(norm.type == "fqua") {
        require("preprocessCore")
        x <- as.matrix(x)
        normalized.data <- as.data.frame(preprocessCore::normalize.quantiles(x, copy=FALSE))
        colnames(normalized.data) <- colnames(x)
        rownames(normalized.data) <- rownames(x)
    }
    
    return(normalized.data)
}


SortDeGenesByPAdj <- function(de.genes) {
    if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
        ordered.de.genes <- de.genes[order(de.genes$padj),]
    } else if("theta" %in% colnames(de.genes)) { ## noiseqbio
        ordered.de.genes <- de.genes[order(de.genes$prob, decreasing = TRUE),]
    } else if("M" %in% colnames(de.genes)) { ## noiseq
        ordered.de.genes <- de.genes[order(de.genes$prob, decreasing = TRUE),]
    } else {
        stop("In SortDeGenesByPAdj not recognized de results table!")
    }
    
    # ordered.de.genes <- de.genes[order(de.genes$padj),]
    return(ordered.de.genes)
}

FilterOutGenesNAPAdj <- function(de.genes) {
    if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
        if(sum(is.na(de.genes$padj))>0) de.genes.not.na <- de.genes[-which(is.na(de.genes$padj)),]
        else de.genes.not.na <- de.genes
    } else {
        de.genes.not.na <- de.genes
    }
    
    # if(sum(is.na(de.genes$padj))>0) de.genes.not.na <- de.genes[-which(is.na(de.genes$padj)),]
    # else de.genes.not.na <- de.genes
    return(de.genes.not.na)
}



TailDeGenesByPAdj <- function(de.genes, n.genes=10, threshold=0.05) {
    ord.de.genes <- SortDeGenesByPAdj(de.genes = de.genes)
    ord.de.genes.not.na <- FilterOutGenesNAPAdj(ord.de.genes)
    if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
        ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$padj >= threshold,]
    } else if("prob" %in% colnames(de.genes)) { ## noiseqbio
        ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$prob <= threshold,]
    }
    # ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$padj >= threshold,]
    tail.genes <- rownames(tail(ord.not.de.genes, n=n.genes))
    return(tail.genes)
}


