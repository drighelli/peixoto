
filterLowCounts <- function(counts.dataframe, is.normalized = c(TRUE, FALSE),
                            design.dataframe, cond.col.name=NULL,  
                            method.type = c("CPM", "Wilcoxon", "Proportion"), 
                            cv.percentage=100, cpm.cutoff=1, seq.depth=NULL,
                            verbose=TRUE) 
{
    
    method.type <- match.arg(method.type)
    stopifnot(is.logical(is.normalized))
    stopifnot(!is.null(cond.col.name))
    stopifnot((cond.col.name %in% colnames(design.dataframe)))
   
    
    design.dataframe <- design.dataframe[order(rownames(design.dataframe)), , 
                                         drop=FALSE]
    counts.dataframe <- counts.dataframe[, order(colnames(counts.dataframe)), 
                                         drop=FALSE]
    
    idx.cols <- which(colnames(counts.dataframe) %in% rownames(design.dataframe))
    sub.counts.dataframe <- counts.dataframe[, idx.cols]
    
    conditions <- design.dataframe[, cond.col.name]
    
    if(method.type == "Proportion") 
    {
        if(is.normalized) {
            ## calcolare seq.depth
            if(is.null(seq.depth)) 
                stop("Proportion test cannot be performed on normalized counts",
                     " without sequencing depth!\nYou need column totals ",
                     " before normalizing the data.")
        } 
        else 
        {
            seq.depth <- NULL
        }
    }
    
    switch( method.type,
            CPM = {
                method.number <- 1
            },
            Wilcoxon = {
                method.number <- 2
            },
            Proportion = {
                method.number <- 3
            }
    )
    
    if(verbose) message("features dimensions before normalization: ",
                        dim(sub.counts.dataframe))
    
    filtered.dataframe <- NOISeq::filtered.data(dataset = sub.counts.dataframe, 
                                                factor = conditions, 
                                                norm = is.normalized, 
                                                depth = seq.depth, 
                                                method = method.number, 
                                                cv.cutoff = cv.percentage, 
                                                cpm = cpm.cutoff, 
                                                p.adj = "BH")
    return(filtered.dataframe)
}
