
#' Title
#'
#' @param counts
#' @param design.matrix
#' @param factors.column
#' @param weight.columns
#' @param contrasts
#' @param useIntercept
#' @param p.threshold use 1 if you want all the results to be returned
#' @param is.normalized 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
applyEdgeR <- function(counts, design.matrix, factors.column=NULL,
                       weight.columns=NULL, contrasts=NULL,
                       useIntercept=FALSE, p.threshold=1,
                       is.normalized=FALSE, verbose=FALSE)
{
    stopifnot(is.character(factors.column))
    stopifnot(!is.null(contrasts))
    stopifnot(is.character(contrasts))
    stopifnot(sum(colnames(design.matrix) %in% factors.column) > 0)

    factors <- design.matrix[[factors.column]]

    if(verbose) message("Setting intercept to: ", useIntercept)
    
    if(!is.null(weight.columns))
    {
        stopifnot(is.character(weight.columns))
        stopifnot(sum(colnames(design.matrix) %in% weight.columns)>0)
        weights <- as.matrix(design.matrix[weight.columns])
        if(useIntercept)
        {
            design <- model.matrix(~ 1 + factors + weights)
        } else {
            design <- model.matrix(~ 0 + factors + weights)
        }
        colnames(design) <- c(as.character(unique(factors)), weight.columns)
    } else {
        if(useIntercept)
        {
            design <- model.matrix(~ 1 + factors)
        } else {
            design <- model.matrix(~ 0 + factors)
        }
        colnames(design) <- c(as.character(unique(factors)))
    }
    
    fit <- applyEdgeRQLFit(counts=counts, factors=factors, design=design,
                        is.normalized=is.normalized, verbose=verbose)

    resClist <- lapply(contrasts, function(c)
    {
        
        resC <- applyEdgeRContrast(contrast=c, design=design,
                                    fit=fit, p.threshold=p.threshold,
                                    verbose=verbose)

        cg <- gsub(pattern=" ", replacement="", x=c)
        cs <- strsplit(cg, split="-")[[1]]
        genes <- rownames(resC)

        ctMeans <- computeMeans(counts=counts, design.matrix=design.matrix,
                                factors.column=factors.column, contrst=cs,
                                genes=genes)
        resCMeans <- cbind(ctMeans, resC)
        if(p.threshold != 1)
        {
            resCMeans <- resCMeans[(resCMeans$FDR < p.threshold),]
        } else {
            resCMeans <- resCMeans[(resCMeans$FDR <= p.threshold),]
        }
            

    })
    names(resClist) <- contrasts

    return(resClist)
}

#' Title
#'
#' @param counts 
#' @param design.matrix 
#' @param factors.column 
#' @param contrst 
#' @param genes 
#'
#' @return
#' @export
#'
#' @examples
computeMeans <- function(counts, design.matrix, factors.column, contrst, genes)
{
    design.factors <- design.matrix[, factors.column, drop=FALSE]
    counts <- counts[match(genes, rownames(counts)),]
    contrMeans <- do.call( cbind,
                           lapply(contrst, function(cc)
                           {
                               idx <- which(design.factors[[factors.column]] 
                                                                        %in% cc)
                               if(length(idx) == 0) stop("contrast ", cc, 
                                                                " not found!")
                               samples <- rownames(design.factors)[idx]
                               subcounts <- counts[,samples]
                               submeans <- apply(X=subcounts, MARGIN=1, mean)
                               return(submeans)
                           })
    )
    colnames(contrMeans) <- contrst
    return(contrMeans)
}

#' Title
#'
#' @param counts 
#' @param factors 
#' @param design 
#' @param verbose 
#' @param is.normalized 
#'
#' @return
#' @export
#'
#' @examples
applyEdgeRQLFit <- function(counts, factors, design, 
                    is.normalized=FALSE, method="TMM", verbose=FALSE)
{    
    if(verbose) message("Fitting edgeR QL model")
    dgel <- edgeR::DGEList(counts=counts, group=factors)
    if(!is.normalized) dgel <- edgeR::calcNormFactors(dgel, method=method)
    edisp <- edgeR::estimateDisp(y=dgel, design=design)
    fit <- edgeR::glmQLFit(edisp, design, robust=TRUE)
    return(fit)
}


#' Title
#'
#' @param counts 
#' @param factors 
#' @param design 
#' @param verbose 
#' @param is.normalized 
#'
#' @return
#' @export
#'
#' @examples
applyEdgeRGLMFit <- function(counts, factors, design, 
                            is.normalized=FALSE, method="TMM", verbose=FALSE)
{    
    if(verbose) message("Fitting edgeR GLM model")
    dgel <- edgeR::DGEList(counts=counts, group=factors)
    if(!is.normalized) dgel <- edgeR::calcNormFactors(dgel, method=method)
    edisp <- edgeR::estimateDisp(y=dgel, design=design)
    fit <- edgeR::glmFit(edisp, design, robust=TRUE)
    return(fit)
}


#' Title
#'
#' @param contrast
#' @param design
#' @param fit
#' @param p.threshold using 1 to be returned all the genes results
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
applyEdgeRContrast <- function(contrast, design, fit, p.threshold=1,
                                verbose=FALSE)
{
    if(verbose) message("contrasting ", contrast)
    contr <- limma::makeContrasts(contrasts=contrast, levels=design)
    qlf <- edgeR::glmQLFTest(fit, contrast=contr)
    res <- edgeR::topTags(qlf, n=Inf, p.value=p.threshold)
    all.gen.res <- res$table
    return(all.gen.res)
}

applyEdgeRLRT <- function(fit, interaction.matrix, interaction.term=NULL,  
                          p.threshold=1, verbose=FALSE)
{
    if(is.null(interaction.term)) 
    {
        warning("interaction.term is missing, using last column of",
                " interaction design!")
        interaction.term <- colnames(interaction.matrix)[[
            ncol(interaction.matrix)]]
    }
    
    if(verbose) message("testing ", interaction.term)
    qlf <- edgeR::glmLRT(fit, coef=interaction.term)
    res <- edgeR::topTags(qlf, n=Inf, p.value=p.threshold)
    all.gen.res <- res$table
    return(all.gen.res)
}


checkColName <- function(colname, colnames) 
{
    if(missing(colname)) stop("Please provide a valid column!")
    if(is.character(colname))
    {
        if(sum(colnames %in% colname)==0)
        {
            stop("Please provide a valid colname!")
        }
    } else {
        if(colname > length(colnames))
        {
            stop("Please provide a valid colname index!")
        }
    }
    return(TRUE)
}


#' constructInteractionMatrix
#'
#' @param design.matrix a design matrix of the experiment where each row 
#' describes a sample
#' @param genotype.col the column name/index describing the genotype in 
#' the design matrix
#' @param genotype.ref an optional character indicating which factor name to 
#' use as genotype reference 
#' @param condition.col the column name/index describing the condition in 
#' the design matrix
#' @param weight.col optional column name/index describing the weigths in 
#' the design matrix (see also weigths parameter)
#' @param weights optional matrix with weights to use in the 
#' interaction matrix (see also weigth.col parameter)
#' @param useIntercept 
#'
#' @return a design matrix with interaction terms
#' @export
#'
#' @examples
constructInteractionMatrix <- function(design.matrix, 
                                       genotype.col, genotype.ref=NULL,
                                       condition.col, weight.col, weights=NULL,
                                       useIntercept=TRUE) 
{
    if(missing(design.matrix)) stop("Please provide a valid Design Matrix!")
    if(missing(genotype.col)) stop("Please provide a valid genotype column!")
    checkColName(genotype.col, colnames(design.matrix))
    checkColName(condition.col, colnames(design.matrix))
    
    if(!is.null(genotype.ref)) 
    {
        geno <- relevel(design.matrix[[genotype.col]], ref=genotype.ref)
    } else {
        geno <- design.matrix[[genotype.col]]
    }
    cond <- design.matrix[[condition.col]]
    
    weigths.flag <- TRUE
    if(is.null(weights))
    {
        if(!missing(weight.col)) 
        {
            checkColName(weight.col, colnames(design.matrix))
            weights <- design.matrix[[weight.col]]
        } else {
            weigths.flag <- FALSE
        }
    } 
    
    if(weigths.flag)
    {
        if(useIntercept) 
        {
            interactionFormula <- "~geno*cond+weights"
        } else {
            interactionFormula <- "~0+geno*cond+weights"
        }
    } else {
        if(useIntercept) 
        {
            interactionFormula <- "~geno*cond"
        } else {
            interactionFormula <- "~0+geno*cond"
        }
    }
    
    interactionMatrix <- model.matrix(as.formula(interactionFormula))
        
    
    
    return(interactionMatrix)
}
