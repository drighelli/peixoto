filteredCountsProp <- filterLowCounts(counts.dataframe=countMatrix, 
                                      is.normalized=FALSE,
                                      design.dataframe=designMatrix,
                                      cond.col.name="gcondition",
                                      method.type="Proportion")


sd.ctrls <- read_excel(path="../data/controls/Additional File 4 full list of BMC genomics SD&RS2.xlsx", sheet=1)
sd.ctrls <- sd.ctrls[order(sd.ctrls$adj.P.Val),]

sd.neg.ctrls <- sd.ctrls[sd.ctrls$adj.P.Val > 0.9, ]

sd.neg.ctrls <- sd.neg.ctrls$`MGI Symbol`
sd.neg.ctrls <- sd.neg.ctrls[-which(is.na(sd.neg.ctrls))]

int.neg.ctrls <- sd.neg.ctrls
neg.map <- convertGenesViaBiomart(specie="mm10", filter="mgi_symbol",
                                  filter.values=int.neg.ctrls, c("external_gene_name",
                                                                 "mgi_symbol", "entrezgene"))

neg.map.nna <- neg.map[-which(is.na(neg.map$entrezgene)),]

neg.ctrls.entrez <- as.character(neg.map.nna$entrezgene)

ind.ctrls <- which(rownames(filteredCountsProp) %in% neg.ctrls.entrez)
counts.neg.ctrls <- filteredCountsProp[ind.ctrls,]

normPropCountsUqua <- NormalizeData(data.to.normalize=filteredCountsProp, 
                                    norm.type="tmm", 
                                    design.matrix=designMatrix)

library(RUVSeq)
neg.ctrl.list <- rownames(counts.neg.ctrls)
groups <- makeGroups(designMatrix$gcondition)
ruvedSExprData <- RUVs(as.matrix(round(normPropCountsUqua)), cIdx=neg.ctrl.list,
                       scIdx=groups, k=5)

normExprData <- ruvedSExprData$normalizedCounts

padj.thr <- 0.05
desMat <- cbind(designMatrix, ruvedSExprData$W)
colnames(desMat) <- c(colnames(designMatrix), colnames(ruvedSExprData$W))

cc <- c("WTSD5 - WTHC5", "KOHC5 - WTHC5",
        "KOSD5 - WTSD5", "KOSD5 - KOHC5")

rescList1 <- applyEdgeR(counts=filteredCountsProp, design.matrix=desMat,
                        factors.column="gcondition", 
                        weight.columns=c("W_1", "W_2", "W_3", "W_4", "W_5"),
                        contrasts=cc, useIntercept=FALSE, p.threshold=1,
                        is.normalized=FALSE, verbose=TRUE)

######## apply edger
counts=filteredCountsProp; design.matrix=desMat;
factors.column="gcondition"; 
weight.columns=c("W_1", "W_2", "W_3", "W_4", "W_5");
contrasts=cc; useIntercept=FALSE; p.threshold=1;
is.normalized=FALSE; verbose=TRUE

factors <- design.matrix[[factors.column]]


    
    weights <- as.matrix(design.matrix[weight.columns])
    if(useIntercept)
    {
        design <- model.matrix(~ 1 + factors + weights)
    } else {
        design <- model.matrix(~ 0 + factors + weights)
    }
    colnames(design) <- c(as.character(unique(factors)), weight.columns)


# fit <- applyEdgeRQLFit(counts=counts, factors=factors, design=design,
#                        is.normalized=is.normalized, verbose=verbose)

### edger fit

dgel <- edgeR::DGEList(counts=counts, group=factors)
if(!is.normalized)
{
    dgel <- edgeR::calcNormFactors(dgel, method="TMM")
}
edisp <- edgeR::estimateDisp(y=dgel, design=design)
fit <- edgeR::glmQLFit(edisp, design, robust=TRUE)

# 
# resClist <- lapply(contrasts, function(c)
# {
    c=contrasts[[1]]
    # resC <- applyEdgeRContrast(contrast=c, design=design,
    #                            fit=fit, p.threshold=p.threshold,
    #                            verbose=verbose)
    # 
    contr <- limma::makeContrasts(contrasts=c, levels=design)
    qlf <- edgeR::glmQLFTest(fit, contrast=contr)
    res <- edgeR::topTags(qlf, n=Inf, p.value=p.threshold)
    
    resC <- res$table
    
    cg <- gsub(pattern=" ", replacement="", x=c)
    cs <- strsplit(cg, split="-")[[1]]
    genes <- rownames(resC)
    
      counts <- qlf$fitted.values ## not working!
    
    # counts <- normExprData ## THIS WORKS on the head!
    ## counts <- normPropCountsUqua
    
        # ctMeans <- computeMeans(counts=counts, design.matrix=design.matrix,
        #                         factors.column=factors.column, contrst=cs,
        #                         genes=genes)
        counts=counts; design.matrix=design.matrix;
                                 factors.column=factors.column; contrst=cs;
                                 genes=genes
        
        design.factors <- design.matrix[, factors.column, drop=FALSE]
        counts <- counts[match(genes, rownames(counts)),]
        ctMeans <- do.call( cbind,
                               lapply(contrst, function(cc)
                               {
                                   idx <- which(design.factors[[factors.column]] 
                                                %in% cc)
                                   if(length(idx) == 0) stop("contrast ", cc, 
                                                             " not found!")
                                   samples <- rownames(design.factors)[idx]
                                   subcounts <- counts[,samples]
                                   print(head(subcounts))
                                   submeans <- apply(X=subcounts, MARGIN=1, mean)
                                   print(head(submeans))
                                   return(submeans)
                               })
        )
        colnames(contrMeans) <- contrst
        
        resCMeans <- cbind(ctMeans, resC)
        
        head(ctMeans)
        ctMeansLogfc <- log2((ctMeans[,1]/ctMeans[,2]) )
        ctMeansLogfcdf <- cbind(ctMeans, ctMeansLogfc)
        head(resC)
        head(ctMeansLogfcdf)
        tail(resC)
        tail(ctMeansLogfcdf)
        
    if(p.threshold != 1)
    {
        resCMeans <- resCMeans[(resCMeans$FDR < p.threshold),]
    } else {
        resCMeans <- resCMeans[(resCMeans$FDR <= p.threshold),]
    }
    
    
# })
# names(resClist) <- contrasts








