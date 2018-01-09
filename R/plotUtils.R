
#########################
ProcessDEResultsForPlot <- function(de.results, threshold, 
                                    counts.dataframe=NULL, design.matrix=NULL) 
{
    de.results.new <- de.results
    de.results.new <- de.results.new[order(rownames(de.results.new)),]
    counts.dataframe.ord <- counts.dataframe[order(rownames(counts.dataframe)),]

    if("baseMean" %in% colnames(de.results.new)) { ## working on DESeq2 results

        de.results.new <- de.results.new[, c("padj", "log2FoldChange", "symbol"), drop=FALSE ]
        conds <- as.character(unique(design.matrix$Conditions))
        sub.des.A <- subset(design.matrix, design.matrix$Conditions %in% conds[1])
        sub.counts.A <- counts.dataframe.ord[, which(colnames(counts.dataframe.ord) %in% rownames(sub.des.A) ), drop=FALSE ]
        sub.des.B <- subset(design.matrix, design.matrix$Conditions %in% conds[2])
        sub.counts.B <- counts.dataframe.ord[, which(colnames(counts.dataframe.ord) %in% rownames(sub.des.B) ), drop=FALSE ]
        de.results.new$meanA <- apply(sub.counts.A, 1, mean)
        de.results.new$meanB <- apply(sub.counts.B, 1, mean)
        #########
        de.results.new$log2Counts <- (1/2) * log2((de.results.new$meanA * de.results.new$meanB))
        de.results.new$log10FoldChange <- log10( (de.results.new$meanA / de.results.new$meanB ) )
        de.results.new$minuslog10PAdj <- (-1) * log10(de.results.new$padj)
        de.results.new$significance <- "not-significative"
        de.results.new$significance[which(de.results.new$padj < threshold)] <- "significative"
        de.results.new <- de.results.new[order(de.results.new$padj, decreasing=FALSE),]
        de.results.new$method <- rep(x="DESeq2", times=dim(de.results.new)[1])

    } else if("theta" %in% colnames(de.results.new)) { ## working on NOISeqbio results

        de.results.new <- de.results.new[, c(1, 2, 6)]
        de.results.new$prob <- de.results$prob
        de.results.new$log2FoldChange <- de.results$log2FC
        de.results.new$log10FoldChange <- log10( (de.results[,1]/de.results[,2]) )
        de.results.new$log2Counts <- (1/2) * log2( (de.results[,1] * de.results[,2]) )
        de.results.new$significance <- "not-significative"
        de.results.new$significance[which(de.results.new$prob > threshold)] <- "significative"
        de.results.new <- de.results.new[order(de.results.new$prob, decreasing=TRUE),]
        de.results.new$minuslog101minuspp <- (-1) * log10( (1 - de.results.new$prob + 0.000001))
        de.results.new$method <- rep(x="NOISeqBio", times=dim(de.results.new)[1])

    } else if("M" %in% colnames(de.results.new)) { ## working on NOISeq results

        de.results.new <- de.results.new[, c(1, 2, 7)]
        de.results.new$prob <- de.results$prob
        de.results.new$log2FoldChange <- de.results$M
        de.results.new$log10FoldChange <- log10( (de.results[,1]/de.results[,2]) )
        de.results.new$log2Counts <- (1/2) * log2( (de.results[,1] * de.results[,2]) )
        de.results.new$significance <- "not-significative"
        de.results.new$significance[which(de.results.new$prob > threshold)] <- "significative"
        de.results.new <- de.results.new[order(de.results.new$prob, decreasing=TRUE),]
        de.results.new$minuslog101minuspp <- (-1) * log10( (1 - de.results.new$prob + 0.000001))
        de.results.new$method <- rep(x="NOISeq", times=dim(de.results.new)[1])

    }else if("F" %in% colnames(de.results.new)) { ## working on edgeR results

        de.results.new <- de.results.new[, c(1:3, 6:7)]
        de.results.new$padj <- de.results.new$FDR
        de.results.new$log2FoldChange <- de.results.new$logFC
        de.results.new$log10FoldChange <- log10( (de.results.new[,1]/de.results.new[,2]))
        de.results.new$minuslog10pval <- -log10(de.results.new$PValue)
        de.results.new$log2Counts <- (1/2) * log2((de.results.new[,1] * de.results.new[,2]))
        de.results.new$significance <- "not-sign"
        idx <- which(de.results.new$FDR < threshold)
        de.results.new$significance[idx] <- "sign"
        de.results.new <- de.results.new[order(de.results.new$padj, decreasing=FALSE),]
        de.results.new$minuslog10PAdj <- (-1) * log10(de.results.new$FDR)
        de.results.new$method <- rep(x="edgeR", times=dim(de.results.new)[1])
        
    }

    de.results.new$gene <- rownames(de.results.new)

    return(de.results.new)
}
