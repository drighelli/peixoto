
 library(plyr)
computeGeneMeansOverGroups <- function(counts, design, groupColumn)
{
    groups <- unique(design[[groupColumn]])
    aa <- ldply(groups, function(x)
    {
        idx <- which(design[[groupColumn]] %in% x)
        apply(counts[,rownames(design)[idx]], 1, mean)
    })
    rownames(aa) <- groups
    aa <- t(aa)
}





