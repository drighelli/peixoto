


gene.name.normalized.counts <- counts
design.matrix <- designMatrix



gene.name <- rownames(res.o)[which(res.o$gene == "Olfr287")]

###### plot counts along times
ProcessCountDataFrameForPlotCountsAcrossTimes <- function(gene.name.normalized.counts, design.matrix, gene.name) {
    gene.name.norm.counts.t <- as.data.frame(t(gene.name.normalized.counts))
    
    gene.name.norm.counts.t$GeneName <- gene.name
    gene.name.norm.counts.t$condition <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "condition"])
    gene.name.norm.counts.t$genotype <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "genotype"])
    gene.name.norm.counts.t$Counts <- gene.name.norm.counts.t[,1]
    return(gene.name.norm.counts.t)
}

ggplot(gene.name.norm.counts.t, mapping=aes(x=condition, y=genotype, color=genotype, group=genotype)) + 
    geom_point() + 
    stat_smooth(se=FALSE, method="loess") + 
    scale_y_log10() + 
    ggtitle(paste( "aaa", gene.name, "gene", sep=" "))
