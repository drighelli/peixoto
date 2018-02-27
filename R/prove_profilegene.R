


gene.name.normalized.counts <- normPropCountsUqua
design.matrix <- designMatrix



gene.name <- rownames(res.o)[which(res.o$gene == "Olfr287")]

###### plot counts along times
ProcessCountDataFrameForPlotCountsAcrossTimes <- function(gene.name.normalized.counts, design.matrix, gene.name) {
    gene.name.normalized.counts <- gene.name.normalized.counts[which(rownames(gene.name.normalized.counts) %in% gene.name),]
    gene.name.norm.counts.t <- as.data.frame(t(gene.name.normalized.counts))
    gene.name.norm.counts.t$GeneName <- gene.name
    gene.name.norm.counts.t$condition <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "condition"])
    gene.name.norm.counts.t$genotype <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "genotype"])
    gene.name.norm.counts.t$Counts <- log(gene.name.norm.counts.t[,1])
    return(gene.name.norm.counts.t)
}

wtcounts <- gene.name.norm.counts.t[gene.name.norm.counts.t$genotype == "WT",]

# gene.name.norm.counts.t$condition[gene.name.norm.counts.t$condition=="HC5"] <- 0
# gene.name.norm.counts.t$condition[gene.name.norm.counts.t$condition=="SD5"] <- 0.01
# gene.name.norm.counts.t$condition <- as.numeric(gene.name.norm.counts.t$condition)

pp <- ggplot(gene.name.norm.counts.t, mapping=aes(x=condition, y=Counts, color=genotype)) + 
    geom_point() + 
    stat_smooth(aes(x=as.numeric(as.factor(condition)), y=Counts, color=genotype),  method = "lm", se=FALSE, fullrange=FALSE) +
    facet_grid(.~genotype) +
    ggtitle(paste( "Profile of", gene.name, "gene", sep=" "))

pp
ggplotly(pp) 
