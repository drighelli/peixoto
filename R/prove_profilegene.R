


gene.name.normalized.counts <- normExprData
design.matrix <- designMatrix


ProcessCountDataFrameForPlotCounts <- function(gene.name.normalized.counts, design.matrix, gene.name) {
    gene.name.normalized.counts <- gene.name.normalized.counts[which(rownames(gene.name.normalized.counts) %in% gene.name), , drop=FALSE]
    gene.name.norm.counts.t <- as.data.frame(t(gene.name.normalized.counts))
    gene.name.norm.counts.t$genename <- gene.name
    gene.name.norm.counts.t$condition <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "condition"])
    gene.name.norm.counts.t$genotype <- as.character(design.matrix[which(rownames(design.matrix) %in% rownames(gene.name.norm.counts.t)), "genotype"])
    gene.name.norm.counts.t$log2counts <- log2(gene.name.norm.counts.t[,1])
    gene.name.norm.counts.t$counts <- gene.name.norm.counts.t[,1]
    return(gene.name.norm.counts.t)
}

# "Olfr287"
geneProfileLucia <- function(normalized.counts, design.matrix, 
                            gene.name, res.o=NULL, show.plot=FALSE, 
                            plotly.flag=FALSE) 
{
    idx <- which(res.o$gene==gene.name)
    if(length(idx) > 0 )
    {
        gene.name.r <- rownames(res.o)[idx]
    } else {
        ## take the gene directly from the counts rownames 
        ## res.o not useful in this case
        gene.name.r <- gene.name
    }
    gn.counts <- ProcessCountDataFrameForPlotCounts(
                    gene.name.normalized.counts=normalized.counts,
                    design.matrix=design.matrix, gene.name=gene.name.r)
    
    pp <- ggplot(gn.counts, 
                aes_string(x="condition", y="counts", color="genotype")) + 
                geom_point() + 
        stat_smooth(data=gn.counts, 
                    mapping=aes(x=as.numeric(as.factor(condition)), 
                        y=counts, 
                        color=genotype), 
                    method = "lm", se=FALSE, fullrange=FALSE) +
        facet_grid(.~genotype) +
        ggtitle(paste( "Profile of", gene.name.r, "gene", sep=" "))
    
    if(show.plot) 
    {
        if(plotly.flag)
        {
            ggplotly(pp)
        } else {
            pp
        }
    } else {
        return(pp)
    }
}




###### plot counts along times



# gene.name.norm.counts.t$condition[gene.name.norm.counts.t$condition=="HC5"] <- 0
# gene.name.norm.counts.t$condition[gene.name.norm.counts.t$condition=="SD5"] <- 0.01
# gene.name.norm.counts.t$condition <- as.numeric(gene.name.norm.counts.t$condition)


pp
ggplotly(pp) 
