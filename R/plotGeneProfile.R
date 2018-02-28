


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
        idx <- which(rownames(res.o)==gene.name)
        if(length(idx) > 0 )
        {
            gene.name.r <- gene.name
        } else {
            stop("gene ", gene.name," not present!")
        }
        
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
        ggtitle(paste( "Profile of", gene.name, "gene", sep=" "))
    
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
