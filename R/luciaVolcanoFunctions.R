createSignLabelsNumbers <- function(proc.df, column.name)
{
    values <- unique(proc.df[[column.name]])
    tot.de <- sum(proc.df[[column.name]] == values[1])
    tot.not.de <- dim(proc.df)[1] - tot.de 
    
    idxsign <- which(proc.df[[column.name]]==values[1])
    proc.df[[column.name]][idxsign] <- paste0(values[1], 
                                            " [", tot.de, "]")
    proc.df[[column.name]][-idxsign] <- paste0(values[2], 
                                            " [", tot.not.de, "]")
    return(proc.df)
}


luciaVolcanoPlot <- function(res.o, positive.controls.df, prefix, 
                             threshold=0.01, plotly.flag=TRUE)
{
    if(plotly.flag)
    {
        xlabl <- paste0("log<sub>2</sub>(FC)")
        ylabl <- paste("-log<sub>10</sub>(PValue)")
    } else {
        xlabl <- bquote(~log[2]~"(FC)")
        ylabl <- bquote(~-log[10]~"(PValue)")
    }
    
    title <- paste(prefix, "Volcano Plot")
    
    new.de <- ProcessDEResultsForPlot(de.results=res.o, threshold=threshold, 
                                      design.matrix=desMat)
    pos.contr <- positive.controls.df
    with.pos.de <- new.de
    
    without.pos.de <- createSignLabelsNumbers(with.pos.de, "significance")
    without.pos.de$Significance <- without.pos.de$significance
    pp <- ggplot(data=without.pos.de) +
        geom_point(aes(x=log2FoldChange, y=minuslog10pval, color=Significance, 
                        text=paste0("padj=", padj, "\nname=", gene)), 
                        size=0.7) +
        scale_color_manual(values=c("red2", "blue2"))
    
    if(!is.null(positive.controls.df) )
    {
        with.pos.de$hit <- NA
        idxpc <- which(tolower(with.pos.de$gene) %in% tolower(pos.contr[,1]))
        if(length(idxpc) > 0) 
        {
            with.pos.de$hit[idxpc] <- "pctr"

            sub.de <- with.pos.de[which(with.pos.de$hit=="pctr"),] 
            sub.de$col <- "pctr"
            sub.de$lit <- "est"
            pos.lit <- pos.contr[pos.contr[,2]=="lit",, drop=FALSE]
            idx.lit <- which(tolower(sub.de$gene) %in% tolower(pos.lit[,1]))
            sub.de$lit[idx.lit] <- "lit"
            
            sign <- which(sub.de$significance==paste("padj <", threshold))
            sub.de$hit[sign] <- paste("pctr <", threshold)
            notsign <- which(sub.de$significance==paste("padj >=", threshold))
            sub.de$hit[notsign] <- paste("pctr >=", threshold)
            
            sub.de <- createSignLabelsNumbers(sub.de, "hit")
            pp <- pp + 
                geom_point(data=sub.de, aes(x=log2FoldChange, y=minuslog10pval, 
                                shape=hit, color=col,
                                text=paste0("padj=", padj, "\nname=", gene)), 
                                size=0.9) +
                scale_color_manual(values=c("red2", "blue2", "green3"))
        }
    }
   
    
    
    pp <- pp + labs(list(title=title, x=xlabl, y=ylabl))
    return(pp)
}

luciaVolcanoPlot1 <- function(res.o, positive.controls.df, prefix, 
                             threshold=0.01, plotly.flag=TRUE)
{
    if(plotly.flag)
    {
        xlabl <- paste0("log<sub>2</sub>(FC)")
        ylabl <- paste("-log<sub>10</sub>(PValue)")
    } else {
        xlabl <- bquote(~log[2]~"(FC)")
        ylabl <- bquote(~-log[10]~"(PValue)")
    }
    
    title <- paste(prefix, "Volcano Plot")
    
    new.de <- ProcessDEResultsForPlot(de.results=res.o, threshold=threshold, 
                                      design.matrix=desMat)
    pos.contr <- positive.controls.df
    with.pos.de <- new.de
    
    without.pos.de <- createSignLabelsNumbers(with.pos.de, "significance")
    without.pos.de$Significance <- without.pos.de$significance
    pp <- ggplot(data=without.pos.de) +
        geom_point(aes(x=log2FoldChange, y=minuslog10pval, color=Significance, 
                       text=paste0("padj=", padj, "\nname=", gene)), 
                   size=0.7) +
        scale_color_manual(values=c("red2", "blue2")) #+
    #     geom_vline(xintercept=-1, color="green4", linetype="dashed", size=0.5) +
    #     geom_vline(xintercept=1, color="green4", linetype="dashed", size=0.5) + 
         # geom_hline(yintercept=1.715, color="green4", linetype="dashed", size=0.5) 
     # pp
    
    if(!is.null(positive.controls.df) )
    {
        with.pos.de$hit <- NA
        idxpc <- which(tolower(with.pos.de$gene) %in% tolower(pos.contr[,1]))
        if(length(idxpc) > 0) 
        {
            with.pos.de$hit[idxpc] <- "Pos. Controls"
            sub.de <- with.pos.de[which(with.pos.de$hit=="Pos. Controls"),] 
            sub.de$col <- "Pos. Controls"
            
            # sub.de$lit <- "est"
            # pos.lit <- pos.contr[pos.contr[,2]=="lit",, drop=FALSE]
            # idx.lit <- which(tolower(sub.de$gene) %in% tolower(pos.lit[,1]))
            # sub.de$lit[idx.lit] <- "lit"
            
            sign <- which(sub.de$significance==paste("P. Adj <", threshold))
            sub.de$hit[sign] <- paste("Pos. Controls <", threshold)
            notsign <- which(sub.de$significance==paste("P. Adj >=", threshold))
            sub.de$hit[notsign] <- paste("Pos. Controls >=", threshold)
            
            sub.de <- createSignLabelsNumbers(sub.de, "hit")
            
            # print(head(sub.de))
            pp <- pp + 
                geom_point(data=sub.de, aes(x=log2FoldChange, y=minuslog10pval, 
                                            color=hit,
                                            text=paste0("padj=", padj, "\nname=", gene)), 
                           size=0.9) +
                # scale_color_manual(values=c("blue2", "gray28", "red2", "gray28"))
                scale_color_manual(values=c("gray56", "gray23", "red2", "red2"))
        }
    }
    
    
    
    pp <- pp + labs(list(title=title, x=xlabl, y=ylabl))
    return(pp)
}

