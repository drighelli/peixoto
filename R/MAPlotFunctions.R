#' Generate an MA-Plot using ggplot
#' 
#' @param processed.de.results dataframe. D-E results processed with ProcessDEResultsForPlot
#' @param strings character list. Returned by GeneratePlotStrings
#' @param plotly.flag logical. Is it a ggplotly object?
#' 
#' @return none
#' 
#' @export 
#' 
#' @importFrom ggplot2 geom_point aes scale_color_manual labs
#' @import ggplot2
GenerateGGMA <- function(processed.de.results, strings, plotly.flag=FALSE) {
    
    if(plotly.flag) {
        xlabl <- paste("log<sub>2</sub>(counts)")
        ylabl <- paste("-log<sub>2</sub>(FC)")
    }else {
        xlabl <- bquote(~log[2]~"(counts)")
        ylabl <- bquote(~log[2]~"(FC)")
    }
    
    switch(processed.de.results$method[1],
           DESeq={
               ggp <- ggplot2::ggplot(processed.de.results) + 
                   geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                  color=significance, ensembl=gene, 
                                  symbol=symbol, 
                                  padj=format(padj, nsmall=10)), 
                              size=0.7)    
               + labs( list(title=strings$title, x=xlabl, y=ylabl))
               + scale_color_manual(values=c("blue2", "red2")) 
               
               if(!plotly.flag) {
                   ggp <- ggp + 
                       geom_point(data=subset(processed.de.results, 
                                              significance=="significative"), 
                                  aes(x=log2Counts, y=log2FoldChange, 
                                      color=significance, ensembl=gene, 
                                      symbol=symbol, 
                                      padj=format(padj, nsmall=10)),
                                  size=0.7 ) 
               }
           },
           NOISeqBio={
               ggp <- ggplot2::ggplot(processed.de.results)
               + geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                color=significance, ensembl=gene, 
                                symbol=symbol, prob=prob), size=0.7)
               + labs( list(title=strings$title, x=xlabl, y=ylabl))
               + scale_color_manual(values=c("blue2", "red2")) 
               if(!plotly.flag) {
                   ggp <- ggp 
                   + geom_point(data=subset(processed.de.results, 
                                            significance=="significative"), 
                                aes(x=log2Counts, y=log2FoldChange,
                                    color=significance, 
                                    ensembl=gene, 
                                    symbol=symbol, 
                                    prob=prob), 
                                size=0.7 ) 
               }     
           },
           NOISeq={
               ggp <- ggplot2::ggplot(processed.de.results) +
                   geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                color=significance, ensembl=gene, 
                                symbol=symbol, prob=prob), 
                            size=0.7)
               + labs( list(title=strings$title, x=xlabl, y=ylabl))
               + scale_color_manual(values=c("blue2", "red2")) 
               if(!plotly.flag) {
                   ggp <- ggp +
                       geom_point(data=subset(processed.de.results, 
                                            significance=="significative"), 
                                aes(x=log2Counts, y=log2FoldChange,
                                    color=significance, 
                                    ensembl=gene, symbol=symbol, 
                                    prob=prob), size=0.7 ) 
               }     
           },
           edgeR={
               ggp <- ggplot2::ggplot(processed.de.results) + 
                   geom_point(aes(x=log2Counts, y=log2FoldChange, 
                                  color=significance, name=gene, 
                                  # symbol=symbol, 
                                  padj=format(padj, nsmall=10)), 
                              size=0.7) +
                   labs(list(title=strings$title, x=xlabl, y=ylabl)) +
                   scale_color_manual(values=c("blue2", "red2")) 
               
               # if(!plotly.flag) {
               #     ggp <- ggp + 
               #         geom_point(data=subset(processed.de.results, 
               #                                significance=="significative"), 
               #                    aes(x=log2Counts, y=log2FoldChange, 
               #                        color=significance, ensembl=gene, 
               #                        symbol=symbol, 
               #                        padj=format(padj, nsmall=10)),
               #                    size=0.7 ) 
               # }
           }
           
    )
    ggp <- ggp +
            geom_hline(yintercept=0) +
            geom_hline(yintercept=1, colour="darkgreen", linetype="dashed") +
            geom_hline(yintercept=-1, colour="darkgreen", linetype="dashed")
    
    return(ggp)
}


#' Title
#'
#' @param de.results 
#' @param counts.dataframe 
#' @param design.matrix 
#' @param show.plot.flag 
#' @param plotly.flag 
#' @param save.plot 
#' @param plot.folder 
#' @param prefix.plot 
#' @param threshold 
#'
#' @return
#' @export
#'
#' @examples
PlotMAPlotCounts <- function(de.results, counts.dataframe, design.matrix, 
                            show.plot.flag=TRUE, plotly.flag=FALSE, 
                            save.plot=FALSE, plot.folder=NULL, 
                            prefix.plot="MA_Plot", threshold=0.05) 
{
    require("plotly")
    
    # title <- paste0(prefix.plot, " MA-Plot")
    
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="MAPlot")
    #de.results=de.results; threshold=threshold; counts.dataframe=counts.dataframe; design.matrix=design.matrix
    processed.de.results <- ProcessDEResultsForPlot(de.results, threshold=threshold, counts.dataframe=counts.dataframe, design.matrix=design.matrix)
    # if(plotly.flag) {
    #     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7)    + ggtitle(strings$title) + scale_color_manual(values=c("blue2", "red2")) 
    # } else {
    #     ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7) + scale_color_manual(values=c("blue2", "red2"))    + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2Counts, y=log2FoldChange, color=significance, name=gene, padj=padj), size=0.7 )    + ggtitle(strings$title) 
    # }
    ggp <- GenerateGGMA(processed.de.results=processed.de.results, strings=strings, plotly.flag=plotly.flag)
    
    if(save.plot) {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the MA-Plot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        }
        
    } 
    
    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
    
    ##return(ggp)
}