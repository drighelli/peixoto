source("R/plotUtils.R")
############ VOLCANO FUNCTIONS


GenerateGGVolcano <- function(processed.de.results, strings, plotly.flag) {
    require("ggplot2")
    switch(processed.de.results$method[1],
        DESeq2={
            if(plotly.flag) {
                xlabl <- paste0("log<sub>2</sub>(FC)")
                ylabl <- paste("-log<sub>10</sub>(padj)")
            }else {
                xlabl <- bquote(~log[2]~"(FC)")
                ylabl <- bquote(~-log[10]~"(padj)")
            }

            ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2FoldChange, y=minuslog10PAdj, color=significance, ensembl=gene, symbol=symbol, padj=format(padj, nsmall=10) ), size=0.7) + labs(list(title=strings$title, x=xlabl, y=ylabl)) + scale_color_manual(values=c("blue2", "red2"))
            if(!plotly.flag) {
                ggp <- ggp + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2FoldChange, y=minuslog10PAdj, color=significance, ensembl=gene, symbol=symbol, padj=format(padj, nsmall=10)), size=0.7 )
            }
        },
        NOISeqBio={
            if(plotly.flag) {
                xlabl <- paste0("log<sub>2</sub>(FC)")
                ylabl <- paste0("-log<sub>10</sub>(1-prob)")
            }else {
                xlabl <- bquote(~log[2]~"(FC)")
                ylabl <- bquote(~-log[10]~"(1-Prob)")
            }

            ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7) + labs(list(title=strings$title, x=xlabl, y=ylabl)) + scale_color_manual(values=c("blue2", "red2"))
            if(!plotly.flag) {
                ggp <- ggp + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7 )
            }
        },
        NOISeq=
        {
            if(plotly.flag) {
                xlabl <- paste0("log<sub>2</sub>(FC)")
                ylabl <- paste0("-log<sub>10</sub>(1-prob)")
            }else {
                xlabl <- bquote(~log[2]~"(FC)")
                ylabl <- bquote(~-log[10]~"(1-Prob)")
            }

            ggp <- ggplot2::ggplot(processed.de.results) + geom_point(aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7) + labs(list(title=strings$title, x=xlabl, y=ylabl)) + scale_color_manual(values=c("blue2", "red2"))
            if(!plotly.flag) {
                ggp <- ggp + geom_point(data= subset(processed.de.results, significance=="significative"), aes(x=log2FoldChange, y=minuslog101minuspp, color=significance, ensembl=gene, symbol=symbol, prob=format(prob, nsmall=10)), size=0.7 )
            }
        },
        edgeR={
            if(plotly.flag) {
                xlabl <- paste0("log<sub>2</sub>(FC)")
                ylabl <- paste("-log<sub>10</sub>(PValue)")
            } else {
                xlabl <- bquote(~log[2]~"(FC)")
                ylabl <- bquote(~-log[10]~"(PValue)")
            }

            ggp <- ggplot2::ggplot(processed.de.results) +
                geom_point(
                    aes(x=log2FoldChange, y=minuslog10pval, color=significance,
                         padj=format(padj, nsmall=10),
                        name=gene), size=0.7) + 
                 labs(list(title=strings$title, x=xlabl, y=ylabl)) + 
                scale_color_manual(values=c("blue2", "red2"))#, "orange2", "orange2"))
            
            # if(!plotly.flag) {
            #     ggp <- ggp + geom_point(data=subset(processed.de.results, 
            #                                     significance=="significative"),
            #                             aes(x=log2FoldChange, y=minuslog10pval,
            #                                 color=significance,
            #                                 padj=format(padj, nsmall=10)),
            #                             size=0.7 )
            # }
        }
    )

    ggp <- ggp + 
        geom_vline(xintercept=0) + 
        geom_vline(xintercept=1, colour="darkgreen", linetype="dashed") + 
        geom_vline(xintercept=-1, colour="darkgreen", linetype="dashed")

    return(ggp)
}

PlotVolcanoPlot <- function(de.results, 
                            counts.dataframe=NULL, 
                            design.matrix=NULL, 
                            show.plot.flag=TRUE, 
                            plotly.flag=FALSE, 
                            save.plot=FALSE, 
                            plot.folder=NULL, 
                            prefix.plot=NULL, 
                            threshold) 
{
    require("plotly")
    # title <- paste0(prefix.plot, " Volcano Plot")
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="VolcanoPlot")

    processed.de.results <- ProcessDEResultsForPlot(de.results=de.results, threshold=threshold, counts.dataframe=counts.dataframe, design.matrix=design.matrix)

    ggp <- GenerateGGVolcano(processed.de.results, strings, plotly.flag)


    if(save.plot) {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the Volcano-Plot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
        } else {
            stop("Please set a name for the Volcano-Plot!")
        }


    }

    if(show.plot.flag) {
        if(plotly.flag) {
            ggplotly(ggp)
        } else {
            plot(ggp)
        }
    }
}
