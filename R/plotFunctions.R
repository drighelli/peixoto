
###### PCA FUNCTIONS

PlotPCAPlotlyFunction <- function(counts.data.frame, design.matrix,
                                  shapeColname=NULL, colorColname=NULL,
                                  scale=FALSE,
                                  xPCA="PC1",
                                  yPCA="PC2",
                                  plot.folder=NULL, prefix.plot=NULL,
                                  show.plot.flag=TRUE, plotly.flag=FALSE,
                                  save.plot=FALSE, ellipse.flag=FALSE,
                                  size=2
                                  )
{
    print(plotly.flag)
    if(is.null(shapeColname) || is.null(colorColname))
    {stop("Please provide the column names for the shape and the color")}

    ## check colnames design matrix
    # require("ggfortify")
    # require("plotly")
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, 
                                plot.type="PCA")
    sub.counts.dataframe <- counts.data.frame[ , 
                which(colnames(counts.data.frame) %in% rownames(design.matrix)), 
                drop=FALSE]
    PCA <- prcomp(t(sub.counts.dataframe), scale.=scale, center=TRUE)
    propvar <- summary(PCA)$importance[2,,drop=FALSE]
    xproppca <- paste0(as.character( propvar[,xPCA]*100 ),"%")
    yproppca <- paste0(as.character( propvar[,yPCA]*100 ),"%")
    sub.pca <- as.data.frame(PCA$x[,c(xPCA, yPCA)])
    sub.pca$samples <- rownames(sub.pca)
    sub.pca <- sub.pca[order(sub.pca$samples), ]
    design.matrix.o <- design.matrix[order(rownames(design.matrix)), , 
                                    drop=FALSE]

    if(dim(sub.pca)[1] == dim(design.matrix.o)[1])
    {
        new.sub.pca <- cbind(sub.pca, design.matrix.o)
    }
    else
    {
        sub.design.matrix <- subset(design.matrix.o,
                            rownames(design.matrix.o) %in% rownames(sub.pca))
        new.sub.pca <- cbind(sub.pca, sub.design.matrix)
    }

    if(ellipse.flag=="both")
    {
        noEllipse <- PlotPCAPlotlyFunction(counts.data.frame=counts.data.frame,
                                  design.matrix=design.matrix,
                                  scale=scale,
                                  shapeColname=shapeColname,
                                  colorColname=colorColname,
                                  plot.folder=plot.folder, 
                                  prefix.plot=prefix.plot,
                                  show.plot.flag=show.plot.flag,
                                  xPCA=xPCA,
                                  yPCA=yPCA,
                                  plotly.flag=plotly.flag, save.plot=save.plot,
                                  ellipse.flag=FALSE,
                                  size=size)
        ellipse <- PlotPCAPlotlyFunction(counts.data.frame=counts.data.frame,
                                  design.matrix=design.matrix,
                                  scale=scale,
                                  shapeColname=shapeColname,
                                  colorColname=colorColname,
                                  plot.folder=plot.folder, 
                                  prefix.plot=prefix.plot,
                                  show.plot.flag=show.plot.flag,
                                  xPCA=xPCA,
                                  yPCA=yPCA,
                                  plotly.flag=plotly.flag, 
                                  save.plot=save.plot,
                                  ellipse.flag=TRUE,
                                  size=size)
        return(list("noEllipse"=noEllipse, "ellipse"=ellipse))
    }
    else if(ellipse.flag)
    {
        if(length(which(colnames(new.sub.pca) %in% colorColname)) >0 )
        {
            setcolname <- colorColname
            # aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorColname,
            #                         shape=shapeColname, name="samples")
            # aesStrObjEll <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname)
        }
        else
        {
            setcolname <- shapeColname
            # aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname,
            #                         shape=shapeColname, name="samples")
            # aesStrObjEll <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname)
        }
        
        aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=setcolname,
                                         shape=shapeColname, name="samples")
        aesStrObjEll <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname)
        
        ggp <- ggplot2::ggplot(new.sub.pca) + 
            ggplot2::geom_point(aesStrObj, size=size) +
            ggplot2::stat_ellipse(aesStrObjEll) + 
            ggplot2::ggtitle(strings$title) +
            ggplot2::xlab(paste(xPCA, xproppca)) + 
            ggplot2::ylab(paste(yPCA, yproppca))
        strings$plot.file.name <- paste0(strings$plot.file.name, "_ellipse")
    }
    else if(!ellipse.flag)
    {
        # if( sum(colnames(new.sub.pca) %in% colorColname) > 0 ) {
        #     aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorColname,
        #                             shape=shapeColname, name="samples")
        # } else {
            aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorColname,
                                    shape=shapeColname, name="samples")
        # }
        ggp <- ggplot2::ggplot(new.sub.pca) + 
            ggplot2::geom_point(aesStrObj, size=size) +
            ggplot2::ggtitle(strings$title) + 
            ggplot2::xlab(paste(xPCA, xproppca)) + 
            ggplot2::ylab(paste(yPCA, yproppca))
    }


    if(save.plot) {
        if(is.null(plot.folder)) {
            stop("Please set a folder where to plot the PCA!")
        }
        if(!is.null(strings$plot.file.name)) {
            # SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder,
                       # plot.file.name=strings$plot.file.name,
                       # plotly.flag=plotly.flag) ## to check and improve ## saw something somewhere on the internet
        }
    }

    if(show.plot.flag) {
        if(plotly.flag) {
            plotly::ggplotly(ggp)
        } else {
            plot(ggp)
        }
    } else {
        return(ggp)
    }
    
}



PlotCountsAlongTimes <- function(normalized.counts, design.matrix, 
                                 gene.name, gene.name.column.name="gene.names", 
                                 show.plot.flag=TRUE, plotly.flag=FALSE, 
                                 save.plot=FALSE, plot.folder=NULL, 
                                 prefix.plot="Counts Plot")
{
    
    strings <- GeneratePlotStrings(path=plot.folder, 
                                   prefix=prefix.plot, 
                                   plot.type="CountsTimePlot")
    
    sub.normalized.counts <- normalized.counts[,
                which(colnames(normalized.counts) %in% rownames(design.matrix))]
    sub.normalized.counts[,gene.name.column.name] <- normalized.counts[,gene.name.column.name]
    
    gene.name.norm.counts <- sub.normalized.counts[which(tolower(sub.normalized.counts[,gene.name.column.name]) %in% tolower(gene.name)), ]
    
    gene.name.norm.counts <- gene.name.norm.counts[, which(colnames(gene.name.norm.counts) %in% rownames(design.matrix))]
    
    if(dim(gene.name.norm.counts)[1] == 1) {
        processed.counts.df <- ProcessCountDataFrameForPlotCountsAcrossTimes(gene.name.norm.counts, design.matrix, gene.name)
    } else if(dim(gene.name.norm.counts)[1] == 0) {
        stop(gene.name, " gene Not Found!")
    } else if(dim(gene.name.norm.counts)[1] > 1) {
        stop(gene.name, " founded in more than one row!")
    }
    require("plotly")
    ggp <- ggplot(processed.counts.df, mapping=aes(x=Times, y=Counts, color=Conditions, group=Conditions)) + geom_point() + stat_smooth(se=FALSE, method="loess") + scale_y_log10() + ggtitle(paste( strings$title, gene.name, "gene", sep=" "))
    
    if(save.plot) {
        if(is.null(strings$plot.folder)) {
            stop("Please set a folder where to plot the boxplot!")
        }
        if(!is.null(strings$plot.file.name)){
            SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=paste(strings$plot.file.name, gene.name, sep="_"), plotly.flag=plotly.flag)
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

