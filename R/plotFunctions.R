
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
    strings <- GeneratePlotStrings(path=plot.folder, prefix=prefix.plot, plot.type="PCA")
    sub.counts.dataframe <- counts.data.frame[ , which(colnames(counts.data.frame) %in% rownames(design.matrix)), drop=FALSE]
    PCA <- prcomp(t(sub.counts.dataframe), scale.=scale, center=TRUE)

    sub.pca <- as.data.frame(PCA$x[,c(xPCA, yPCA)])
    sub.pca$samples <- rownames(sub.pca)
    sub.pca <- sub.pca[order(sub.pca$samples), ]
    design.matrix.o <- design.matrix[order(rownames(design.matrix)), , drop=FALSE]

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
                                  plot.folder=plot.folder, prefix.plot=prefix.plot,
                                  show.plot.flag=show.plot.flag,
                                  plotly.flag=plotly.flag, save.plot=save.plot,
                                  ellipse.flag=FALSE,
                                  size=size)
        ellipse <- PlotPCAPlotlyFunction(counts.data.frame=counts.data.frame,
                                  design.matrix=design.matrix,
                                  scale=scale,
                                  shapeColname=shapeColname,
                                  colorColname=colorColname,
                                  plot.folder=plot.folder, prefix.plot=prefix.plot,
                                  show.plot.flag=show.plot.flag,
                                  plotly.flag=plotly.flag, save.plot=save.plot,
                                  ellipse.flag=TRUE,
                                  size=size)
        return(list("noEllipse"=noEllipse, "ellipse"=ellipse))
    }
    else if(ellipse.flag)
    {
        if(length(which(colnames(new.sub.pca) %in% colorColname)) >0 )
        {
            aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorColname,
                                    shape=shapeColname, name="samples")
            aesStrObjEll <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname)
        }
        else
        {
            aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname,
                                    shape=shapeColname, name="samples")
            aesStrObjEll <- ggplot2::aes_string(x=xPCA, y=yPCA, color=shapeColname)
        }
        ggp <- ggplot2::ggplot(new.sub.pca) + ggplot2::geom_point(aesStrObj, size=size) +
            ggplot2::stat_ellipse(aesStrObjEll) + ggplot2::ggtitle(strings$title) +
            ggplot2::xlab(xPCA) + ggplot2::ylab(yPCA)
        strings$plot.file.name <- paste0(strings$plot.file.name, "_ellipse")
    }
    else if(!ellipse.flag)
    {
        if( sum(colnames(new.sub.pca) %in% colorColname) > 0 ) {
            aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorColname,
                                    shape=shapeColname, name="samples")
        } else {
            aesStrObj <- ggplot2::aes_string(x=xPCA, y=yPCA, color=colorColname,
                                    shape=shapeColname, name="samples")
        }
        ggp <- ggplot2::ggplot(new.sub.pca) + ggplot2::geom_point(aesStrObj, size=size) +
            ggplot2::ggtitle(strings$title) + ggplot2::xlab(xPCA) + ggplot2::ylab(yPCA)
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
    }
    return(ggp)
}


GeneratePlotStrings <- function(path=NULL, prefix, plot.type) {
    title <- gsub(pattern = "_", replacement = " ", x = UpdatePrefix(prefix, plot.type))

    plot.folder <- gsub(pattern = " ", replacement = "_", x = file.path(path, plot.type))

    plot.file.name <- gsub(pattern = " ", replacement = "_", x = UpdatePrefix(prefix, plot.type))
    if(!is.null(path)) dir.create(plot.folder, showWarnings = FALSE, recursive = TRUE)

    return(list("title"= title, "plot.folder"=plot.folder, "plot.file.name"=plot.file.name))
}


UpdatePrefix <- function(prefix, ...) {
    # new.prefix <- paste(prefix, postix, sep=sep)
    dots <- list(...)
    if( length(dots) != 0 ) {
        for (str in dots) {
            # str <- gsub(pattern = ".", replacement = "_", str)
            prefix <- paste(prefix, str, sep = " " )
        }

    } else {
        stop("provide a string to append to ", new.prefix)
    }
    return(prefix)
}
