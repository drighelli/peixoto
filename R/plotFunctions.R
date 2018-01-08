
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

#'
#'
#' ############## BOXPLOT FUNCTIONS
#' PrepareDataFrameForGGBoxplot <- function(data.frame.to.plot, design.matrix, to.log=TRUE, groupColName) {
#'     ## control if design matrix contains times and conditions
#'
#'     if(to.log){
#'         new.df <- as.data.frame(stack(log(data.frame.to.plot+1)))
#'     } else{
#'         new.df <- as.data.frame(stack(data.frame.to.plot))
#'     }
#'     new.df <- new.df[,c(1:2,4)]
#'     group.p.s <- c(rep(NA, dim(new.df)[1]))
#'     conditions.p.s <- c(rep(NA, dim(new.df)[1]))
#'     new.df <- cbind(new.df, group.p.s, conditions.p.s)
#'
#'     colnames(new.df) <- c("names", "samples", "values", groupColName, "conditions")
#'     head(new.df)
#'
#'     groups <- unique(design.matrix[[groupColName]])
#'     for(group in groups) {
#'         samples.t <- rownames(design.matrix)[which(design.matrix[[groupColName]] %in% group)]
#'         new.df[[groupColName]][which(new.df$samples %in% samples.t)] <- as.character(group)
#'     }
#'     new.df[[groupColName]] <- as.factor(new.df[[groupColName]])
#'     rm(groups)
#'
#'     conditions <- unique(design.matrix$Conditions)
#'     for(condition in conditions) {
#'         samples.t <- rownames(design.matrix)[which(design.matrix$Conditions %in% condition)]
#'         new.df$conditions[which(new.df$samples %in% samples.t)] <- as.character(condition)
#'     }
#'     # new.df$conditions <- as.integer(new.df$conditions)
#'     rm(conditions)
#'     return(new.df)
#' }
#'
#'
#'
#' PlotGroupsBoxplot <- function(data.frame.to.plot, design.matrix, output.path=NULL, prefix.plot=NULL, show.plot.flag=TRUE, save.plot=FALSE, plotly.flag=FALSE, groupColName=NULL, toLog=FALSE) {
#'     ## control if design matrix contains times and conditions
#'     require("plotly")
#'     strings <- GeneratePlotStrings(path=output.path, prefix=prefix.plot, plot.type="Boxplot")
#'
#'     new.df <- PrepareDataFrameForGGBoxplot(data.frame.to.plot=data.frame.to.plot, design.matrix=design.matrix, groupColName=groupColName, to.log=toLog)
#'
#'     ggp <- ggplot(new.df, aes_string(x="samples", y="names", fill=groupColName)) + geom_boxplot(position=position_dodge(2)) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle(strings$title)
#'
#'
#'     if(save.plot) {
#'         if(is.null(strings$plot.folder)) {
#'             stop("Please set a folder where to plot the boxplot!")
#'         }
#'         if(!is.null(strings$plot.file.name)){
#'             SaveGGplot(ggplot.to.save=ggp, plot.folder=strings$plot.folder, plot.file.name=strings$plot.file.name, plotly.flag=plotly.flag)
#'         }
#'
#'     }
#'
#'     if(show.plot.flag) {
#'         if(plotly.flag) {
#'             ggplotly(ggp)
#'         } else {
#'             plot(ggp)
#'         }
#'     }
#'
#' }
#'
#'
#' SaveGGplot <- function(ggplot.to.save, plot.folder,
#'                        plot.file.name, plotly.flag=FALSE) {
#'     #' Save a plot generated with GGPlot on disk.
#'     #'
#'     #' @param ggplot.to.save a ggplot or ggplotly object.
#'     #' @param plot.folder Character. the folder where to save the plot.
#'     #' @param plot.file.name Character. The file.name for the plot file.
#'     #' @param plotly.flag logical. Is it a ggplotly object?
#'     #'
#'     #' @return none
#'     #'
#'     #' @export
#'     #'
#'     #' @importFrom grDevices pdf dev
#'     #' @import htmlwidgets
#'
#'     require("htmlwidgets")
#'     plotname <- file.path(plot.folder, paste0(plot.file.name, ".pdf"))
#'     i<-1
#'     while(file.exists(plotname)) {
#'         plotname <- file.path(plot.folder, paste0(plot.file.name, i,".pdf"))
#'         i<-i+1
#'     }
#'
#'     if(plotly.flag) {
#'         htmlwidgets::saveWidget(ggplotly(ggplot.to.save),
#'                                 paste0(getwd(),"/", plotname, ".html"))
#'         #paste0(getwd(), "/", plotname, ".html"))
#'         message("Plot saved as ", paste0(plotname, ".html"))
#'     } else {
#'         pdf(plotname)
#'         # plot(ggplot.to.save)
#'         print(ggplot.to.save)
#'         dev.off()
#'         message("Plot saved as ", plotname)
#'     }
#' }
