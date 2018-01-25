
FilterLowCounts <- function(counts.dataframe, design.dataframe, is.normalized=c(TRUE, FALSE), method.type=c("CPM", "Wilcoxon", "Proportion"), cv.percentage, cpm.cutoff, seq.depth=NULL) {
  ## TO BE TESTED
  #try(
    if( !("Conditions" %in% colnames(design.dataframe))) stop("Design dataframe must contain Conditions on colnames")#,
    if( sum(rownames(design.dataframe) %in% colnames(counts.dataframe)) != length(rownames(design.dataframe)) ) stop("Rownames of Design dataframe must be equal to Colnames of Counts dataframe ")
  #)
  # design.dataframe <- design.dataframe[rand,]
  # design.dataframe[which( rownames(design.dataframe) %in% colnames(counts.dataframe)), "Conditions"]

  design.dataframe <- design.dataframe[order(rownames(design.dataframe)), , drop=FALSE]
  counts.dataframe <- counts.dataframe[, order(colnames(counts.dataframe)), drop=FALSE]

  sub.counts.dataframe <- counts.dataframe[, which( colnames(counts.dataframe) %in% rownames(design.dataframe)) ]

  conditions <- design.dataframe[,"Conditions"]

  if(method.type == "Proportion") {
    if(is.normalized) {
      ## calcolare seq.depth
      if(is.null(seq.depth)) stop("Proportion test cannot be performed on normalized counts without sequencing depth!\nYou need column totals before normalizing the data.")
    } else {
      seq.depth <- NULL
    }
  }

  switch( method.type,
          CPM={
            method.number <- 1
          },
          Wilcoxon={
            method.number <- 2
          },
          Proportion={
            method.number <- 3
          }
  )

  require("NOISeq")
  print(dim(sub.counts.dataframe))
  filtered.dataframe <- NOISeq::filtered.data(dataset=sub.counts.dataframe, factor=conditions, norm=is.normalized, depth=seq.depth, method=method.number, cv.cutoff=cv.percentage, cpm=cpm.cutoff, p.adj="BH")
  # filtered.data(x, factor=conditions, norm=norm, depth=depth, method=3, cpm=cpm)
  return(filtered.dataframe)
}

DeSeqTimeCourse <- function(counts.dataframe, design.dataframe, test.type="LRT") { ## DA PERFEZIONARE CON ULTERIORI ARGOMENTI
  ## This function computes differential expressed genes in a time course experiment with DeSeq2 package
  ## input:
  ## @counts.dataframe: a dataframe of counts with colnames equal to rownames of design.dataframe
  ## @design.dataframe: a dataframe with colnames named as Times and Conditions and rownames equal to colnames of counts.dataframe
  ## output:
  ## @readable.results: a dataframe within results for DE genes

  try(
    if( (!("Times" %in% colnames(design.dataframe)) ||  !("Conditions" %in% colnames(design.dataframe))))
      stop("Design data frame must contain Times and Conditions colnames")
  )

  require("DESeq2")


  if(test.type != "LRT") {
    deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Conditions)
    # deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
    de.results <- DESeq2::DESeq(deseq.dataset, test="Wald")
  } else {
    deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Times + Conditions + Times:Conditions)
    deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
    de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times + Conditions)
  }

  readable.results <- DESeq2::results(de.results)
  #readable.results.df <- as.data.frame(readable.results)
  # return(readable.results.df)
  return(readable.results)
}

ApplyDeSeqBackup <- function(counts.dataframe, design.dataframe, test.type="LRT") { ## DA PERFEZIONARE CON ULTERIORI ARGOMENTI
  ## This function computes differential expressed genes in a time course experiment with DeSeq2 package
  ## input:
  ## @counts.dataframe: a dataframe of counts with colnames equal to rownames of design.dataframe
  ## @design.dataframe: a dataframe with colnames named as Times and Conditions and rownames equal to colnames of counts.dataframe
  ## output:
  ## @readable.results: a dataframe within results for DE genes

  try(
    if( (!("Times" %in% colnames(design.dataframe)) ||  !("Conditions" %in% colnames(design.dataframe))))
      stop("Design data frame must contain Times and Conditions colnames")
  )

  require("DESeq2")

  if(test.type == "Wald") {
    deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Conditions)
    deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
    de.results <- DESeq2::DESeq(deseq.dataset)
  } else if (test.type == "LRT") {
    deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Times + Conditions + Times:Conditions)
    deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
    # de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times + Conditions)
    de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times)
  }


  # else {
  #   deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Times + Conditions + Times:Conditions)
  #   deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
  #   de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times + Conditions)
  #   # de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Conditions)
  # }

  readable.results <- DESeq2::results(de.results)
  #readable.results.df <- as.data.frame(readable.results)
  # return(readable.results.df)
  return(as.data.frame(readable.results))
}


ApplyDeSeq <- function(counts.dataframe, design.dataframe, test.type=c("LRT_T", "LRT_TC", "LRT_NoInteraction", "Wald")) { ## DA PERFEZIONARE CON ULTERIORI ARGOMENTI
    ## This function computes differential expressed genes in a time course experiment with DeSeq2 package
    ## input:
    ## @counts.dataframe: a dataframe of counts with colnames equal to rownames of design.dataframe
    ## @design.dataframe: a dataframe with colnames named as Times and Conditions and rownames equal to colnames of counts.dataframe
    ## output:
    ## @readable.results: a dataframe within results for DE genes

    try(
        if( (!("Times" %in% colnames(design.dataframe)) ||  !("Conditions" %in% colnames(design.dataframe))))
            stop("Design data frame must contain Times and Conditions colnames")
    )

    require("DESeq2")

    if(test.type == "Wald") {
        deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Conditions)
        deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])

        de.results <- DESeq2::DESeq(deseq.dataset)
        readable.results <- as.data.frame(DESeq2::results(de.results))
    } else if (test.type == "LRT_T") {
        deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Times + Conditions + Times:Conditions)
        deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
        # de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times + Conditions)
        de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times)
        wald.conditions.colname <- colnames(de.results@rowRanges@elementMetadata)[grep("Conditions_" , colnames(de.results@rowRanges@elementMetadata))][1]

        wald.results <- as.data.frame(DESeq2::results(de.results, name=wald.conditions.colname, test="Wald"))
        lrt.results <- as.data.frame(DESeq2::results(de.results))
        readable.results <- list("LRT"=lrt.results, "waldLRT"=wald.results)
    } else if(test.type == "LRT_TC") {
        deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Times + Conditions + Times:Conditions)
        deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
        de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times + Conditions)

        wald.conditions.colname <- colnames(de.results@rowRanges@elementMetadata)[grep("Conditions_" , colnames(de.results@rowRanges@elementMetadata))][1]
        wald.results <- as.data.frame(DESeq2::results(de.results, name= wald.conditions.colname, test="Wald"))
        lrt.results <- as.data.frame(DESeq2::results(de.results))
        readable.results <- list("LRT"=lrt.results, "waldLRT"=wald.results)
    } else if(test.type == "LRT_NoInteraction") {
        deseq.dataset <- DESeq2::DESeqDataSetFromMatrix(countData=counts.dataframe, colData=design.dataframe, design=~ Times + Conditions)
        deseq.dataset$Conditions <- relevel(deseq.dataset$Conditions, ref=levels(deseq.dataset$Conditions)[2])
        de.results <- DESeq2::DESeq(deseq.dataset, test="LRT", reduced=~ Times)
        wald.conditions.colname <- colnames(de.results@rowRanges@elementMetadata)[grep("Conditions_" , colnames(de.results@rowRanges@elementMetadata))][1]

        wald.results <- as.data.frame(DESeq2::results(de.results, name= wald.conditions.colname, test="Wald"))
        lrt.results <- as.data.frame(DESeq2::results(de.results))
        readable.results <- list("LRT"=lrt.results, "waldLRT"=wald.results)
    }

    return(readable.results)
}


### NOISEQ
ApplyNOISeq <- function(counts.dataframe, design.dataframe, factors.column=NULL,
                        conditions.vector=NULL,
                        test.type=c("NOISeqBio", "NOISeq")) { ## DA PERFEZIONARE CON ULTERIORI ARGOMENTI
    ## This function computes differential expressed genes with NOISeq package
    ## input:
    ## @counts.dataframe: a dataframe of counts with colnames equal to rownames of design.dataframe
    ## @design.dataframe: a dataframe with colnames named as Times and Conditions and rownames equal to colnames of counts.dataframe
    ## output:
    ## @readable.results: a dataframe within results for DE genes

    try(
        if( !(factors.column %in% colnames(design.dataframe)))
            stop("Design data frame must contain ", factors.column, "in colnames")
    )

    require("NOISeq")
    idx.sam <- which(colnames(counts.dataframe) %in% rownames(design.dataframe))
    sub.counts.dataframe <- counts.dataframe[ , idx.sam, drop=FALSE ]

    if(test.type == "NOISeqBio") {
        noiseq.dataset <- NOISeq::readData(data=sub.counts.dataframe, factors=design.dataframe)
        de.results <- noiseqbio(input=noiseq.dataset, norm="n", factor=factors.column, conditions=conditions.vector, filter=0, r=20)
    }

    if(test.type == "NOISeq") {
        # print(head(sub.counts.dataframe))
        # print(design.dataframe)
        noiseq.dataset <- NOISeq::readData(data=sub.counts.dataframe,
                                           factors=design.dataframe)
        de.results <- noiseq(input=noiseq.dataset, norm="n", factor=factors.column,
                             conditions=conditions.vector,
                             k=0.5, pnr=0.2, nss=5, v=0.02, lc=0,
                             replicates="technical")
    }

    # length(which(de.results@results[[1]]$prob > 0.99))
    # readable.results <- DESeq2::results(de.results)
    #readable.results.df <- as.data.frame(readable.results)
    # return(readable.results.df)
    return(as.data.frame(de.results@results[[1]]))
}



SelectGenesFromHistBreaks <- function(data.counts, breakss=100, n.genes.per.break=1) {
  data.counts.mean <- apply(data.counts, 1, mean)

  data.counts.mean.filtered <- data.counts.mean[which(data.counts.mean <= summary(data.counts.mean)["3rd Qu."])]
  h <- hist(data.counts.mean.filtered, breaks=breakss, plot=FALSE)
  data.counts.mean.sorted <- sort(data.counts.mean.filtered)
  gene.names <- c()
  j<-1
  for(i in 1:length(h$counts)) {
    names <- names(data.counts.mean.sorted[ j : (j+h$counts[i]-1) ])
    if(length(names) >= n.genes.per.break) {
      n.g.names <- names[c(1:n.genes.per.break)]
    } else {
      n.g.names <- names[c(1:length(names))]
    }
    gene.names <- c(gene.names, n.g.names)

    j <- j+h$counts[i]
  }
  return(gene.names)
}

EstimateNegativeControlGenesForRUV <- function(de.genes, n.tail.genes=2000, counts.dataset, n.genes.per.hist.break=1, threshold=0.05) {
  de.genes.tail <- TailDeGenesByPAdj(de.genes=de.genes, n.genes=n.tail.genes, threshold=threshold)
  counts.tail <- counts.dataset[which(rownames(counts.dataset) %in% de.genes.tail),]
  estimated.genes <- SelectGenesFromHistBreaks(counts.tail, n.genes.per.break=n.genes.per.hist.break)
  return(estimated.genes)
}


RUVgNormalizationFunction <- function(data.to.normalize, design.matrix, estimated.gene.names) {

  if( length( which(colnames(design.matrix) == "Conditions")) == 0 ) {
    stop("Please provide a design matrix with a named \"Conditions\" column")
  }
  require("RUVSeq")
  dataset.set <- newSeqExpressionSet( as.matrix(data.to.normalize),
                                               phenoData= data.frame(
                                                 design.matrix$Conditions, row.names=rownames(design.matrix)
                                               ))

  ruved.set <- RUVg(dataset.set, estimated.gene.names, k=1)

  return(as.data.frame(ruved.set@assayData$normalizedCounts))
}

NormalizeData <- function(data.to.normalize, norm.type=c("fqua", "uqua", "tmm", "ruvg"), design.matrix=NULL, estimated.genes=NULL) {
  ## @ norm.type can be uqua, tmm, fqua or ruvg

  x <- data.to.normalize

  if(norm.type != "fqua") {
    if(norm.type == "ruvg") {
      if(!(is.null(design.matrix) && is.null(estimated.genes))){ ## test
        normalized.data <- RUVgNormalizationFunction(data.to.normalize, design.matrix, estimated.genes)
      } else {
        stop("Please select a design matrix and a list of negative control genes for RUVSeq normalization")
      }
    } else {
      require("edgeR")
      x <- edgeR::DGEList(counts=x)
      if(norm.type=="uqua") {
        x <- edgeR::calcNormFactors(x, method='upperquartile')
      }

      if(norm.type=="tmm") {
        x <- edgeR::calcNormFactors(x, method='TMM')
      }

      x <- edgeR::estimateCommonDisp(x, verbose=FALSE)
      x <- edgeR::estimateTagwiseDisp(x)
      normalized.data <- as.data.frame(x$pseudo.counts)
    }

  } else if(norm.type == "fqua") {
    require("preprocessCore")
    x <- as.matrix(x)
    normalized.data <- as.data.frame(preprocessCore::normalize.quantiles(x, copy=FALSE))
    colnames(normalized.data) <- colnames(x)
    rownames(normalized.data) <- rownames(x)
  }

  return(normalized.data)
}



ConvertGeneNamesToEnsembl <- function(map.dataframe, named.dataset) {

  datanames <- data.frame(rownames= rownames(named.dataset), ensembl_ids=NA, gene_names=NA, ind_map=NA, ens_present=FALSE)

  ind.genes.to.map <- which(tolower(datanames$rownames) %in% tolower(map.dataframe$external_gene_name))
  dataset.ens <- named.dataset
  i <- 0
  for (ind in ind.genes.to.map) {
    gene.name <- rownames(dataset.ens)[ind]

    ind.map <- which(tolower(map.dataframe$external_gene_name) == tolower(gene.name))

    if(length(ind.map) == 0) {

      cat("gene name, ", gene.name," not present! STRANGE!\n")
      #readline()

    } else if(length(ind.map) == 1) {

      ind.check <- which( tolower(rownames(dataset.ens)) == tolower(map.dataframe$ensembl_gene_id[ind.map]))

      if(length(ind.check) == 0) {
        # rownames(counts.o.ens)[ind] <- interesting.map.o$ensembl_gene_id[ind.map]
        datanames[ind,"ensembl_ids"] <- map.dataframe$ensembl_gene_id[ind.map]
        datanames[ind,"ind_map"] <- ind.map

      } else {
        datanames[ind.check,"gene_names"] <- gene.name
        datanames[ind,"ensembl_ids"] <- map.dataframe$ensembl_gene_id[ind.map]
        datanames[ind,"ens_present"] <- TRUE
        # repeated.ensembl <- c(repeated.ensembl, map.dataframe$ensembl_gene_id[ind.map])
        # repreated.name <- c(repeated.name, gene.name)
        # cat("[", ind, "] The ", map.dataframe$ensembl_gene_id[ind.map], " is already present, check ", gene.name, "gene!\n")
        # readline()
      }

    } else if(length(ind.map) > 1) {
      datanames[ind, "ensembl_ids"] <- paste(map.dataframe$ensembl_gene_id[ind.map], collapse="; ")
      datanames[ind,"ind_map"] <- paste(ind.map, collapse="; ")
      ## do nothing, for now
    }
    i <- i+1
    if( i %% 100 == 0 ) cat(i, " of ", length(ind.genes.to.map), "\n")
  }
  write.table(datanames, file="./rownames_map.txt", sep="\t", row.names=FALSE, quote=FALSE)
}

CountBamFilesFeatureCounts <- function(bam.files, annotation.file, output.folder=NULL, prefix.output.file, gtf.attr.type="gene_id", n.thread=8, strand.specific=0, is.paired.end=FALSE) {
  # annotation.file <- "/media/dario/dati/time_course/Cervical/cuffmerge/merged_assembly/merged.gtf"
  # counts.folder <- file.path("/media/dario/dati/time_course/Cervical", "results", "counts", "FeatureCounts")
  fc_SE <- Rsubread::featureCounts(files=bam.files,
                                   annot.ext=annotation.file,
                                   isGTFAnnotationFile=TRUE,
                                   GTF.featureType='exon', ## def
                                   GTF.attrType=gtf.attr.type, ## gene_id def
                                   useMetaFeatures=TRUE, ## def
                                   allowMultiOverlap=FALSE, ## def
                                   nthreads=n.thread,
                                   strandSpecific=strand.specific, ## 0 def
                                   countMultiMappingReads=FALSE, ## def
                                   isPairedEnd=is.paired.end) ## FALSE def
  # print(head(fc_SE[1]))
  #fc.bam.file.name <- paste0( tail(unlist(strsplit(bam.file, "/")), n=1), "_FeatureCounts.txt")
  fc.count.file.name <- paste0(prefix.output.file, "_Counts_FeatureCounts.txt")
  counts.folder <- UpdateFolderPath(output.folder, "counts/FeatureCounts")
  fc.count.file.path <- file.path(counts.folder, fc.count.file.name)
  # dir.create(counts.folder, recursive=TRUE, showWarnings=FALSE)
  write.table(fc_SE[1]$counts, file=fc.count.file.path , quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
  message("file ", fc.count.file.path, " written on disk!\n")
  return(fc_SE[1]$counts)
}

SortDeGenesByPAdj <- function(de.genes) {
  if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
    ordered.de.genes <- de.genes[order(de.genes$padj),]
  } else if("theta" %in% colnames(de.genes)) { ## noiseqbio
    ordered.de.genes <- de.genes[order(de.genes$prob, decreasing=TRUE),]
  } else if("M" %in% colnames(de.genes)) { ## noiseq
    ordered.de.genes <- de.genes[order(de.genes$prob, decreasing=TRUE),]
  } else {
    stop("In SortDeGenesByPAdj not recognized de results table!")
  }

  # ordered.de.genes <- de.genes[order(de.genes$padj),]
  return(ordered.de.genes)
}

FilterOutGenesNAPAdj <- function(de.genes) {
  if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
    if(sum(is.na(de.genes$padj))>0) de.genes.not.na <- de.genes[-which(is.na(de.genes$padj)),]
    else de.genes.not.na <- de.genes
  } else {
    de.genes.not.na <- de.genes
  }

  # if(sum(is.na(de.genes$padj))>0) de.genes.not.na <- de.genes[-which(is.na(de.genes$padj)),]
  # else de.genes.not.na <- de.genes
  return(de.genes.not.na)
}

SignificantDeGenesPAdj <- function(de.genes, threshold=0.05) {

  if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
    de.genes.not.na <- FilterOutGenesNAPAdj(de.genes)
    sign.de.genes <- subset(de.genes, padj < threshold)
    de.up <- subset(sign.de.genes, log2FoldChange >= 0)
    de.down <- subset(sign.de.genes, log2FoldChange < 0)
    # sign.de.genes <- de.genes[ de.genes.not.na$padj < threshold,]
    # de.up <- sign.de.genes[sign.de.genes$log2FoldChange >= 0,]
    # de.down <- sign.de.genes[sign.de.genes$log2FoldChange < 0,]
  } else if("theta" %in% colnames(de.genes)) { ## noiseqbio
    print(head(de.genes))
    de.genes.not.na <- de.genes
    sign.de.genes <- subset(de.genes, prob > threshold)
    de.up <- subset(sign.de.genes, log2FC >= 0)
    de.down <- subset(sign.de.genes, log2FC < 0)
    # sign.de.genes <- de.genes[ de.genes$prob > threshold, ]
    # de.up <- sign.de.genes[ sign.de.genes$log2FC >= 0, ]
    # de.down <- sign.de.genes[ sign.de.genes$log2FC < 0, ]
  } else if("M" %in% colnames(de.genes)) { ## noiseq
    print(head(de.genes))
    de.genes.not.na <- de.genes
    sign.de.genes <- subset(de.genes, prob > threshold)
    de.up <- subset(sign.de.genes, M >= 0)
    de.down <- subset(sign.de.genes, M < 0)
    # sign.de.genes <- de.genes[ de.genes$prob > threshold, ]
    # de.up <- sign.de.genes[ sign.de.genes$log2FC >= 0, ]
    # de.down <- sign.de.genes[ sign.de.genes$log2FC < 0, ]
  } else {
    stop("In SignificantDeGenesPAdj not recognized de results table!")
  }
  return( list("de.total"=de.genes, "de.not.na"=de.genes.not.na, "SIGN"=sign.de.genes, "UP"= de.up, "DOWN"=de.down))
}

TailDeGenesByPAdj <- function(de.genes, n.genes=10, threshold=0.05) {
  ord.de.genes <- SortDeGenesByPAdj(de.genes=de.genes)
  ord.de.genes.not.na <- FilterOutGenesNAPAdj(ord.de.genes)
  if("baseMean" %in% colnames(de.genes)) { ## working on DESeq results
    ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$padj >= threshold,]
  } else if("prob" %in% colnames(de.genes)) { ## noiseqbio
    ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$prob <= threshold,]
  }
  # ord.not.de.genes <- ord.de.genes.not.na[ ord.de.genes.not.na$padj >= threshold,]
  tail.genes <- rownames(tail(ord.not.de.genes, n=n.genes))
  return(tail.genes)
}

