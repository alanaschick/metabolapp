## Read multiple files in from csv
## requires rNMR
readBatch <- function (inPath ){
  
  if(missing(inPath))
    inPath <- myOpen()
  if( length(inPath) < 1 )
    return(invisible())
  
  outList <- list()
  for( i in 1:length(inPath)){
    outList[[i]] <- readPeak(inPath = inPath[i] ) 
  }
  
  names(outList) <- inPath
  
  return(outList)
  
}


## General functions
myMax <- function(x){max(x, na.rm = T)}
myWhich.max <- function(x) { out <- myMax(x); return(which(x == out)[1])}

## Log trasform function
## Bioinformatics: Rocke and Durbin 2003
durbin = function( y, c = 0 )
{
  if( !is.null(c))
    y = log2(( y + sqrt( y^2 + c^2 )) / 2 )
  tmp = !is.finite(as.matrix(y))
  if( sum( tmp ))
    y[tmp] = NA
  y
}

## Standard error function that removes NAs
stdErr <- function( x ){
  std.err <- function( y ){ 
    if( any( is.na(y)) )
      y <- y[-which(is.na(y)) ]
    return(sqrt( var(y) / length(y) ))
  }
  if( is.null(nrow(x)) || nrow(x) < 2 )
    return( std.err(x) )
  return(apply( x, 2, std.err)) 
}


## Function for generating barplots
plotBar <- function( sLab, y, cpdNames, stdError = T, byGroup = T, ylim, col, 
                     ... ){
  
  ## Check sLab format
  sam <- unique(sLab)	
  if( !is.factor(sLab) ){
    sLab <- factor(sLab)
  }
  
  if(!is.data.frame(y))
    y <- data.frame(y)
  
  ## Check y format
  if( is.null(nrow(y)) ){
    if( !length(y) )
      stop( 'y must be a dataframe or vector')
    y <- data.frame( matrix(y, nrow = 1 ))
  }
  
  ## Check cpdNames
  if(missing(cpdNames) ){
    cpdNames <- 1:nrow(y)
  }
  
  ## Find sample mean and sd 
  if( !stdError ){
    inSd <- aggregate(t(y), FUN = function(x) sd(x, na.rm = T), 
                      by = list(sLab))		
  }else
    inSd <- aggregate(t(y), FUN = function(x) stdErr(x), by = list(sLab))	
  inDat <- aggregate(t(y), FUN = function(x) mean(x, na.rm = T), 
                     by = list(sLab) )
  
  ## Fix sample names
  uName <- levels(inDat[,1]) 
  inDat <- data.frame(t(inDat[,-1]))
  inSd <- data.frame(t(inSd[,-1]))
  names(inDat) <- names(inSd) <- uName
  
  ## Reorder colums to match input list
  idx <- match(sam, names(inDat))
  if( ncol(inDat) > 1){
    inDat <- inDat[,idx]
    inSd <- inSd[,idx]
  }
  
  ## set ylim
  if( missing(ylim)){
    ylim <- range(inDat)
    if( ylim[1] > 0 )
      ylim[1] <- 0 
    ylim[2] <- ylim[2] + max(unlist(inSd))
  }
  
  ## Set dim
  if( nrow(inDat) == 1)
    byGroup <- TRUE
  if( ncol(inDat) == 1)
    byGroup <- FALSE
  
  ## Format data frame	
  if( !byGroup ){
    inDat <- t(inDat)
    inSd <- t(inSd)
    colnames(inDat) <- colnames(inSd) <-  cpdNames
  }else
    row.names(inDat) <- row.names(inSd) <- cpdNames
  
  ## Set colors
  if(missing(col)){
    col <- 1:length(cpdNames)
    col[1] <- 'gold' 
  }
  
  ## Setup for plot
  out <- 	unlist(barplot( as.matrix(inDat), beside = T, 
                          legend = rownames(inDat), args.legend = list(bty = 'n'), ylim = ylim, 
                          col = col, ...))
  if( nrow(out) < 2 )
    int <- diff(as.vector(unlist(out)))[1]/8
  else
    int <- diff(out)[1]/4
  
  for( i in 1:length(out) ){
    lines( rep(out[i], 2), c(unlist(inDat)[i] - unlist(inSd)[i], 
                             unlist(inDat)[i] + unlist(inSd)[i] ))
    lines( c(out[i]-int, out[i]+int) , 
           rep(unlist(inDat)[i] - unlist(inSd)[i], 2))
    lines( c(out[i]-int, out[i]+int) , 
           rep(unlist(inDat)[i] + unlist(inSd)[i], 2))
  }
  
  
  return(invisible(list(mean = inDat, sd = inSd)))
  
  
}


## make fancy boxplot like figures
## sLab - char vector with the sample group name for each data column
##         eg sLab = c( 'A', 'A', 'A', 'B', 'B', 'B')
## y - numeric table with one column for each sample, all rows will be plotted
## main - plot title
## xlab - x plot label
## ylab - y plot label
## sub - sub x plot label
## add - logical argument - if TRUE, data will be plotted on top of existing plot
## pch - symbol used for plotting
## col - the plot color
## ... - other plotting related arguments, see the points command

## example application
## splat <- readPeak()
## dat <- splat[, -(1:13)] ## remove any non-numeric data columns
## groups <- c(rep('A', 3), rep('B', 3))
## plotXY( sLab = groups, y = dat[1,]) 
plotXY <- function( sLab, y, cpdNames, main = '', xlab = '', ylab = '', sub = '', 
                    add = F, pch = 1, cex = 1, col, xlim, ylim, stack = T, 
                    stdError = T, zPlot = F, vioPlot = F, fixedYlim = F, zeroRemove = T, ... ){
  
  if( !all(is.numeric(sLab)) ){
    idx <- match(sLab, unique(sLab))
    if(!zPlot)
      x <- jitter(idx, 1)
    else
      x <- idx
  }else{
    idx <- sLab
    if(!zPlot)
      x <- jitter(sLab, .25)
    else
      x <- idx
  }
  
  if( is.null(nrow(y)) ){
    if( !length(y) )
      stop( 'y must be a dataframe or vector')
    y <- data.frame( matrix(y, nrow = 1 ))
  }
  
  
  #	if( any(unlist(y) == Inf )  )
  #		y[ y == Inf ] <- max(y[is.finite(unlist(y))], na.rm = T) * 100 	
  #	if( any(unlist(y) == -Inf )  )
  #		y[ y == -Inf ] <- min(y[is.finite(unlist(y))], na.rm = T) / 100
  
  myPad <- function(vec){
    
    delt <- abs(diff(vec))*.1
    return(vec +c(-delt, delt))
  }
  
  if( missing( xlim) )
    xlim <- myPad(range( x, na.rm = T, finite = T))
  if( missing(ylim))
    ylim <- range(as.vector(unlist(y)), na.rm = T, finite = T)
  
  if(!zPlot & !add & !stack & !vioPlot){		
    plot(x = range(x), y = range(as.vector(unlist(y)), na.rm = T, finite = T), 
         axes = F, main = main, xlab = xlab, ylab = ylab, sub = sub, 
         xlim = xlim, ylim = ylim, type = 'n')
    box()
    axis(2)
    axis(1, at = idx, labels = sLab)
    
  }
  
  
  if(zPlot){
    col <- as.vector(matrix(rep(x, nrow(y)), nrow = nrow(y), byrow = T))
    
    if(!add ){
      mar <- par("mar")
      par(mar = c(7, mar[2:4]))			
    }
    
    
    plot(rep(1:nrow(y), ncol(y)), as.vector(unlist(y)), axes = F, main = main, 
         xlab = xlab, ylab = ylab, sub = sub, 
         pch = "-", col = col, cex = 2, ..., ylim = ylim)
    
    box()
    axis(2)
    axis(1, at = 1:nrow(y), labels = row.names(y),las = 2)
    legend( 'topright', unique(sLab), text.col = unique(col), bty = 'n')
    
    
    return(invisible(y))
    
  }
  
  
  if( missing (col) || length(col) != nrow(y) ){
    if( vioPlot )
      col <- c("gold", 2:nrow(y))
    else
      col <- 1:nrow(y)
  }
  
  
  
  if( stack || vioPlot){
    if( nrow(y) > 1 )
      par( mfrow = c(nrow( y ), 1), mar = c(.5, 3, .5, 2), tck = 0)
  }
  
  for (i in 1:nrow(y) ){
    
    if( vioPlot ){
      myVioPlot(sLab = sLab, y = y[i,], col = col[i], padj = -1.5)
      legend( 'topleft', cpdNames[i], text.col = col[i], bty = 'n')
      next
    }
    
    
    if( stack ){
      
      if( !fixedYlim )
        ylim <- NULL
      plot(x, y[i,], axes = F, main = main, xlab = xlab, ylab = ylab, sub = sub, 
           pch = pch, col = col[i], cex = cex, ylim = ylim, ...) 
      
      if( any(is.na(y[i,])) || any( y[i,] == 0)  ){
        tY <- y[i,]
        idxS <- which(is.na(tY))
        tY[ idxS ] <- 0
        idxS <- which(tY == 0)
        points(x[idxS], tY[idxS], pch = pch, col = 'grey', cex = cex)				
      }
      
      box()
      axis(2)
      axis(1, at = idx, labels = sLab, padj = -1.5, tick = F)
      #abline( h = 1, pch = 2)
      
    }else{
      tY <- as.vector(unlist(y))
      points(x, y[i,], pch = pch, col = col[i], cex = cex, ...)	
      if( any(is.na(tY)) ){
        idxS <- which(is.na(tY))
        tY[ idxS ] <- 0
        points(x[idxS], tY[idxS], pch = pch, col = 'grey', cex = cex)				
      }
    }
    
    for( j in unique(idx) ){
      if( zeroRemove == T && any(y[i, which(idx == j)] == 0)){
        y[i, which(idx == j)][which(y[i, which(idx == j)] == 0)] <- NA
      }
      
      
      ySub <- mean(unlist(y[i, which(idx == j)]), na.rm = TRUE)
      if( stdError )
        y2Sub <- stdErr(unlist(y[i, which(idx == j)]))
      else
        y2Sub <- sd(unlist(y[i, which(idx == j)]), na.rm = TRUE)
      xSub <- c( j - .1, j + .1 )
      lines( xSub, rep(ySub, 2), lwd = 2)
      lines( rep(j, 2), c(ySub - y2Sub, ySub + y2Sub)) 
      lines( xSub, rep(ySub - y2Sub, 2), lwd = 2)
      lines( xSub, rep(ySub + y2Sub, 2), lwd = 2)
      
      if(all(is.na(ySub)) ){
        lines( xSub, rep(0, length(xSub)), lwd = 2)
      }
      
    }
    
    
    if( stack && !missing(cpdNames) && length(cpdNames) == nrow(y))
      legend( 'topleft', cpdNames[i], text.col = par('fg'), bty = 'n')	
    
    
  }
  
  if( !stack && !add && !vioPlot && !missing(cpdNames) && length(cpdNames) == nrow(y))
    legend( 'topleft', cpdNames, text.col = par('fg'), bty = 'n')
  
  
}

## Internal function for finding a mode of a dataset; for multimode data
## this function only returns the first index
myMode <- function(x) {
  ux <- sort(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}



## Internal function for compressing multiple comound records
## excludeName - character, names of columns to be excluded from analysis
findMaxRecord <- function( inDat, excludeName ){
  
  ## Get compound name
  com <- inDat$compound
  uCpd <- unique(com)
  
  if( length(com) != nrow(inDat) )
    stop( "compound column required")
  
  tDat <- metaSep(inDat, excludeName = excludeName)
  outDat <- tDat[[2]]
  
  
  ## Remove duplicate metabolite records and find max for each
  usrWarn <- F
  for( i in 1:length(uCpd) ){
    idx <- which( com == uCpd[i] )
    
    if( nrow(tDat[[2]][idx,]) < 2  )
      next
    
    subRow <- apply(tDat[[2]][idx,], 2, myMax)
    subIdx <- apply(tDat[[2]][idx,], 2, myWhich.max)
    if( length(unique(subIdx)) > 1 )		
      usrWarn <- T		
    subIdx <- myMode(apply(tDat[[2]][idx,], 2, myWhich.max))
    tRow <- tDat[[1]][idx,][subIdx,]
    
    
    for( j in 1:nrow(tDat[[2]][idx,]) ){
      tDat[[2]][idx,][j,] <- subRow
      tDat[[1]][idx,][j,] <- tRow
    }
  }
  
  
  outDat <- unique( cbind(tDat[[1]], tDat[[2]]) )
  
  ## reshuffle names to match input list
  outDat <- outDat[,match(names(inDat), names(outDat))]
  
  ## Warn user if mutiple max records were found
  if(usrWarn)
    warning("Some compounds had maxima, meta data report the most common record")
  
  return(outDat)
  
}



## Make realistic noise data
randomize <- function (inDat){
  
  for(i in 1:nrow(inDat)){
    tmp <- unlist(inDat[i,])
    
    inDat[i,] <- rnorm(length(tmp), mean(tmp, na.rm = TRUE), sd(tmp, na.rm = T))
  }
  
  
  return(inDat)
  
}


## Internal function for seperating MAVEN metadata from datasets
## excludeName - character, names of columns to be excluded from analysis
## mavNam - character string - default names of Maven format metadata
metaSep <- function(inDat, 
                    mavNam = c('label', 'metaGroupId', 'groupId', 'goodPeakCount', 
                               'medMz', 'medRt', 'maxQuality', 'note', 'compound', 'compoundId', 
                               'expectedRtDiff', 'ppmDiff', 'parent', 'pVal'), excludeName){
  
  ## Check input dataset
  if( missing(inDat) || !is.data.frame(inDat) || (nrow(inDat) < 1) )
    stop('A MAVEN format dataset is required')
  
  ## Combine default and user-provided names
  if( !missing(excludeName) )
    mavNam <- c(mavNam, excludeName)	
  
  ## Remove all standard Maven meta data
  mIdx <- na.omit(match(mavNam, names(inDat)))
  
  ## Find metadaFta
  if( length(mIdx) == 1 ){
    metaDat <- data.frame(compound = inDat[,mIdx], stringsAsFactors = F)
  }else
    metaDat <- inDat[,mIdx]
  
  ## remove the metaData
  outDat <- inDat[,-na.omit(match(mavNam, names(inDat)))]
  
  return(list( metaDat, outDat))
  
}




## Default maven plotter
## This function takes a maven-format list and makes several types of plots
## inDat - maven format table
## rCst - Logical - TRUE clusters by row
## cCst - Logical - TRUE clusters by column
## pCst - Logical - TRUE ranks rows by P value 
## nCst - Logical - TRUE sorts samples by name
## grid - Logical - TRUE plots the a sample grid over the heatmap
## pca - Logical - TRUE makes the PCA plots
## linePlot - Logical - TRUE makes a line plot
## stackPlot - Logical - TRUE makes a staked line plot
## heatMap - Logical - TRUE makes a classic heatmap
## barPlot - Logical - Plots bar figure, note this needs replicates see "by"
## zPlot - Logical - this plots data relative to a reference sample set
## scale - controls how data are normalized
##			must be one of the following: "log", "durbin", "row", "column", or "none"
##			'log' applies base 2 log transform
##          'durbin' applies durbin transform
##          'row' normalizes by row
##          'column' noralizes by column
##			'fold'...?
##			'ratio'...?
## 			'zScore' normalizes relative to a control, sets equal variance
##          'none' does not apply any normalization
## c - integer - constant for durbin trasform (0 is ln)
## n - number of levels used in heatmap
## zlim - numeric vector of length 2 for the c(min, max) plotted value
##      if zlim is missing or NULL, then all data are plotted
## excludeName - character, names of columns to be excluded from analysis
## na.rm - Logical - TRUE converts NA to 0, INF to positive max, 
##                   -INF to negative max
## maxDup - logical, finds the maximum record in cases of duplicated compounds
## inverse - logical, inverts the incoming data (1/data)
## pThresh - logical - excludes data above the pValue threshold
## bonferroni - logical - calculates the bonferroni-corrected pThresh
## returnPCA - logical - TRUE returnes PCA scores
## returnMeta - Logical, if T, will return all input meta data
##             note: this cannot be used with column/p clustering or thresholding
## rmZero - logical - removes rows with all zero intensities
## permutation - logical, if TRUE it will radomize sample names
## random - Logical, TRUE will generate random number distribution that match
##          the mean/sd for each compound
## returns the dataframe used to plot the data
mavPlot <- function(	inDat, rCst = F, cCst = F, pCst = F, nCst = F, 
                     grid = FALSE, individuals = T, variables = F, pcacolour = F,
                     scale = c("log", "durbin", "row", "column", "fold", "ratio", "zScore", "none"), 
                     inverse = F, c = 1, n = 10, zlim, maxDup = F, avgRep = F, 
                     by = "\\.", excludeName = NULL, scaleby = NULL, compound_list = NULL,
                     byClass = F, na.rm  = T, sLab, dThresh, group_list = NULL, pheatmap = F,
                     pca = F, linePlot = F, stackPlot = F, heatMap = F, barPlot = F,  zPlot = F, 
                     vioPlot = F, pCalc = F, pThresh = 1, bonferroni = F, ylim, xlim, 
                     returnMeta = T, returnPCA = F, rmZero = F, permutation = F, random = F, ...){
  
  ## check format of peaklist
  if( !is.data.frame(inDat) || nrow(inDat) < 1 || ncol(inDat) < 2 
      || !any(names(inDat) == 'compound'))
    stop("A dataframe with a 'compound' and at least one data column is required")
  
  ## Set Scale argument
  scale <- if (missing(scale)) 
    "none"
  else match.arg(scale)
  
  ## Check pThresh
  if( avgRep && pThresh < 1)
    stop("avgRep and pThresh can not be used in combination")
  if( pThresh < 1 && bonferroni )
    stop("pThresh and bonferroni can not be used in combination")
  
  ## Make all compound names unique
  #######BUG to fix: numberic data crash the compound list
  com <- inDat$compound
  if(any(is.na(com)) &&  any(names(inDat) == 'medMz') 
     && any(names(inDat) == 'medRt')){
    com[which(is.na(com))] <- paste(inDat$medMz[which(is.na(com))], 
                                    inDat$medRt[which(is.na(com))], sep = "@")
  }
  
  for( i in unique(com) ){
    sub <- com[ com == i ]
    if( length(sub) == 1 )
      next
    com[com == i] <- paste(sub, 1:length(sub), sep = '_')
  }
  inDat$compound <- com
  
  
  ## Remove duplicate metabolite records and find max for each
  if(maxDup)
    inDat <- findMaxRecord( inDat, excludeName = excludeName )
  
  ## Combine any metadata fields to be excluded
  tDat <- metaSep(inDat, excludeName = excludeName)
  metaDat <- tDat[[1]]
  inDat <- tDat[[2]]
  
  ## Order samples by name
  if (nCst){
    inDat <- inDat[,order(names(inDat))]
  }
  
  ## hold sample names
  sam <- names(inDat) ## read tab header to get sample names
  rownames(inDat) <- metaDat$compound
  
  ## Zero any data below the data threshold
  if( !missing(dThresh) )
    inDat[ inDat < dThresh ] <- 0
  
  ## Remove all zero data
  if( any(!apply(inDat, 1, function(x) all(x == 0))) && rmZero ){
    inDat <- inDat[- which(apply(inDat, 1, function(x) all(x == 0))),]		
  }
  
  
  ## Randomize data labels
  if( permutation ){
    tNam <- names(inDat)
    inDat <- inDat[sample(1:ncol(inDat), ncol(inDat))]
    names(inDat) <- tNam
  }
  
  if( random )
    inDat <- randomize(inDat)
  
  
  ## Condenses metabolite signals according to compound classes
  if (byClass ){
    
    if( any(is.na(metaDat$note)) )
      stop("The note column must contain compound classes")
    
    cpdClass <- metaDat$note			
    out <- list()
    for(i in 1:ncol(inDat))	
      out[[i]] <- aggregate(inDat[,i], by = list(cpdClass), "sum")
    inDat <- data.frame(sapply(out, function(x) x[,2]))
    metaDat <- data.frame( compound = out[[1]][,1], stringsAsFactors = F)
    names(inDat) <- sam
    rownames(inDat) <- metaDat$compound
  }
  
  
  
  ## Get sample labels
  if( (barPlot || scale == 'zScore' || pCst || stackPlot || vioPlot ||
       linePlot || zPlot || pCalc || pThresh < 1) && missing (sLab))
    sLab <- sapply( strsplit(sam, split = by), function (x) x[1])		
  
  
  ## Average sample replicates
  if( avgRep ){
    
    tLab <- list(factor(
      sapply(strsplit(names(inDat), split = by), function (x) x[1])))
    
    inDat <- aggregate(t(inDat), FUN = function(x) mean(x, na.rm = T), by =tLab)
    
    
    uName <- levels(inDat[,1]) 
    inDat <- data.frame(t(inDat[,-1]))
    names(inDat) <- uName
    
    ## Reorder colums to match input list
    sam <- unique(sapply(strsplit(sam, split = by), function (x) x[1]))
    inDat <- inDat[,c(which(names(inDat) == 'compound'), 
                      match(sam, names(inDat)))]
  }
  
  
  ## Scale data
  if( scale == 'log'){
    inDat <- log2( inDat )	
  }else if( scale == 'durbin' ){
    inDat <- durbin(inDat, c = c )	
  }else if (scale == "row") {
    inDat <- sweep(inDat, 1L, rowMeans(inDat, na.rm = TRUE), check.margin = FALSE)
    sx <- apply(inDat, 1L, sd, na.rm = TRUE)
    inDat <- sweep(inDat, 1L, sx, "/", check.margin = FALSE)
  }else if (scale == "column") {
    inDat <- sweep(inDat, 2L, colMeans(inDat, na.rm = TRUE), check.margin = FALSE)
    sx <- apply(inDat, 2L, sd, na.rm = TRUE)
    inDat <- sweep(inDat, 2L, sx, "/", check.margin = FALSE)
  }else if (scale == "fold" || scale == "ratio" || scale == 'zScore'){
    userSel <- scaleby
    normCol <- match(userSel, names(inDat))
    inDat <- foldChange( inTable = inDat, normVec = normCol, scale = scale)	
  }
  
  ## Remove NA and Inf data
  if( na.rm  ){
    printList <- newList <- NULL
    
    if( any(is.na(inDat)) ){
      inDat[is.na(inDat)] <- 0
      printList <- c(printList, "NA")
      newList <- c(newList, 0)
    }
    
    if( any(inDat == Inf)){
      inDat[ inDat == Inf ] <- range(unlist(inDat), finite = T)[2]	
      newList <- c(printList, "Inf")
      newList <- c(newList, max(inDat, na.rm = T))
    }
    
    if( any(inDat == -Inf)){
      inDat[ inDat == -Inf ] <- range(unlist(inDat), finite = T)[1]	
      printList <- c(printList, "-Inf")
      newList <- c(newList, min(inDat, na.rm = T))
    }
    if(!is.null(printList))
      print( paste(c(paste(printList, collapse = ", "), 
                     "data values were changed to:",	
                     paste(newList, collapse = ", ")), collapse = " "))
  }	
  
  
  if( inverse )
    inDat <- 1/inDat
  
  
  ## Rank values according to ANOVA p value
  if( pCst || pCalc || pThresh < 1){
    
    metaDat$pVal <- findPval(sLab = sLab, inDat, type = 'p')
    
    if( pCst ){
      sidx <- order(metaDat$pVal)
      inDat <- inDat[sidx,]
      if(ncol( metaDat) > 1 )
        metaDat <- metaDat[sidx,]
      else
        metaDat <- data.frame( compound = metaDat[sidx,], stringsAsFactors = F)
    }		
    
    
    ## Calculate pThreshold
    if( bonferroni ){
      
      ## this needs some research, how to correct ANOVA?
      #pThresh <- .05 / (nrow(inDat) * length(unique(sLab)))
      pThresh <- .05 / nrow(inDat) 
      
    }
    
    
    if( pThresh < 1 ){
      idxP <- which(metaDat$pVal < pThresh)
      inDat <- inDat[idxP,]
      if(ncol( metaDat) > 1 )
        metaDat <- metaDat[idxP,]
      else
        metaDat <- data.frame( compound = metaDat[idxP,], stringsAsFactors = F)	
    }
  }
  
  
  
  ## Generate PCA plots
  if( pca ){
    outPCA <- pcaPlot( dat = inDat, samples = sam, compounds = metaDat$compound)
    
    ## Get eigenvalues to scale axes appropriately
    eigenvalues <- get_eigenvalue(outPCA)
    eig <- eigenvalues[,1]
    
    ## Set colours to sample types
    sLab <- sapply( strsplit(sam, split = "\\."), function (x) x[1])
    num_groups <- length(unique(sLab))
    set.seed(102)
    cc20 <- sample(rainbow(20, v = 0.8))
    cols <- cc20[1:num_groups]
    SampleType <- sLab
    
    ## Plot
    if (individuals & !variables){
      if (!pcacolour){
        pp <- fviz_pca_ind(outPCA) + coord_fixed(sqrt(eig[2]/eig[1]))
      }
      if (pcacolour){
        pp <- fviz_pca_ind(outPCA, palette = cols) + coord_fixed(sqrt(eig[2]/eig[1])) + geom_point(aes(color = SampleType))
      }
    }
    if (!individuals & variables){
      pp <- fviz_pca_var(outPCA) + coord_fixed(sqrt(eig[2]/eig[1]))
    }		
    if (individuals & variables){
      if (!pcacolour){
        pp <- fviz_pca_biplot(outPCA) + coord_fixed(sqrt(eig[2]/eig[1]))
      }
      if (pcacolour){
        pp <- fviz_pca_biplot(outPCA) + coord_fixed(sqrt(eig[2]/eig[1])) + geom_point(aes(color = SampleType)) + scale_color_manual(values = cols)
      }	
    }
      print(pp)
  }			
    
    ## Get user selections
    if( barPlot || linePlot || zPlot || stackPlot || vioPlot){
      userSel <- compound_list
      usrRow <- match(userSel, metaDat$compound )
    }
    
    ## Generate bar plot	
    if( barPlot ){
      
      userSel <- group_list
      usrCol <- myMatch(userSel, sLab )
      plotBar(sLab = sLab[usrCol], y = inDat[usrRow,usrCol], 
              cpdNames = row.names(inDat)[usrRow], ylim = ylim, ...)
      
      invisible( inDat[usrRow,usrCol] ) 
      
    }
    
    ## Generate XY plot
    if(linePlot){
      plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
             stack = F, xlim = xlim, ylim = ylim, ...)			
    } 
    
    
    ## Generate stacked XY plot
    if(stackPlot){
      plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
             stack = T, xlim = xlim, ylim = ylim, ...)
      
    } 
    
    ## Generate violin plot
    if( vioPlot ){
      
      plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
             stack = F, vioPlot = T, xlim = xlim, ylim = ylim, ...)
      
    }
    
    ## Generate Z plot
    if(zPlot){
      plotXY(sLab = sLab, y = inDat[usrRow,], cpdNames = row.names(inDat)[usrRow],
             zPlot = T, xlim = xlim, ylim = ylim, ...)
      if( scale == "zScore" )
        abline( h=c(-1.644852, 1.644852), lty = 2)
      
    } 
    
    
    ## Generate heatmap
    if( heatMap ){
      lDat <- levelFun(inDat, n = n, zlim = zlim)
      metImage( as.matrix(inDat[nrow(inDat):1,]), col = lDat$col, breaks = lDat$breaks, 
                rCst = rCst, cCst = cCst, grid = grid )		
    }
    
    ## Pheatmap
    if (pheatmap){
      sLab <- sapply( strsplit(sam, split = "\\."), function (x) x[1])
      group_org <- data.frame(row.names = sam, SampleType = sLab)
      num_groups <- length(unique(sLab))
      set.seed(102)
      cc20 <- sample(rainbow(20, v = 0.8))
      cols <- cc20[1:num_groups]
      names(cols) <- unique(sLab)
      my_colour <- list(SampleType = cols)
      
      pheatmap(inDat,
               color = colorRampPalette(c("white", "#66CCFF", "#000033"))(100),
               cluster_cols = FALSE,
               annotation_col = group_org,
               annotation_colors = my_colour,
               cellwidth = 9,
               cellheight = 12)
    }		
    
    ## Return PCA scores
    if(returnPCA)
      return(outPCA)
    
    ## Add metadata
    if( returnMeta ){
      inDat <- cbind( metaDat, inDat)
    }else{
      inDat <- data.frame( compound = row.names(inDat), inDat, stringsAsFactors = F)
      names(inDat) <- c('compound', sam)
    }
    
    
    invisible(inDat)
  }
  
  
  
  
  ###############################################################################
  ## Internal functions to support mavPlot
  ###############################################################################
  
  
  foldChange <- function( inTable, normVec, zeroOffset = F, normMin = .01,
                          scale = 'fold' ){
    
    ## Internal functions
    myMin <- function( x ){ min( x, na.rm = TRUE ) * normMin }
    mySum <- function( x ){ sum(x, na.rm = TRUE) }
    mySd <- function( x ){ sd(x, na.rm = TRUE) }
    myMean <- function( x ){ mean(x, na.rm = TRUE) }
    myNaN <- function( x ){
      if( any(is.nan(x)) )
        x[ is.nan(x) ] <- 0
      return(x)
    } 
    
    ## Calculate standard deviation for normalization vector
    if( scale == "zScore"){
      if( length(normVec) < 1  )
        stop('Z scales require a standard deviation from multiple replicates')
      sdTab <- apply(inTable[, normVec], 1, mySd )
      normTab <- apply(inTable[, normVec], 1, myMean )
    }else{
      if( length(normVec) > 1 )
        normTab <- apply(inTable[, normVec], 1, myMean)
      else
        normTab <- inTable[, normVec ]		
    }
    
    
    
    ## Make norm matrix
    normTab <- as.vector(unlist(normTab))
    normTab <- data.frame(matrix(rep(normTab, ncol(inTable)), 
                                 ncol = ncol(inTable)))
    
    ## Make sd matrix
    if( scale == "zScore"){
      sdTab <- as.vector(unlist(sdTab))
      sdTab <- data.frame(matrix(rep(sdTab, ncol(inTable)), 
                                 ncol = ncol(inTable)))
    }
    
    ## Set all 0 values to 1% of min observed value
    if( zeroOffset ){
      zSet <- inTable == 0 
      if( mySum(zSet) ){
        tTab <- inTable 
        tTab[ zSet ] <- NA
        mSet <- as.vector(unlist(apply(tTab, 1, myMin )))
        for( i in 1:length(mSet) ){
          if( any(inTable[i,] == 0, na.rm =T) ){
            idx <- which( inTable[i,] == 0 )
            inTable[i, idx] <- mSet[i]				
          }
        }
      }
    }
    
    
    
    ## Calculate z score
    if( scale == 'zScore' ){
      out <- (inTable - normTab)/sdTab 
      #out <- apply( out, 2, myNaN )		
      return(out)		
    }
    
    
    
    ## Normalize
    out <- inTable /  normTab
    
    
    
    if( scale == 'ratio' ){
      out <- apply( out, 2, myNaN )		
      return(out)		
    }
    
    
    
    if( any(out < 1 ) ){
      foldFun <- function(x) { 
        tIdx <- which( x < 1)
        x[ tIdx ] <- abs(1/ x[tIdx]) * -1
        return(x)
      }
      out <- apply(out, 2, foldFun)	
      out <- apply( out, 2, myNaN )
    }
    
    return(out)
    
  }
  
  ## Function for plotting metabolite heatmaps
  metImage <- function( x, col, breaks, cCst = FALSE, rCst = FALSE, grid = TRUE ){
    
    ## Set col/row dendrogram
    Colv <- Rowv <- NA
    if( cCst )
      Colv <- NULL
    if( rCst )
      Rowv <- NULL
    
    w <- (3 + par('mar')[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par( mar = c( 5.1, 1, 4.1, 4.1), las = 1)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(breaks, na.rm = T), xaxs = "i", yaxs = "i")
    rect(0, breaks[-length(breaks)], 1, breaks[-1L], col = col)
    axis(4, at = breaks)
    
    
    #######################################	
    myHM(as.matrix(x), Colv = Colv, Rowv = Rowv, na.rm = T,
         col = col, breaks = breaks, grid = grid )
    
    
  }
  
  
  
  ## Make pca plots of maven data
  ## samples - vector of sample names
  ## compounds - vector of compound names
  pcaPlot <- function( dat, samples, compounds){	
    
    ## Prep datafame for pca
    dat <- data.frame(t(dat), stringsAsFactors = FALSE)
    
    names(dat) <- compounds
    row.names(dat) <- samples
    
    ## Generate the pca object
    pcaDat <- prcomp(dat)
    
    ## return the pca data invisibly
    invisible(pcaDat)
  }
  
  
  levelFun <- function(x, n = 10, r = c(0, 1), g = c(0,0), b = c(1,0), zlim, returnInterval = FALSE){
    
    x <- suppressWarnings(as.numeric(as.vector(unlist(x))))
    
    
    ## Remove NA and Inf data
    if( any(is.na(x)) )
      x[is.na(x)] <- 0
    if( any(abs(x) == Inf))
      x[ abs(x) == Inf ] <- 0
    
    ## Threshold data
    if( !missing(zlim) && length(zlim) == 2 && all(is.numeric(zlim)) ){
      zlim <- sort(abs(zlim) )
    }else
      zlim <- c(  min(abs(x)),  max(abs( x )) )
    
    ## Set breaks
    xSeq <- c(-rev(seq( zlim[1], zlim[2], length.out = n)), seq( zlim[1], zlim[2], length.out = n))
    
    
    ## Set colors	 
    col <- c(
      rgb(rep(r[1],(n-1)),
          rep(g[1],(n-1)), 
          rep(b[1],(n-1)), alpha = rev((1:n)[-1] / n)), 
      'transparent',
      rgb(rep(r[2],(n-1)),
          rep(g[2],(n-1)), 
          rep(b[2],(n-1)),  alpha = (1:n)[-1] / n) 		
    )		 
    
    if( !returnInterval )
      return( list(breaks = xSeq, col = col))
    return(list(breaks = xSeq, 
                col = col, 
                interval = findInterval(x, xSeq))) 
  }
  
  
  
  ## Internal plotting function for generating violin plots
  ## This function was modified from the R package vioplot
  myVioPlot <- function (y, sLab, range = 1.5, h = NULL, ylim = NULL, 
                         names = NULL, 
                         horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
                         lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                         at, add = FALSE, wex = 1, drawRect = TRUE, padj = 0) 
  {
    
    
    ## Load the required sm package
    suppressWarnings(msg <- require( sm) )
    if( !msg ){
      install.packages(sm)
      require(sm)
    }
    
    uLab <- unique(sLab)
    
    datas <- list()
    for( i in 1:length(uLab)){
      datas[[ i ]] <- as.numeric(y[ which(sLab == uLab[i]) ])		
    }
    
    names(datas) <- uLab	
    
    n <- length(datas)
    if (missing(at)) 
      at <- 1:n
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    if (!(is.null(h))) 
      args <- c(args, h = h)
    for (i in 1:n) {
      data <- datas[[i]]
      data.min <- min(data)
      data.max <- max(data)
      q1[i] <- quantile(data, 0.25)
      q3[i] <- quantile(data, 0.75)
      med[i] <- median(data)
      iqd <- q3[i] - q1[i]
      upper[i] <- min(q3[i] + range * iqd, data.max)
      lower[i] <- max(q1[i] - range * iqd, data.min)
      est.xlim <- c(min(lower[i], data.min), max(upper[i], 
                                                 data.max))
      smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
                                       args))
      hscale <- 0.4/max(smout$estimate) * wex
      base[[i]] <- smout$eval.points
      height[[i]] <- smout$estimate * hscale
      t <- range(base[[i]])
      baserange[1] <- min(baserange[1], t[1])
      baserange[2] <- max(baserange[2], t[2])
    }
    if (!add) {
      xlim <- if (n == 1) 
        at + c(-0.5, 0.5)
      else range(at) + min(diff(at))/2 * c(-1, 1)
      if (is.null(ylim)) {
        ylim <- baserange
      }
    }
    if (is.null(names)) {
      label <- names(datas)
    }
    else {
      label <- names
    }
    boxwidth <- 0.05 * wex
    if (!add) 
      plot.new()
    if (!horizontal) {
      if (!add) {
        plot.window(xlim = xlim, ylim = ylim)
        axis(2)
        axis(1, at = at, label = label, padj = padj, tick = F)
      }
      box()
      for (i in 1:n) {
        polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                c(base[[i]], rev(base[[i]])), col = col, border = border, 
                lty = lty, lwd = lwd)
        
        tcol = rep(gray(4/8),  length(datas[[i]]))
        if( any( datas[[i]] == 0 ) )
          tcol[datas[[i]] == 0] <- grey(7/8)
        
        points(x = at[i]+jitter(rep(0, length(datas[[i]])), 10), datas[[i]], 
               col = tcol, cex = .75, pch = 1)
        
        
        if (drawRect) {
          lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
                lty = lty)
          rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
               q3[i], col = rectCol)
          points(at[i], med[i], pch = pchMed, col = colMed)
          
        }		
        
      }
    }
    else {
      if (!add) {
        plot.window(xlim = ylim, ylim = xlim)
        axis(1)
        axis(2, at = at, label = label, padj = padj, tick = F)
      }
      box()
      for (i in 1:n) {
        polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                                rev(at[i] + height[[i]])), col = col, border = border, 
                lty = lty, lwd = lwd)
        if (drawRect) {
          lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
                lty = lty)
          rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
                 boxwidth/2, col = rectCol)
          points(med[i], at[i], pch = pchMed, col = colMed)
        }
      }
    }
    invisible(list(upper = upper, lower = lower, median = med, 
                   q1 = q1, q3 = q3))
  }
  
  
  
  ## Internal bridge function to R package sm
  sm.density <- function (x, h, model = "none", weights = NA, group = NA, ...) 
  {
    x.name <- deparse(substitute(x))
    data <- sm:::sm.check.data(x, NA, weights, group, ...)
    x <- data$x
    weights <- data$weights
    group <- data$group
    nobs <- data$nobs
    ndim <- data$ndim
    opt <- data$options
    if (ndim > 3) 
      stop("data with >3 dimensions are not allowed.")
    if (!all(is.na(group))) {
      if (!all(weights == 1) & opt$verbose > 0) 
        cat("Warning: weights ignored in sm.density.compare\n")
      return(sm:::sm.density.compare(x, group, h, model, ...))
    }
    sm:::replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs)/ndim))
    rawdata <- list(nbins = opt$nbins, x = x, nobs = nobs, ndim = ndim)
    if (opt$nbins > 0) {
      if (!all(weights == 1) & opt$verbose > 0) 
        cat("Warning: weights overwritten by binning\n")
      bins <- binning(x, nbins = opt$nbins)
      x <- bins$x
      weights <- bins$x.freq
      nx <- length(bins$x.freq)
      if (!all(is.na(opt$h.weights))) 
        stop("use of h.weights is incompatible with binning - set nbins = 0")
    }
    else nx <- nobs
    if (opt$positive) {
      sm:::replace.na(opt, delta, apply(as.matrix(x), 2, min))
      if ((ndim == 3) & (opt$verbose > 0)) 
        cat("the positive estimation is not available for 3 variables.\n")
    }
    if (missing(h)) {
      if (opt$positive) {
        xlog <- log(as.matrix(x) + outer(rep(1, nx), opt$delta))
        if (ndim == 1) 
          xlog <- as.vector(xlog)
        h <- h.select(xlog, y = NA, weights = weights, ...)
      }
      else h <- h.select(x = x, y = NA, weights = weights, 
                         ...)
    }
    
    if (is.na(opt$band)) {
      if (model == "none") 
        opt$band <- FALSE
      else opt$band <- TRUE
    }
    if ((model == "none") && opt$band) 
      opt$band <- FALSE
    if (ndim == 1) {
      if (length(h) != 1) 
        stop("length(h) != 1")
      sm:::replace.na(opt, xlab, x.name)
      sm:::replace.na(opt, ylab, "Probability density function")
      if (!opt$panel) 
        est <- sm:::sm.density.1d(x, h, model, weights, rawdata, 
                                  options = opt)
      else rp.density1(x, h, model, weights, rawdata, opt)
    }
    if (ndim == 2) {
      if (length(h) != 2) 
        stop("length(h) != 2")
      dimn <- dimnames(x)[[2]]
      name.comp <- if (length(dimn) == 2) 
        dimn
      else outer(x.name, c("[1]", "[2]"), paste, sep = "")
      sm:::replace.na(opt, xlab, name.comp[1])
      sm:::replace.na(opt, ylab, name.comp[2])
      if (!opt$panel) 
        est <- sm:::sm.density.2d(x, h, weights, rawdata, options = opt)
      else rp.density2(x, h, model, weights, rawdata, opt)
    }
    if (ndim == 3) {
      dimn <- dimnames(x)[[2]]
      name.comp <- if (length(dimn) == 3) 
        dimn
      else outer(x.name, c("[1]", "[2]", "[3]"), paste, sep = "")
      sm:::replace.na(opt, xlab, name.comp[1])
      sm:::replace.na(opt, ylab, name.comp[2])
      sm:::replace.na(opt, zlab, name.comp[3])
      opt$nbins <- 0
      if (!opt$panel) 
        est <- sm:::sm.density.3d(x, h = h, weights, rawdata, 
                                  options = opt)
      else rp.density3(x, h, model, weights, rawdata, opt)
    }
    if (!opt$panel) {
      est$data <- list(x = x, nbins = opt$nbins, freq = weights)
      est$call <- match.call()
      invisible(est)
    }
    else invisible()
  }
  
  
  myHM <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
                    reorderfun = function(d,w) reorder(d, w), add.expr, 
                    symm = FALSE, revC = identical(Colv, "Rowv"), na.rm = TRUE, 
                    grid = TRUE, breaks = NULL, margins = c(5, 5), 
                    ColSideColors, RowSideColors, cexRow = 0.2 + 
                      1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
                    verbose = getOption("verbose"), ...) 
  {
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
      stop("'x' must be a numeric matrix")
    
    x[ x == Inf ] <- max(x[is.finite(unlist(x))], na.rm = T) * 100 	
    x[ x == -Inf ] <- min(x[is.finite(unlist(x))], na.rm = T) / 100 
    
    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1) 
      stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L) 
      stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv")) 
      doCdend <- FALSE
    if (is.null(Rowv)) 
      Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
      Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
      if (inherits(Rowv, "dendrogram")) 
        ddr <- Rowv
      else {
        hcr <- hclust(dist(x))
        ddr <- as.dendrogram(hcr)
        if (!is.logical(Rowv) || Rowv) 
          ddr <- reorderfun(ddr, Rowv)
      }
      if (nr != length(rowInd <- order.dendrogram(ddr))) 
        stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
      if (inherits(Colv, "dendrogram")) 
        ddc <- Colv
      else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
          stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        ddc <- ddr
      }
      else {
        hcc <- hclust(dist(if (symm) 
          x
          else t(x)))
        ddc <- as.dendrogram(hcc)
        if (!is.logical(Colv) || Colv) 
          ddc <- reorderfun(ddc, Colv)
      }
      if (nc != length(colInd <- order.dendrogram(ddc))) 
        stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
      if (is.null(rownames(x))) 
        (1L:nr)[rowInd]
    else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
      if (is.null(colnames(x))) 
        (1L:nc)[colInd]
    else colnames(x)
    else labCol[colInd]
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
              4)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) != 
          nc) 
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1L], 0.2, lhei[2L])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) != 
          nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                                     1), lmat[, 2] + 1)
      lwid <- c(lwid[1L], 0.2, lwid[2L])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
      cat("layout: widths = ", lwid, ", heights = ", lhei, 
          "; lmat=\n")
      print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
      par(mar = c(margins[1L], 0, 0, 0.5))
      image(rbind(1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2L]))
      image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1L], 0, 0, margins[2L]))
    x <- t(x)
    if (revC) {
      iy <- nr:1
      if (doRdend) 
        ddr <- rev(ddr)
      x <- x[, iy]
    }
    else iy <- 1L:nr
    
    
    if( !is.null(breaks) && length(breaks) > 1 ){
      x[ x > max(breaks) ] <- max(breaks)
      x[ x < min(breaks) ] <- min(breaks)
    }
    
    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
            c(0, nr), axes = FALSE, xlab = "", ylab = "", breaks = breaks, ...)
    axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
         cex.axis = cexCol, col.axis = par('fg'))
    box()
    if( grid )
      abline(v = 1L:nc + .5, h = iy + .5)
    if (!is.null(xlab)) 
      mtext(xlab, side = 1, line = margins[1L] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
         cex.axis = cexRow, col.axis = par('fg'))
    if (!is.null(ylab)) 
      mtext(ylab, side = 4, line = margins[2L] - 1.25)
    if (!missing(add.expr)) 
      eval(substitute(add.expr))
    par(mar = c(margins[1L], 0, 0, 0))
    if (doRdend) 
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
    if (doCdend) 
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
      frame()
    if (!is.null(main)) {
      par(xpd = NA)
      title(main, cex.main = 1.5 * op[["cex.main"]])
    }
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                                doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
  }
  
  
  pepFinder <- function( pSeq, frag = 2:10, charge = -1 ){
    pSeq <- unlist(strsplit(pSeq, ''))
    pRes <- 1:length(pSeq)
    
    out <- outIdx <- NULL
    for( i in frag ){
      
      for( j in 1:length(pSeq) ){
        if( any (is.na(pSeq[ j:(j+i-1) ])))
          next
        out <- c(out, paste(pSeq[ j:(j+i-1) ], collapse = ''))
        outIdx <- c(outIdx, range(pRes[ j:(j+i-1) ]) )
      }	
    }
    outIdx <- data.frame(matrix(outIdx, ncol = 2, byrow = T))
    names(outIdx) <- c('start', 'end')
    return(cbind(outIdx, aaAdd(out, charge = charge)))
  }
  
  
  
  fragFinder <- function( pSeq, charge = -1 ){
    
    frag = 1:nchar(pSeq)
    out <- pepFinder( pSeq = pSeq, frag = frag, charge = charge)
    out <- unique(out[,match(c('compound', 'mz', 'nRes'), names(out))])
    out2 <- out
    out2$compound <- paste(out2$compound, '-H2O', sep = '')
    out2$mz <- out2$mz - mzCalc('H2O', charge = 0)$mz
    out <- rbind(out, out2)
    
    out2 <- out
    out2$compound <- paste(out2$compound, '-CO2', sep = '')
    out2$mz <- out2$mz - mzCalc('CO2', charge = 0)$mz
    out <- rbind(out, out2)
    
    return(out[order(out$mz),])
    
    
  }
  
  aaAdd <- function( inPep, charge = 0 ){
    
    aa <- data.frame(
      Name = c('Alanine', 'Arginine', 'Asparagine', 'AsparticAcid', 
               'Cystine', 'GlutamicAcid', 'Glutamine', 'Glycine', 'Histidine', 
               'Isoleucine', 'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 
               'Proline', 'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'), 
      Code = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 
               'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
      C = c(3, 6, 4, 4, 6, 5, 5, 2, 6, 6, 6, 6, 5, 9, 5, 3, 4, 11, 9, 5),
      H = c(7, 14, 8, 7, 12, 9, 10, 5, 9, 13, 13, 14, 11, 11, 9, 7, 9, 12, 11, 11),
      N = c(1, 4, 2, 1, 2, 1, 2, 1, 3, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1),
      O = c(2, 2, 3, 4, 4, 4, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2),
      S = c(0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
      stringsAsFactors = FALSE
    )
    
    eMass <- list( C = 12.000000, H = 1.007825, N = 14.003074, O = 15.994915, 
                   S = 31.972072, P = 30.97376163, I = 126.904473)
    
    
    ## Sort peptide list
    aa <- aa[order(aa$Code),]
    
    if( missing (inPep) || length(inPep) < 1 ){
      diPep <- triPep <- tetPep <- NULL
      
      ## Find unique dipeptides
      diPep <- list()
      for( i in 1:20)
        diPep[[i]] <- paste(aa$Code[i], aa$Code, sep = '')
      diPep <- unique(sort(sapply(strsplit(unlist(diPep),''), 
                                  function(x) paste(sort(x), collapse = ''))))
      ## Find unique tripeptides
      triPep <- list()
      for( i in 1:20)
        triPep[[i]] <- paste(aa$Code[i], diPep, sep = '')	
      triPep <- unique(sort(sapply(strsplit(unlist(triPep),''), 
                                   function(x) paste(sort(x), collapse = ''))))
      
      ## Find unique tetra
      tetraPep <- list()
      for( i in 1:20)
        tetraPep[[i]] <- paste(aa$Code[i], triPep, sep = '')
      tetraPep <- unique(sort(sapply(strsplit(unlist(tetraPep),''), 
                                     function(x) paste(sort(x), collapse = ''))))
      
      ## Find unique penta
      pentaPep <- list()
      for( i in 1:20)
        pentaPep[[i]] <- paste(aa$Code[i], tetraPep, sep = '')
      pentaPep <- unique(sort(sapply(strsplit(unlist(pentaPep),''), 
                                     function(x) paste(sort(x), collapse = ''))))
      
      
      pepList <- sort(c(diPep, triPep, tetraPep, pentaPep))
    }else
      pepList <- inPep
    
    ## Calculate Formulae
    tNames <- pepList
    if( length(inPep) > 1 ){
      pepList <- sapply(sapply(pepList, function(x) strsplit(x, '')), 
                        function(x) match(x, aa$Code))				
      if( is.matrix(pepList) ){
        pepList <- as.list(data.frame(pepList))
        names(pepList) <- tNames			
      }
    }else{
      pepList <- list(as.vector(sapply(strsplit(pepList, ''), 
                                       function(x) match(x, aa$Code))))
      names(pepList) <- tNames
    }
    
    outList <- sapply(pepList, function(x){
      paste('C', sum(aa[x,]$C), 'H', sum(aa[x,]$H) - (2 * (length(x) - 1) ),
            'N', sum(aa[x,]$N), 'O', sum(aa[x,]$O) - (length(x) - 1), 
            'S', sum(aa[x,]$S),sep = '', collapse = '' )
    })
    
    
    ## Calculate mass to charge 
    mzDat <- mzCalc(as.vector(outList), charge = charge)
    outList <- data.frame( compound = names(pepList), 
                           formula = mzDat$formula,
                           mz = mzDat$mz, nRes = sapply(strsplit(names(pepList), ''), length),
                           z = mzDat$z, 
                           stringsAsFactors = FALSE)
    return(outList)
    
  }
  
  findPval <- function( sLab, mzData, type = 'sd'  ){
    
    nr <- nrow( mzData )
    pVal <- list()
    sLab <- factor(sLab)
    
    for( i in 1:nr ){	
      
      tVal <- tryCatch( pVal[[i]] <- aov(as.vector(unlist(mzData[ i, ])) ~ sLab), 
                        error=function(er) return(NA))
      
      if( is.na(tVal[1]) ){
        pVal[[i]] <- NA
        next
      }
      
      if(type == 'p')
        pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$Pr[1]
      
      if(type =='f')
        pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$"F value"[1]
    }
    
    
    #	if(type == 'p')
    #		return(unlist(
    #						lapply(pVal, function(x) unclass(summary.aov(x))[[1]]$Pr[1])))
    #	if(type =='f')
    #		return(unlist(
    #						lapply(pVal, function(x) unclass(summary.aov(x))[[1]]$"F value"[1])))
    
    return(as.vector(unlist(pVal)))
  }
  
  compProt <- function( fasta ){
    
    ## Find protein headers
    sidx <- which(substring(fasta, 1, 1) == ">")
    eidx <- c(sidx[-1] -1, length(fasta))
    
    out <- list()
    for( i in 1:length(sidx) ){
      out[[i]] <- paste(fasta[(sidx[i] + 1):eidx[i]], collapse = '')
    }
    
    ## Find unique genes
    gene <- sapply(strsplit(fasta[sidx], 'ENSG'), function(x) substr(x[2], 1, 11))
    uGene <- unique(gene)
    pOut <- list()
    
    for( i in 1:length(uGene) ){
      pOut[[i]] <- paste( out[which(gene == uGene[i])], collapse = '#' )
    }
    
    return(list(uGene, pOut))
  }
  
  matchSeq <- function( proList, inPep, preCal, hbIdx, thresh ){
    
    midx <- nchar(inPep)
    out <- rep(FALSE, length(proList))
    
    ## Find hb matches first
    if( !missing(hbIdx) ){
      for( i in hbIdx ){
        tPro <- proList[[i]]
        for( j in 1:nchar(tPro) ){
          if(substr(tPro, j, j+midx-1) == inPep){
            out[i] <- TRUE
            break
          }
        }
      }
    }
    outHb <- out
    
    out <- rep(FALSE, length(proList))
    for( i in 1:length(proList) )	{
      
      if( !missing (preCal) ){
        if( i < preCal )
          next
        if( i == preCal ){
          out[i] <- TRUE
          next
        }			
      }
      
      if( !missing(thresh) ){
        if( sum(out) >= thresh )
          break
      }
      
      
      tPro <- proList[[i]]
      for( j in 1:nchar(tPro) ){
        if(substr(tPro, j, j+midx-1) == inPep){
          out[i] <- TRUE
          break
        }
      }
    }
    
    out[which(outHb)] <- TRUE
    return(out)
    
  }
  
  ## Extracts the ms1 mz and RT from a mascot .mgf file
  mascotRead <- function( inPath, formula, charge ){
    
    ## Read mascot .mgf peak file
    print('Reading file...'); flush.console()
    if(missing(inPath))
      mzDat <- readLines(file.choose())
    else
      mzDat <- readLines( inPath )
    
    ## Parse file
    print('Parsing file...'); flush.console()
    idx1 <- grep("BEGIN IONS", mzDat)
    idx2 <- c(idx1[-1] -3, length(mzDat)-2)
    out <- mzDat[sort(c(idx1 + 2, idx1+4))]
    out <- sapply(strsplit(
      sapply(strsplit(out, '='), function(x) x[2]), " "), function(x) x[1])
    out <- (matrix(as.numeric(out), ncol = 2, byrow = T))
    out[,2] <- out[,2]/60
    z <- as.numeric(sapply(strsplit(sapply(strsplit(mzDat[idx1+3], 'CHARGE='), 
                                           function(x) x[2]), '+'), function(x) x[1]))
    specScan <- matrix(as.numeric(unlist(
      strsplit(sub(',', '', 
                   sub(' scans: ', ' ', 
                       sub('TITLE=Spectrum', '', 
                           mzDat[idx1+1]))), ' '))), ncol = 2, byrow = T)
    
    ## Generate ms1 list
    out <- data.frame( mz = out[,1], 
                       rt = out[,2], 
                       z = z, 
                       spectrum = specScan[,1],
                       scan = specScan[,2] )
    
    if( missing( formula) || missing(charge) )
      return(out)
    
    ## Calculate expected observed masses
    print('Calculating masses...'); flush.console()
    mzList <- NULL
    for( i in unique(charge) ){
      mzList <- rbind(mzList, 
                      data.frame( formula = unique(formula), z = i, mz = mzCalc(formula, i)$mz, 
                                  stringsAsFactors = F ))
    }
    
    ## Match MS1 signals with expected masses	
    print('Matching masses...'); flush.console()
    ms1Match <- NULL
    for( i in unique(mzList$z) ){
      sOut <- unique(out[out$z == i,]$mz)
      sKno <- unique(mzList[mzList$z == i,]$mz)
      ms1Match <- sort(unique(c(ms1Match, 
                                sOut[mzMatch( mzObs = sOut, mzKno = sKno)$mzObs])))
    }
    
    ## Filter data for matching input formula and find ms2 data
    idx3 <- which(!is.na(match(out$mz, ms1Match)))
    out <- out[idx3,]
    out$ms2 <- NA
    for( i in 1:length(idx3) )
      out[i,]$ms2 <- paste( mzDat[(idx1[idx3[i]]+6):idx2[idx3[i]]], collapse = ' ')		
    
    return(out)		
    
    
    ## This code can be used for ms2 matching, but it is slow		
    #	if(!missing(inPep)){
    #		zKno <- data.frame(		z1 = mzCalc( inPep$formula, 1)$mz,
    #													z2 = mzCalc( inPep$formula, 2)$mz,
    #													z3 = mzCalc( inPep$formula, 3)$mz,
    #													z4 = mzCalc( inPep$formula, 4)$mz )
    #		out1 <- out[out$z == 1,][mzMatch( out[out$z == 1,]$mz, zKno$z1)$mzObs,]
    #		out2 <- out[out$z == 2,][mzMatch( out[out$z == 2,]$mz, zKno$z2)$mzObs,]
    #		out3 <- out[out$z == 3,][mzMatch( out[out$z == 3,]$mz, zKno$z3)$mzObs,]
    #		out4 <- out[out$z == 4,][mzMatch( out[out$z == 4,]$mz, zKno$z4)$mzObs,]
    #	}
    #	
    #	## Find ms2 signals that are in the top 2/3 of max signal
    #	out$ms2 <- rep(NA, nrow(out))
    # for( i in 1:nrow(out) ){
    #		
    #		if( i < nrow(out) )
    #			j <- i + 1
    #		else
    #			j <- length(mzDat) + 1
    #		
    #		tmp <- matrix(as.numeric(unlist(
    #								(strsplit(mzDat[(idx1[i]+6) : (idx1[j]-3)], ' ')))), 
    #				         ncol = 2, byrow = T)
    #		out[i,]$ms2 <- paste(tmp[,1][which(tmp[,2] > max(tmp[,2], na.rm = T) / 3)], 
    #				collapse = ' ')
    #	}
    
    return( out )
  }
  
  
  ms2Match <- function( ms2List, inPep, ms2Charge = 1, ms2Res = 500, intThr = .25){
    
    ## Error checking
    if(any(is.na(match(c('sample', 'mz', 'rt', 'z', 'ms2'), names(obsMz)))))
      stop('ms2List must have: sample, mz, rt, z, and ms2 column names')
    if(any(is.na(match(c('formula', 'compound'), names(inPep)))))
      stop('inPep must have: formula and compound column names')	
    
    ## Calculate masses for all observed charge states
    print('Calculating masses...'); flush.console()
    inPep <- unique(data.frame(compound = inPep$compound, 
                               formula = inPep$formula, z = NA, stringsAsFactors = F))	
    outPep <- NULL
    for( i in sort(unique(ms2List$z)) ){
      outPep <- rbind(outPep, 
                      data.frame( 
                        compound = inPep$compound,
                        formula = inPep$formula, 
                        z = i, 
                        mz = mzCalc(inPep$formula, i)$mz, 
                        stringsAsFactors = F ))
    }
    
    ms2Parse <- function(x){ 
      out <- data.frame(matrix(as.numeric(unlist(strsplit(x, ' '))), 
                               ncol = 2, byrow = T), stringsAsFactors = F)
      names(out) <- c('mz', 'int' )
      return(out)
    }
    
    ## Match ms2 fragments - this is basic matching
    print('Matching fragments...'); flush.console()
    H2O <- mzCalc( 'H2O', charge = 0)$mz
    ms1 <- ms2 <- putId <- tPep <- tFrag <- tList <- NULL
    outList <- list()
    for( i in 1:nrow(ms2List) ){
      ms1 <- ms2List[i,]
      ms2 <- ms2Parse( ms1$ms2 )
      ms1$ms2I <- max(ms2$int)
      ms1$ms2 <- NULL
      tPep <- outPep[outPep$z == ms1$z,]
      putId <- tPep[mzMatch( mzObs = ms1$mz, mzKno = tPep$mz)$mzKno,]
      subList <- list()
      subList[[1]] <- ms1
      if( !nrow(putId) ){
        subList[[2]] <- NA
        outList[[i]] <- subList
        next
      }
      
      ## Match expected fragments based MS1 IDs		
      for( j in 1:nrow(putId) ){
        
        ## create synthetic fragment list
        tFrag <- tFragDH <- pepFinder(putId$compound[j], 
                                      frag = 1:10, charge = ms2Charge)
        tFragDH$compound <- paste( tFragDH$compound, 'H2O', sep = '-')
        tFragDH$mz <- tFragDH$mz - H2O
        tFragDH$formula <- NA
        tFrag <- rbind(tFrag, tFragDH)
        
        ## match the fragments to observed ms2
        tList <- mzMatch( mzObs = ms2$mz, mzKno = tFrag$mz, ppm = ms2Res)
        if( !length(tList$mzKno) ){
          subList[[j + 1]] <- NA
          next				
        }
        
        tList$mzObs <- ms2[tList$mzObs,]
        tList$mzKno <- tFrag[tList$mzKno,]
        subList[[j + 1]] <- tList
      }
      
      names(subList) <- c('ms1', putId$compound)
      outList[[i]] <- subList
      
      if( (i %% 100) == 0 ){
        print( i )
        flush.console()
      }
    }
    
    
    ## Score each ms2 frag spectrum - this returns pairwise matches
    print( 'Scoring fragmentation spectra' )
    flush.console()	
    for( i in 1:length(outList) ){
      maxInt <- outList[[i]][[1]]$ms2I 
      
      ## Filter mass list
      for( j in length(outList[[i]]):2){
        ## Make sure that max matched frag is within x% of max observed frag
        if( is.na(outList[[i]][[j]]) || 
            max(outList[[i]][[j]]$mzObs$int / maxInt) < intThr){
          outList[[i]][[j]] <- NULL
          next
        }
        
        ## Find ppm match for each fragment
        obs <- outList[[i]][[j]]$mzObs$mz
        kno <- outList[[i]][[j]]$mzKno$mz
        ppm <- mzMatch(mzObs = obs, mzKno = kno, ms2Res, probMat = T)
        frag <- list()
        for( k in 1:length(ppm) ){
          frag[[k]] <- outList[[i]][[j]]$mzKno$compound[ppm[[k]][
            which.min(abs(obs[[k]] - kno[ppm[[k]]]) / kno[ppm[[k]]])]]
          ppm[[k]] <- min(abs(obs[[k]] - kno[ppm[[k]]]) / kno[ppm[[k]]] * 10^6)
        }
        
        outList[[i]][[j]]$mzObs$ppm <- unlist(ppm)
        outList[[i]][[j]]$mzObs$frag <- unlist(frag)
      }	
    }
    outList <- outList[ sapply(outList, length) > 1 ]
    
    return(outList)	
  }
  
  #scoreFrag <- function( inList, cpdList  ){
  #	ms1 <- inList[[1]]
  #	ms2 <- inList[-1]
  #	N <- length(ms2)
  #	cpd <- names(ms2)
  #	
  #	cpd <- cbind(cpd, 
  #			mzCalc(cpdList[match(cpd, Hb$compound),]$formula, charge = ms1$z), 
  #		 stringsAsFactors = F)
  #	cpd$ppm <- abs(ms1$mz - cpd$mz )/cpd$mz *10^6
  #	
  #	for( i in 1:length(N) ){
  #	x <- abs(ms2[[i]]$mzObs$ppm - 500)/500
  #	y <- ms2[[i]]$mzObs$int /ms1$ms2I		
  #
  #
  #	}
  #	
  #}
  
  ## Plots a list of peptides relative to their parent protein
  seqPlot <- function( parentSeq, inList, add = F, cex = .5, vScale = 1, y, group, 
                       col, zCol = FALSE, main = '', xlab = 'Sequence', ylab = '', plotHash = F, 
                       n, zlim, xlim = NULL, ylim = NULL, ... ){
    
    parentSeq <- strsplit(parentSeq, '')[[1]]
    if( !add ){
      
      if( missing(y) || zCol)
        yOut <- c(0, nrow(inList))
      else
        yOut <- range(y)
      if(yOut[1] > 0 )
        yOut[1] <- 0
      
      if( missing(group) ){
        plot( c(0, length(parentSeq)), y = yOut, type = 'n', 
              xlab = xlab, ylab = ylab, main = main )
        text(x = 1:length(parentSeq), y = 0, labels = parentSeq, cex = cex  )
      }else{
        plot( c(0, length(parentSeq)), y = max(summary(factor(group))) * c(-1,1), 
              type = 'n', xlab = xlab, ylab = ylab, main = main ) #xlim = c(0,150) )
        #if( !plotHash )
        text(x = 1:length(parentSeq), y = 0, labels = parentSeq, cex = cex  )
      }
    }
    
    if( missing(col) )
      col <- 'black'
    if( missing(zlim) )
      zlim <- NULL
    if( missing(n) )
      n <- NULL
    
    ## Plot scale bar
    if( zCol && !missing(y) ){
      levelDat <- levelFun( y, returnInterval = T, n = n, zlim = zlim)
      zIdx <- which(levelDat$col == "transparent" )	
      
      ## Remove next line of code for final verasion
      levelDat$col[zIdx:(zIdx+4)] <- c('black', "tan", 'orange', 'darkorange1', 'red')
      
      levelDat$col[zIdx] <- 'black'
      col <- levelDat$col[ levelDat$interval ]
      
      
      lX <- .85 * par()$usr[2]
      lY <- .95 * par()$usr[4]
      for( i in 1:length(levelDat$col) ){
        tY <- lY - (.05*lY*i)
        if( i > zIdx)
          j <- i+1 
        else
          j <- i
        
        #### Remove next line of code for final version
        if( i != j )
          next
        
        if( i == zIdx ){
          points( lX, tY, col = 'black', pch = 15, cex = 2.5)
          text( lX + (.05*lX), tY, 0)
          next
        }
        
        points( lX, tY, col = rev(levelDat$col)[i], pch = 15, cex = 2.5)
        text( lX + (.05*lX), tY, rev(levelDat$breaks)[j])
      }
    }
    
    rSum <- function(x, idx ){
      out <- x
      N <- length(x)
      for( i in 1:length(x) ){
        out[i] <- sum(x[i:N])
      }
      return(out[idx])
    }
    
    
    for ( i in 1:nrow( inList) ){
      
      ## Set y axis
      yOut <- i * vScale
      if( !missing(y) && !zCol )
        yOut <- y[i]
      if( !missing(group) ){
        tGrp <- group[i]
        if( tGrp %% 2 )
          yOut <- rSum(group == tGrp, i)
        else
          yOut <- -rSum(group == tGrp, i)
      }
      
      
      ## Set v scale
      #		if( missing(y) || zCol )
      #			yOut <- i * vScale
      #		else if( missing(group) )
      #			yOut <- y[i]
      #		else{
      #			if( group[i] %% 2 )
      #				yOut <- i 
      #			else
      #				yOut <- i * -1
      
      ##############work here to do #################		
      #}
      
      ## Set color
      if( length(col) != nrow(inList) )
        colOut <- col[1]
      else
        colOut <- col[i]
      
      xRange <- inList[i,]$start : inList[i,]$end
      
      if( !plotHash )
        text(x = xRange, y = yOut, labels = parentSeq[xRange], 
             col = colOut, cex = cex, ... )
      else{
        points(x = xRange, y = rep(yOut, length(xRange)), pch = 15, 
               col = colOut, cex = cex, ... )			
      }
      
    }	
  }
  
  
  ## General function for finding the intersection of two mass lists
  ## mzObs - numeric vector - input observed masses
  ## mzKno - numeric vecor - known masses to be matched
  ## probMat - logical - If TRUE will return the all of the mzKno masses matching 
  ##                     each of the mzObs entries
  ##                   - if FALSE will return a list of both obsrved and known
  ##                     masses with matching signals
  mzMatch <- function( mzObs, mzKno, ppm = 10, probMat = FALSE ){
    
    ## Find mz diff matrix 
    mzMat <- matrix(NA, nrow = length(mzObs), ncol = length(mzKno) )
    for( i in 1:length(mzKno) ){
      mzMat[,i] <- (abs(mzObs - mzKno[i]) / mzObs * 10^6 <= ppm )
    }
    
    if( probMat ){
      out <- apply( mzMat, 1, which) 
      if( is.list(out) )
        return(out)
      if( length(mzObs) == 1 )	
        return( list(unlist( out )))
      if( is.vector(out) )
        return( as.list(out) )			
      
      outAlt <- list()
      for( i in 1:ncol(out))
        outAlt[[i]] <- out[,i]
      return(outAlt)
    }
    
    
    
    mzMat <- list( mzObs = which(apply(mzMat, 1, any)), 
                   mzKno = which(apply(mzMat, 2, any)))
    
    return(mzMat)
  }
  
  
  
  ## mzList - maven format input list
  ## name - char - name used to designate blank
  ## thresh - numeric - minimum threshold for matching
  blankFilter <- function( mzList, name = 'blank', thresh = .1 ){
    
    ## Remove all standard Maven meta data
    mavNam <- c('label', 'metaGroupId', 'groupId', 'goodPeakCount', 'medMz', 'medRt', 
                'maxQuality', 'note', 'compound', 'compoundId', 'expectedRtDiff', 
                'ppmDiff', 'parent')	
    
    ## Find the blank to sample ratio
    tDat <- mzList[,-na.omit(match(mavNam, names(mzList)))] ## remove the non numeric data
    if( length(grep(name, names(tDat))) < 1 )
      stop('Sample name not found')
    if( length(grep(name, names(tDat))) == 1 )
      blk <- tDat[,grep(name, names(tDat))]
    else
      blk <- as.vector(apply(tDat[,grep(name, names(tDat))], 1, mean))
    datMax <- as.vector(apply(tDat[,-grep(name, names(tDat))], 1, mean))
    return( mzList[!((blk/datMax ) > thresh),] )
    
  }
  
  isotopeCheck <- function( mzList, mz, ppm = 10 ){
    
    
  }
  
  
  ## Function for clustering untargeted mz data
  ## mzTable - dataframe - Maven format dataframe from untargeted csv maven output
  ## mzThresh - numeric - Mass tolerance in ppm
  ## rtThresh - numeric - Renention time tolerance in minutes
  ## pW - numeric vector, 0 to 1 - weight given to mass, 
  ##                      retention time, and correlation
  ## pThresh - numeric, 0 to 1 - overall probablility threshold
  ## adList - numeric vector - list of known adducts
  ## collapse - logical - T collapses groups to highest intensity signal
  ## clusterOnly - logical - T returns collapsed groups with multiple signals 
  mzDiff <- function( 
    mzTable, 
    mzThresh = 10, 
    rtThresh = .5,
    pW = c(1,1,1),
    pThresh = 0.75,
    adList = c(0, 1.003355, 1.007825, 22.9897692809, 18.010565, 43.99038), #, 34.968853, 6.016936)
    collapse = FALSE, clusterOnly = FALSE
  ){
    
    mzList <- mzTable$medMz
    rtList <- mzTable$medRt
    pW <- pW/sum(pW) 
    
    
    ## Remove all standard Maven meta data
    mavNam <- c('label', 'metaGroupId', 'groupId', 'goodPeakCount', 'medMz', 'medRt', 
                'maxQuality', 'note', 'compound', 'compoundId', 'expectedRtDiff', 
                'ppmDiff', 'parent', 'pVal')	
    mzData <- mzTable[,-na.omit(match(mavNam, names(mzTable)))] ## remove the non numeric data
    
    ## Find mz diff matrix 
    mzMat <- matrix(NA, ncol = length(mzList), nrow = length(mzList) )
    for( i in 1:length(mzList) ){
      mzMat[i,] <- abs(mzList[i] - mzList)
    }
    
    ## Calculate mass prob matrix
    myMin <- function( x ){ return( min(x, na.rm = TRUE))}
    mzProb <- matrix(NA, ncol = ncol(mzMat), nrow = nrow(mzMat) )	
    for( i in adList){
      mzTmp <- abs(mzMat - i) 
      mzTmp <- (mzTmp / mzList) * 10^6
      
      ## Keep the best mz match
      mzIdx <- which(mzTmp < mzThresh) 
      tmpDat <- cbind(mzTmp[mzIdx], mzProb[mzIdx])
      mzProb[mzIdx] <- apply( tmpDat, 1, myMin )	
    }
    
    mzProb <- (abs(mzProb - mzThresh)/mzThresh)
    mzProb[which(is.na(mzProb))] <- 0
    
    
    ## Find rt diff matix
    rtMat <- matrix(NA, ncol = length(rtList), nrow = length(rtList) )
    for( i in 1:length(rtList) ){
      rtMat[i,] <- abs(rtList[i] - rtList)
    }
    
    ## Create linear probability matix with a threshold
    rtMat[rtMat > rtThresh ] <- NA
    rtMat <- abs((rtMat - rtThresh)/rtThresh)
    rtMat[which(is.na(rtMat))] <- 0
    
    
    ## Find covarance
    corMat <- abs(cor( t(mzData) ))
    
    ## Make weighted prob matrix
    probMat <- (mzProb * pW[1]) + (rtMat * pW[2]) + (corMat * pW[3])
    
    
    ## Find compound grouping
    group <- apply(probMat, 1, function(x) which(x > pThresh))
    oLen <- sapply(group, length)
    nLen <- 1
    while( any(oLen != nLen) ){
      oLen <- sapply(group, length)
      for( i in 1:length(group) ){
        group[[i]] <- sort(unique(unlist(
          group[sapply( group, 
                        function(x) any(match(x, group[[i]]), na.rm = T) )])))
      }
      nLen <- sapply(group, length)
    }
    
    
    group <- unique( group )
    group <- unique(group[order(sapply(group, length), decreasing = T)])
    grp <- nGrp <-  rep( NA, length(unlist(group)) )
    for( i in 1:length(group) ){
      grp[ group[[i]] ] <- i		
      nGrp[ group[[i]] ]  <- length(group[[i]])
    }
    
    ## Reorder list
    mzTable <- cbind( grp, nGrp, gPar = NA, mzTable)[order(grp),]
    grp <- grp[order(grp)]
    
    ## Find parent isotope for each group
    uGrp <- unique(grp)
    datMax <- as.vector(apply(mzTable[,grep('_', names(mzTable))], 1, sum))
    for( i in 1:length(uGrp) ){
      idx <- which( grp == uGrp[i] )
      if( any(mzTable[idx,]$nGrp == 1)) 
        next
      idx1 <- idx[ which.max(datMax[idx]) ]
      idx2 <- idx[ which.min(mzTable$medMz[idx]) ]
      
      if( idx1 == idx2 )
        mzTable[idx1,]$gPar <- 'Parent'
      else{
        mzTable[idx1,]$gPar <- 'MaxI'
      }
    }
    
    mzTable <- unique(mzTable)
    
    ## Collapse groups to parent masses
    if( collapse ){
      sTab <- mzTable[ mzTable$nGrp > 1, ]
      sTab <- sTab[ !is.na(sTab$gPar), ]
      if( clusterOnly )
        mzTable <- sTab
      else
        mzTable <- rbind( sTab, mzTable[ -(which(mzTable$nGrp > 1)), ])		
      mzTable$nGrp <- 1			
    }
    
    return(mzTable)
  }
  
  ## Calculates exact mass from a formula
  ## formula - character vector - must use standard atomic formula nomenclature
  ## charge - integer - subtracts or adds protons and divides by specified charge
  ## eMass - numeric list - exact masses and atom names
  mzCalc <- function( formula, charge = -1, 
                      eMass = list( C = 12.000000, H = 1.007825, N = 14.003074, O = 15.994915, 
                                    S = 31.972072, Na = 22.9897692809, K = 38.96370668, Fe = 55.9349375, 
                                    P = 30.97376163, I = 126.904473, Cl = 34.96885268)
  ){
    
    ## Generate list of all input atoms in formula list
    alpha <- c(LETTERS, letters)
    atomList <- match(strsplit(paste(formula, collapse = ''), '')[[1]], alpha)
    if( length(which(is.na(atomList))) )
      atomList <- as.numeric(atomList[-which(is.na(atomList))])
    tlAtom <- NULL
    if( any(atomList > 26 )){
      lc <- which( atomList > 26 )
      tlAtom <- paste( alpha[atomList[ lc - 1 ]], alpha[atomList[ lc ]], 
                       sep = '')
      atomList <- atomList[ -c(lc - 1, lc) ]
    }
    atomList <- unique(c(alpha[atomList], tlAtom))
    
    ## Check to make sure that input atoms are in the eMass list
    if( any(is.na(match(atomList, names(eMass))))){
      stop(c('The following atoms are not in the lookup table: \n', 
             paste(atomList[ is.na(match(atomList, names(eMass)))], ' ')))		
    }
    
    ## Sub formula for parsing formula
    parseForm <- 	function(x, atom = 'C'){ 
      
      if( is.na(match( atom, names(eMass)) ))
        stop('Atom is not in lookup table')
      
      if( !length(grep(atom, x)) )
        return(0)
      
      ## Set all other atoms to 'z'
      falseMatch <- names(eMass)[- match( atom, names(eMass) )]
      for( i in falseMatch ){
        x <- sub( i, 'z', x)
      }
      
      ## Find atom
      subForm <- unlist( strsplit(x, atom) ) 
      if( length( subForm ) == 1 )
        return(1)
      atomNum <- suppressWarnings( as.numeric( strsplit(subForm[2], 'z')[[1]][1]))
      if( is.na(atomNum) )
        return(1)
      
      return(as.numeric(strsplit(subForm[2], 'z')[[1]][1]))
    }
    
    ## Parse all incoming formulas into an atom table
    atomList <- atomList[ rev(order(nchar(atomList)))]
    atomTable <- data.frame( matrix(NA, ncol = length(atomList), 
                                    nrow = length(formula)))
    names(atomTable) <- atomList
    for( i in 1:length(atomList) ){
      atomTable[,i] <- as.vector(sapply(formula, 
                                        function(x, y = atomList[i]) parseForm(x, y)))
    }
    
    atomTable <- atomTable[, order(match( names(atomTable), names(eMass)))]
    formulaList <- apply( atomTable, 1, 
                          function(x){ 
                            paste(names(x)[ x > 0], x[x > 0], sep = '', collapse = '')
                          } ) 
    
    ## Make atomic mass table
    mzTable <- data.frame( matrix(0, ncol = length(atomList), 
                                  nrow = length(formula)))
    for( i in 1:ncol(atomTable) )
      mzTable[,i] <- atomTable[,i] * 
      eMass[[ which(names(eMass) == names(atomTable)[i] )]]
    
    ## Calculate mass
    mzTable <- apply( mzTable, 1, sum )
    if( length(charge) == length(mzTable) ){
      idx <- which(charge != 0 )
      if( length(idx) ){
        mzTable[idx] <- abs((mzTable[idx] + 
                               (charge[idx] * (eMass$H - .00054857990946)))/charge[idx]) # mass of proton
      }
    }else{
      if( length(charge) != 1 )
        warning( 'The lengths of the charge and formula vectors are not equal')
      charge = charge[1]
      if( charge != 0 ){
        mzTable <- abs((mzTable + 
                          (charge * (eMass$H - .00054857990946)))/charge) # mass of proton 
      }		
    }
    
    
    return(data.frame(formula = formulaList, mz = mzTable, z = charge, 
                      stringsAsFactors = F))	
  }
  
  rtCal <- function( calData, outData, knownList = knowns, charge = -1 ){
    
    if( missing(calData) )
      calData <- readPeak()
    if(any( is.na( match( outList$compound, knownList$compound)) ))
      stop( 'Compound names must match entries in the knowns list')
    
    ## Remove all standard Maven meta data
    mavNam <- c('label', 'metaGroupId', 'groupId', 'goodPeakCount', 'medMz', 
                'medRt', 'maxQuality', 'note', 'compound', 'compoundId', 'expectedRtDiff', 
                'ppmDiff', 'parent')	
    inDat <- calData[,-na.omit(match(mavNam, names(calData)))] 
    
    ## remove blanks
    blk <- grep('blank', names(inDat), ignore.case = TRUE)
    if( length( blk ) )
      inDat <- inDat[,-blk ]
    
    ## Find data max
    datMax <- apply( inDat, 1, sum )
    
    ## Find max intensity for each compound ID
    uCpd <- unique( calData$compound )
    obs <- data.frame(matrix(NA, ncol = ncol(calData), nrow = length(uCpd)))
    for( i in 1:length(uCpd) ){
      idx <- which( calData$compound == uCpd[i] )
      obs[i, ] <- calData[ idx[which.max( datMax[idx] )], ]			
    }
    names(obs) <- names(calData)
    
    ## Find RT and mz matches
    exp <- knownList[match( obs$compound, knownList$compound) ,]
    exp$mz <- mzCalc(exp$formula, charge = charge )$mz
    
    print('Mass calibration in ppm:') 
    print(summary((obs$medMz - exp$mz)/exp$mz *10^6))
    bm <- coef(lm(obs$medRt ~ exp$rt) ) 
    print('Retention calibration:')
    print(summary((obs$medRt - exp$rt)))
    print( bm )
    
    knownList$rt <- knownList$rt * bm[2] + bm[1]
    plot( exp$rt, obs$medRt )
    abline( coef = bm, lty = 2)	
    
    if( missing(outData) )
      writePeak( knownList )
    
    return(invisible(knownList))	
  } 
  
  ## Utility function for matching multiples with multiples
  ## x input list
  ## y index list
  ## returns match index in units of Y
  myMatch <- function(x, y){ 
    outI <- which(!is.na(match(y, x))) 
    out <- rep(NA, length(y))
    uX <- unique(x)
    for( i in 1:length(uX) ){
      if( is.na(uX[i]) )
        next
      if( any(y == uX[i]) )
        out[ y == uX[i] ] <- i
    }
    
    if( any(is.na(out)))
      out <- out[ -which(is.na(out)) ]	
    return(outI[order(out)])
  }
  
  ## Model for competative growth 
  ## pop1 - numeric - fraction of population 1 at time 0
  ## pop2 - numeric - fraction of population 2 at time 0
  ## lc1  - numeric - life cycle length pop1 (hours)
  ## lc2  - numeric - life cycle length pop2 (hours)
  ## r    - numeric - differential growth rate 
  ## g    - numeric - numer of generations (assumes 48h/g)
  compGrowth <- function(pop1 = .5, pop2 = .5, lc1 = 47, lc2 = 49, r1 = .065, 
                         r2 = -r1, g = 35, t0 = 43, sc = 8, ri = 1/sc, pSync = 1, sample48h = T, 
                         syncC = T){
    
    ## Define experimental setup
    t0 = t0 - 1
    Tte <- (-(48*g):(48*g*2)) - 1
    Ge <- Tte %/% 48
    Te <- Tte %% 48
    
    ## Calculate number of generations for each population
    g1 =(Tte+t0) %/% lc1
    g2 =(Tte+t0) %/% lc2 
    
    ## Calculate population sizes
    out <- data.frame(
      Tte = Tte, 
      Ge = Ge,
      Te = Te,
      g1 = g1,
      t1 =((Tte+t0) %% lc1) + 1,
      p1 =compG(p1 = pop1, g = g1, r1 = 1+r1 )$p1,
      DNA1 = rep(NA, length(g1)),
      g2 = g2, 
      t2 =((Tte+t0) %% lc2) + 1,
      p2 = compG(p2 = pop2, g = g2, r2 = 1+r2 )$p2,
      DNA2 = rep(NA, length(g1))
    )
    
    ## Calculate amount of DNA for each population
    D1 <- DNAr(t = 1:lc1, lc = lc1, sc = sc, ri = ri) 
    D2 <- DNAr(t = 1:lc2, lc = lc2, sc = sc, ri = ri)
    out$DNA1 <- D1[match(out$t1, 1:lc1)] 
    out$DNA2 <- D2[match(out$t2, 1:lc2)]
    out$DNA1 <- (out$DNA1 * out$p1) 
    out$DNA2 <- (out$DNA2 * out$p2) 
    
    ## Adjust for synchronicity
    if( syncC ){
      idx1 <- which( out$Tte == 0)
      idx2 <- rev(which( out$Ge <= g ))[1]
      out$DNAS1 <- out$DNAS2 <- NA
      syncT <- sync(c = 2)$p
      N <- (length(syncT) - 1)/2
      for( i in  idx1:idx2 ){
        syncT <- sync( c = out$Ge[i]/g * 6 + pSync )$p
        out[i,]$DNAS1 <- sum(out$DNA1[ (i - N) : (i + N)] * syncT)
        out[i,]$DNAS2 <- sum(out$DNA2[ (i - N) : (i + N)] * syncT)
      }		
      
      ## Calculate fractional pool
      TOT <- out$DNAS1  + out$DNAS2 
      out$DNAS1 <- out$DNAS1 / TOT 
      out$DNAS2 <- out$DNAS2 / TOT
    }
    
    ## Calculate the fractional pool
    TOT <- out$DNA1 + out$DNA2
    out$DNA1 <- out$DNA1 / TOT 
    out$DNA2 <- out$DNA2 / TOT
    
    out <- out[ out$Tte >= 0 & out$Ge <= g,]
    if( !sample48h )
      return(out)
    return(out[ out$Te == 0,])
  }
  
  compG <- function( p1 = 0.5, p2 = 1 - p1, r = .08, r1 = 1 + r, r2 = 1 - r, 
                     g = 1:35){
    
    P1T <- p1 * exp( r1 * g)
    P2T <- p2 * exp( r2 * g)
    tot <- P1T + P2T
    
    return( data.frame(g = g, p1 = P1T, p2 = P2T, tot ))
  }
  
  
  DNAr <- function(lc = 48, t = 1:lc, sc = 4, ri = .25){ 
    return( 1/(1+exp(-t/4 + lc/8)) * (sc - ri) + ri )
  }
  
  
  sync <- function(x = -240:240, t = 0, c = 4 ){
    gc <- exp(-(x-t)^2/(2*c^2))
    tot <- sum(exp(-(x-t)^2/(2*c^2)))
    return( data.frame(t = x, p = gc/tot) )
  }
  
  ## Mesure the RSQ for multiple regression curves
  fitComp <- function(obs, pop1, lc1, lc2, r1, t0, sc, pSync, 
                      N = 10000, rSample = T, syncC = F ){
    
    if( missing(obs) )
      stop( 'A vector of observed populations must be provided')
    
    ## Set sampling parameters	
    if( missing(pop1) )
      pop1 <- ((1:19) * 5) / 100
    pop2 <- 1 - pop1
    if( missing(lc1) )
      lc1 <- 38:58 #lc1 <- 45:51
    if( missing(lc2) )
      lc2 <- 38:58 #lc2 <- 45:51
    if( missing(r1) )
      r1 <- (-5:5) *2 #r1 <- .050+(0:10)/1000
    if( missing(t0) )
      t0 <- 38:48  #t0 <- 41:43
    if( missing(sc) )
      sc <- 1:10 #sc <- (-4:4)/4 + 4
    if( missing(pSync) )
      pSync <- 1:10 #pSync <- (1:4)*2
    if( !syncC )
      pSync <- 0
    
    fitT <- unique(expand.grid( pop1 = pop1, pop2 = pop2, lc1 = lc1, 
                                lc2=lc2, r1 = r1, t0=t0, sc=sc, pSync=pSync ))	
    if( rSample && nrow(fitT) > (N*2) )
      fitT <- fitT[ sort(sample(1:nrow(fitT), N)), ]
    
    if( !syncC )
      print( paste('Processing time:', round((nrow(fitT)/36)/60,2), 'minutes' ))
    else
      print( paste('Processing time:', round((nrow(fitT)*18)/60,2), 'minutes' ))
    flush.console()
    
    ## Calculate regression curves 
    print('Calculating regression curves ...')
    flush.console()
    fitData <- list()
    for( i in 1:nrow(fitT) ){
      fitData[[i]] <- compGrowth( 
        pop1 = fitT$pop1[i], 
        pop2 = fitT$pop2[i], 
        lc1 = fitT$lc1[i], 
        lc2 = fitT$lc2[i], 
        r1 = fitT$r1[i], 
        t0 = fitT$t0[i], 
        sc = fitT$sc[i], 
        pSync = fitT$pSync[i], 
        g = length(obs) - 1, syncC = syncC )
    }
    
    ## Measure the fit 
    print('Calculating RSQ ...')
    flush.console()
    RSQ1 <- list()
    if( syncC ){
      for( i in 1:length(fitData) )
        RSQ1[[i]] <- RSQ( obs = obs, fit = fitData[[i]]$DNAS1 )
    }else{
      for( i in 1:length(fitData) )
        RSQ1[[i]] <- RSQ( obs = obs, fit = fitData[[i]]$DNA1 )
    }
    fitT$RSQ <- unlist(RSQ1)
    
    return( fitT )
  }
  
  ## Calculate RSQ for each fit 
  RSQ <- function( obs, fit, na.omit = F ){
    
    if( !na.omit ){
      if( any(is.na(fit)) )
        return(NA)		
    }
    
    mod <- summary( lm(fit ~ obs) )
    return(mod$r.squared * (abs(mod$coefficients[2,1]) / mod$coefficients[2,1]) )
  }
  
  ## Table Mode
  tMode <- function(x){
    out <- apply(x, 2, table)
    out <- lapply(out, function(x) (rev(sort(x))[1:3]))
    out <- lapply( out, function(x) as.numeric(names(x[!is.na(x)])))
    return(out)
  }
  
  ## Increment vector
  inc <- function( x, y){
    
    out <- list()
    for( i in 1:length(y) )
      out[[i]] <- x + y[i]
    return( sort(unique(as.vector( unlist( out) ))))
  }
  
  ## predict peptide frag spectrum
  pepFrag <- function( pSeq, inMass, charge = 1, ppm = 1000 ){
    #	pepF <- pepFinder(pSeq = pSeq)
    #	pepF$mz <- mzCalc( pepF$formula, charge = charge)$mz 
    #	pepF$mz_H2O <- pepF$mz - (mzCalc('H2O', charge = 0)$mz / charge) 
    #	
    #	out <- data.frame( 
    #			compound = c(pepF$compound, paste(pepF$compound, '-H2O', sep ='')),
    #			mz = c( pepF$mz, pepF$mz_H2O), stringsAsFactors = FALSE	)
    #	out <- out[order(out$mz),]
    out <- fragFinder( pSeq = pSeq, charge = charge)
    
    if( missing(inMass))
      return( out )
    
    out <- out[mzMatch(mzObs = inMass, mzKno =out$mz, ppm = ppm)$mzKno,]
    return(out)
  }
  
  ## Internal function for processing big lists of knowns
  cpdId <- function(inList, knownList, blankName, mThresh, adList ){
    
    nLevel <- nrow(knownList) %/%  1000
    if( nrow(knownList) == (rev(nLevel)[1]*1000))
      nLevel <- c(nLevel, max(nLevel)+1 )
    
    out <- tmp <- NULL
    for( i in 1:nLevel ){
      idx <- ((i -1) * 1000 + 1) : ( i * 1000)
      tmp <- cpdMatch(inList = inList, knownList = knownList[idx,], 
                      blankName = blankName, mThresh = mThresh, adList = adList )
      out <- rbind(out, tmp)
    }
    
    if( is.null(out) || !nzchar(out) || nrow(out) < 1 )
      return(out)
    if( nrow(out) > 1000  ){
      print('WARNING: input list too large for final filter step')
      return(out)
    }
    uCompound <- unique(unlist(strsplit(out$compound, '_')))
    knownList <- knownList[na.omit(match(uCompound, knownList$compound)),]
    
    return(cpdMatch(inList = inList, knownList = knownList, blankName = blankName, 
                    mThresh = mThresh, adList = adList ))
    
  }
  
  ## function for matching untargeted lists with putative Ids
  cpdMatch <- function( inList, knownList, blankName = c('blank', 'RBC'),
                        mThresh = .8, 
                        adList = c(0, 1.003355, 1.007825, 22.9897692809-1.007825, 22.9897692809,
                                   17.00274, 18.01056, 43.98983, 35.03711) ){
    
    if( missing( inList) )
      inList <- readPeak() ## untargeted list	
    if( missing (knownList))
      inList <- readPeak() ## read known list
    
    if( nrow(knownList) > 1000 )
      return( cpdId(inList = inList, knownList = knownList, 
                    blankName = blankName, mThresh = mThresh, adList = adList ))
    
    nullList <- data.frame( compound = 'tmp', formula = 'tmp', 
                            mz = 1, rt = 1, note = NA, stringsAsFactors = F)[-1,]
    
    ## Match compound list with untargeted masses and allow for isotope overlap
    uMz <- sort(unique(knownList$mz) )
    uMz <- unique(c(uMz - 1.003355, uMz, uMz + 1.003355))
    dat2 <- inList[mzMatch(inList$medMz, uMz)$mzObs,]
    
    ## Find parent isotopomer and restrict data to samples with isotopes observed
    dat3 <- mzDiff(dat2, pThresh = mThresh, pW = c(1,1,1), adList = adList, 
                   collapse = T, clusterOnly = T)
    if( nrow(dat3) < 1 )
      return(nullList)
    
    ## Sequentally remove signals with high blank or user-related categories
    dat4 <- blankFilter(dat3[,-c(1:3)], name = blankName[1], thresh = .1)
    if( length(blankName) > 1 ){
      for( i in 2:length(blankName) ){
        if(!length(grep( blankName[i], names(dat4))))
          next
        dat4 <- blankFilter(dat4, name = blankName[i], thresh = .1)	
      }
    }
    if( nrow(dat4) < 1 )
      return(nullList)
    
    ## Return data to original inList format
    dat4 <- inList[match(dat4$groupId, inList$groupId),]
    
    ## calculate predicted masses from formula list
    uForm <- mzCalc(unique(knownList$formula))
    
    ## find all predicted formula that match observed masses
    dat5 <- sapply(mzMatch(dat4$medMz, uForm$mz, probMat = T), 
                   function(x) uForm$formula[x])
    if( length(dat5) < 1 )
      return(nullList)
    dat4$formula <- sapply(dat5, function(x) x[1])
    
    ## Put extra possilble formula matches in the note column	
    mzCount <- sapply(dat5, length)
    dat4[ mzCount > 1, ]$note <- sapply(dat5[ mzCount > 1 ], 
                                        function(x) paste(x, collapse = '_'))
    dat4 <- dat4[mzCount > 0,]
    
    ## Remove dimers 
    datN <- mzDiff( dat4, pThresh = .9, rtThresh = .1, pW = c(0,1,1)) 
    rmList <- NULL
    for( i in unique(datN$grp)){
      subL <- mzMatch( datN[datN$grp == i,]$medMz, 
                       datN[datN$grp == i,]$medMz * 2 + 1.007825)$mzObs
      if( !length(subL) )
        next
      rmList <- rbind(rmList, data.frame( grp = i, obs = subL)) 
    }
    tList <- myMatch(datN[datN$grp == rmList$grp,][rmList$obs,]$groupId, 
                     dat4$groupId)
    if( length(tList))
      dat4 <- dat4[-tList,]
    if( nrow(dat4) < 1 )
      return(nullList)
    
    ## Compile output table
    out <- data.frame( compound = dat4$compound, formula = dat4$formula, 
                       mz = dat4$medMz, rt = dat4$medRt, stringsAsFactors = F)
    out$note <- dat4$note
    out <- out[order(out$rt),]
    out <- out[order(out$mz),]
    for( i in 1:nrow(out)){
      out$compound[i] <- paste(
        knownList$compound[myMatch( out$formula[i], knownList$formula)], 
        collapse = '_')
    }
    out <- unique(out)
    
    ## Identify duplicate formulas 
    dupForm <- unique(out$formula[duplicated(out$formula)])	
    for( i in dupForm ){
      idx <- which(out$formula == i)
      out[idx,]$compound <- paste(out[idx,]$compound[1], 1:length(idx), sep='_') 
    }
    
    return(out)
  }
  
  
  
  
  ## Amino acid hydrophobicity
  calcHI <- function( inPep ){
    hiTab <- list(
      D  = 1, E  = 1, P  = 1, F  = 1, L  = 1, V  = 1, I  = 1,
      R  = -1, K  = -1, Q  = -1, N  = -1, H  = -1, 
      Y  = 0, W  = 0, S  = 0, T  = 0, G  = 0, A  = 0, M  = 0, C  = 0)
    hiTab <- data.frame( compound = names(hiTab), hi = unlist(hiTab), stringsAsFactors = F)
    hiTab2 <- hiTab
    hiTab2[ -match(c('D', 'E', 'R', 'K', 'Q', 'H'), hiTab$compound),]$hi <- 0
    
    return(
      data.frame(
        compound = inPep, 
        HI = sapply(strsplit(inPep, ''), 
                    function(x){ 
                      sum(hiTab$hi[match(x, hiTab$compound)]) + 
                        hiTab2$hi[match(rev(x)[1], hiTab$compound)]  +
                        hiTab2$hi[match((x)[1], hiTab$compound)]/2 }) 
      )
    )
  }
  
  
  
  ## Observed retention times for compounds that have been run on the exactive
  knowns <- data.frame(
    compound = c("cholesteryl sulfate", "Deoxycholic acid", 
                 "trans_trans-farnesyl diphosphate", "Cholic acid", "Geranyl-PP", 
                 "Taurodeoxycholic acid", "lipoate", "3-hydroxy-3-methylglutaryl-CoA-nega", 
                 "butyryl-CoA", "malonyl-CoA", "propionyl-CoA", "acetoacetyl-CoA", 
                 "acetyl-CoA", "dGTP", "succinyl-CoA/methylmalonyl-CoA", "malonyl-CoA", 
                 "3-methylphenylacetic acid", "coenzyme A", "acetoacetyl-CoA", 
                 "Phenylpropiolic acid", "3-hydroxybutyryl-CoA", 
                 "guanosine 5'-diphosphate-3'-diphosphate", 
                 "Indoleacrylic acid", "5-phosphoribosyl-1-pyrophosphate", "dephospho-CoA", 
                 "dTTP", "GTP", "ATP", "dATP", "NADPH", "phenylpyruvate", "UTP", 
                 "dUTP", "dephospho-CoA", "Phenyllactic acid", "2_3-Diphosphoglyceric acid", 
                 "CTP", "dCTP", "FMN", "UDP-D-glucuronate", "dGDP", 
                 "Hydroxyphenylacetic acid", 
                 "IDP", "Indole-3-carboxylic acid", "Kynurenic acid",
                 "cyclic bis(3'->5') dimeric GMP", 
                 "4-Pyridoxic acid", "Xanthurenic acid", "NADH", "ADP-D-glucose", 
                 "GDP", "ADP", "Hydroxyisocaproic acid", "5-methyl-THF", "aconitate", 
                 "2-Isopropylmalic acid", "Pyrophosphate", "citrate/isocitrate", 
                 "isocitrate", "citrate", "1_3-diphopshateglycerate", "phosphoenolpyruvate", 
                 "NADP+", "prephenate", "UDP", "2-oxo-4-methylthiobutanoate", 
                 "oxaloacetate", "Sedoheptoluse bisphosphate", "Octoluse Bisphosphate", 
                 "7_8-dihydrofolate", "folate", "hydroxyphenylpyruvate", "dCDP", 
                 "2_3-dihydroxybenzoic acid", "quinolinate", "fructose-1-6-bisphosphate", 
                 "3-phosphoglycerate", "UDP-N-acetyl-glucosamine", "dTDP", 
                 "adenosine 5'-phosphosulfate", 
                 "orotidine-5'-phosphate", "CDP", "fumarate", "anthranilate", 
                 "shikimate-3-phosphate", "UDP-D-glucose", "cyclic-AMP", 
                 "6-phospho-D-gluconate", 
                 "Citraconic acid", "a-ketoglutarate", "2-keto-isovalerate",
                 "N-acetyl-glutamate", 
                 "glutathione disulfide", "xanthosine-5-phosphate", "malate", 
                 "acetylphosphate", "N-carbamoyl-L-aspartate", "riboflavin", 
                 "2-Hydroxy-2-methylbutanedioic acid", 
                 "Carbamoyl phosphate", "dAMP", "3-S-methylthiopropionate", 
                 "Methylmalonic acid", 
                 "AMP", "dTMP", "nicotinate", "succinate", "thiamine-phosphate", 
                 "dGMP", "2-oxobutanoate", "shikimate", "p-hydroxybenzoate", "GMP", 
                 "IMP", "pantothenate", "dUMP", "UMP", "dCMP", "Thiamine pyrophosphate", 
                 "p-aminobenzoate", "NAD+", "CMP", "dihydroxy-acetone-phosphate", 
                 "pyruvate", "orotate", "3-phospho-serine", "deoxyadenosine", 
                 "N-Acetyl-L-alanine", "glutathione", "deoxyribose-phosphate", 
                 "acetoacetate", "xanthosine", "Pyroglutamic acid", "indole", 
                 "D-erythrose-4-phosphate", "glyoxylate", "D-glyceraldehdye-3-phosphate", 
                 "N-acetyl-glucosamine-1/6-phosphate", "sn-glycerol-3-phosphate", 
                 "myo-inositol", "tryptophan", "glucose-1-phosphate", "dihydroorotate", 
                 "S-adenosyl-L-homocysteine", "fructose-6-phosphate", "hexose-phosphate", 
                 "N-acetyl-glutamine", "S-adenosyl-L-methionine", 
                 "D-sedoheptulose-1/7-phosphate", 
                 "lactate", "glycolate", "trehalose-6-Phosphate", "glucose-6-phosphate", 
                 "glycerate", "Ascorbic acid", "glucono-δ-lactone", "CDP-ethanolamine", 
                 "CDP-choline", "D-glucono-δ-lactone-6-phosphate", "Uric acid", 
                 "ribose-phosphate", "Octoluse 8/1P", "S-adenosyl-L-homoCysteine", 
                 "allantoate", "thymidine", "2-dehydro-D-gluconate", "D-gluconate", 
                 "homocysteic acid", "deoxyguanosine", "adenosine", "aspartate", 
                 "glutamate", "guanosine", "deoxyinosine", "Cellobiose", "phenylalanine", 
                 "deoxyuridine", "inosine", "adenine", "guanine", "Aminoadipic acid", 
                 "thymine", "tyrosine", "leucine/isoleucine", "Guanidoacetic acid", 
                 "pyridoxine", "uridine", "methionine", "biotin", "cytidine", 
                 "xanthine", "hypoxanthine", "uracil", "trehalose/sucrose", "proline", 
                 "glutamine", "O-acetyl-L-serine", "cysteine", "cystathionine", 
                 "N-acetyl-L-ornithine", "D-glucosamine-6-phosphate", 
                 "D-glucosamine-1-phosphate", 
                 "histidine", "ornithine", "valine", "asparagine", "cytosine", 
                 "threonine", "homoserine", "homocysteine", "serine", "alanine", 
                 "citrulline", "allantoin", "taurine", "lysine", "histidinol", 
                 "thiamine", "glucosamine", "S-ribosyl-L-homocysteine-nega", "arginine", 
                 "Creatinine", "Imidazoleacetic acid", "DL-Pipecolic acid", 
                 "N-Acetylputrescine", 
                 "Methylcysteine", "Pyridoxamine", "1-Methyl-Histidine", "Acetyllysine", 
                 "Acetylcarnitine", "Kynurenine", "Flavone", "Glycerophosphocholine", 
                 "1-Methyladenosine", "Nicotinamide ribotide", "Diiodothyronine", 
                 "Atrolactic acid", "D-glucarate", 
                 "aminoimidazole carboxamide ribonucleotide", 
                 "acadesine", "carnitine", "nicotinamide", "choline", "4-aminobutyrate", 
                 "glycine", "FAD", "S-adenosyl-L-methioninamine", 
                 "S-methyl-5'-thioadenosine", 
                 "L-arginino-succinate", "ethanolamine", "betaine aldehyde", 
                 "dimethylglycine", 
                 "betaine", "purine", "creatine", "hydroxyproline", "methylnicotinamide", 
                 "2-Aminooctanoic acid", "Methionine sulfoxide", "Phosphorylcholine", 
                 "NG-dimethyl-L-arginine", "N-acetyl-glucosamine", "5-methoxytryptophan", 
                 "Cystine", "7-methylguanosine", "sarcosine"),
    formula = c("C27H46O4S1", "C26H43N1O5", "C15H28O7P2", 
                "C24H40O5", "C10H20O7P2", "C26H45N1O6S1", "C8H14O2S2", "C27H44N7O20S1P3", 
                "C25H42N7O17S1P3", "C24H38N7O19S1P3", "C24H40N7O17S1P3", "C25H40N7O18S1P3", 
                "C23H38N7O17S1P3", "C10H16N5O13P3", "C25H40N7O19S1P3", "C24H38N7O19S1P3", 
                "C9H10O2", "C21H36N7O16S1P3", "C25H40N7O18S1P3", "C9H6O2", "C25H42N7O18S1P3", 
                "C10H11N5O17P4", "C11H9N1O2", "C5H13O14P3", "C21H35N7O13S1P2", 
                "C10H17N2O14P3", "C10H16N5O14P3", "C10H16N5O13P3", "C10H16N5O12P3", 
                "C21H30N7O17P3", "C9H8O3", "C9H15N2O15P3", "C9H15N2O14P3", "C21H35N7O13S1P2", 
                "C9H10O3", "C3H8O10P2", "C9H16N3O14P3", "C9H16N3O13P3", "C17H21N4O9P1", 
                "C15H22N2O18P2", "C10H15N5O10P2", "C8H8O3", "C10H14N4O11P2", 
                "C9H7N1O2", "C10H7N1O3", "C20H24N10O14P2", "C8H9N1O4", "C10H7N1O4", 
                "C21H29N7O14P2", "C16H25N5O15P2", "C10H15N5O11P2", "C10H15N5O10P2", 
                "C6H12O3", "C20H25N7O6", "C6H6O6", "C7H12O5", "H4O7P2", "C6H8O7", 
                "C6H8O7", "C6H8O7", "C3H8O10P2", "C3H5O6P1", "C21H28N7O17P3", 
                "C10H10O6", "C9H14N2O12P2", "C5H8O3S1", "C4H4O5", "C7H16O13P2", 
                "C8H18O14P2", "C19H21N7O6", "C19H19N7O6", "C9H8O4", "C9H15N3O10P2", 
                "C7H6O4", "C7H5N1O4", "C6H14O12P2", "C3H7O7P1", "C17H27N3O17P2", 
                "C10H16N2O11P2", "C10H14N5O10S1P1", "C10H13N2O11P1", "C9H15N3O11P2", 
                "C4H4O4", "C7H7N1O2", "C7H11O8P1", "C15H24N2O17P2", "C10H12N5O6P1", 
                "C6H13O10P1", "C5H6O4", "C5H6O5", "C5H8O3", "C7H11N1O5", "C20H32N6O12S2", 
                "C10H13N4O9P1", "C4H6O5", "C2H5O5P1", "C5H8N2O5", "C17H20N4O6", 
                "C5H8O5", "C1H4N1O5P1", "C10H14N5O6P1", "C4H8O2S1", "C4H6O4", 
                "C10H14N5O7P1", "C10H15N2O8P1", "C6H5N1O2", "C4H6O4", "C12H18N4O4S1P1", 
                "C10H14N5O7P1", "C4H6O3", "C7H10O5", "C7H6O3", "C10H14N5O8P1", 
                "C10H13N4O8P1", "C9H17N1O5", "C9H13N2O8P1", "C9H13N2O9P1", "C9H14N3O7P1", 
                "C12H19N4O7S1P2", "C7H7N1O2", "C21H27N7O14P2", "C9H14N3O8P1", 
                "C3H7O6P1", "C3H4O3", "C5H4N2O4", "C3H8N1O6P1", "C10H13N5O3", 
                "C5H9N1O3", "C10H17N3O6S1", "C5H11O7P1", "C4H6O3", "C10H12N4O6", 
                "C5H7N1O3", "C8H7N1", "C4H9O7P1", "C2H2O3", "C3H7O6P1", "C8H16N1O9P1", 
                "C3H9O6P1", "C6H12O6", "C11H12N2O2", "C6H13O9P1", "C5H6N2O4", 
                "C14H20N6O5S1", "C6H13O9P1", "C6H13O9P1", "C7H12N2O4", "C15H23N6O5S1", 
                "C7H15O10P1", "C3H6O3", "C2H4O3", "C12H23O14P1", "C6H13O9P1", 
                "C3H6O4", "C6H8O6", "C6H10O6", "C11H20N4O11P2", "C14H27N4O11P2", 
                "C6H11O9P1", "C5H4N4O3", "C5H11O8P1", "C8H17O11P1", "C14H20N6O5S1", 
                "C4H8N4O4", "C10H14N2O5", "C6H10O7", "C6H12O7", "C4H9N1O5S1", 
                "C10H13N5O4", "C10H13N5O4", "C4H7N1O4", "C5H9N1O4", "C10H13N5O5", 
                "C10H12N4O4", "C12H22O11", "C9H11N1O2", "C9H12N2O5", "C10H12N4O5", 
                "C5H5N5", "C5H5N5O1", "C6H11N1O4", "C5H6N2O2", "C9H11N1O3", "C6H13N1O2", 
                "C3H7N3O2", "C8H11N1O3", "C9H12N2O6", "C5H11N1O2S1", "C10H16N2O3S1", 
                "C9H13N3O5", "C5H4N4O2", "C5H4N4O1", "C4H4N2O2", "C12H22O11", 
                "C5H9N1O2", "C5H10N2O3", "C5H9N1O4", "C3H7N1O2S1", "C7H14N2O4S1", 
                "C7H14N2O3", "C6H14N1O8P1", "C6H14N1O8P1", "C6H9N3O2", "C5H12N2O2", 
                "C5H11N1O2", "C4H8N2O3", "C4H5N3O1", "C4H9N1O3", "C4H9N1O3", 
                "C4H9N1O2S1", "C3H7N1O3", "C3H7N1O2", "C6H13N3O3", "C4H6N4O3", 
                "C2H7N1O3S1", "C6H14N2O2", "C6H11N3O1", "C12H17N4O1S1", "C6H13N1O5", 
                "C9H17N1O6S1", "C6H14N4O2", "C4H7N3O1", "C6H8N2O2", "C6H11N1O2", 
                "C6H14N2O1", "C4H9N1O2S1", "C8H12N2O2", "C7H11N3O2", "C8H16N2O3", 
                "C9H18N1O4", "C10H12N2O3", "C15H10O2", "C8H21N1O6P1", "C11H15N5O4", 
                "C11H15N2O8P1", "C15H13N1O4I2", "C9H10O3", "C6H10O8", "C9H15N4O8P1", 
                "C9H14N4O5", "C7H15N1O3", "C6H6N2O1", "C5H14N1O1", "C4H9N1O2", 
                "C2H5N1O2", "C27H33N9O15P2", "C14H23N6O3S1", "C11H15N5O3S1", 
                "C10H18N4O6", "C2H7N1O1", "C5H12N1O1", "C4H9N1O2", "C5H11N1O2", 
                "C5H4N4", "C4H9N3O2", "C5H9N1O3", "C7H9N2O1", "C8H17N1O2", "C5H11N1O3S1", 
                "C5H15N1O4P1", "C8H18N4O2", "C8H16N1O9P1", "C12H14N2O3", "C6H12N2O4S2", 
                "C11H16N5O5", "C3H7N1O2"),
    rt = c( 17.25, 16.79, 16.74, 16.69, 16.46, 16.23, 15.97, 15.72, 15.72, 15.7,
            15.64, 15.5, 15.46, 15.46, 15.41, 15.41, 15.41, 15.4, 15.36, 15.2, 15.2,
            15.18, 14.94, 14.94, 14.94, 14.94, 14.94, 14.94, 14.94, 14.87, 15.1, 
            14.81, 14.81, 14.7, 14.68, 14.68, 14.9, 14.8, 14.4, 14.68, 14, 14.68, 
            13.9, 14.31, 14.22, 14.2, 14.14, 14.14, 14.14, 14.14, 13.9, 13.8, 14.14,
            14.1, 14, 13.87, 13.87, 13.6, 13.6, 13.6, NA, 13.87, 13.87, 13.87, 13.81, 
            13.75, 13.64, 13.64, 13.64, NA, 13.6, 13.58, 13.58, 13.58, 13.58, 13.5, 
            13.58, 13.58, 13.58, 13.58, 13.58, 13.53, 13.5, 13.47, 13.44, 13.44, 
            13.44, 13.38, 13.29, 13.15, 13.06, 13, 12.94, 12.7, 12.7, 12.38, 12.25,
            12.1, 12.06, 12.06, 12, 11.92, 11.72, 11.3, 11.27, 11.2, 11.6, 11, 11, 
            10.95, 10.92, 10.88, 10.5, 10.5, 11, 10.28, 10, 9.8, 9.56, 8.7, 8.62, 
            8.5, 9, 8.39, 8.1, 8.2, 8, 7.97, 7.7, 7.74, 7.74, 7.74, 7.69, 7.5, 7.45, 
            7.45, 7.35, 7.3, 7.3, 7.2, 7.5, 7.3, 6.9, 7.1, 7.2, 7, 7, 7, 6.95, 6.95, 
            6.95, 6.95, 6.85, 6.69, 6.59, 6.59, 6.38, 6.33, 6.17, 6.17, 6, 6, 6, 
            5.62, 5.4, 5, 5, 4.8, 4.8, 6.7, 4.05, 4.2, 4, 3.9, 3.86, 3.8, 3.54, 3.5, 
            3.5, 3.4, 3.22, 2.6, 2, 1.9, 1.89, 1.8, 1.7, 1.5, 12.6, 1.2, 1.18, 1.18, 
            1.18, 1.18, 1.12, 1.1, 1.1, 1.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
            1, 0.9, 0.82, 0.82, 0.8, 0.8, 0.8, 0.7, 0.68, 0.6, NA, NA, NA, 1, 1, 
            7.6, 2, 1.1, 2.4, 3.6, 7.7, NA, 3, NA, NA, NA, 13, NA, NA, 1, NA, NA, 4, 
            1, 15, 1, 11.6, 7.7, NA, NA, NA, NA, NA, NA, 1, 3.5, NA, NA, NA, NA, 2, 
            NA, NA, 12, 1), 
    mz = c(465.304405579909, 448.306847579909, 
           381.123751839909, 407.280298579909, 313.061151839909, 498.289484579909, 
           205.036247579909, 910.150198469909, 836.149803469909, 852.108333469909, 
           822.13415346991, 850.129068469909, 808.118503469909, 505.988473469909, 
           866.123983469909, 852.108333469909, 149.060803579909, 766.10793846991, 
           850.129068469909, 145.029503579909, 852.144718469909, 595.90277009991, 
           186.056052579909, 388.944543469909, 686.14160683991, 480.981991469909, 
           521.983388469909, 505.988473469909, 489.993558469909, 744.083831469909, 
           163.040068579909, 482.961256469909, 466.966341469909, 686.14160683991, 
           165.055718579909, 264.951996839909, 481.977240469909, 465.982325469909, 
           455.097341209909, 579.027014839909, 426.022141839909, 151.040068579909, 
           427.006157839909, 160.040402579909, 188.035317579909, 689.087596839909, 
           182.045882579909, 204.030232579909, 664.117499839909, 588.07496683991, 
           442.017056839909, 426.022141839909, 131.071368579909, 458.179356579909, 
           173.009163579909, 175.061198579909, 176.935951839909, 191.019728579909, 
           191.019728579909, 191.019728579909, 264.951996839909, 166.975100209909, 
           742.06818146991, 225.040463579909, 402.994924839909, 147.012140579909, 
           130.998598579909, 368.999341839909, 399.009906839909, 442.148056579909, 
           440.132406579909, 179.034983579909, 386.015993839909, 153.019333579909, 
           166.014582579909, 338.988776839909, 184.985665209909, 606.074298839909, 
           401.015659839909, 426.012627209909, 367.018423209909, 402.010908839909, 
           115.003683579909, 136.040402579909, 253.011880209909, 565.047749839909, 
           328.045245209909, 275.017360209909, 129.019333579909, 145.014248579909, 
           115.040068579909, 188.056447579909, 611.144691579909, 363.034741209909, 
           133.014248579909, 138.980185209909, 175.036046579909, 375.131009579909, 
           147.029898579909, 139.975434209909, 330.060895209909, 119.017225579909, 
           117.019333579909, 346.055810209909, 321.049328209909, 122.024752579909, 
           117.019333579909, 344.071363209909, 346.055810209909, 101.024418579909, 
           173.045548579909, 137.024418579909, 362.05072520991, 347.039826209909, 
           218.103397579909, 307.033678209909, 323.02859320991, 306.049662209909, 
           424.037694839909, 136.040402579909, 662.10184983991, 322.044577209909, 
           168.990750209909, 87.0087685799095, 155.009831579909, 184.001649209909, 
           250.094563579909, 130.050967579909, 306.076532579909, 213.016965209909, 
           101.024418579909, 283.068409579909, 128.035317579909, 116.050572579909, 
           199.001315209909, 72.9931185799095, 168.990750209909, 300.048994209909, 
           171.006400209909, 179.056113579909, 203.082601579909, 259.022445209909, 
           157.025481579909, 383.114314579909, 259.022445209909, 259.022445209909, 
           187.072431579909, 398.137789579909, 289.033010209909, 89.0244185799095, 
           75.0087685799095, 421.075270209909, 259.022445209909, 105.019333579909, 
           175.024813579909, 177.040463579909, 445.053107839909, 488.107882839909, 
           257.006795209909, 167.021064579909, 229.011880209909, 319.043575209909, 
           383.114314579909, 175.047279579909, 241.082996579909, 193.035378579909, 
           195.051028579909, 182.012869579909, 266.089478579909, 266.089478579909, 
           132.030232579909, 146.045882579909, 282.084393579909, 251.078579579909, 
           341.108938579909, 164.071702579909, 227.067346579909, 267.073494579909, 
           134.047218579909, 150.042133579909, 160.061532579909, 125.035651579909, 
           180.066617579909, 130.087352579909, 116.046550579909, 168.066617579909, 
           243.062261579909, 148.043774579909, 243.080888579909, 242.078245579909, 
           151.026149579909, 135.031234579909, 111.020001579909, 341.108938579909, 
           114.056052579909, 145.061866579909, 146.045882579909, 120.012474579909, 
           221.060153579909, 173.093166579909, 258.038429209909, 258.038429209909, 
           154.062200579909, 131.082601579909, 116.071702579909, 131.046216579909, 
           110.035985579909, 118.050967579909, 118.050967579909, 134.028124579909, 
           104.035317579909, 88.0404025799095, 174.088415579909, 157.036714579909, 
           124.007389579909, 145.098251579909, 140.082935579909, 264.105031579909, 
           178.072097579909, 266.07038457991, 173.104399579909, 112.051635579909, 
           139.051301579909, 128.071702579909, 129.103336579909, 134.028124579909, 
           167.082601579909, 168.077850579909, 187.108816579909, 203.116307579909, 
           207.077516579909, 221.060803579909, 257.103374209909, 280.105128579909, 
           333.049328209909, 523.886128579909, 165.055718579909, 209.030293579909, 
           337.055476209909, 257.089144579909, 160.097917579909, 121.040736579909, 
           103.100262579909, 102.056052579909, 74.0247525799095, 784.149862839909, 
           354.147959579909, 296.082285579909, 289.115359579909, 60.0454875799095, 
           101.084612579909, 102.056052579909, 116.071702579909, 119.036319579909, 
           130.062200579909, 130.050967579909, 136.064211579909, 158.118652579909, 
           164.038689579909, 183.066594209909, 201.135699579909, 300.048994209909, 
           233.093166579909, 239.016575579909, 297.107868579909, 88.0404025799095), 
    stringsAsFactors = FALSE )
  
  
  ## Read in EupathDB protein FASTA files
  readEuPath <- function(inPath = file.choose() ){
    inList <- scan(inPath, sep="\n", what = 'A', allowEscapes = F,  
                   strip.white = T, blank.lines.skip = T)
    inList <- paste(inList, collapse = "")
    inList <- strsplit(inList, ">")[[1]][-1]
    inList <- strsplit(inList, " | ", fixed = T)
    return(inList)
  }
  
  ## Read FASTA
  readFASTA <- function(inPath = file.choose() ){
    inList <- scan(inPath, sep="\n", what = 'A', allowEscapes = F,  
                   strip.white = T, blank.lines.skip = T)
    inList <- paste(inList, collapse = "")
    inList <- gsub(' OS=', "~", inList, fixed = T)
    inList <- gsub(' PE=', "~", inList, fixed = T)
    inList <- gsub('>', "", inList, fixed = T)
    inList <- gsub(' SV=', "~", inList, fixed = T)
    inList <- gsub('tr|', "sp|", inList, fixed = T)
    inList <- strsplit(inList, "sp|", fixed = T)[[1]][-1]
    
    ## remove proteins that aren't tied to genes
    if( any( !grepl('GN=', inList) ))
      inList <- inList[-which(!grepl('GN=', inList))]
    
    inList <- gsub(' GN=', "~", inList, fixed = T)
    inList <- sapply(inList, function(x) strsplit(x, "~", fixed = T))
    names(inList) <- NULL
    
    ## Remove SV from start of protein sequence
    for(i in 1:length(inList ))
      inList[[i]][5] <- substr(inList[[i]][5], start = 2, stop = nchar(inList[[i]][5]))
    
    return(inList)
  }
  
  readBrenda <- function( inDir = choose.dir(), printSum = T ){
    inList <- list.files( inDir, full.names = T )
    outList <- list()
    for(i in 1:length(inList))
      outList[[i]] <- scan(inList[i], sep="\n", what = 'A', allowEscapes = T,  
                           strip.white = F, blank.lines.skip = F, quiet = T)
    outList <- as.vector(unlist(outList))
    outList <- sapply(outList, function(x) strsplit(x, "\t", fixed = T) )
    outList <- matrix(as.vector(unlist(outList)), ncol = 3, byrow = T)
    
    if(printSum)
      print(paste(names(summary(factor(outList[,3]))[1]), ": ", 
                  length(unique(outList[,1])), sep = ""), quote = F)
    
    return(outList)
  }
  
  ## Read human gene map
  readDat <- function( inPath = file.choose(), uIdn = c("PF3D7", "API", "mito"), 
                       outPath ){
    
    ## Read input file
    print(paste('reading file:', basename(inPath)), quote = F); flush.console()
    inList <- scan(inPath, sep="\n", what = 'A', allowEscapes = F,  
                   strip.white = F, blank.lines.skip = F)
    
    ## setup outgoing folder
    print('   parsing file ...', quote = F); flush.console()
    out <- NULL
    for( i in 1:length(uIdn) ){
      out <- cbind(out, sapply(inList[grep(uIdn[i], inList)],	function(x){
        tList <- unlist(strsplit(x, ',', fixed = T))
        return(c(tList[1:10], paste(tList[-(1:10)], collapse='~')))
      }))
      dimnames(out) <- NULL		
      
    }
    out <- t(out)
    
    ## Find the ions and gene records and parse the data
    idx <- substring(out[,1], 0, 1) == 'q'
    out <- list( dat = out[idx,], gene = 
                   apply(out[!idx,], 1, function(x) paste(na.omit(x), collapse = '')))
    names(out$gene) <- NULL
    out$gene <- strsplit(out$gene, "|", fixed = T)
    out$gene <- data.frame(matrix(unlist(out$gene), ncol = 5, byrow = T),
                           stringsAsFactors = F)
    out$gene$gene <- sapply(strsplit(out$gene[,1], "\"", fixed = T), 
                            function(x) x[2])
    out$gene$chr <- sapply(strsplit(out$gene$X2, "_", fixed = T), function(x) x[2])
    out$gene$loc <- sapply(
      strsplit(
        sub("(", "-", sapply(
          strsplit(out$gene$X2, ":", fixed = T), 
          function(x)x[2]), fixed = T), "-", fixed = T), 
      function(x) mean(as.numeric(x[1:2])))
    out$gene$len <- as.numeric(sapply(strsplit(out$gene$X3, "=", 
                                               fixed = T), function(x) x[2]))
    
    print('   making output file ...', quote = F); flush.console()
    out$match <- strsplit(out$dat[,11], "~")
    for( i in 1:length(out$match) )
      out$match[[i]] <-  paste(out$match[[i]], i, length(out$match[[i]]), sep = "\"")
    out$match <- t(sapply(
      strsplit(sub("0;", "", unlist(out$match), fixed = T), "\""),
      function(x) x[2:5]))
    out$match <- cbind(out$match, t(sapply(strsplit(out$match[,2], ":", 
                                                    fixed = T), function(x) x[3:4])))
    dimnames(out$match) <- NULL	
    
    out$match <- data.frame(out$match, stringsAsFactors = F)
    names(out$match) <- c('gene', 'matchInfo', 'datIdx', 'nGen', 'start', 'end')
    
    out$out <- merge(out$gene[, -(1:5)], out$match, by.x = 'gene', by.y = 'gene') 
    out$out$datIdx <- as.numeric(out$out$datIdx)
    out$out <- cbind(out$out, out$dat[out$out$datIdx,1:10], stringsAsFactors = F)	
    names(out$out) <- c('gene', 'chr', 'loc', 'len', 'matchInfo', 'datIdx', 
                        'nGen', 'start', 'end', 'ionInfo', 'mz', 'ppm', 'x1', 'pepSeq', 'x2', 
                        'x3', 'Mscore', 'x4', 'x5')
    
    if(!missing(outPath))
      write.csv( out$out, file = outPath, quote = F, row.names = F )
    else
      return(out$out)
  }
  
  ## Read Panther
  readPanth <- function( inPath = file.choose() ){
    inTable <- scan(file.choose(), what = 'A', sep = "\t")
    inTable <- gsub("ENSEMBL=", "|", inTable, fixed = T)	
    inTable <- as.vector(unlist(strsplit( inTable, "|", fixed = T)))
    inTable <- inTable[grep('ENSG', inTable)]
    return(inTable)
  }
  
  ## Parse batch NCBI gene queery
  readNCBI <- function( inPath = file.choose() ){
    inList <-  scan(inPath, sep="\n", what = 'A', blank.lines.skip = F )
    inList <- paste(inList[-1], collapse = "~")
    inList <- strsplit(inList, "~~")[[1]]
    xT <- length(inList)
    outList <- data.frame(
      geneID = rep(NA, xT), 
      chr = rep(NA, xT), 
      loc = rep(NA, xT)) 
    
    for(i in 1:length(inList)){
      
      subL <- unlist(strsplit( inList[i], "~"))
      if( subL[1] == "" && (length(subL) > 1) )
        subL <- subL[ -1 ]
      
      outList$gene[i] <- unlist(strsplit(subL[1], " "))[2]
      if(any(grepl("Chromosome:", subL )))
        outList$chr[i] <- sub(",", "", sub(";", "", 
                                           unlist(strsplit(subL[ grepl("Chromosome:", subL )], " "))[2]), fixed = T) 
      
      if( is.na(outList$chr[i]) & any(grepl("Genomic context: Apicoplast", subL )))
        outList$chr[i] <- "Api"
      
      if(any(grepl("Annotation:", subL ))){
        tmpList <- strsplit(subL[ grepl("Annotation:", subL )], " ")[[1]]
        tmpList <- tmpList[grep( "..", tmpList, fixed = T)]
        outList$loc[i] <- mean(as.numeric(unlist(strsplit(sub(",", "", sub(")", "", sub("(", "", 
                                                                                        tmpList, 
                                                                                        , fixed = T), fixed = T) , fixed = T), "..", fixed = T))))  			
      }
      
    }
    
    outList$chr[ is.na(outList$chr) ] <- 'Ukn'
    
    outList$nLoc <- NA
    uChr <- sort(unique(outList$chr))
    for( i in 1:length(uChr) ){
      chrR <- range(outList[outList$chr == uChr[i],]$loc, na.rm = T)	
      outList[outList$chr == uChr[i],]$nLoc <- 
        (outList[outList$chr == uChr[i],]$loc - chrR[1])/(chrR[2] - chrR[1])
    }
    
    return(outList)
  }
  
  
  ## Make gene/chromosome table
  genLocPf <- function( inList ){
    if( missing(inList) )
      inList <- readEuPath()
    
    chr <- sapply(inList, function(x) x[1] )
    chr <- sapply(chr, function(x) substring (strsplit(x, "_"   )[[1]][2],1,2))
    loc <- sapply(inList, function(x) 
      mean(as.numeric(strsplit(strsplit(
        strsplit(x[4], ":")[[1]][2], "(", fixed = T)[[1]][1], "-")[[1]]))
    )
    
    seq <- sapply(inList, function(x) gsub('SO=protein_coding', "", x[7] ))
    
    out <- data.frame( gene = names(chr), chr = as.vector(chr), loc = loc, 
                       nLoc = loc, seq = seq, stringsAsFactors = F)
    
    uChr <- sort(unique(out$chr))
    for( i in 1:length(uChr) ){
      
      chrR <- range(out[out$chr == uChr[i],]$loc)	
      out[out$chr == uChr[i],]$nLoc <- (out[out$chr == uChr[i],]$loc - chrR[1])/
        (chrR[2] - chrR[1])
    }
    
    return(out)
  }
  
  genLocHs <- function( inList, inMap ){
    if( missing(inList) )
      inList <- readFASTA()
    
    out <- data.frame(matrix(as.vector(unlist(inList)), ncol = 5, byrow = T), 
                      stringsAsFactors = F)
    names(out) <- c('name', 'species', 'gene', 'ukn', 'seq')
    out <- merge(out, inMap, by.x = 'gene', by.y = 'Input' )
    
    names(out)[ names(out) == 'Chr'] <- 'chr'
    names(out)[ names(out) == 'Start'] <- 'start'
    names(out)[ names(out) == 'End'] <- 'end'
    out$loc <- apply(out[,match(c('start', 'end'), names(out))], 1, mean)
    
    return(out)
  }
  
  
  ## Plot genome-wide hits
  pepPlot <- function( chr, obs, gene, topQuant = F, type = 'map', geneName = T, 
                       logS = T, ...){
    
    obs <- obs[order(obs$start),]
    
    uX <- unique(chr$chr)
    uX <- suppressWarnings(uX[ order(as.numeric(uX)) ])
    
    x <- match( obs$chr, uX ) + obs$nLoc *.8 - .4
    y <- (obs$pepOL)
    ylab = 'Matches'
    
    if( logS ){
      y <- log(y)
      ylab = "ln(Matches)"
    }
    
    
    out <- unique(data.frame(x = x, y = y, gene = obs$gene, 
                             stringsAsFactors = F))
    x <- out$x
    y <- out$y
    
    if (type == 'map' ){
      
      plot( 1:length(uX), rep(0, length(uX)), type = 'n', 
            ylim = c(0, max(y)), axes = F, xlab = 'Chromosome', ylab= ylab, ...)
      
      axis( 1, at = 1:length(uX), labels = uX, tick = F)
      axis(2)
      
      fN <- fivenum( y )
      
      for( i in 1:length(uX))
        lines( c(i - .4, i + .4), rep(0, 2), lwd = 2, col = 'black')
      points(x[y < fN[2]],y[y < fN[2]], pch = 16, col = 'grey')
      points(x[(y >= fN[2]) & (y < fN[3])],y[(y >= fN[2]) & (y < fN[3])], pch = 16, col = 'blue')
      points(x[(y >= fN[3]) & (y < fN[4])],y[(y >= fN[3]) & (y < fN[4])], pch = 16, col = 'orange')
      points(x[(y >= fN[4])],y[(y >= fN[4])], pch = 16, col = 'red')
      
      if( geneName ){
        text(x,y, labels = out$gene, pos = 4)
        
      }
    }
    
    if( type == 'gene' ){
      
      if( missing (gene))
        gene <- obs[obs$pepOL == max( obs$pepOL),]$gene[1]
      
      if( topQuant ){
        fN <- fivenum( y )
        gene <- unique(obs[ obs$pepOL >= fN[4], ]$gene)			
      }
      
      
      for( i in 1:length(gene)){
        sub <- obs[obs$gene == gene[i], ]
        seqPlot( sub$seq[1], sub, main = gene[i], ...)
      }		
    }
  }
  
  plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                         col = "black", sep = 0.5, ...) {
    height <- 1
    if (is(xlim, "Ranges"))
      xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins) * (height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom +
           height, col = col, ...)
    title(main)
    axis(1)
  }
  
  
  ## Read mascot output and create a Maven inpupt file
  mascot2Maven <- function( inTable, Modifications = F, PPM = 10, nMatch = 5 ){
    
    inName <- c('Sequence', 'IonScore', 'z', 
                'mz_obs', 'mz_z1', 'ppm', 'rt', 'Protein') 
    
    ## Enforce specific names as input
    if(any(is.na(match(inName, names(inTable)))))
      stop( paste(c('input files must have the following names:', 
                    paste(" ",inName))))
    
    ## Don't allow modified peptides 
    if( !Modifications )
      inTable <- inTable[ inTable$Modifications == "", ]
    
    ## PPM filter
    idx <- grep('ppm', names(inTable))
    if( length(idx) != 1 )
      stop('check the ppm names in input file')
    inTable <- inTable[abs(inTable[ , idx]) < PPM,]
    
    ## Create output table
    out <- data.frame( compound = sort(unique(inTable$Sequence)), 
                       mz = NA, rt = NA, z = NA, mz_z1 = NA, score = NA, ppm = NA, 
                       nMatch = NA, note = NA, stringsAsFactors = F)
    
    ## find the high probability peptides
    for( i in 1:nrow(out) ){
      tCpd <- out$compound[i]
      tScore <- which.max( inTable[ inTable$Sequence == tCpd,]$IonScore )
      if( length(tScore) == 0)
        tScore <- 1
      tTab <- inTable[ inTable$Sequence == tCpd,][tScore,]
      out$mz[i] <- tTab$mz_obs
      out$rt[i] <- tTab$rt
      out$z[i] <- tTab$z
      out$mz_z1[i] <- tTab$mz_z1
      out$score[i] <- tTab$IonScore
      out$ppm[i] <- tTab$ppm
      out$nMatch[i] <- nrow(inTable[ inTable$Sequence == tCpd,] )
      out$note[i] <- tTab$Protein
    }
    
    ## Find number of AA in each peptide
    keep <- (nchar(out$compound) < 6 & !(out$score < 7 & out$nMatch < 5))
    
    ## Filter outgoing list
    out <- out[ keep | (out$score >= 15 & nMatch > 3),]
    
    ## Hemoglobin specific sorting
    out$note <- sub("Hemoglobin subunit alpha OS=Homo sapiens GN=HBA1 PE=1 SV=2 - [HBA_HUMAN]", 
                    "HBA", out$note, fixed = T)
    out$note <- sub("Hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2 - [HBB_HUMAN]", 
                    "HBB", out$note, fixed = T)
    out$note <- sub("Hemoglobin subunit zeta OS=Homo sapiens GN=HBZ PE=1 SV=2 - [HBAZ_HUMAN]", 
                    "HBAZ", out$note, fixed = T)	
    out$note <- sub("Hemoglobin subunit delta OS=Homo sapiens GN=HBD PE=1 SV=2 - [HBD_HUMAN]", 
                    "HBD", out$note, fixed = T)
    out$note <- sub("Sickle mutant hemoglobin subunit beta OS=Homo sapiens GN=HBB PE=1 SV=2 - [HBS_HUMAN]", 
                    "HBS", out$note, fixed = T)
    out$note <- sub("Hemoglobin subunit epsilon OS=Homo sapiens GN=HBE1 PE=1 SV=2 - [HBE_HUMAN]", 
                    "HBE", out$note, fixed = T)
    out$note <- sub("Hemoglobin subunit gamma-1 OS=Homo sapiens GN=HBG1 PE=1 SV=2 - [HBG1_HUMAN]", 
                    "HBG1", out$note, fixed = T)
    out$note <- sub("Hemoglobin subunit gamma-2 OS=Homo sapiens GN=HBG2 PE=1 SV=2 - [HBG2_HUMAN]", 
                    "HBG2", out$note, fixed = T)	
    out <- out[order(out$note),]	
    out <- out[order(grepl('HBB', out$note)),]	
    out[nrow(out):1,] <- out[order(grepl('HBA', out$note)),]
    
    
    return(out)
    
  }
  
  
  ## Parse mascot result files, produce nested list 
  ## Current parsing function
  ## Return sample-by-sample record
  findScore <- function(x) { 
    seqIdx <- which(names(x)  == "Sequence" )
    idx <- grep( "Score", names(x) )	
    mScore <- unlist(apply(x[,idx],1, myMax))
    mRec <- names(x)[idx[unlist(apply(x[,idx],1, myWhich.max))]]
    mRec <- unlist(lapply(strsplit(mRec,"\\."), function(x) x[2]))
    nam <- paste(".", unlist(lapply(strsplit(names(x)[idx], "\\."), 
                                    function(x) x[2])), sep = "")
    tmpList <- list()
    for( i in 1:length(nam) ){
      tmpList[[i]] <- x[,c(seqIdx, grep(nam[i], names(x), fixed = T))]	
      tmpList[[i]]$Sample <- nam[i]		
      
      idxC <- c(grep("Sample", names(tmpList[[i]])), 
                grep("Sequence", names(tmpList[[i]])),
                grep("Score", names(tmpList[[i]])), 
                grep("RT", names(tmpList[[i]])), 
                grep("Charge", names(tmpList[[i]])), 
                grep("m.z", names(tmpList[[i]])))		
      tmpList[[i]] <- tmpList[[i]][, idxC]
      
      names(tmpList[[i]]) <- c("sample", "compound", "score", "rt", "z", "mz")
      
      tmpList[[i]]$mScore <- mScore
      tmpList[[i]]$mRec <- mRec	
      tmpList[[i]]$modifications <- x[,grep("Modifications", names(x))]
      tmpList[[i]]$parent <- x[,grep("MH...Da.", names(x))]
      
      
      #grep("X..PSMs", names(x))
      #grep("X..Proteins", names(x))
      #grep("X..Protein.Groups", names(x))
      #grep("Protein.Group.Accessions", names(x))
      
    }
    
    out <- NULL
    for( i in 1:length(tmpList)){
      if( i == 3 ){		
      }
      out <- rbind(out, tmpList[[i]])
    }
    
    out$sample <- unlist(lapply(strsplit(out$sample, "\\."), function(x) x[2]))
    
    return(out)
  }
  
  
  ## Finds the best RT and Mz for a parsed Mascot input and creates MAVEN output
  makeMavenRTMZ <- function( inTable){
    
    
    outTab <- NULL
    
    subTab <- data.frame( 	
      compound = 'compound', 
      rt = 0.0,
      mz = 0.0,
      file = 'file',
      score = 0, 
      note = 'mod', stringsAsFactors = F)
    
    ## Find best record for each file
    for(i in unique(inTable$file)){
      fIdx <- which(inTable$file == i)
      tTab <- inTable[fIdx,] 
      
      ## Find best score for each file
      for (j in unique(tTab$compound)){
        cIdx <- which(tTab$compound == j)
        ttTab <- tTab[cIdx,]
        ttTab <- ttTab[which( ttTab$score == ttTab$mScore[1]),]
        
        subTab$compound <- ttTab$compound[1]
        subTab$rt <- mean(ttTab$rt)
        subTab$mz <- mean(ttTab$mz)
        subTab$file <- ttTab$file[1]
        subTab$score <- max(ttTab$mScore)
        subTab$note <- ttTab$modifications[1]
        
        outTab <- rbind(outTab, subTab)
        
      }
      
    }
    
    return( outTab )
  }
  
  
  ################################################################################
  ## Digestomics functions
  
  ## Function for matching peptides to proteome
  seqMap <- function(seq, pep){
    subPep <- paste(rep("#", nchar(pep)), sep = "", collapse = "")
    return(which(unlist(strsplit((gsub(pep, subPep, seq)), '')) == "#"))
  }
  
  
  
  ###############################################################################
  ## date calculator
  dateCalc <- function(start, end){
    require(chron)
    start <- strsplit(start, " ")
    end <- strsplit(end, " ")
    
    dtsS <- suppressWarnings(dates(sapply( start, 
                                           function(x) paste((unlist(strsplit(x[1], "-"))[c(2,3,1)]), 
                                                             collapse = "/"))))
    timeS <- suppressWarnings(times(sapply( start, function(x) x[2] )))
    dtsE <- suppressWarnings(dates(sapply( end, 
                                           function(x) paste((unlist(strsplit(x[1], "-"))[c(2,3,1)]), 
                                                             collapse = "/"))))
    timeE <- suppressWarnings(times(sapply( end, function(x) x[2] )))
    elapsed <- suppressWarnings(
      apply(cbind(dtsS, dtsE), 1, diff) + 
        (apply(cbind(timeS, timeE), 1, diff) ))	
    return(elapsed)
    
  }
  
  
  ## Find time to positive culture versus total stay
  clinicalStats <- function(inDat){
    myMin <- function( x ){ min( x, na.rm = TRUE )}
    
    uPID <- unique(inDat$PID)
    outList <- data.frame(
      PID = uPID, 
      Duration = rep(NA, length(uPID)), 
      CultureTime = rep(NA, length(uPID)))
    
    for( i in 1:length(uPID)){
      
      idx <- which(inDat$PID == uPID[i])
      subDat <- inDat[idx,]
      
      
      
      tmp <- try(c(uPID[i], 
                   myMin(dateCalc( subDat$ENCNTR_ADMIT_DTM, subDat$DSCHG_DTM )), 
                   myMin(dateCalc( subDat$CULT_START_DTM, subDat$VERIFY_DTM ))),
                 silent = T)
      if(class(tmp) == "try-error")
        next
      
      outList[i,] <- tmp
      
    }
    return(outList)
    
    
  }
  
  
  ## Aggregate stats
  
  summaryStats <- function( inDat1, inDat2, breaks = 500, pad = T, ...){
    
    
    ## Make a remove list
    idx <- rep(F, length(inDat1))
    
    if( any (is.na(inDat1) ))
      idx[ which(is.na(inDat1)) ] <- T
    if( any (inDat1 <= 0 ) )
      idx[ which(inDat1 <= 0) ] <- T
    if( any (is.na(inDat2) ))
      idx[ which(is.na(inDat2)) ] <- T
    if( any (inDat2 <= 0 ) )
      idx[ which(inDat2 <= 0 ) ] <- T
    
    inDat1 <- inDat1[!idx]
    inDat2 <- inDat2[!idx]
    
    idx2 <- order(inDat1)
    inDat1 <- inDat1[idx2]
    inDat2 <- inDat2[idx2]
    
    if( pad ){
      idx3 <- 1:length(inDat1)
      padF <- ceiling(.01 * length(idx3) )
      idx3 <- padF:(length(idx3)-padF)
      
      inDat1 <- inDat1[idx3]
      inDat2 <- inDat2[idx3]
      
    }
    
    
    if(breaks == 0)
      return(invisible(data.frame(inDat1, inDat2)))
    
    
    histDat <- hist(inDat1, breaks = breaks)
    hb <- histDat$breaks
    
    sLab <- rep(NA, length(inDat1))
    for(i in 1:length(inDat1)){
      
      if(is.na(inDat1[i]))
        sLab[i] <- NA
      else
        sLab[i] <- hb[which.min(abs(inDat1[i] - hb))]
      
    }
    
    
    
    outDat <- aggregate(inDat2, FUN = function(x) mean(x, na.rm = T), 
                        by = list(sLab) )
    
    plot(outDat, ...)
    
    return(invisible(outDat))
    
    
  }
  
  
  ## File Open GUI
  
  
  
  ## User function for ordering datasets
  #### UNDER CONSTRUCTION
  #groupSort <- function( inDat, listInput = FALSE ){
  #	
  #	
  #	
  #	if( !listInput ){
  #		tDat <- metaSep(inDat)	
  #		nameList <- names(tDat[[2]])
  #		outList <- data.frame( index = 1:length(nameList),
  #				sample = nameList, group = nameList, stringsAsFactors=FALSE)		
  #	}else{
  #		tDat <- inDat
  #		outList <- inDat[[2]]
  #		outList$index <- 1:nrow(outList)
  #	}
  #
  #	
  #
  #	editVec <- rep(TRUE, ncol(outList))
  #	editVec[1] <- FALSE
  #	
  #
  #	
  #	outList <- MZtableEdit( outList, editVec )
  #	
  #	return(tDat)
  #	
  #	outList <- list( 
  #			data = cbind(tDat[[1]], tDat[[2]][,outList$index], stringsAsFactors=FALSE), 
  #			groups = outList, stringsAsFactors=FALSE )
  #	
  #	return(outList)
  #	
  #}
  
  
  
  
  
  ## Interactive GUI for editing tables
  ## data - data.frame; the table to edit
  ## editable - logical vector; indicates whether entries in each column in the 
  ##	table should be editable, should be the same length as the number of columns
  ## title - character string; title for the GUI
  ## colVer - function list; functions, one for each column, used to verify the 
  ##	entries in each column.  Functions should return TRUE or FALSE.
  ## errMsgs - character vector; error messages to display if a function in colVer
  ##	returns FALSE, should be the same length as as the number of columns.  If 
  ##	NULL, no error messages are displayed
  MZtableEdit <- function(data, editable=rep(TRUE, ncol(data)),	title='rNMR', 
                          colVer=NULL, errMsgs=rep(paste('Data type for new entry must',  
                                                         'match the data type for a given column.'), ncol(data))){
    
    ##check colVer argument
    if (is.null(colVer)){
      verFun <- function(x) return(TRUE)
      colVer <- rep(list(verFun), ncol(data))
    }
    if (length(colVer) != ncol(data))
      stop('length of colVer argument must equal the number of columns in data')
    
    ##check errMsgs argument
    if (!is.null(errMsgs) && length(errMsgs) != ncol(data))
      stop('length of errMsgs argument must equal the number of columns in data')
    
    ##creates main window
    dlg <- tktoplevel()
    tcl('wm', 'attributes', dlg, topmost=TRUE)
    returnVal <- NULL
    tkwm.title(dlg, title)
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    colNames <- colnames(data)
    
    ##create tablelist widget
    tableFrame <- ttklabelframe(dlg, text='Data Table:')
    xscr <- ttkscrollbar(tableFrame, orient='horizontal', command=function(...) 
      tkxview(tableList, ...))
    yscr <- ttkscrollbar(tableFrame, orient='vertical', command=function(...) 
      tkyview(tableList, ...))
    colVals <- NULL
    for (i in colNames)
      colVals <- c(colVals, '0', i, 'center')
    tableList <- tkwidget(tableFrame, 'tablelist::tablelist', columns=colVals, 
                          activestyle='underline', height=11, width=110, exportselection=FALSE,
                          labelcommand='tablelist::sortByColumn', selectmode='extended', bg='white', 
                          spacing=3, stretch='all', editselectedonly=TRUE, selecttype='cell',
                          xscrollcommand=function(...) tkset(xscr, ...),
                          yscrollcommand=function(...) tkset(yscr, ...))
    
    ##add data to tablelist widget
    for (i in 1:nrow(data))
      tkinsert(tableList, 'end', unlist(data[i, ]))
    for (i in 1:ncol(data))
      tcl(tableList, 'columnconfigure', i - 1, sortmode='dictionary', 
          editable=editable[i], width=0, align='left')
    
    ##get the data types for each column
    colTypes <- NULL
    for (i in 1:ncol(data))
      colTypes <- c(colTypes, storage.mode(data[, i]))
    
    ##selects all rows Ctrl+A is pressed
    tkbind(tableList, '<Control-a>', function(...) 
      tkselect(tableList, 'set', 0, 'end'))
    
    ##rewrites data after GUI is updated
    writeData <- function(){
      
      ##get the data from the GUI
      newData <- NULL
      numRows <- as.numeric(tcl(tableList, 'index', 'end'))
      if (numRows == 0){
        data <<- newData
        return(invisible())
      }
      for (i in 0:numRows)
        newData <- rbind(newData, as.character(tcl(tableList, 'get', i)))
      
      ##format data
      colnames(newData) <- colNames
      newData <- as.data.frame(newData, stringsAsFactors=FALSE)
      for (i in 1:ncol(newData))
        suppressWarnings(storage.mode(newData[, i]) <- colTypes[i])
      data <<- newData
    }
    
    ##save tableList data after table is sorted
    tkbind(dlg, '<<TablelistColumnSorted>>', writeData)
    
    ##create top button
    optionFrame <- ttkframe(tableFrame)
    moveFrame <- ttklabelframe(optionFrame, text='Move selected rows', padding=6)
    onTop <- function(){
      usrSel <- as.numeric(tcl(tableList, 'curselection'))
      if (!length(usrSel) || usrSel == 0)
        return(invisible())
      tkselection.set(tableList, usrSel)
      for (i in seq_along(usrSel))
        tkmove(tableList, usrSel[i], 0 + i - 1)
      tcl(tableList, 'see', 0)
      writeData()
    }
    topButton <- ttkbutton(moveFrame, text='Top', width=11, command=onTop)
    
    ##create up button
    onUp <- function(){
      usrSel <- as.numeric(tcl(tableList, 'curselection'))
      if (!length(usrSel) || usrSel == 0)
        return(invisible())
      tkselection.set(tableList, usrSel)
      for (selItem in usrSel)
        tkmove(tableList, selItem, selItem - 1)
      tcl(tableList, 'see', min(usrSel) - 1)
      writeData()
    }
    upButton <- ttkbutton(moveFrame, text='^', width=9, command=onUp)
    
    ##create down button
    onDown <- function(){
      usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
      if (!length(usrSel) || usrSel == nrow(data) - 1)
        return(invisible())
      tkselection.set(tableList, usrSel)
      for (selItem in usrSel)
        tkmove(tableList, selItem, selItem + 2)
      tcl(tableList, 'see', max(usrSel) + 1)
      writeData()
    }
    downButton <- ttkbutton(moveFrame, text='v', width=9, command=onDown)
    
    ##create bottom button
    onBottom <- function(){
      usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
      if (!length(usrSel) || usrSel == nrow(data) - 1)
        return(invisible())
      tkselection.set(tableList, usrSel)
      for (i in seq_along(usrSel))
        tkmove(tableList, usrSel[i], nrow(data) - i + 1)
      tcl(tableList, 'see', nrow(data) - 1)
      writeData()
    }
    bottomButton <- ttkbutton(moveFrame, text='Bottom', width=11, 
                              command=onBottom)
    
    ##create sig. fig. spinbox
    sigFigFrame <- ttklabelframe(optionFrame, text='Display', padding=6)
    onSigFig <- function(){
      if (tclvalue(sigFigVal) == 'max'){
        for (i in 1:nrow(data))
          tcl(tableList, 'rowconfigure', i - 1, text=unlist(data[i, ]))
        return(invisible())
      }
      sigFig <- as.numeric(tclvalue(sigFigVal))
      for (i in seq_along(data[1, ])){
        if (any(is.logical(data[, i])))
          next
        newData <- tryCatch(signif(data[, i], sigFig), 
                            error=function(er) return(data[, i]))
        newData[is.na(newData)] <- 'NA'
        tcl(tableList, 'columnconfigure', i - 1, text=newData)
      }
    }
    sigFigVal <- tclVar('max')
    sigFigBox <- tkwidget(sigFigFrame, 'spinbox', width=6, wrap=TRUE,
                          textvariable=sigFigVal, values=c('max', 1:9), command=onSigFig)
    sigFigLabel <- ttklabel(sigFigFrame, text='significant figures')
    
    ##create table edit widgets
    if (any(editable)){
      
      ##check interactively edited cells using functions provided in colVer
      onEdit <- function(widget, rowNum, colNum, newVal, tclReturn=TRUE){
        rowNum <- as.numeric(rowNum) + 1
        colNum <- as.numeric(colNum) + 1
        if (newVal == 'NA'){
          if (tclReturn)
            return(tclVar(as.character(newVal)))
          else
            return(TRUE)
        }
        suppressWarnings(storage.mode(newVal) <- colTypes[colNum])
        if (!colVer[[colNum]](newVal)){
          if (!is.null(errMsgs))
            myMsg(errMsgs[colNum], icon='error', parent=dlg)
          if (tclReturn)
            tcl(tableList, 'cancelediting')
          else
            return(FALSE)
        }else
          data[rowNum, colNum] <<- newVal
        if (tclReturn)
          return(tclVar(as.character(newVal)))
        else
          return(TRUE)
      }
      tkconfigure(tableList, editendcommand=function(...) onEdit(...))
      
      ##create cell editing textbox
      ceditFrame <- ttklabelframe(optionFrame, text='Edit selected cells', 
                                  padding=6)
      usrEntry <- tclVar(character(0))
      textEntry <- ttkentry(ceditFrame, width=13, justify='center', 
                            textvariable=usrEntry)
      
      ##update cell editing textbox with current cell selection value
      onCellSel <- function(){
        usrSel <- as.character(tcl(tableList, 'curcellselection'))
        if (!length(usrSel))
          tclObj(usrEntry) <- character(0)
        selVals <- as.character(tcl(tableList, 'getcells', usrSel))
        if (length(grep(selVals[1], selVals, fixed=TRUE)) == length(selVals))
          tclObj(usrEntry) <- selVals[1]
        else
          tclObj(usrEntry) <- character(0)
      }
      tkbind(tableList, '<<TablelistSelect>>', onCellSel)
      
      ##create apply button
      onApply <- function(){
        tcl(tableList, 'finishediting')
        newVal <- tclvalue(usrEntry)
        usrSel <- as.character(tcl(tableList, 'curcellselection'))
        for (i in usrSel){
          rowNum <- unlist(strsplit(i, ','))[1]
          colNum <- unlist(strsplit(i, ','))[2]
          isValid <- onEdit(rowNum=rowNum, colNum=colNum, newVal=newVal, 
                            tclReturn=FALSE)
          if (isValid)
            tcl(tableList, 'cellconfigure', i, text=newVal)
          else
            return(invisible())
        }
        writeData()
      }
      applyButton <- ttkbutton(ceditFrame, text='Apply', width=8, command=onApply)
      
      ##create copy button
      reditFrame <- ttklabelframe(optionFrame, text='Edit rows', padding=6)
      clipboard <- NULL
      onCopy <- function(){
        usrSel <- as.numeric(tcl(tableList, 'curselection'))
        if (!length(usrSel) || usrSel == 0)
          return(invisible())
        tkselection.set(tableList, usrSel)
        selVals <- NULL
        for (i in usrSel)
          selVals <- rbind(selVals, as.character(tcl(tableList, 'get', i)))
        clipboard <<- selVals
      }
      copyButton <- ttkbutton(reditFrame, text='Copy', width=10, command=onCopy)
      
      ##create paste button
      onPaste <- function(){
        if (is.null(clipboard))
          return(invisible())
        for (i in 1:nrow(clipboard))
          tkinsert(tableList, 'end', unlist(clipboard[i, ]))
        writeData()
        tcl(tableList, 'see', nrow(data) - 1)
      }
      pasteButton <- ttkbutton(reditFrame, text='Paste', width=10, 
                               command=onPaste)
      
      ##create insert button
      onInsert <- function(){
        tkinsert(tableList, 'end', as.character(rep(NA, ncol(data))))
        writeData()
        tcl(tableList, 'see', nrow(data) - 1)
      }
      insertButton <- ttkbutton(reditFrame, text='Insert', width=10, 
                                command=onInsert)
      
      ##create delete button
      onDelete <- function(){
        usrSel <- as.numeric(tcl(tableList, 'curselection'))
        if (!length(usrSel))
          return(invisible())
        tkselection.set(tableList, usrSel - 1)
        tkdelete(tableList, usrSel)
        writeData()
      }
      deleteButton <- ttkbutton(reditFrame, text='Delete', width=10, 
                                command=onDelete)
    }
    
    ##create ok button
    bottomFrame <- ttkframe(dlg)
    onOk <- function(){
      
      ##verify data
      tcl(tableList, 'finishediting')
      if (!is.null(data)){
        for (i in 1:ncol(data)){
          suppressWarnings(storage.mode(data[, i]) <- colTypes[i])
          if (!colVer[[i]](data[, i])){
            if (!is.null(errMsgs))
              myMsg(errMsgs[i], icon='error', parent=dlg)
            return(invisible())
          }
        }
      }
      
      ##return the data and close the GUI
      returnVal <<- data
      tkgrab.release(dlg)
      tkdestroy(dlg)
      return(returnVal)
    }
    okButton <- ttkbutton(bottomFrame, text='OK', width=10, command=onOk)
    
    ##create cancel button
    onCancel <- function(){
      tkgrab.release(dlg)
      tkdestroy(dlg)
      return(returnVal)
    }
    cancelButton <- ttkbutton(bottomFrame, text='Cancel', width=10, 
                              command=onCancel)
    
    ##create export button
    onExport <- function(){
      tkwm.iconify(dlg)
      if ('ACTIVE' %in% names(data))
        initFile <- 'roiTable'
      else if ('Index' %in% names(data))
        initFile <- 'peakList'
      else
        initFile <- 'roiSummary'
      fileName <- mySave(initialfile=initFile, defaultextension='txt', 
                         title='Export', filetypes=list('xls'='Excel Files', 'txt'='Text Files'))
      if (length(fileName) == 0 || !nzchar(fileName)){
        tkwm.deiconify(dlg)
        return(invisible())
      }
      write.table(data, file=fileName, quote=FALSE, sep='\t', row.names=FALSE, 
                  col.names=TRUE)
      tkwm.deiconify(dlg)
    }
    exportButton <- ttkbutton(bottomFrame, text='Export', width=10, 
                              command=onExport)
    
    ##add widgets to treeFrame
    tkgrid(tableFrame, column=1, row=1, sticky='nswe', pady=6, padx=6)
    tkgrid(tableList, column=1, row=1, sticky='nswe')
    tkgrid(xscr, column=1, row=2, sticky='we')
    tkgrid(yscr, column=2, row=1, sticky='ns')
    
    ##make treeFrame stretch when window is resized
    tkgrid.columnconfigure(dlg, 1, weight=1)
    tkgrid.rowconfigure(dlg, 1, weight=1)
    tkgrid.columnconfigure(tableFrame, 1, weight=1)
    tkgrid.rowconfigure(tableFrame, 1, weight=1)
    
    ##add widgets to moveFrame
    tkgrid(optionFrame, column=1, columnspan=2, row=3, pady=8)
    tkgrid(moveFrame, column=1, row=1, padx=8)
    tkgrid(topButton, column=1, row=1, pady=2, padx=c(0, 4))
    tkgrid(upButton, column=2, row=1, pady=2, padx=1)
    tkgrid(downButton, column=3, row=1, padx=1, pady=2)
    tkgrid(bottomButton, column=4, row=1, pady=2, padx=c(4, 0))
    
    ##add widgets to sigFigFrame
    tkgrid(sigFigFrame, column=2, row=1, padx=8)
    tkgrid(sigFigBox, column=1, row=1, padx=c(4, 2), pady=c(2, 4))
    tkgrid(sigFigLabel, column=2, row=1, padx=c(0, 4), pady=c(2, 4))
    
    ##add editing widgets
    if (any(editable)){
      
      ##add widgets to rowFrame
      tkgrid(reditFrame, column=1, row=2, pady=4, padx=8)
      tkgrid(copyButton, column=1, row=1, padx=c(0, 2))
      tkgrid(pasteButton, column=2, row=1, padx=c(0, 8))
      tkgrid(insertButton, column=3, row=1, padx=c(0, 2))
      tkgrid(deleteButton, column=4, row=1, padx=c(0, 0))
      
      ##add widgets to ceditFrame
      tkgrid(ceditFrame, column=2, row=2, pady=4, padx=8)
      tkgrid(textEntry, column=1, row=1, padx=2)
      tkgrid(applyButton, column=3, row=1, padx=2)
    }
    
    ##add widgets to bottomFrame
    tkgrid(bottomFrame, column=1, row=2, pady=c(6, 0))
    tkgrid(okButton, column=1, row=1, padx=4)
    tkgrid(cancelButton, column=2, row=1, padx=4)
    tkgrid(exportButton, column=3, row=1, padx=c(20, 4))
    tkgrid(ttksizegrip(dlg), column=1, row=3, sticky='se')
    
    ##Allows users to press the 'Enter' key to make selections
    onEnter <- function(){
      focus <- as.character(tkfocus())
      if (length(grep('.2.2.3.1$', focus)))
        onApply()
      else
        tryCatch(tkinvoke(focus), error=function(er){})
    }
    tkbind(dlg, '<Return>', onEnter)
    
    ## configure the toplevel
    tkfocus(tableList)
    if (as.logical(tkwinfo('viewable', dlg)))
      tkgrab.set(dlg)
    tkwait.window(dlg)
    return(returnVal)
    
    invisible()
  }
  
  
  
  ## Junk temp code for finding isotope ratio
  isotopeRatio <- function (inDat, IR = TRUE){
    
    uCom <- unique (inDat$compound)
    
    outDat <- NULL
    for (i in 1:length(uCom) ){
      tRow <- inDat[inDat$compound == uCom[i], ]
      oRow <- tRow[1,]
      oRow$note <- "IsotopeRatio"
      
      if( IR ){
        oRow[,-(1:13)] <- 1-(tRow[tRow$note == "C12 PARENT",-(1:13)] / 
                               apply(tRow[,-(1:13)], 2, sum))
      }else{
        oRow[,-(1:13)] <- apply(tRow[,-(1:13)], 2, sum)
      }
      
      
      outDat <- rbind(outDat, oRow)
      
    }
    
    return(outDat)
    
    
  }