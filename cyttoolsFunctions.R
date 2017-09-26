libraryList <- c("reshape2","parallel", "devtools", "stringr", "pvclust",
                 "flowCore", "limma", "edgeR", "FlowSOM", "flowType",
                 "cytofCore", "tidyverse")

lapply(libraryList, require, quietly = T, character.only = TRUE)

# performs arcsinH transform on flowFrame, borrowed from FlowVS
transFlowVS <- function (fs, channels, cofactors){   
    nmatch = which(channels %in% colnames(fs))
    if (length(nmatch) != length(channels)) 
      stop(" At least one channel name is not present in the flowSet.")
    cat("Transforming flowSet by asinh function with supplied cofactors.\n")
    fs1 = fs[1:length(fs)]
    for (i in 1:length(fs)) {
      mat = as.matrix(exprs(fs[[i]]))
      for (j in 1:length(channels)) {
        mat[, channels[j]] = asinh(mat[, channels[j]]/cofactors[j])
      }
    fs1[[i]] = flowFrame(mat)
    }
  return(fs1)
}

checkDesignCols <- function(targets, colsToCheck){
                      checkCols <- grep(paste(colsToCheck, collapse = "|"),
                                        colnames(targets))
                      if(length(checkCols) != length(colsToCheck)){
                        return(TRUE)
                      }else{
                        return(FALSE)
                      }
}

NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

read.flowSet.transVS <- function(targets, file){
  parameterIndex <- which(targets$Lineage == 1 | targets$Functional == 1)
  channels <- as.character(targets$name[parameterIndex])
  cofactors <- targets$TransformCofactor[parameterIndex]

  flowSet <- read.flowSet(file, transformation = F) # reads in files as flowSet, required for flowType

  flowSet.trans <- transFlowVS(flowSet,
                               channels,
                               cofactors)
  return(flowSet.trans)
}

fs.as.list <- function(flowSet.trans){
  flowSet.list <- fsApply(flowSet.trans, function(x){return(x)}, simplify = F)
  return(flowSet.list)
}

logitTransform <- function(p){ 
  log(p/(1-p)) 
}

checkLmFit <- function(fit, propData){
  numUnEst <- rowSums(is.na(fit$coefficients)) 
  numUnEst <- sum(numUnEst > 0 & numUnEst < NCOL(fit$coefficients))
  if(numUnEst == nrow(propData)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

recoderFunc <- function(data, oldvalue, newvalue) {
  
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  return(newvec)
  
}

# order one vector by the contents of another vector
orderVectorByVector <- function(x,y){
  orderOrder <- vector()
  for (i in 1:length(y)){
    index <- which(x == as.character(y[i]))
    orderOrder <- c(orderOrder, index)
  }
  return(orderOrder)
}

orderVectorByListOfTerms <- function(x,y){
  y <- rev(y)
  finalVector <- c(1:length(x))
  for (i in 1:length(y)){
    index <- grep(as.character(y[i]), x[finalVector])
    other <- grep(as.character(y[i]), x[finalVector], invert = T)
    nextIndex <- c(index, other)
    finalVector <- finalVector[nextIndex]
  }
  return(finalVector)
}