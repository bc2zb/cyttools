libraryList <- c("flowCore",
                 "limma",
                 "edgeR",
                 "FlowSOM",
                 "ncdfFlow",
                 "openCyto",
                 "tidyverse")

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

  flowSet <- as.flowSet(read.ncdfFlowSet(file, truncate_max_range = F)) # reads in files as flowSet, required for flowType
  
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

estParamFlowVS <- function (fs, channels) 
{
  nmatch = which(channels %in% colnames(fs))
  if (length(nmatch) != length(channels)) 
    stop(" At least one channel name is not present in the flowSet.")
  cofactors = NULL
  for (col in channels) {
    cat("====================================================================\n")
    cat("Channel ", col, " : ")
    cat("Finding optimum cofactor for asinh transformation\n")
    cat("====================================================================\n")
    fs1D = fs[, col]
    cf = optimStat(fs1D)
    cofactors = c(cofactors, cf)
    cat("\n Optimum cofactor for ", col, " : ", cf, "\n")
    cat("====================================================================\n\n")
  }
  return(cofactors)
}

optimstat <- function (fs1D, cfLow = -10, cfHigh = 3, MAX_BT = 10^9)
{
    if (cfLow >= cfHigh) {
        print("Warning: cfLow>=cfHigh, using default values")
        cfLow = -1
        cfHigh = 10
    }
    cf = cfLow:cfHigh
    ncf = length(cf)
    cfopt = rep(0, ncf - 1)
    btopt = rep(0, ncf - 1)
    cat(sprintf("%18s       %10s    %15s  %8s \n", "cf range", 
        "opt cf", "Bartlett's stat", "time"))
    cat("====================================================================\n")
    for (i in 1:(ncf - 1)) {
        ptm <- proc.time()
        tol = (exp(cf[i + 1]) - exp(cf[i]))/10
        opt = suppressWarnings(optimize(f = flowVS1D, interval = c(exp(cf[i]), 
            exp(cf[i + 1])), fs1D, tol = tol, plot = FALSE, MAX_BT = MAX_BT))
        btopt[i] = opt$objective
        cfopt[i] = opt$minimum
        cat(sprintf("[%9.2f, %9.2f ] %10.2f ", exp(cf[i]), exp(cf[i + 
            1]), opt$minimum))
        if (opt$objective == MAX_BT) {
            cat(sprintf("%10s ", " MAX (10^9)"))
        }
        else cat(sprintf("%15.2f ", opt$objective))
        cat(sprintf("%13.2f \n", (proc.time() - ptm)[1]))
    }
    minIdx = which.min(btopt)
    del = cfopt[minIdx]/10
    btLocal = rep(0, 11)
    btLocal[6] = btopt[minIdx]
    cfLocal = c(5:1, cfopt[minIdx], 1:5)
    cfLocal[1:5] = cfopt[minIdx] - 5:1 * del
    cfLocal[7:11] = cfopt[minIdx] + 1:5 * del
    for (i in c(1:5, 7:11)) {
        btLocal[i] = flowVS1D(cfLocal[i], fs1D)
    }
    minIdx = which.min(btLocal)
    plot(cfLocal, btLocal, type = "o", pch = 16, xlab = "Cofactors", 
        ylab = "Bartlett's statistics", main = paste("Optimum cofactor for ", 
            colnames(fs1D), " : ", format(round(cfLocal[minIdx], 
                2), nsmall = 2), sep = ""))
    points(cfLocal[minIdx], btLocal[minIdx], pch = 16, col = "red")
    return(cfLocal[minIdx])
}


optimStat <- function (fs1D, cfLow = -10, cfHigh = 3, MAX_BT = 10^9) 
{
  if (cfLow >= cfHigh) {
    print("Warning: cfLow>=cfHigh, using default values")
    cfLow = -1
    cfHigh = 10
  }
  cf = cfLow:cfHigh
  ncf = length(cf)
  cfopt = rep(0, ncf - 1)
  btopt = rep(0, ncf - 1)
  cat(sprintf("%18s       %10s    %15s  %8s \n", "cf range", 
              "opt cf", "Bartlett's stat", "time"))
  cat("====================================================================\n")
  for (i in 1:(ncf - 1)) {
    ptm <- proc.time()
    tol = (exp(cf[i + 1]) - exp(cf[i]))/10
    opt = suppressWarnings(optimize(f = flowVS1D,
                                    interval = c(exp(cf[i]), 
                                                 exp(cf[i + 1])),
                                    fs1D,
                                    tol = tol,
                                    plot = FALSE,
                                    MAX_BT = MAX_BT))
    btopt[i] = opt$objective
    cfopt[i] = opt$minimum
    cat(sprintf("[%9.2f, %9.2f ] %10.2f ", exp(cf[i]), exp(cf[i + 
                                                                1]), opt$minimum))
    if (opt$objective == MAX_BT) {
      cat(sprintf("%10s ", " MAX (10^9)"))
    }
    else cat(sprintf("%15.2f ", opt$objective))
    cat(sprintf("%13.2f \n", (proc.time() - ptm)[1]))
  }
  minIdx = which.min(btopt)
  del = cfopt[minIdx]/10
  btLocal = rep(0, 11)
  btLocal[6] = btopt[minIdx]
  cfLocal = c(5:1, cfopt[minIdx], 1:5)
  cfLocal[1:5] = cfopt[minIdx] - 5:1 * del
  cfLocal[7:11] = cfopt[minIdx] + 1:5 * del
  for (i in c(1:5, 7:11)) {
    btLocal[i] = flowVS1D(cfLocal[i], fs1D)
  }
  minIdx = which.min(btLocal)
  plot(cfLocal, btLocal, type = "o", pch = 16, xlab = "Cofactors", 
       ylab = "Bartlett's statistics", main = paste("Optimum cofactor for ", 
                                                    colnames(fs1D), " : ", format(round(cfLocal[minIdx], 
                                                                                        2), nsmall = 2), sep = ""))
  points(cfLocal[minIdx], btLocal[minIdx], pch = 16, col = "red")
  return(cfLocal[minIdx])
}

flowVS1D <- function (cofactor, fs1D, bwFac = 2, plot = FALSE, MAX_BT = 10^9, 
                      save.plot = FALSE, populationQuant = 0.01, borderQuant = 0.001, 
                      annotation = "", save.folder = "peaks", ...) 
{
  peakStats = matrix(nrow = 0, ncol = 5)
  colnames(peakStats) = c("n", "mean", "median", "variance", 
                          "sample")
  stain = colnames(fs1D)[1]
  for (i in 1:length(fs1D)) {
    y = exprs(fs1D[[i]])[, 1]
    y = asinh(y/cofactor)
    peakStat = densityPeaks(y, stain, bwFac = bwFac, plot = plot, 
                            save.plot = save.plot, populationQuant = populationQuant, 
                            borderQuant = borderQuant, annotation = paste(i), 
                            save.folder = "", ...)
    peakStat = cbind(peakStat, rep(i, nrow(peakStat)))
    peakStats = rbind(peakStats, peakStat)
  }
  sampleIdx = peakStats[, 5]
  nsamples = length(unique(sampleIdx))
  npeaks = table(sampleIdx)
  sampleOrder = order(abs(npeaks - median(npeaks)), decreasing = TRUE)
  rm = round(nsamples/10)
  if (rm > 0 && nsamples > (rm + 1)) {
    rmIdx = sampleOrder[1:rm]
    rmSamples = as.integer(names(npeaks[rmIdx]))
    selcted = !(sampleIdx %in% rmSamples)
    peakStats = peakStats[selcted, ]
  }
  pkVar = peakStats[, 4]
  npeaks = length(pkVar)
  varOrder = order(abs(pkVar - median(pkVar)), decreasing = TRUE)
  rm = round(npeaks/10)
  if (rm > 0 && npeaks > (rm + 1)) {
    rmIdx = varOrder[1:rm]
    peakStats = peakStats[-rmIdx, ]
  }
  if (nrow(peakStats) <= 1) 
    bt = MAX_BT
  else bt = bartlettTest(peakStats)
  return(bt)
}

densityPeaks <- function (y, stain, bwFac = 2, plot = FALSE, save.plot = FALSE, 
                          populationQuant = 0.01, borderQuant = 0.001, peak.density.thr = 0.05, 
                          peak.distance.thr = 0.05, annotation = "", save.folder = "", 
                          ...) 
{
  y = y[y > min(y) & y < max(y)]
  min.range = quantile(y, borderQuant)
  max.range = quantile(y, 1 - borderQuant)
  y = y[y > min.range & y < max.range]
  mat = as.matrix(y)
  colnames(mat) = stain
  ff = flowFrame(mat)
  curv1.filter <- filter(ff, curv1Filter(stain, bwFac = bwFac))
  bnds <- curvPeaks(curv1.filter, y)
  peaks = bnds$peaks
  peak.rm.idx = vector()
  mrange = max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
  npeaks = nrow(peaks)
  if (npeaks >= 2) {
    p1 = 1
    p2 = 2
    while (p2 <= npeaks) {
      peak1 = peaks[p1, ]
      peak2 = peaks[p2, ]
      if (abs(peak1[1] - peak2[1]) < mrange * peak.distance.thr) {
        peak.rm.idx = ifelse(peak1[2] < peak2[2], c(peak.rm.idx, 
                                                    p1), c(peak.rm.idx, p2))
        p1 = ifelse(peak.rm.idx == p2, p1, p2)
        p2 = p2 + 1
      }
      else {
        p1 = p2
        p2 = p2 + 1
      }
    }
  }
  if (length(peak.rm.idx) > 0) {
    peaks = peaks[-peak.rm.idx, ]
    if (!is(peaks, "matrix")) {
      peaks = matrix(peaks, ncol = 2)
      colnames(peaks) = colnames(bnds$peaks)
    }
  }
  dens <- density(y)
  peaksStats = matrix(nrow = 0, ncol = 4)
  colnames(peaksStats) = c("n", "mean", "median", "variance")
  if (nrow(peaks) == 0) 
    return(peaksStats)
  left = min(y, na.rm = TRUE)
  bounds = vector()
  separators = vector()
  populations = list()
  for (i in 1:nrow(peaks)) {
    p = peaks[i, ]
    if (i == nrow(peaks)) {
      right = max(y, na.rm = TRUE)
    }
    else {
      sel <- (dens$x > peaks[i, 1] & dens$x < peaks[i + 
                                                      1, 1])
      right <- dens$x[sel][which.min(dens$y[sel])]
    }
    y.sel = y[y >= left & y <= right]
    min.range = quantile(y.sel, populationQuant)
    max.range = quantile(y.sel, 1 - populationQuant)
    if (!is.na(p[1])) {
      r = (max.range - p[1])/(p[1] - min.range)
      if (r >= 3 && i == nrow(peaks)) {
        right = min(max.range, p[1] + 3 * (p[1] - min.range))
        y.sel = y[y >= left & y <= right]
        min.range = quantile(y.sel, populationQuant)
        max.range = quantile(y.sel, 1 - populationQuant)
      }
      else if (r <= 1/3 && i == 1) {
        left = max(min.range, p[1] - 3 * (max.range - 
                                            p[1]))
        y.sel = y[y >= left & y <= right]
        min.range = quantile(y.sel, populationQuant)
        max.range = quantile(y.sel, 1 - populationQuant)
      }
      population = y.sel[y.sel >= min.range & y.sel <= 
                           max.range]
      if (length(population) != 0) {
        peaksStats = rbind(peaksStats, c(length(population), 
                                         mean(population), median(population), var(population)))
        populations[[length(populations) + 1]] = population
        bounds = c(bounds, min.range, max.range)
        if (i < nrow(peaks)) 
          separators = c(separators, right)
      }
      left = right
      if (r >= 5 && i == nrow(peaks) && length(separators) == 
          0) {
        rare = y[y >= left]
        min.range = quantile(rare, populationQuant)
        max.range = quantile(rare, 1 - populationQuant)
        rare = rare[rare >= min.range & rare <= max.range]
        if (length(rare) != 0) {
          populations[[length(populations) + 1]] = rare
          separators = c(separators, left)
          bounds = c(bounds, min.range, max.range)
        }
      }
    }
  }
  if (plot && nrow(peaksStats) > 0) {
    main.text = paste("Peaks for", stain)
    if (annotation != "") 
      main.text = paste("Sample=", annotation, " : ", main.text, 
                        sep = "")
    plot(dens, main = main.text, cex.main = 1, ...)
    ramp <- colorRamp(c("blue", "green", "red"))
    cols = rgb(ramp(seq(0, 1, length = length(bounds)/2)), 
               alpha = 50, maxColorValue = 255)
    for (i in seq(1, length(bounds), 2)) {
      sel <- dens$x >= bounds[i] & dens$x <= bounds[i + 
                                                      1]
      polygon(c(dens$x[min(which(sel))], dens$x[sel], dens$x[max(which(sel))]), 
              c(0, dens$y[sel], 0), col = cols[(i + 1)/2], 
              border = NA)
    }
    lines(dens, ...)
    for (sep in separators) abline(v = sep, col = 2, lwd = 2)
    if (save.plot == TRUE && annotation != "") {
      file.name = paste(stain, "_", annotation, ".pdf", 
                        sep = "")
      if (save.folder != "") {
        dir.create(save.folder, showWarnings = FALSE)
        file.name = paste(save.folder, file.name, sep = "/")
      }
      dev.copy2pdf(file = file.name)
    }
  }
  return(peaksStats)
}

curvPeaks <- function (x, dat, borderQuant = 0.01, n = 201, from, to, densities = NULL) 
{
  if (missing(from)) 
    from <- min(dat)
  if (missing(to)) 
    to <- max(dat)
  bound <- attr(x@subSet, "boundaries")
  dens <- if (!is.null(densities)) 
    densities
  else density(dat, n = n, from = from, to = to, na.rm = TRUE)$y
  regPoints <- list()
  peaks <- midpoints <- regions <- densFuns <- NULL
  i <- 1
  if (!all(is.na(bound[[1]]))) {
    for (b in bound) {
      if (b[2] > quantile(c(from, to), borderQuant) && 
          b[1] < quantile(c(from, to), 1 - borderQuant)) {
        afun <- approxfun(seq(from, to, len = n), dens)
        sel <- seq(b[1], b[2], len = 50)
        regPoints[[i]] <- cbind(x = sel, y = afun(sel))
        m <- optimize(afun, b, maximum = TRUE)
        peaks <- rbind(peaks, cbind(x = m$maximum, y = m$objective))
        regions <- rbind(regions, cbind(left = b[1], 
                                        right = b[2]))
        midpoints <- c(midpoints, mean(b))
        densFuns <- c(densFuns, afun)
        i <- i + 1
      }
    }
  }
  if (i == 1) 
    return(list(peaks = cbind(x = NA, y = NA), regions = cbind(left = NA, 
                                                               right = NA), midpoints = NA, regPoints = list(cbind(x = NA, 
                                                                                                                   y = NA)), densFuns = NA))
  return(list(peaks = peaks, regions = regions, midpoints = midpoints, 
              regPoints = regPoints, densFuns = densFuns))
}

bartlettTest <- function (peakStats) 
{
  peakStats = peakStats[peakStats[, 1] > 1, , drop = FALSE]
  size = peakStats[, 1]
  variance = peakStats[, 4]
  sum1 = sum((size - 1) * log(variance))
  sum2 = sum((size - 1) * variance)
  N = sum(size)
  sum3 = sum(1/(size - 1))
  k = nrow(peakStats)
  bt = ((N - k) * log(sum2/(N - k)) - sum1)/(1 + (sum3 - 1/(N - 
                                                              k))/(3 * (k - 1)))
  return(bt)
}

preprocess_frame <- function(flow_frame){
  # gate out debris
  nondebris_gate <- gate_mindensity(flow_frame, "FSC-A")
  # split flow_frame into positive and negative frames
  nondebris_set <- split(flow_frame,
                         flowCore::filter(flow_frame, nondebris_gate))
  # beads and cells have different responses to this gate, grab whichever has the higher number of events
  num_pos_events <- nrow(nondebris_set$`+`)
  num_neg_events <- nrow(nondebris_set$`-`)
  if(num_pos_events > num_neg_events){
    nondebris_frame <- nondebris_set$`+`
  }else{
    nondebris_frame <- nondebris_set$`-`
  }
  # gate singlets
  singlets_gate <- gate_singlet(nondebris_frame, maxit = 20, wider_gate = T)
  singlets_set <- split(nondebris_frame,
                        flowCore::filter(nondebris_frame, singlets_gate))
  # retrive singlets only
  singlets_frame <- singlets_set$`singlet+`
  # perform ellipsoid gating to tighten up population
  ellipsoid_gate <- gate_flowClust_2d(singlets_frame,
                                      xChannel = "FSC-A",
                                      yChannel = "SSC-A")
  ellipsoid_set <- split(singlets_frame,
                         flowCore::filter(singlets_frame, ellipsoid_gate))
  # retrieve events from within ellipsoid
  preprocessed_frame <- ellipsoid_set$`defaultEllipsoidGate+`
  return(preprocessed_frame)
}
safe_preprocess_frame <- safely(preprocess_frame)

find_peak_emission_detector <- function(preprocess_frame){
  channel_info <- data.frame()
  for( j in colnames(preprocess_frame[,grep("SSC|FSC|Time|\\-H", colnames(preprocess_frame), invert = T)]
  )){
    g <- gate_mindensity(preprocess_frame[,grep("SSC|FSC|Time|\\-H", colnames(preprocess_frame), invert = T)],
                         channel = j)
    peak_frame <- split(preprocess_frame[,grep("SSC|FSC|Time|\\-H", colnames(preprocess_frame), invert = T)],
                        flowCore::filter(preprocess_frame[,grep("SSC|FSC|Time|\\-H", colnames(preprocess_frame), invert = T)],
                                         g))[[1]]
    num_events <- nrow(peak_frame)
    median_intensity <- median(exprs(peak_frame)[,colnames(peak_frame) %in% j])
    channel_info <- channel_info %>% 
      bind_rows(data.frame("channel" = j,
                           "num_events" = num_events,
                           "median_intensity" = median_intensity))
  }

  channel_info <- channel_info %>%
    dplyr::filter(num_events > 100)
  peak_channel <- channel_info$channel[which.max(channel_info$median_intensity)]
  g <- gate_mindensity(preprocess_frame, channel = peak_channel)
  peak_frame <- split(preprocess_frame, flowCore::filter(preprocess_frame, g))[[1]]
  # extract median intesities for all other channels
  return(apply(exprs(peak_frame), 2, median))
}
