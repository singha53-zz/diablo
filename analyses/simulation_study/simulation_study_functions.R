################################################################################
#
# simulation_study_functions.R
# August 22, 2017
#
################################################################################

library(mixOmics)
## parameters to tune (fold-change, noise)
## types of variables (p_relevant-corDis, p_relevant-NonCorDis, p_irrelevant-corNonDis, p_irrelevant-NonCorNonDis)

simData = function(fc, noise, J, n, p_relevant, p_irrelevant){
  ## load libraries
  library(tidyverse)
  library(amritr)
  library(mvtnorm)
  
  ## Generate variates that contribute to the full design
  rho = 1
  sigma = matrix(rho, J, J); diag(sigma) = 1
  ## Variate for Group 1
  bComp1 <- as.data.frame(rmvnorm(n, rep(-fc/2, J), sigma))
  ## Variate for Group 2
  bComp2 <- as.data.frame(rmvnorm(n, rep(fc/2, J), sigma))
  bComp.full <- rbind(bComp1, bComp2) 
  bComp_full_irrelevant <- as.data.frame(rmvnorm(2*n, rep(0, J), sigma))
  colnames(bComp.full) <- colnames(bComp_full_irrelevant) <- paste("Dataset", 1:J, sep="_")
  
  ## Generate variates that contribute to the null design
  rho = 0
  sigma = matrix(rho, J, J); diag(sigma) = 1
  ## Variate for Group 1
  bComp1 <- as.data.frame(rmvnorm(n, rep(-fc/2, J), sigma))
  ## Variate for Group 2
  bComp2 <- as.data.frame(rmvnorm(n, rep(fc/2, J), sigma))
  bComp.null <- rbind(bComp1, bComp2)  #+ matrix(rnorm(2*n*J, mean = 0, sd = noise), nc = J)
  bComp_null_irrelevant <- as.data.frame(rmvnorm(2*n, rep(0, J), sigma))
  colnames(bComp.null) <- colnames(bComp_null_irrelevant) <- paste("Dataset", 1:J, sep="_")
  
  ## Generate 3 datasets (first p_relevant variables contribute to the full design, 
  ## the next to the null design)
  allLoadings <- c(runif(1000, min = -0.3, max = -0.2), runif(1000, min = 0.2, max = 0.3))  # loadings 
  X_relevant <- lapply(1 : J, function(i){
    w.full <- sample(allLoadings, p_relevant); w.full <- w.full/sqrt(sum(w.full^2));
    fullDat <- as.matrix(bComp.full[, i, drop = FALSE]) %*% matrix(w.full, nrow = 1)
    w.null <- sample(allLoadings, p_relevant); w.null <- w.null/sqrt(sum(w.null^2));
    nullDat <- as.matrix(bComp.null[, i, drop = FALSE]) %*% matrix(w.null, nrow = 1)
    cbind(fullDat, nullDat)
  })
  
  ## Add p_irrelevant correlative and non-correlative variables that are not discriminatory
  ## correlative variables
  allLoadings <- c(runif(1000, min = -0.3, max = -0.2), runif(1000, min = 0.2, max = 0.3))  # loadings 
  X_irrelevant <- lapply(1 : J, function(i){
    w.full <- sample(allLoadings, p_irrelevant); w.full <- w.full/sqrt(sum(w.full^2));
    corDat <- as.matrix(bComp_full_irrelevant[, i, drop = FALSE]) %*% matrix(w.full, nrow = 1)
    w.null <- sample(allLoadings, p_irrelevant); w.null <- w.null/sqrt(sum(w.null^2));
    noncorDat <- as.matrix(bComp_null_irrelevant[, i, drop = FALSE]) %*% matrix(w.null, nrow = 1)
    cbind(corDat, noncorDat)
  })
  
  ## final dataset
  X2 <- lapply(1 : J, function(i){
    cbind(X_relevant[[i]], X_irrelevant[[i]])
  })
  names(X2) <- c("Dataset1", "Dataset2", "Dataset3")
  
  ## add noise
  data <- lapply(1 : J, function(i){
    i2 <- X2[[i]] + matrix(rnorm(length(X2[[i]]), mean = 0, sd = noise), nr = nrow(X2[[i]]))
    colnames(i2) <- paste(paste(rep(c("corDis", "unCorDis", "corNonDis", "unCorNonDis"), 
      c(p_relevant,p_relevant,p_irrelevant,p_irrelevant)), 
      c(1:p_relevant, 1:p_relevant, 1:p_irrelevant, 1:p_irrelevant), sep = "."), "Dat", i, sep = "_")
    rownames(i2) <- paste0("Subj", 1:nrow(i2))
    i2
  })
  names(data) <- paste0("Dataset", 1:J)
  
  ## response variables
  Y <- rep(c("group1", "group2"), each = n)
  names(Y) <- paste0("Subj", 1:nrow(data[[1]]))
  
  return(list(data=data, Y=Y, 
    bComp.full=bComp.full, bComp_full_irrelevant=bComp_full_irrelevant,
    bComp.null=bComp.null, bComp_null_irrelevant=bComp_null_irrelevant))
}

pairPlot = function(mat, group){
  opar = par()[!names(par()) %in% c("cin", "cra", "csi", "cxy", 
    "din", "page")]
  VarX = mat
  datNames = colnames(VarX)
  if (ncol(VarX) <= 2) 
    stop("This function is only available when there are more than 3 blocks")
  Y = group
  legend.ncol = min(5, nlevels(Y))
  numberOfCols = ncol(VarX)
  numberOfRows = numberOfCols
  mat = matrix(0, nrow = numberOfRows, ncol = numberOfRows)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) mat[i, j] = paste(i, j, sep = "_")
  }
  plotType = list(cor = mat[lower.tri(mat)], scatter = mat[upper.tri(mat)], 
    lab = diag(mat))
  par(mfrow = c(numberOfRows + 1, numberOfCols), mar = rep.int(1/2, 
    4), oma = c(2, 2, 2, 2))
  layout(matrix(c(1:(numberOfCols)^2, rep((numberOfCols)^2 + 
      1, numberOfCols)), numberOfRows + 1, numberOfCols, byrow = TRUE), 
    heights = c(rep(1, numberOfRows), 0.25 * floor(nlevels(Y)/legend.ncol)))
  for (i in 1:numberOfRows) {
    for (j in 1:numberOfCols) {
      ptype = unlist(lapply(plotType, function(x) {
        intersect(paste(i, j, sep = "_"), x)
      }))
      splotMatPlot(x = VarX[, i], y = VarX[, j], datNames, 
        Y, ptype)
      if (i == 1 & j %in% seq(2, numberOfRows, 1)) 
        Axis(side = 3, x = VarX[, i])
      if (j == numberOfRows & i %in% seq(1, numberOfRows - 
          1, 1)) 
        Axis(side = 4, x = VarX[, i])
    }
  }
  plot(1:3, 1:3, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = levels(Y), col = color.mixo(1:nlevels(Y)), 
    pch = 19, ncol = legend.ncol, cex = 1.5)
  par(opar)
}

splotMatPlot = function (x, y, datNames, Y, ptype){
  if (names(ptype) == "cor") {
    plot(1, type = "n", axes = FALSE)
    r = round(cor(x, y), 2)
    text(1, 1, labels = r, cex = 0.6/strwidth(abs(r)) * abs(r))
    box()
  }
  if (names(ptype) == "scatter") 
    panel.ellipses(x = x, y = y, Y = Y)
  if (names(ptype) == "lab") {
    plot(1, type = "n", axes = FALSE)
    ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), 
      unique)))
    text(x = 1, y = 1, labels = datNames[ind], cex = 2)
    box()
  }
  if (FALSE) {
    if (names(ptype) == "bar") {
      Y2 = factor(as.character(Y), levels = groupOrder)
      boxplot(x ~ Y2, horizontal = TRUE, axes = FALSE, 
        ylim = c(min(x) - 3, max(x)), col = color.mixo(match(levels(Y2), 
          levels(Y))))
      axis(4, at = 1:nlevels(Y2), labels = levels(Y2))
    }
    if (names(ptype) == "stackedbar") {
      Y2 = factor(as.character(Y), levels = groupOrder)
      bars = table(Y2)
      barplot(bars, col = color.mixo(match(levels(Y2), 
        levels(Y))), axes = FALSE)
      axis(4, at = seq(0, max(bars), length.out = 5), labels = seq(0, 
        max(bars), length.out = 5))
    }
  }
}

panel.ellipses = function (x, y, Y = Y, pch = par("pch"), col.lm = "red", axes = FALSE, ...){
  library(ellipse)
  ind.gp = matrice = cdg = variance = list()
  for (i in 1:nlevels(Y)) ind.gp[[i]] = which(as.numeric(Y) == 
      i)
  matrice = lapply(ind.gp, function(z) {
    matrix(c(x[z], y[z]), ncol = 2)
  })
  cdg = lapply(matrice, colMeans)
  variance = lapply(matrice, var)
  coord.ellipse = lapply(1:nlevels(Y), function(x) {
    ellipse(variance[[x]], centre = cdg[[x]], ellipse.level = 0.95)
  })
  max.ellipse = sapply(coord.ellipse, function(x) {
    apply(x, 2, max)
  })
  min.ellipse = sapply(coord.ellipse, function(x) {
    apply(x, 2, min)
  })
  ind.names = names(Y)
  cex = 0.5
  plot(x, y, xlab = "X.label", ylab = "Y.label", col = color.mixo(as.numeric(Y)), 
    pch = 20, axes = axes, xlim = c(min(x, min.ellipse[1, 
      ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, 
        ]), max(y, max.ellipse[2, ])))
  box()
  for (z in 1:nlevels(Y)) points(coord.ellipse[[z]], type = "l", 
    col = color.mixo(c(1:nlevels(Y))[z]))
}


## Cross-validation function for ensemble classifier
perf.ensemble.splsda = function(data, Y, keepX, ncomp, nrepeat, M){
  BER <- list()
  for(x in 1:nrepeat){
    #stratified subsamplings of Y
    folds = amritr::createFolds(Y, k = M) 
    
    predLabel <- list()
    for(j in 1: length(folds)){
      omit = folds[[j]]
      data.train <- lapply(data, function(i){
        i[-omit, , drop = FALSE]
      })
      data.test <- lapply(data, function(i){
        i[omit, , drop = FALSE]
      })
      Y.train = Y[-omit]
      Y.test = Y[omit]
      
      pred <- lapply(1 : J, function(i){
        result <- splsda(X = data.train[[i]], Y = Y.train, keepX = keepX, ncomp = ncomp)
        predict(result, data.test[[i]])$class$max.dist
      })
      predLabel[[j]] <- do.call(cbind, pred)
    }
    
    preds0 <- do.call(rbind, predLabel)
    preds <- apply(preds0, 1, function(x){a=table(x); if (length(which(a==max(a)))==1) {b=names(which.max(a))}else{b=NA}; b})
    mat0 <- table(pred = preds, truth = Y[names(preds)]); mat <- mat0; diag(mat0) <- 0;  
    BER[[x]] <- mean(colSums(mat0)/colSums(mat))
  }
  data.frame(Mean = mean(unlist(BER)), SD = sd(unlist(BER)))
}
