###############################################################################
# DIABLO fit
diabloFit = function(keepX, ncomp, X_train, X_test, y_train, y_test, design, folds, iter, modName){
  diabloMod <- block.splsda(X = X_train, Y = y_train, ncomp = ncomp, keepX = keepX, scheme = "centroid", design = design)  ## model fitting
  cv <- perf(diabloMod, validation = "Mfold", folds = folds, nrepeat = iter)  ## cross-validation
  
  ## Training error
  train <- data.frame(mean_err = cv$MajorityVote.error.rate$centroids.dist["Overall.BER", ncomp],
             sd_err = cv$MajorityVote.error.rate.sd$centroids.dist["Overall.BER", ncomp]) %>% 
    mutate(cohort = "train", modName = modName)
  ## Test error
  pred <- predict(diabloMod, X_test, method="centroids")
  mat <- table(factor(pred$WeightedVote$centroids.dist[, ncomp], levels(y_test)), y_test)
  test <- data.frame(mean_err = mean((colSums(mat)-diag(mat))/colSums(mat)), sd_err = NA) %>% 
    mutate(cohort = "test", modName = modName)
  
  ## summary
  error <- rbind(train, test)
  panels <- data.frame(panelLength = sapply(keepX, sum),
    dataset = names(sapply(keepX, sum)),
    modName = modName)
  
  return(list(error=error, panels=panels))
}

###############################################################################
# Concatenation_sPLSDA
concatSplsdaFit = function(keepX, ncomp, X_train, X_test, y_train, y_test, folds, iter, modName){
  splsdaMod <- splsda(X = X_train, Y = y_train, ncomp = ncomp, keepX = keepX)  ## model fitting
  cv <- perf(splsdaMod, validation = "Mfold", folds = folds, nrepeat = iter)  ## cross-validation

  ## Training error
  train <- data.frame(mean_err = cv$error.rate$BER[ncomp, "centroids.dist"],
    sd_err = cv$error.rate.sd$BER[ncomp, "centroids.dist"]) %>% 
    mutate(cohort = "train", modName = modName)
  ## Test error
  pred <- predict(splsdaMod, X_test, method="centroids")
  mat <- table(factor(pred$class$centroids.dist[, ncomp], levels(y_test)), y_test)
  test <- data.frame(mean_err = mean((colSums(mat)-diag(mat))/colSums(mat)), sd_err = NA) %>% 
    mutate(cohort = "test", modName = modName)
  
  ## summary
  error <- rbind(train, test)
  features <- table(sapply(strsplit(unlist(lapply(1 : ncomp, function(i){selectVar(splsdaMod, ncomp = i)$name})), "_"), function(i) i[1]))
  panels <- data.frame(panelLength = sum(features),
    dataset = names(features),
    modName = modName)
  
  return(list(error=error, panels=panels))
}

###############################################################################
## Ensemble_sPLSDA
ensembleSplsdaFit = function(keepX, ncomp, X_train, X_test, y_train, y_test, folds, iter, modName){
  ## Estimate training error
  classes_ensemble_splsda <- mapply(function(keepX, Xtrain, dataset){
    result <- splsda(X = Xtrain, Y = y_train, keepX = keepX, ncomp = ncomp)
    cv <- perf(result, validation = "Mfold", folds = folds, nrepeat = iter, progressBar = FALSE)
    do.call(rbind, lapply(1 : ncomp, function(i) cv$predict[[ncomp]][,, i])) %>%
      as.data.frame() %>%
      mutate(sample = rownames(.),
             Dataset = dataset,
             Iteration = rep(1:ncomp, each = length(y_train)))
  }, keepX = keepX, Xtrain = X_train[names(X.test)], dataset = names(X.test), SIMPLIFY = FALSE)
  ### Average predictions across datasets
  train <- do.call(rbind, classes_ensemble_splsda) %>%
    group_by(sample, Iteration) %>%
    dplyr::select(-c(Dataset)) %>%
    summarise_all(funs(mean)) %>%
    gather(pam50, pb, -c(sample:Iteration)) %>%
    group_by(sample, Iteration) %>%
    filter(pb == max(pb)) %>%
    dplyr::select(Iteration, sample, pam50) %>%
    group_by(Iteration) %>%
    nest() %>%
    mutate(err = purrr::map(data, ~{
      mat <- table(.$pam50, y_train[.$sample])
      mean((colSums(mat) - diag(mat))/colSums(mat))
    })) %>%
    unnest(err) %>%
    dplyr::select(err) %>%
    summarise(mean_err = mean(err), sd_err = sd(err), cohort = "train", modName = modName)

  ## Estimate training error
  test_classes_ensemble_splsda <- mapply(function(keepX, Xtrain, Xtest, dataset){
    result <- splsda(X = Xtrain, Y = y_train, keepX = keepX, ncomp = 3)
    pred <- predict(result, Xtest, "all")
    pred$predict[,, paste("dim", ncomp)] %>%
      as.data.frame() %>%
      mutate(sample = rownames(.), Dataset = dataset) %>%
      gather(pam50, pb, -c(sample, Dataset))
  }, keepX = keepX, Xtrain = X_train[names(X.test)], Xtest = X_test, dataset = names(X.test), SIMPLIFY = FALSE)
  ### Average predictions across datasets
  test_classes_ensemble_splsda <- do.call(rbind, test_classes_ensemble_splsda) %>%
    group_by(sample, pam50) %>%
    summarise(pb = mean(pb)) %>%
    group_by(sample) %>%
    filter(pb == max(pb))
  mat <- table(test_classes_ensemble_splsda$pam50, y_test[test_classes_ensemble_splsda$sample])
  test <- data.frame(mean_err = mean((colSums(mat) - diag(mat))/colSums(mat)),
    sd_err = NA, cohort = "test", modName = modName)
  
  error <- rbind(train, test)
  panels = data.frame(panelLength = sapply(keepX, sum), dataset = names(sapply(keepX, sum)), modName = modName)

  return(list(error=error, panels=panels))
}

###############################################################################
## Concatenation_Enet
concatEnetFit  = function(X_train, X_test, y_train, y_test, penalties, folds, iter, modName){
  computeError = function(pred, Y.test){
    mat <- table(factor(pred, levels(Y.test)), Y.test)
    error <- mean((colSums(mat) - diag(mat))/colSums(mat))
    names(error) <- c("MajVote")
    error
  }
  
  ## Run caret (glmnet) to predicted probablities for each dataset
  library(caret); library(glmnet);
  ctrl <- trainControl(method = "repeatedcv",
    number = folds,
    repeats = iter,
    classProbs = TRUE,
    savePredictions = TRUE)
  set.seed(123)
  ctrl$index <- createMultiFolds(y_train, k = folds, times = iter)
  
  ## fit training data
  enetGrid = expand.grid(alpha = 1, lambda = penalties)
  mod <- train(x=X_train, y=y_train,
    # preProc=c("center", "scale"),
    method = "glmnet",
    type.multinomial = "grouped",
    tuneGrid = enetGrid,
    trControl = ctrl)
  
  ## extract probabilities 
  prob <- mod$pred %>% 
    mutate(Resample = sapply(strsplit(Resample, "\\."), function(i) i[2])) %>% 
    dplyr::select(-c(Basal:LumB))
  
  lambdas <- prob %>% 
    group_by(lambda, Resample) %>% 
    nest() %>% 
    mutate(error = purrr::map(data, ~{
      confMat <- table(factor(.$pred, levels(.$obs)), .$obs)
      mean((colSums(confMat)-diag(confMat))/colSums(confMat))
    })) %>% 
    unnest(error) %>% 
    group_by(lambda) %>% 
    dplyr::summarise(mean = mean(error), sd = sd(error)) %>% 
    filter(mean == min(mean)) %>% 
    rename(dat=lambda)
  
  ## test performance
  alpha = 1
  lambda = as.numeric(lambdas[, "dat"])
  mod = glmnet(X_train, y_train, alpha = alpha, lambda = lambda, family="multinomial", type.multinomial = "grouped")
  predictors <- names(which(coef(mod)[[1]][,1] != 0))[-1]
  pred <- predict(mod, X_test, type = "class")
  
  # compute error and panel lengths
  testErr <- computeError(pred, y_test)
  panelLength <- table(sapply(strsplit(predictors, "_"), function(i) i[1]))
  
  error  <-  data.frame(mean_err = c(as.numeric(lambdas[, "mean"]), testErr),
    sd_err = c(as.numeric(lambdas[, "sd"]), NA),
    cohort = c("train", "test"),
    modName = modName)
  panels <- data.frame(panelLength = as.numeric(panelLength),
    dataset = names(panelLength)) %>% 
    gather(dataset, panelLength) %>% 
    mutate(modName = modName)
  return(list(error=error, panels=panels))
}

###############################################################################
## Ensemble_Enet
ensembleEnetFit = function(X_train, X_test, y_train, y_test, penalties, folds, iter, modName){
  ## Run caret (glmnet) to predicted probablities for each dataset
  library(caret); library(glmnet);
  ctrl <- trainControl(method = "repeatedcv",
    number = folds,
    repeats = iter,
    classProbs = TRUE,
    savePredictions = TRUE)
  set.seed(123)
  ctrl$index <- createMultiFolds(y_train, k = folds, times = iter)
  
  mods <- vector("list", length(X_train))
  names(mods) <- names(X_train)
  for(dat in names(mods)){
    enetGrid = expand.grid(alpha = 1, lambda = penalties)
    mods[[dat]] <- train(x=X_train[[dat]], y=y_train,
      preProc=c("center", "scale"),
      method = "glmnet",
      type.multinomial = "grouped",
      tuneGrid = enetGrid,
      trControl = ctrl)
  }
  
  pb <- lapply(mods, function(i){
    i$pred %>% 
      filter(lambda == i$bestTune$lambda) %>% 
      mutate(Resample = sapply(strsplit(Resample, "\\."), function(i) i[2])) %>% 
      gather(pam50, pb, Basal:LumB) %>% 
      dplyr::select(rowIndex, obs, Resample, pam50, pb)
  })
  
  ## Training error
  trainErr_Ensemble_enet <- do.call(rbind, pb) %>% 
    group_by(rowIndex, obs, Resample, pam50) %>% 
    summarise(pb = mean(pb)) %>% 
    filter(pb == max(pb)) %>% 
    group_by(Resample) %>% 
    nest() %>% 
    mutate(avgErr = purrr::map(data, ~{
      mat <- table(.$obs, .$pam50)
      mean((colSums(mat) - diag(mat))/colSums(mat))
    })) %>% 
    unnest(avgErr) %>% 
    dplyr::summarise(mean_err = mean(avgErr), sd_err = sd(avgErr), cohort = "train", modName = modName)
  
  ## Test error
  testEnet <- lapply(names(X_train), function(j){
    alpha = 1
    lambda = mods[[j]]$bestTune$lambda
    mod = glmnet(X_train[[j]], y_train, alpha = 1, lambda = lambda, family="multinomial", type.multinomial = "grouped")
    predictors <- names(which(coef(mod)[[1]][,1] != 0))[-1]
    pred <- predict(mod, X_test[[j]], type = "response")[,,1] %>%
      as.data.frame() %>% 
      mutate(sample = rownames(.)) %>% 
      gather(pam50, pb, -sample)
    list(predictors=predictors, pred=pred)
  }) %>% zip_nPure() 
  test_classes_ensemble_enet <- do.call(rbind, testEnet$pred) %>% 
    group_by(sample, pam50) %>% 
    summarise(pb = mean(pb)) %>% 
    group_by(sample) %>% 
    filter(pb == max(pb))
  mat <- table(test_classes_ensemble_enet$pam50, y_test[test_classes_ensemble_enet$sample])
  testErr_Ensemble_enet <- data.frame(mean_err = mean((colSums(mat) - diag(mat))/colSums(mat)),
    sd_err = NA, cohort = "test", modName = modName)
  
  error = rbind(trainErr_Ensemble_enet, testErr_Ensemble_enet)
  panels = data.frame(panelLength = sapply(testEnet$predictors, length), dataset = names(X_test), modName = modName)
  
  return(list(error=error, panels=panels))
}
