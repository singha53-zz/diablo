##################################################################
#
# asthma_functions.R
#
##################################################################

leave.one.out.cv = function(X, Y, ncomp, keepX, design, sample){
  
  averageVotes_unpaired <- averageVotes_paired <- list()
  
  for(rm in 1 : length(unique(sample))){
    sequence <- which(rm == sample)
    trainX <- lapply(X, function(i) i[-sequence, ])
    trainY <- Y[-sequence]
    xtest <- lapply(X, function(i) i[sequence, ])
    xwtest <- lapply(X, function(i) i[sequence, ]-matrix(colMeans(i[sequence,]),ncol=ncol(i),nrow=length(sequence),byrow=TRUE))
    samplew <- sample[-sequence]
    
    Cov <- data.frame(sample=samplew, stimu=trainY)
    trainXw <-  lapply(trainX, function(i){
      withinVariation(i, Cov)
    })
    
    ## unpaired diablo
    diablo_unpaired = block.splsda(X = trainX, Y = trainY, ncomp = ncomp, keepX = keepX, design = design)
    test_unpaired = predict(diablo_unpaired, xtest, method = "all")
    averageVotes_unpaired[[rm]] <- test_unpaired$AveragedPredict[,,2][,1]
    
    ## paired diablo
    diablo_paired = block.splsda(X = trainXw, Y = trainY, ncomp = ncomp, keepX = keepX, design = design)
    test_paired = predict(diablo_paired, xwtest, method = "all")
    averageVotes_paired[[rm]] <- test_paired$AveragedPredict[,,2][,1]
  }
  avgVote <-  list(unpaired = unlist(averageVotes_unpaired), paired=unlist(averageVotes_paired))
  avgVote <- lapply(avgVote, function(i){
    roc.score = roc(response = sapply(strsplit(names(unlist(averageVotes_unpaired)), "\\."), function(i) i[3]), predictor = i, plot = FALSE, percent = TRUE, na.rm =TRUE, ci = TRUE)
    roc.res = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
    roc.res$fp = 100 - as.numeric(roc.res$Specificity)
    roc.res$tp = as.numeric(roc.res$Sensitivity)
    auc= paste(c(paste0("AUC = ", round(roc.score$auc, 1), "%"), paste("95% CI:", paste(paste0(round(as.numeric(roc.score$ci),1)[-2], "%"), collapse = "-"))), collapse = "\n")
    return(list(roc.res=roc.res, auc=auc))
  })
  avgVote
}