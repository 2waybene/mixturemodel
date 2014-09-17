
library (Metrics) 
library (caret)

#CLEAR WORKSPACE
rm(list = ls(all = TRUE))

require(compiler)
multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL)
{

    #Load Libraries
 #   require(Metrics)
 #   require(caret)
    #Check data
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
      stop("levels of observed and predicted data do not match")
      #Calculate custom one-vs-all stats for each class
     prob_stats <- lapply(levels(data[, "pred"]), function(class)
       {
          #Grab one-vs-all data for the class
           pred <- ifelse(data[, "pred"] == class, 1, 0)
           obs <- ifelse(data[, "obs"] == class, 1, 0)
           prob <- data[,class]

         #Calculate one-vs-all AUC and logLoss and return

           cap_prob <- pmin(pmax(prob, .000001), .999999)
           prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
           names(prob_stats) <- c("ROC", "logLoss")
           return(prob_stats)
      })

      prob_stats <- do.call(rbind, prob_stats)
      rownames(prob_stats) <- paste( "Class:" , levels(data[, "pred"]))
    
      #Calculate confusion matrix-based statistics
      CM <- confusionMatrix(data[, "pred"], data[, "obs"])
                            #Aggregate and average class-wise stats
                            #Todo: add weights
      class_stats <- cbind(CM$byClass, prob_stats)
      class_stats <- colMeans(class_stats)
                            #Aggregate overall stats
      overall_stats <- c(CM$overall)
    
      #Combine overall with class-wise stats and remove some stats we don't want
      stats <- c(overall_stats, class_stats)
      stats <- stats[! names(stats) %in% c("AccuracyNull","Prevalence", "Detection Prevalence")]
                            
    #Clean names and return
                  
    names(stats) <- gsub('[[:blank:]] +', '_' , names(stats))
    return(stats)
})



gc(reset=TRUE)

#Setup parallel cluster
#If running on the command line of linux, use method='fork'
library(doParallel)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
 
#Fit model
library(caret)
set.seed(19556)
model <- train(
Species~.,
data=iris,
method='knn',
tuneGrid=expand.grid(.k=1:30),
metric='Accuracy',
trControl=trainControl(
method='repeatedcv',
number=10,
repeats=15,
classProbs=TRUE,
summaryFunction=multiClassSummary))
  
#Stop parallel cluster
stopCluster(cl)
 
#Save pdf of plots
dev.off()
pdf('plots.pdf')
for(stat in c('Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyPValue',
'Sensitivity', 'Specificity', 'Pos_Pred_Value',
'Neg_Pred_Value', 'Detection_Rate', 'ROC', 'logLoss')) {
print(plot(model, metric=stat))
}
dev.off()
