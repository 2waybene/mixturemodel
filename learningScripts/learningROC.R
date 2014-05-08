library(pROC)

data(aSAH)

# Syntax (response, predictor):
plot.roc(aSAH$outcome, aSAH$s100b)

# With a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)

# identical:
plot(rocobj)
plot.roc(rocobj)
# Add a smoothed ROC:
plot.roc(smooth(rocobj), add=TRUE, col="blue")
legend("bottomright", legend=c("Empirical", "Smoothed"),
       col=c(par("fg"), "blue"), lwd=2)
# With more options:
plot(rocobj, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="blue", print.thres=TRUE)
# To plot a different partial AUC, we need to ignore the existing value
# with reuse.auc=FALSE:
plot(rocobj, print.auc=TRUE, auc.polygon=TRUE, partial.auc=c(1, 0.8),
     partial.auc.focus="se", grid=c(0.1, 0.2), grid.col=c("green", "red"),
     max.auc.polygon=TRUE, auc.polygon.col="blue", print.thres=TRUE,
     reuse.auc=FALSE)
# Add a line to the previous plot:
plot.roc(aSAH$outcome, aSAH$wfns, add=TRUE)
# Alternatively, you can get the plot directly





data(aSAH)
# Syntax (response, predictor):
auc(aSAH$outcome, aSAH$s100b)
# With a roc object:
rocobj <- roc(aSAH$outcome, aSAH$s100b)
# Full AUC:
auc(rocobj)
# Partial AUC:
auc(rocobj, partial.auc=c(1, .8), partial.auc.focus="se", partial.auc.correct=TRUE)
# Alternatively, you can get the AUC directly from roc():
roc(aSAH$outcome, aSAH$s100b)$auc
roc(aSAH$outcome, aSAH$s100b,
partial.auc=c(1, .8), partial.auc.focus="se",
partial.auc.correct=TRUE)$auc



multiclass.roc(aSAH$gos6, aSAH$s100b)
# Produces an innocuous warning because one level has no observation
# Select only 3 of the aSAH$gos6 levels:
multiclass.roc(aSAH$gos6, aSAH$s100b, levels=c(3, 4, 5))
# Give the result in percent
multiclass.roc(aSAH$gos6, aSAH$s100b, percent=TRUE)

aSAH$gos6




# Basic example
multiclass.roc(aSAH$gos6, aSAH$s100b)
# Produces an innocuous warning because one level has no observation
 
# Select only 3 of the aSAH$gos6 levels:
multiclass.roc(aSAH$gos6, aSAH$s100b, levels=c(3, 4, 5))
 
# Give the result in percent
multiclass.roc(aSAH$gos6, aSAH$s100b, percent=TRUE)


# randomForest & pROC packages should be installed:
# install.packages(c('randomForest', 'pROC'))
data(iris)
library(randomForest)
library(pROC)
set.seed(1000)
# 3-class in response variable
rf = randomForest(Species~., data = iris, ntree = 100)

predictions <- as.numeric(predict(rf, iris, type = 'response'))
multiclass.roc(iris$Species, predictions)