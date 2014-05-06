<<<<<<< HEAD
##	File:   preModHong_02.R
##	Author: Hong Xu



library(caret)
library(pROC)
library(Metrics)


##### REF: http://stats.stackexchange.com/questions/31579/what-is-the-optimal-k-for-the-k-nearest-neighbour-classifier-on-the-iris-dat
# https://gist.github.com/zachmayer/3061272
#Multi-Class Summary Function
#Based on caret:::twoClassSummary

#====================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

#root <- windows
root <- mac.os




require(compiler)
## Loading required package: compiler

multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL)
{

#Load Libraries

	require(Metrics)
	require(caret)

	#Check data
	if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
		stop("levels of observed and predicted data do not match")
	
	#Calculate custom one-vs-all stats for each class
	prob_stats <- lapply(levels(data[, "pred"]), function(class){


		#Grab one-vs-all data for the class
		pred <- ifelse(data[, "pred"] == class, 1, 0)
		obs <- ifelse(data[, "obs"] == class, 1, 0)
		prob <- data[,class]

		#Calculate one-vs-all AUC and logLoss and return
		cap_prob <- pmin(pmax(prob, .000001), .999999)
		prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
		names(prob_stats) <- c('ROC', 'logLoss')
		return(prob_stats)
	})



	prob_stats <- do.call(rbind, prob_stats)
	rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))

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
	stats <- stats[! names(stats) %in% c('AccuracyNull', 'Prevalence', 'Detection Prevalence')]

	#Clean names and return
	names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
	return(stats)
})


	## Note: no visible binding for global variable Metrics 
	## Note: no visible binding for global variable caret 
	## set up working directory

	setwd(paste (root, "/myGit/mixturemodel/reconData/para1/", sep=""))

	##	param1
	data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")
	
	# remove columns with zero variance
	# http://stackoverflow.com/questions/8805298/quickly-remove-zero-variance-variables-from-a-data-frame
	var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
	dataN0 <- data[,-which(var0)]
	# drop the first column of ID?
	dataN0[,1] <- NULL

	##### BEGIN: data partition >>>>>
	## set random seed
	set.seed(12345)
	## create data partition
	inTrainingSet <- createDataPartition(dataN0$label, p=.7, list=FALSE)
	labelTrain <- dataN0[ inTrainingSet,]
	labelTest <- dataN0[-inTrainingSet,]
	nrow(labelTrain)

	nrow(labelTest)

	##### END: data partition <<<<<
	##### BEGIN: tune the parameters >>>>>
	## control:
	# resampling technique: 5-repeat 10-fold cross-validation
	# performance metrics: ROC AUC curve
	ctrl <- trainControl(method = "repeatedcv",
				repeats = 5,
				summaryFunction = multiClassSummary,
				classProbs = TRUE)
	
	##### END: tune the parameters <<<<<
	##### BEGIN: train model - svm >>>>>
	set.seed(1024)
	svmFit <- train(label ~ ., data = labelTrain,
		## training model: svm >>>

		method = "svmRadial",
		metric = "ROC",
		# pre-process
		# preProc = c("center", "scale", "YeoJohnson"),
		tuneLength = 10,
		trControl = ctrl)

	## prediction
	svmPred <- predict(svmFit, labelTest)
	str(svmPred)
	## Factor w/ 3 levels "c","k","n": 1 1 1 1 1 1 1 1 1 1 ...
	## predicted probabilities
	svmProbs <- predict(svmFit, labelTest, type = "prob")
	str(svmProbs)
	## 
	## 
	##data.frame: 81 obs. of 3 variables:
	## $ c: num 0.909 0.9 0.94 0.907 0.918 ...
	## $ k: num 0.058 0.0644 0.0371 0.0599 0.0519 ...
	## $ n: num 0.0326 0.0355 0.0227 0.0334 0.0297 ...
	## evaluation
	confusionMatrix(svmPred, labelTest$label)

	## plot
	for(stat in c(
		'Accuracy',
		'Kappa', 
		'AccuracyLower', 
		'AccuracyUpper', 
		'AccuracyPValue',
		'Sensitivity', 
		'Specificity', 
		'Pos_Pred_Value',
		'Neg_Pred_Value', 
		'Detection_Rate', 
		'ROC', 
		'logLoss'))
	{
		print(plot(svmFit, metric=stat))
	}

##### END: train model - svm <<<<<
	
	##### BEGIN: train model - nb >>>>>
	#	https://stat.ethz.ch/pipermail/r-help/2007-October/144592.html
	library(klaR)
	nbFit <- train(label ~ .,
		method = "nb",
		data = labelTrain,
		# pre-process
		# preProc = c("center", "scale", "YeoJohnson"),
		rControl = trainControl(method = "repeatedcv", repeats = 5))
	nbPred <- predict(nbFit, labelTest)
	confusionMatrix(nbPred, labelTest$label)


	pcannFit <- train(label ~ .,
		method = "pcaNNet",
		data = labelTrain,
		# pre-process
		# preProc = c("center", "scale", "YeoJohnson"),
		trace = FALSE)
	
	pcannPred <- predict(pcannFit, labelTest)
	confusionMatrix(pcannPred, labelTest$label)


	##### END: train model - pcaNNet <<<<<
	##### BEGIN: train model - rf >>>>>
	rfFit <- train(label ~ .,
			method = "rf",
			# pre-process
			# preProc = c("center", "scale", "YeoJohnson"),
			data = labelTrain)
	rfPred <- predict(rfFit, labelTest)
	confusionMatrix(rfPred, labelTest$label)



	##### END: train model - rf <<<<<
	##### BEGIN: train model - RRF >>>>>
	rrfFit <- train(label ~ .,
			method = "RRF",
			# pre-process
			# preProc = c("center", "scale", "YeoJohnson"),
			data = labelTrain)
	rrfPred <- predict(rrfFit, labelTest)
	confusionMatrix(rrfPred, labelTest$label)








=======
##	File:   preModHong_02.R
##	Author: Hong Xu

library(caret)
library(pROC)
library(Metrics)
##### REF: http://stats.stackexchange.com/questions/31579/what-is-the-optimal-k-for-the-k-nearest-neighbour-classifier-on-the-iris-dat
# https://gist.github.com/zachmayer/3061272
#Multi-Class Summary Function
#Based on caret:::twoClassSummary


require(compiler)
## Loading required package: compiler

#====================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

#root <- windows
root <- mac.os


multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL)
{

#Load Libraries

	require(Metrics)
	require(caret)

	#Check data
	if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
		stop("levels of observed and predicted data do not match")
	
	#Calculate custom one-vs-all stats for each class
	prob_stats <- lapply(levels(data[, "pred"]), function(class){


		#Grab one-vs-all data for the class
		pred <- ifelse(data[, "pred"] == class, 1, 0)
		obs <- ifelse(data[, "obs"] == class, 1, 0)
		prob <- data[,class]

		#Calculate one-vs-all AUC and logLoss and return
		cap_prob <- pmin(pmax(prob, .000001), .999999)
		prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
		names(prob_stats) <- c('ROC', 'logLoss')
		return(prob_stats)
	})



	prob_stats <- do.call(rbind, prob_stats)
	rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))

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
	stats <- stats[! names(stats) %in% c('AccuracyNull', 'Prevalence', 'Detection Prevalence')]

	#Clean names and return
	names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
	return(stats)
})


	## Note: no visible binding for global variable Metrics 
	## Note: no visible binding for global variable caret 
	## set up working directory

	setwd(paste (root, "/myGit/mixturemodel/reconData/para1/", sep=""))

	##	param1
	data <- read.table("recon_3classes_para1.txt", header=TRUE, sep = "\t")
	
	# remove columns with zero variance
	# http://stackoverflow.com/questions/8805298/quickly-remove-zero-variance-variables-from-a-data-frame
	#var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
  var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.numeric(x) else x)))
	dataN0 <- data[,-which(var0)]
	# drop the first column of ID?
	dataN0[,1] <- NULL

	##### BEGIN: data partition >>>>>
	## set random seed
	set.seed(12345)
	## create data partition
	inTrainingSet <- createDataPartition(dataN0$label, p=.7, list=FALSE)
	labelTrain <- dataN0[ inTrainingSet,]
	labelTest <- dataN0[-inTrainingSet,]
	nrow(labelTrain)

	nrow(labelTest)

	##### END: data partition <<<<<
	##### BEGIN: tune the parameters >>>>>
	## control:
	# resampling technique: 5-repeat 10-fold cross-validation
	# performance metrics: ROC AUC curve
	ctrl <- trainControl(method = "repeatedcv",
				repeats = 5,
				summaryFunction = multiClassSummary,
				classProbs = TRUE)
	
	##### END: tune the parameters <<<<<
	##### BEGIN: train model - svm >>>>>
	set.seed(1024)
	svmFit <- train(label ~ ., data = labelTrain,
		## training model: svm >>>

		method = "svmRadial",
		metric = "ROC",
		# pre-process
		# preProc = c("center", "scale", "YeoJohnson"),
		tuneLength = 10,
		trControl = ctrl)

	## prediction
	svmPred <- predict(svmFit, labelTest)
	str(svmPred)
	## Factor w/ 3 levels "c","k","n": 1 1 1 1 1 1 1 1 1 1 ...
	## predicted probabilities
	svmProbs <- predict(svmFit, labelTest, type = "prob")
	str(svmProbs)
	## 
	## 
	##data.frame: 81 obs. of 3 variables:
	## $ c: num 0.909 0.9 0.94 0.907 0.918 ...
	## $ k: num 0.058 0.0644 0.0371 0.0599 0.0519 ...
	## $ n: num 0.0326 0.0355 0.0227 0.0334 0.0297 ...
	## evaluation
	confusionMatrix(svmPred, labelTest$label)

	## plot
	for(stat in c(
		'Accuracy',
		'Kappa', 
		'AccuracyLower', 
		'AccuracyUpper', 
		'AccuracyPValue',
		'Sensitivity', 
		'Specificity', 
		'Pos_Pred_Value',
		'Neg_Pred_Value', 
		'Detection_Rate', 
		'ROC', 
		'logLoss'))
	{
		print(plot(svmFit, metric=stat))
	}

##### END: train model - svm <<<<<
	
	##### BEGIN: train model - nb >>>>>
	#	https://stat.ethz.ch/pipermail/r-help/2007-October/144592.html
	library(klaR)
	nbFit <- train(label ~ .,
		method = "nb",
		data = labelTrain,
		# pre-process
		# preProc = c("center", "scale", "YeoJohnson"),
		rControl = trainControl(method = "repeatedcv", repeats = 5))
	nbPred <- predict(nbFit, labelTest)
	confusionMatrix(nbPred, labelTest$label)


	pcannFit <- train(label ~ .,
		method = "pcaNNet",
		data = labelTrain,
		# pre-process
		# preProc = c("center", "scale", "YeoJohnson"),
		trace = FALSE)
	
	pcannPred <- predict(pcannFit, labelTest)
	confusionMatrix(pcannPred, labelTest$label)


	##### END: train model - pcaNNet <<<<<
	##### BEGIN: train model - rf >>>>>
	rfFit <- train(label ~ .,
			method = "rf",
			# pre-process
			# preProc = c("center", "scale", "YeoJohnson"),
			data = labelTrain)
	rfPred <- predict(rfFit, labelTest)
	confusionMatrix(rfPred, labelTest$label)



	##### END: train model - rf <<<<<
	##### BEGIN: train model - RRF >>>>>
	rrfFit <- train(label ~ .,
			method = "RRF",
			# pre-process
			# preProc = c("center", "scale", "YeoJohnson"),
			data = labelTrain)
	rrfPred <- predict(rrfFit, labelTest)
	confusionMatrix(rrfPred, labelTest$label)








>>>>>>> 91a98f03dae9a38b1a7c69f2b14bd565a6f8a7ea
