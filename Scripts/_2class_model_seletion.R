##=========================================================
##  File: _2class_model_selection_02.R
##  Author: Jianying Li
##  Comment: compare different data processing parameter 
##	       settings and compare seven models (ROC - metrics)
##		 produce model comparison box plots
##
##  Extended from _2class_prediction_model_comparison_02.R
##=========================================================


#===============================================
#	Set up os paths
#==============================================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

root <- windows
#root <- linux 
#root <- mac.os
##===============================================


library(caret)
##===================================
## set up working directory
## and getting the dataset
##===================================

setwd(paste (root, "/myGit/mixturemodel/reconData/para2/", sep=""))

reconFiles <- list.files (pattern = "recon_*")

for (i in length(reconFiles))
{

	parSet <- sub (".txt", "", sub ("recon_3classes_", "", reconFiles[i]))
	comFigName <- paste ("Model_comparison_", parSet, ".png", sep = "")
	title <- paste ("Model comparison with ", parSet,sep = "")
	comFigName <- paste (root, "/myGit/mixturemodel/modeling/model_selection/", comFigName, sep="")
	data <- read.table(reconFiles[i], header=TRUE, sep = "\t")


	##	data cleaning

	var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
	dataN0 <- data[,-which(var0)]
	#dataN0[,1] <- paste ("ID_", dataN0[,1], sep = "")
	rownames(dataN0) <- paste ("ID_", dataN0[,1], sep = "")

 dataN0[,1] <- null


 

  	##	Retain data ONLY with two classes

	data.2.classes <- dataN0[-which (dataN0$label == "k"),]
	data.k <- dataN0[which (dataN0$label == "k"),]



	labels <- as.vector(data.2.classes$label)
	data.2.classes <- data.2.classes[,-17]
	data.2.classes <- cbind (data.2.classes, label=labels)
	
	file2classes  <- data.2.classes
	file.olk <- data.k 

	## create data partition
	inTrainingSet <- createDataPartition(file2classes$label, p=.7, list=FALSE)
	labelTrain <- file2classes[ inTrainingSet,]
	labelTest <- file2classes[-inTrainingSet,]
	nrow(labelTrain)
	nrow(labelTest)


	##=======================
	##	rpart
	##=======================
	require(rpart)
	require(partykit)
	
	rpartFull <- rpart(label ~ ., data = labelTrain)



	require(partykit)
	rpartFulla <- as.party(rpartFull)


	rpartPred <- predict(rpartFull, labelTest, type = "class")
	confusionMatrix(rpartPred, labelTest$label)   # requires 2 factor vectors


	##################################################################

	cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)

	cvCtrl.2 <- trainControl(method = "LOOCV", 
		      summaryFunction = twoClassSummary,
            	           classProbs = TRUE)

	set.seed(12345)
	rpartTune <- train(label ~ ., data = labelTrain, 
                   method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
                   trControl = cvCtrl)

	rpartPred2 <- predict(rpartTune, labelTest)
	confusionMatrix(rpartPred2, labelTest$label)



	##====================================
	##	SVM Example 
	##====================================
	set.seed(12345)
	svmTune <- train(label ~ ., 
		    data = labelTrain,

                 method = "svmRadial",
                 # The default grid of cost parameters go from 2^-2,
                 # 0.5 to 1, 
                 # We'll fit 9 values in that sequence via the tuneLength
                 # argument.
                 tuneLength = 9,
                 preProc = c("center", "scale"),
                 metric = "ROC",   
                 trControl = cvCtrl)

	svmPred <- predict(svmTune, labelTest)
	confusionMatrix(svmPred, labelTest$label)


	##======================================
	##	partial logistic regression
	##========================================
	plrTune <- train(label ~ ., data = labelTrain,  
                 method = "multinom",
                 preProc = c("center", "scale"),
                 tuneGrid = data.frame(decay = c(0.1, 1, 10, 20, 40)),
                 trace = FALSE, maxit = 1000,
                 metric = "ROC",
                 trControl = cvCtrl)

	plrPred <- predict(plrTune, labelTest)
	confusionMatrix(plrPred, labelTest$label)

	##========================
	##	randomForest
	##========================

	require(randomForest)

	RFTune <- train(label ~.,data=labelTrain, 
				method = "rf",
		     tuneLength = 12,
                 metric = "ROC",
                 trControl = cvCtrl)

	RFPred2 <- predict(RFTune, labelTest)
	confusionMatrix(RFPred2, labelTest$label)


RFtrain <- randomForest(label ~.,data=labelTrain )
RFPred <- predict(RFtrain, labelTest)
confusionMatrix(RFPred, labelTest$label)


rrfFit <- train(label ~ ., method = "RRF", data = labelTrain)
rrfPred <- predict(rrfFit, labelTest)
confusionMatrix(rrfPred, labelTest$label)

	##========================
	##	knn
	##========================

	knnFit1 <- train(label ~.,data=labelTrain, 
                 method = "knn",
                 preProcess = c("center", "scale"),
 		     metric = "ROC",
                 tuneLength = 10,
                 trControl = cvCtrl)


	knnPred <- predict(knnFit1, labelTest)
	confusionMatrix(knnPred , labelTest$label)


	require(MASS)
	nnetFit <- train(label ~.,data=labelTrain, 
                 method = "nnet",
                 preProcess = "range", 
                 tuneLength = 2,
 		     metric = "ROC",
                 trace = FALSE,
                 maxit = 100,
			trControl = cvCtrl)



	nnetPred <- predict(nnetFit, labelTest)
	confusionMatrix(nnetPred , labelTest$label)


	##===================================
	## Collecting Results With resamples
	##===================================
	cvValues <- resamples(list(CART = rpartTune, SVM = svmTune, 
                           plr = plrTune, nnet = nnetFit,
					knn = knnFit1, rrf = RFTune))
            

	##=================
	#	box plot
	##=================
	png (comFigName)
	trellis.par.set()
	bwplot(cvValues, layout = c(3, 1), main =  title)
	dev.off()
}



##====================================
##	SVM Example 
##====================================
set.seed(12345)
svmTune <- train(label ~ ., 
		    data = labelTrain,

                 method = "svmRadial",
                 # The default grid of cost parameters go from 2^-2,
                 # 0.5 to 1, 
                 # We'll fit 9 values in that sequence via the tuneLength
                 # argument.
                 tuneLength = 9,
                 preProc = c("center", "scale"),
                 metric = "ROC",   
                 trControl = cvCtrl.2)

svmPred <- predict(svmTune, labelTest)

confusionMatrix(svmPred, labelTest$label)

svmPredProb <- predict(svmTune, labelTest , type = "prob")

boxplot(svmPredProb)
#points(svmPredProb)

str(svmPredProb)
label.c.as.c <- svmPredProb$c[labelTest$label=="c"]
label.n.as.c <- svmPredProb$c[labelTest$label=="n"]

svmPredProb$c =="c"


boxplot(svmPredProb, outpch = NA)
stripchart(svmPredProb,
            vertical = TRUE, method = "jitter",
            pch = 21, col = "maroon", bg = "bisque",
            add = TRUE)



##=============================
##	Test on olk sample
##=============================
dim(file.olk)
svmPred.k.prob  <- predict(svmTune, file.olk, type = "prob")



den.c <- density(svmPred.k.prob$c)
label.k.as.c <- svmPred.k.prob$c
den.n <- density(svmPred.k.prob$n)


length(which(svmPred.k.prob$c > 0.5))
length(which(svmPred.k.prob$n > 0.5))


plot( den.c)
lines( den.n, col = "red")

file.olk[-which(svmPred.k.prob$c >0.5),]

plot(density(svmPred.k.prob$c))
plot(density(svmPred.k.prob$n))
pairs(svmPred.k.prob)

##==============================
##	To get the figure 6
##==============================

rep ("c", length(label.c.as.c))
rep ("n", length(label.n.as.c))
rep ("k", length(label.k.as.c))

predicted.c <- list (label = as.vector(c(rep ("n", length(label.n.as.c)), rep ("k", length(label.k.as.c)), rep ("c", length(label.c.as.c)))), 
prob = as.vector(c( svmPredProb$c[labelTest$label=="n"], svmPred.k.prob$c, svmPredProb$c[labelTest$label=="c"])))
str(predicted.c)

boxplot(prob ~ label, data = as.data.frame(predicted.c), main = "Samples (by label) predicted as OSCC", ylab = "Probability", outpch = NA)
stripchart(prob ~ label, data = as.data.frame(predicted.c),
            vertical = TRUE, method = "jitter",
            pch = 21, col = "maroon", bg = "bisque",
            add = TRUE)
mtext ("Prediction probability")
