##=========================================================
##  File: build_olk_index.R
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

i = 2
reconFiles[2]

	parSet <- sub (".txt", "", sub ("recon_3classes_", "", reconFiles[i]))
	comFigName <- paste ("Model_comparison_", parSet, ".png", sep = "")
	title <- paste ("Model comparison with ", parSet,sep = "")
	comFigName <- paste (root, "/myGit/mixturemodel/modeling/model_selection/", comFigName, sep="")
	data <- read.table(reconFiles[i], header=TRUE, sep = "\t")


	##	data cleaning

	var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
	dataN0 <- data[,-which(var0)]
	rownames(dataN0) <- paste ("ID_", dataN0[,1], sep = "")

 dataN0[,1] <- NULL

  	##	Retain data ONLY with two classes

	data.2.classes <- dataN0[-which (dataN0$label == "k"),]
	data.k <- dataN0[which (dataN0$label == "k"),]

dim(data.2.classes)

	labels <- as.vector(data.2.classes$label)
	data.2.classes <- data.2.classes[,-16]
	data.2.classes <- cbind (data.2.classes, label=labels)

	
	file2classes  <- data.2.classes
	file.olk <- data.k 

	## create data partition


	## create data partition
set.seed(12345)
	inTrainingSet <- createDataPartition(file2classes$label, p=.7, list=FALSE)
	labelTrain <- file2classes[ inTrainingSet,]
	labelTest <- file2classes[-inTrainingSet,]
	nrow(labelTrain)
	nrow(labelTest)


	

	cvCtrl <- trainControl(method = "repeatedcv", repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)

	cvCtrl.2 <- trainControl(method = "LOOCV", 
		      summaryFunction = twoClassSummary,
            	           classProbs = TRUE)

	



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
rownames(file.olk)
rownames(svmPred.k.prob) <- rownames(file.olk)
#write.table (svmPred.k.prob, file = "prediction_on_olk.txt", row.name = TRUE,sep="\t")

den.c <- density(svmPred.k.prob$c)
label.k.as.c <- svmPred.k.prob$c
den.n <- density(svmPred.k.prob$n)


length(which(svmPred.k.prob$c > 0.5))
length(which(svmPred.k.prob$n > 0.5))


plot( den.c)
lines( den.n, col = "red")

file.olk[-which(svmPred.k.prob$c >0.5),]
file.olk[which(svmPred.k.prob$c >0.5),]

##==========================
##	Get three categories
##==========================

lowRisk <- file.olk[which(svmPred.k.prob$c <=0.3),]
highRisk <- file.olk[which(svmPred.k.prob$c >=0.7),]
midRisk <- file.olk[-which((svmPred.k.prob$c <=0.3) | (svmPred.k.prob$c >=0.7)),]

rownames(lowRisk)
rownames(highRisk)
rownames(midRisk)


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
colnames(as.data.frame(predicted.c))
predicted.c$label[predicted.c$label == "n"] <- "Normal"
predicted.c$label[predicted.c$label == "k"] <- "OLK"
predicted.c$label[predicted.c$label == "c"] <- "OSCC"
#boxplot(prob ~ label, data = as.data.frame(predicted.c), main = "Samples (by label) predicted as OSCC", ylab = "Probability", outpch = NA)
boxplot(prob ~ label, data = as.data.frame(predicted.c),  ylab = "Oral Cancer Risk Index (OCRI)", outpch = NA)
stripchart(prob ~ label, data = as.data.frame(predicted.c),
            vertical = TRUE, method = "jitter",
            pch = 21, col = "maroon", bg = "bisque",
            add = TRUE)
#mtext ("Prediction probability")



