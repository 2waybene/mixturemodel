##==========================================================================
##	FileName: NB_summaryDT.R
##	Author: Jianying Li
##==========================================================================

source("x:/myGit/mixturemodel/Scripts/mixtureModelFunctions.R")
library(Rlab)

##=====Sample data
macFileDir <- "/Users/li11/myGit/mixturemodel/data/"
winFileDir <- "x:/myGit/mixturemodel/data/"
f <- "summary_dt.txt";
f_IN <- paste (winFileDir, f , sep="")
dt <- read.table(f_IN, header= T, sep = "\t")
str(dt)

which(is.na(dt$Mean1))	#No missing
which(is.na(dt$Mean2))	#No missing
which(is.na(dt$Mean3))

which(is.na(dt$SD1))	#No missing
which(is.na(dt$SD2))	#No missing
which(is.na(dt$SD3))

##==============================================
##	Use median to impute the missing values
##==============================================
dt4model <- dt[c(1:7)]
dt4model$Mean3[which(is.na(dt4model$Mean3))] <- stats(dt4model$Mean3)[6] #median
dt4model$SD3[which(is.na(dt4model$SD3))] <- stats(dt4model$SD3)[2] #mean
boxplot(dt4model[-1])
title("Box plot on pilot data")
mtext("Missing was imputed with median or mean")
str(dt4model)
setwd("X:/myGit/mixturemodel/firstDraft")
save(dt4model, file ="pilot_data_noMissing.Rda")
##=====================================================================

##=============================================================
##	Split datat into training and testing
##=============================================================
dt4model [,1] <- as.factor(c( rep ("norm", 9), rep("olk", 9), rep ("oscc", 9)))

dt4model2 <- dt4model[c(1:6)]


set.seed(1234)
testID <-c ()
for (i in 1:3)
{
	testID[c((i*2-1), i*2)] <-  (i-1)*9 + (sample(c(1:9),2))
}
testID 	


trainDT <- dt4model2[-testID,]
testDT  <- dt4model2[testID,]

str(trainDT)
str(testDT)

ans <- testDT[,1]
##=================================================
#	Try some prediction models
##=================================================
#	Naive Bayes
##=================================================
library(e1071)
trainNB =  trainDT 
testNB  =  testDT
modelNB2 <- naiveBayes(Group ~ ., data = trainNB)
predNB2 <- predict(modelNB2, testNB[,-1])
predNB2
table(predNB2, ans)

##=================================================
##	svm in e1071 package
##=================================================
trainSVM =  trainDT 
testSVM  =  testDT
modelSVM <- svm (Group ~ ., data = trainSVM)
predSVM  <- predict (modelSVM, testSVM[,-1])
predSVM

ans <- testSVM[,1]
table(predSVM, ans)


##=======================================================================
library(caret)


# train LDA model
ldaFit <- train( x = trainDT[,-1], y = trainDT[,1],  method = "lda")
predLDA <- predict (ldaFit, testDT[,-1])
predLDA


##======================================================================

TrainData <- dt4model[,-1]
TrainClasses <- dt4model[,1]
trt=trainControl(method='cv',number=10)
nbFit <- train(TrainData, TrainClasses, method = "nb", trControl = trt)







#ainControl(method = "boot"))
library(klaR)
nbFit2 <- NaiveBayes(TrainData, TrainClasses)
, usekernel=TRUE)





# train LDA model
ldaFit <- train( x = trainX, y = trainY,  method = "lda")
predLDA <- predict (ldaFit, testX)
predLDA
