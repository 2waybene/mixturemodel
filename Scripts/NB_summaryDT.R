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
##=====================================================================

library(caret)
TrainData <- dt4model[,-1]
TrainClasses <- dt4model[,1]
trt=trainControl(method='cv',number=10)
nbFit <- train(TrainData, TrainClasses, method = "nb", trControl = trt)



##=========================================================================
library(e1071)
trainNB =  dt4model

modelNB2 <- naiveBayes(Group ~ ., data = trainNB)
predNB2 <- predict(modelNB2, trainNB)
predNB2



#ainControl(method = "boot"))
library(klaR)
nbFit2 <- NaiveBayes(TrainData, TrainClasses)
, usekernel=TRUE)



##=============================================================
##	Split datat into training and testing
##=============================================================

# prepare data for caret training: separated labels and features
Ygtm <- dt4model[1] # convert sub data frame to vector
Xcur <- dt4model[-1] # get features without genotype labels




# split train and test sets (half/half)
library(caret)
set.seed(1)
testID <-c ()
for (i in 1:3)
{
	testID[c((i*2-1), i*2)] <-  (i-1)*9 + (sample(c(1:9),2))
}
testID 	


trainX <- Xcur[-testID,]
testX <- Xcur[testID,]

trainY <- Ygtm[-testID,]
testY <- Ygtm[testID,]



# train LDA model
ldaFit <- train( x = trainX, y = trainY,  method = "lda")
predLDA <- predict (ldaFit, testX)
predLDA
