library(earth)
data(etitanic)
head(model.matrix(survived ~ ., data = etitanic))


dummies <- dummyVars(survived ~ ., data = etitanic)
head(predict(dummies, newdata = etitanic))

data(mdrr)
data.frame(table(mdrrDescr$nR11))


nzv <- nearZeroVar(mdrrDescr, saveMetrics= TRUE)
nzv[nzv$nzv,][1:10,]


dim(mdrrDescr)



nzv <- nearZeroVar(mdrrDescr)
filteredDescr <- mdrrDescr[, -nzv]
dim(filteredDescr)

descrCor <-  cor(filteredDescr)
highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .999)

descrCor <- cor(filteredDescr)
summary(descrCor[upper.tri(descrCor)])

highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
filteredDescr <- filteredDescr[,-highlyCorDescr]
descrCor2 <- cor(filteredDescr)
summary(descrCor2[upper.tri(descrCor2)])

