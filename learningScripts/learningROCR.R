## computing a simple ROC curve (x-axis: fpr, y-axis: tpr)
library(ROCR)
data(ROCR.simple)
pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf)
## precision/recall curve (x-axis: recall, y-axis: precision)
perf1 <- performance(pred, "prec", "rec")
plot(perf1)
## sensitivity/specificity curve (x-axis: specificity,
## y-axis: sensitivity)
perf1 <-performance(pred, "sens", "spec")
plot(perf1)



data(ROCR.simple)
pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels )
perf <- performance( pred, "tpr", "fpr" )
plot( perf )
# To entertain your children, make your plots nicer
# using ROCR’s flexible parameter passing mechanisms
# (much cheaper than a finger painting set)
par(bg="lightblue", mai=c(1.2,1.5,1,1))
plot(perf, main="ROCR fingerpainting toolkit", colorize=TRUE,
     xlab="Mary axis", ylab="", box.lty=7, box.lwd=5,
     box.col="gold", lwd=17, colorkey.relwidth=0.5, xaxis.cex.axis=2,
     xaxis.col='blue', xaxis.col.axis="blue", yaxis.col='green', yaxis.cex.axis=2,
     yaxis.at=c(0,0.5,0.8,0.85,0.9,1), yaxis.las=1,xaxis.lwd=2, yaxis.lwd=3,
     yaxis.col.axis="orange", cex.lab=2, cex.main=2)


##============================

data(ROCR.hiv)
attach(ROCR.hiv)
str(ROCR.hiv)


pred.svm <- prediction(hiv.svm$predictions, hiv.svm$labels)
perf.svm <- performance(pred.svm, 'tpr', 'fpr')
pred.nn <- prediction(hiv.nn$predictions, hiv.svm$labels)
perf.nn <- performance(pred.nn, 'tpr', 'fpr')


plot(perf.svm, lty=3, col="red", main="SVMs and NNs for prediction of HIV-1 coreceptor usage")
plot(perf.nn, lty=3, col="blue",add=TRUE)
plot(perf.svm, avg="vertical", lwd=3, col="red", spread.estimate="stderror",plotCI.lwd=2,add=TRUE)
plot(perf.nn, avg="vertical", lwd=3, col="blue", spread.estimate="stderror",plotCI.lwd=2,add=TRUE)


legend(0.6,0.6,c('SVM','NN'),col=c('red','blue'),lwd=3)








