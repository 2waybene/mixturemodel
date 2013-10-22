library(mixdist)
data(cassie)
plot(cassie)
title("Original mixture data")
cassie.mix <- as.mixdata(cassie)
str(cassie)
str(cassie.mix)

plot(cassie.mix)
title("Original mixture data")
mtext("In histogram format")
fitcass1 <- mix(cassie,mixparam(c(3,5,7,9,11),.5),"gamma",mixconstr(consigma="CCV"))
summary(fitcass1)
plot(fitcass1)
title("Gamma fit to mixture data")


fitcass2 <- mix(cassie,coef(fitcass1),"gamma")
fitcass2 <- mix(cassie,coef(fitcass1),"gamma",iterlim=150)
summary(fitcass2)

par(mfrow=c(1,1))
plot(fitcass2)


##==============================
#	Pike lake fishering data
##==============================

data(pikeraw) # load raw data `pikeraw'
pikeraw # display the data set `pikeraw'
mixgroup(pikeraw) # group raw data
pikemd <- mixgroup(pikeraw, breaks = c(0, seq(19.75, 65.75, 2), 80))
plot(pikemd)
mixgroup(pikeraw, breaks = c(0, seq(19.75, 65.75, 2), 80), usecondit = TRUE, k = 5)
# construct grouped data associated with conditional data
mixgroup(pikeraw, usecondit = TRUE)
mixgroup(pikeraw, usecondit = TRUE, k = 3) # grouping data with a warning message
mixgroup(pikeraw, usecondit = TRUE, k = 8)


#fit data

fitPike1 <- mix (pikemd, mixparam(c(23.32727, 33.09053, 41.27244, 51.24468, 61.31818), c(2.436638, 2.997594, 4.274585, 5.077521, 7.070303)), "norm",mixconstr(consigma="CCV"))
title ("Fit with MIX software by Peter McDonald")






