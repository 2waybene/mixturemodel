first.mean
first.sd 


second.man <- first.mean
second.mean <- second.man
second.sd <- first.sd


third.man <- first.mean
third.mean <- third.man
third.sd <- first.sd

save (first.mean, second.mean, third.mean, first.sd , second.sd, third.sd, file = "sample_128110_stats.rda")


weight <- c(250, 80, 5)  # Number of fish
weight <- weight / sum(weight)

mean   <- c(first.mean, second.mean, third.mean)
sigma  <- c(first.sd ,  second.sd, third.sd)

Delta <- 0.1
x <- seq(0,5.5, by=Delta)
y1 <- weight[1]*P(x, mean[1], sigma[1])
y2 <- weight[2]*P(x, mean[2], sigma[2])
y3 <- weight[3]*P(x, mean[3], sigma[3])
y <- y1 + y2 + y3 

par(mfrow=c(1,1))

plot(x,y, type="l", lwd=3,
     main="Mixture of three population at: 250:80:5  - sample 128110",
     xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")



