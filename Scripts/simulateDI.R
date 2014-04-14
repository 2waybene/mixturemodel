#	simulating DNA D.I. values

#	normal population

mean.norm <- c()
for (i in 1:20)
{
	mean.norm[i] <- runif (1, 0.9,1.2)
}
mean.normalP <- mean.norm[sample(1:20,1)]

#	mitotic population
mean.mitotic <- c()
for (i in 1:20)
{
	mean.mitotic[i] <- runif (1, 1.7,2.2)
}
mean.mitoticP <- mean.mitotic[sample(1:20,1)]


#	aneuploidy  population

mean.aneu <- 3.3
st.aneu <-3.5
n = 40

sample.aneu <- rnorm(n, mean = mean.aneu, sd = st.aneu)
sample.aneuP <- sample.aneu[sample.aneu >=0]
plot(density(sample.aneuP))
mean(sample.aneuP)
sd(sample.aneuP)


##=======================================
##	sample D.I. range
##=======================================
Delta <- 0.001
x <- seq(0,11, by=(11/512))
x <- x[-1]
length(x)

y1 <- P(x, mean.normalP, 0.3)
y2 <- P(x, mean.mitoticP, 0.3)



plot(x, y1)
plot(x, y2)

y3 <- P(x, mean(sample.aneuP), sd(sample.aneuP))
plot(x, y3)



##=============================
##	Different weight
##==============================

##==============================
normWeight <- (0.55,0.41,0.04)
mitoticWeight <- (0.41, 0.55,0.04)
OSCCWeight <- (0.3, 0.2,0.5)


weight <- normWeight


y1 <- weight[1]*P(x, mean.normalP, sigmaP)
y2 <- weight[2]*P(x, mean.mitoticP, sigmaM)
y3 <- weight[3]*P(x, mean.aneuP, sigmaA)




##=========================================
par(mfrow=c(1,1))

plot(x,y, type="l", lwd=3,
     main="Heming Lake Pike: Distribution by Age Groups",
     xlab="Length [cm]", ylab="Probability Density")
abline(h=0, lty="dotted")
lines(x,y1, col="red")
lines(x,y2, col="green")
lines(x,y3, col="blue")



x <- 1:12
# a random permutation
sample(x)
