##================================================================#
#	File name: peek_data.R		  					#
#	Author:    Jianying Li							#
#	History:   Initially coded July, 15th 2013	 		#
#	         				      #
##================================================================#
winDir <- "C:/Users/li11/Documents/"
winProjDir <- "x:/"
linuxDir <- "~/"

library(gplots)
library(Rlab)
library(discrimARTs)

lib1 <- paste (winProjDir, "R-project/customPackages/dataManipTools.R", sep ="")
source(lib1)

lib2 <- paste (winProjDir,"R-project/customPackages/arraySeqTools.R", sep ="")
source(lib2)

lib3 <- paste (winProjDir,"R-project/customPackages/plotTools.R", sep ="")
source(lib3)



working.dir <-paste (winProjDir,"/myGit/mixturemodel/workingDir/", sep="")
setwd(working.dir)

file.dir <-paste (winProjDir,"/myGit/mixturemodel/data/", sep="")

#out.dir <-paste (winProjDir,"/project2013/comparativeExpression/fifthTrial/", sep="")

##================================
##	functions
##================================

 
##========================================================
## Stat analysis	
##========================================================
f <- "fmd_DT_raw.txt"
f <- "oscc-olk1_parsed.txt"
f_IN <- paste (file.dir, f , sep="")
dt <- read.table(f_IN, header= F, sep = "\t")
plotDensity(dt)

mix.prob = 0.005
mix.prob = 0.95

fit.default.oscc.olk <- mix.mle(dt$V1)
fit.default.oscc.olk <- mix.mle(dt$V1, method = "normal", mix.prob=mix.prob)

fit.default.oscc.olk
plot(fit.default.oscc.olk)


fit.oscc.olk <- mix.mle(dt$V1, method="normal",
	mix.prob=mix.prob, dist1.par1=1, dist1.par2=0.24, dist2.par1=2.3, dist2.par2=0.24)

fit.oscc.olk <- mix.mle(dt$V1, method="normal",
	mix.prob=mix.prob, dist1.par1=1, dist1.par2 = 0.2, dist2.par1=2, dist2.par2= 0.2)
## Default printing and plotting methods
print(fit.oscc.olk)

layout(1:2)
plot(fit.oscc.olk)
plot(fit.default.oscc.olk)




