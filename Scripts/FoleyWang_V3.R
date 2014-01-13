##======================================================
#	File:   FoleyWang_V3.R
#	Author: Jianying Li
#	Initial coded for "liver" sample in May 2012
#	Get the heatmap, overlap gene list, incoporating 
#	house keeping gene list, tested from linux
#	Date:   01/06/2014
##======================================================
source("~/R-project/customPackages/dataManipTools.R")
source("~/R-project/customPackages/arraySeqTools.R")
source("~/R-project/customPackages/plotTools.R")
source("~/project2012/ASE/R-scripts/functions-ASE.R")
library(gplots)
library(Rlab)
library(affy)

##====================================
##	Functions...
##====================================
stripFileNames <- function (superList)
{
	#grep ("[0-9]+_Liver_", superList, value=TRUE)
	name1 <- sub("[0-9]+_Liver_", "", superList)
	#grep ("_\\(Mouse430_2)", name1, value=TRUE)
	name2 <- sub ("_\\(Mouse430_2)", "", name1)
	return(name2)
}

##====================================

setwd("X:/project2014/Foley_Yu/workingDir/")
pheno.dir <- "X:/project2014/Foley_Yu/doc/"
celpath = "X:/project2014/Foley_Yu/rawData/"



##==================================================================#
#	Read in experimental design information, read raw .cel file	  #
#	with experimental info incooporated					  #
##==================================================================#

expInfo = paste (pheno.dir, "expInfo.txt", sep="")
pheno = read.table(expInfo,sep="\t",header=TRUE, row.names=1)
str(pheno)
##	Unless I need to -- JYL 
# dat <- ReadAffy(filenames=paste(celpath, pheno$FileName,".CEL",sep=""), sampleNames=rownames(pheno), phenoData=pheno)

#	str(dat)
##  save (dat, file="raw_Abatch.rda")

##================#
##	Basic QC	#
##================#

#	library(affyQCReport)
##====Once is enough==========================
#	QCReport(dat, file="DataQC.pdf"); #JYL
##=============================================

#	eset = rma(dat)  # run RMA algorithm
#	pData(eset)
#	save(eset, file=paste('eset_prior_pvac','rda',sep="."))

load ("eset_prior_pvac.rda")

##==================================#
#	pvac-based Gene Filtering	#
##==================================#
library(pvac)
#	ft = pvacFilter(dat)
#	e2 = eset[ft$aset,]
#	save(eset,ft,file=paste('eset_post_pvac','rda',sep="."))
load ("eset_post_pvac.rda")

##========Get PVAC plot=======================================

plot(density(ft$pvac[ft$nullset]), xlab = "PVAC score", main = "", col = "gray",  cex.lab = 0.5, xlim = c(0, 1))
lines(density(ft$pvac), col = 1)
abline(v = ft$cutoff, lty = 2, col = "gray")
print(sessionInfo())
title ("PVAC filtering")
mtext("The verticle dotted line is the cut-off, 0.5")




##========================================#
#	Or, draw correlation plot only	#
##========================================#

toDrawCorrPlot = TRUE # set to FALSE if one wants to skip
if(toDrawCorrPlot){

##Unless needed, ignore this and go straight to correlation plot,
# This sph is missing, so, need to run e2 = eset[ft$aset,] first!!!

  sph = pData(dat)	
  pData(dat) <- sph$BioRep


  correlationPlot(dat)
  dev.off()
  pData(dat) = sph

}

##============================#
#	RMA normalization		#
##============================#
eset = rma(dat)  # run RMA algorithm
pData(eset)
save(eset, file=paste('eset_prior_pvac','rda',sep="."))


##==================================#
#	pvac-based Gene Filtering	#
##==================================#
library(pvac)
ft = pvacFilter(dat)
e2 = eset[ft$aset,]
save(eset,ft,file=paste('eset_post_pvac','rda',sep="."))

num.of.color <- length(dimnames(exprs(e2))[[2]])
plot.colors <- sample(1:200,num.of.color,replace=T) 
plot.colors

fig.title <- "Profile plot after RMA followed by PVAC"
sub.title <- ""

profile.plot(exprs(e2), fig.title, sub.title, plot.colors)


profile.plot <- function (dt, fig.title, fig.subtitle = "", plot.colors)
{
	sampleNames <- stripFileNames(colnames(dt))
	lgd <- sampleNames
	plotDensity (dt, main = fig.title, xlab = "Log(2) Intensity", ylab= "Density", type="n")
	mtext (fig.subtitle)
	for (i in 1:length(plot.colors))
	{
		lines(density(dt[,i]), col=colors()[plot.colors[i]], lty=2)
	}
	legend ("topleft",lgd, text.col = colors() [plot.colors], lty = c(rep(1,12)), col =  colors() [plot.colors])

} 

substr(colnames(exprs(e2))


correlationPlot(e2)
correlationPlot(dat)
cor(exprs(e2))


str(exprs(e2))
# num [1:18619, 1:36] 10.82 11.25 10.83 10.2 9.28 ...
str(exprs(dat))
 #num [1:1004004, 1:36] 
##========Get PVAC plot=======================================

plot(density(ft$pvac[ft$nullset]), xlab = "PVAC score", main = "", col = "gray",  cex.lab = 0.5, xlim = c(0, 1))
lines(density(ft$pvac), col = 1)
abline(v = ft$cutoff, lty = 2, col = "gray")
print(sessionInfo())
title ("PVAC filtering")
mtext("The verticle dotted line is the cut-off, 0.4285")

##==========================================================
#	pca, bio3d is not working with R3.0 and above
#	!!!
##==========================================================
library(bio3d)
data.pca <- pca.xyz(t(exprs(e2)))
save(data.pca,file=paste('pca_post_pvac','rda',sep="."))
plot.pca(data.pca, pch = as.vector(dimnames(data.pca$z)[[1]]), col=rainbow(length(dimnames(data.pca$z)[[1]])),  cex=0.8)


##==========================================
# Or, draw pca plots [samples] (Lu, Jun)
##=========================================
load("eset_post_pvac.rda")
e2 = eset[ft$aset,]
#load("eset_prior_pvac.rda")
#str(eset)

toDrawPCA = TRUE
if(toDrawPCA){
  pca <- prcomp(t(exprs(e2)),scale=TRUE,center=TRUE) # summary(pca)

  pc1 = pca$x[,"PC1"]; pc2 = pca$x[,"PC2"];  pc3 = pca$x[,"PC3"]

  #outf = paste("pca.BoneMarrow",d,"pdf",sep=".")
 # pdf(outf,height=3.5,width=3.5);  par(mar=c(3.5,3.3,1.2,1),mgp=c(2,0.5,0))


  labs = pData(e2)$BioRep
  col = rep('black',nrow(pData(e2)))
  col[pData(e2)$Condition=="L"] = "green"; 
  col[pData(e2)$Condition=="O"] = "red";
  col[pData(e2)$Condition=="PL"] = "blue";
  pa = pc1; pb = pc2 ; cex=0.8
  pa = pc1; pb = pc3 ; cex=0.8
  pa = pc2; pb = pc3 ; cex=0.8
	
  plot(pa,pb,main= "PCA on the PC1 & 2", xlab="PC1",ylab="PC2",pch="",cex=cex,cex.lab=cex)
  plot(pa,pb,main= "PCA on the PC1 & 3", xlab="PC1",ylab="PC3",pch="",cex=cex,cex.lab=cex)
  plot(pa,pb,main= "PCA on the PC2 & 3", xlab="PC2",ylab="PC3",pch="",cex=cex,cex.lab=cex)
  text(pa,pb,labels=labs,col=col,cex=cex)
  legend("bottomright",legend=c("PFPE_tissue","LN","OCT_KCM", "PFPE_LCM"),
         text.col=c("black","green","red", "blue"),bty="n",cex=cex)

 # dev.off()
}



##======================================#
##    limma-t test                      #
##======================================#
library(affy); library(limma);library(affycoretools)
library(mouse4302.db); library(annotate)


##===========Customized functions=============================================


get.glist <- function(fit2,cmp,pcut,fcut,fo)  ## Pull the DEG 
{
  	res = decideTests(fit2,p.value=pcut,adjust.method="BH",lfc=log2(fcut))
  	res = as.matrix(res)
  	af.ids = names(res[abs(res[,cmp])==1,cmp])
  	if(length(af.ids)==0)
	{
		print ("This is contrast   ")
		print (cmp)	
    		print("No genes meet the threshold! \n")
    		return(0)
  	}
  	sym = unlist(lookUp(af.ids,"mouse4302",'SYMBOL')); 
	af.ids = af.ids[!is.na(sym)]
  	if(length(af.ids)==0)
	{
		print ("This is contrast   ")
		print (cmp)
    		print("No genes meet the threshold! \n")
    		return(0)
  	}

  	aa = topTable(fit2,coef=cmp,number=Inf,adjust="BH"); rownames(aa) = aa$ID; aa$ID = NULL
  	aa = aa[rownames(aa)%in%af.ids,]
  	up.down = sign(aa$logFC)
  	fold = with(aa, 2^(abs(logFC)) * sign(logFC))
  	bb = cbind(aa[,c(4,5)],fold,up.down,
            symbol=unlist(lookUp(rownames(aa),annotation(eset),'SYMBOL')),
            description=unlist(lookUp(rownames(aa),annotation(eset),'GENENAME'))
         )
  	names(bb)[3] = paste('Fold',cmp,sep='.')
  	write.csv(bb,fo)
}




##====================End of functions=========================================
##	With pvac filtering
setwd("X:/project2014/Foley_Yu/workingDir")
load("eset_post_pvac.rda")
e2 = eset[ft$aset,]


pData(e2)
grp = factor(pData(e2)$Condition)
design = model.matrix(~0 + grp)
colnames(design) = levels(grp);  
design

cont.matrix = makeContrasts(

		  LN.vs.OCT     = L  - O,
		  LN.vs.PL      = L  - PL,
 		  LN.vs.Pt      = L  - Pt,
		  
		  OCT.vs.PL      = O  - PL,
 		  OCT.vs.Pt      = O  - Pt,

		  PL.vs.Pt       = PL - Pt,
              levels=design)
cont.matrix

fcut = 2;
pcut = 0.05;


fit = lmFit(exprs(e2),design)
fit2 = contrasts.fit(fit,cont.matrix);  
fit2 = eBayes(fit2) ; 
save(fit2,file="lmFitWcntr.rda")

cmps = c(	
		"LN.vs.OCT",    
		"LN.vs.PL",    
 		"LN.vs.Pt",   
		"OCT.vs.PL", 
		"OCT.vs.Pt", 
		"PL.vs.Pt"		
)   

for (cmp in cmps)
{
	if(pcut==1 && fcut==0)
	{
      	fo.glist = paste('glist',cmp,'noCUTOFF','csv',sep='.')
    	}else{
      	fo.glist = paste('glist',cmp,pcut,fcut,'csv',sep='.')
   	}
	get.glist(fit2,cmp,pcut,fcut,fo.glist)
}


##========================================
##	Investigating the gene lists
##========================================

DEGs.files <- list.files(pattern =  '\\.csv')
DEGs = c(list())

for (i  in (1: length(DEGs.files)))
{
	cat (DEGs.files[i]); cat ("\n")
	DEGs.name <- gsub ("glist.", "", DEGs.files[i])
	DEGs.name <- gsub (".0.05.2.csv", "", DEGs.name)
	dt <- read.csv(DEGs.files[i])
	DEGs[[i]] = list (name = DEGs.name, DEG = as.character(dt$symbol))
}

for (i  in (1: length(DEGs.files)))
{
	cat (i); cat ("\n")
	str(DEGs[i]); cat ("\n")
}

VennDiagram <- draw.three.list.venndigram(DEGs[[1]][2]$DEG, DEGs[[2]][2]$DEG, DEGs[[3]][2]$DEG,  DEGs[[1]][1], DEGs[[2]][1], DEGs[[3]][1])
grid.draw(VennDiagram$figure)



VennDiagram <- draw.four.list.venndigram(DEGs[[1]][2]$DEG, DEGs[[2]][2]$DEG, DEGs[[3]][2]$DEG,  DEGs[[4]][2]$DEG, DEGs[[1]][1], DEGs[[2]][1], DEGs[[3]][1], DEGs[[4]][1])
grid.draw(VennDiagram$figure)



VennDiagram <- draw.four.list.venndigram(DEGs[[2]][2]$DEG, DEGs[[4]][2]$DEG, DEGs[[5]][2]$DEG,  DEGs[[6]][2]$DEG, DEGs[[2]][1], DEGs[[4]][1], DEGs[[5]][1], DEGs[[6]][1])
grid.draw(VennDiagram$figure)
