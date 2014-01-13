##======================================================
#	File:   comparativeExpression.R
#		
#	Author: Jianying Li
#	
#	Date:   04/08/2013
#	Comment: Modified based on ASE_gene_FPKM.R
##======================================================
source("C:/Users/li11/Documents/R-project/customPackages/dataManipTools.R")
source("C:/Users/li11/Documents/R-project/customPackages/arraySeqTools.R")
source("X:/project2012/ASE/R-scripts/functions-ASE.R")
source ("x:/project2013/miRNAseq/scripts/miRNA_tools.R")
library(Rlab)
library(gplots)

load("X:/project2012/ASE/FPKM_RNAseq/ASE_w_FPKM.rda")
str(ASE.dt)
rownames(ASE.dt)

##====================================================
#	Getting directory and file set up
##====================================================
setwd ("X:/project2013/comparativeExpression/workingDir/")
dataDir <- "X:/project2013/comparativeExpression/dataDir/"


##===============possibly change======================
firstTrialDir<- "X:/project2013/comparativeExpression/firstTrial/"
filePath <- dataDir
outDir <- firstTrialDir
matrix.file <- "combined_FPKM_genes_mm9.txt"
expInfo <- "expInfo_mm9.txt"
##=======================================================================


##=====================================
##	Read in expInfo and data file
##=====================================

data.dm <- paste (filePath, matrix.file, sep="")
dm <- readarray(data.dm)
dt <- dm$xd
file.names <- colnames(dt)
file.names

str(dt)
raw.FPKM <- dt
plotDensity(log2(dt))
bplot(log2(dt))


expInfo.dm <- paste (filePath, expInfo, sep="")
p.dt <- read.table(expInfo.dm, header = TRUE, sep="\t", row.names=as.vector(file.names))
str(p.dt)
rownames(p.dt)

##====================================================#
#	Processing RNAseq FPKM data				#
##====================================================#
dim(dt)
no.zero.dt <- filtering.zero.dm(dt)
dim(no.zero.dt)
log2.dt <- log2(no.zero.dt)


##	or working with raw dt
#	log2.dt <- log2(dt)
#	length(which (!is.finite(log2.dt)))



log2.dt[which (!is.finite(log2.dt))] <- 0.00000000001  #replacing missings with 10e-11
str(log2.dt)

#processed.dt <- row.filtering(log2.dt, -10, 40) # Overloaded method has been rewritten as row.filtering.dm 

processed.dt <- row.filtering.dm (log2.dt, -10, 40)
processed.dt <- row.flt.boundary(processed.dt, -0.001, 0.001, 40)
str(processed.dt)
stats(processed.dt)
plotDensity(processed.dt)

##========================================================

##================================================
##	Getting an annotated profile plot
##==========================================

num.of.color <- length(dimnames(dt)[[2]])
plot.colors <- sample(1:200,num.of.color,replace=T) 
plot.colors

fig.title <- "Profile plot for log2 FPKM on cuffline with RefSeq"
fig.title <- "Profile plot for log2 FPKM on cuffline with RefSeq 0 filtered"
sub.title <- ""

profile.plot.w.legend(processed.dt, fig.title, sub.title, plot.colors, as.vector(p.dt$Header))



##=========================================================
##	What about correlation?
##	MA-plots?
##
##==========================================================

##=========================================================
##	Get a correlation half heatmap plot on genes
##=========================================================
#pairs(processed.dt)
#pearsons <-  cor(dt)
pearsons <-  cor(processed.dt)

str(pearsons)
pearsons.log2 <- paste (outDir,"pairwise-corr.txt" , sep="")
pearsons.log2 <- paste (outDir,"pairwise-corr-noZero.txt" , sep="")
write.table(pearsons, pearsons.log2, sep="\t")

LDheatmap(pearsons, color = rainbow(75))
heatmap.2(pearsons)
rownames(pearsons)
##====================================================
#	Linear model fitting
##=====================================================
str(p.dt)
grp = factor( paste(p.dt$Gender,p.dt$Specie, p.dt$Generation, sep=".")); 
design = model.matrix(~0 + grp)
colnames(design) = levels(grp);  
fit <-  lmFit(processed.dt,design)


##=====================
#	ANOVA test
##=====================
fit2 = eBayes(fit) ;
topTable(fit2)

#=================================
#	Look for contrasts
#=================================

##===========================
#	Gender difference
##===========================

cont.matrix = makeContrasts(F.vs.M = (F.B.O + F.B.P + F.C.O + F.C.P) - (M.B.O + M.B.P + M.C.O + M.C.P), levels =design)
#cont.matrix
fit2 = contrasts.fit(fit,cont.matrix);  
fit2 = eBayes(fit2) ;
topTable(fit2)
aa <- topTable(fit2,number=Inf,adjust="BH")
str(aa[which(aa$adj.P.Val < 0.01),])

outFile <- paste (outDir, "overall_gender_diff_BH01.txt", sep ="")
outResults <- aa[which(aa$adj.P.Val < 0.01),]
write.table (outResults, outFile,  sep = "\t",  row.names= FALSE)
outResults
zoom.in.dt<-processed.dt[(rownames(processed.dt)%in%outResults$ID),]
ASE.Heatmap(zoom.in.dt, "Between gender difference")
str(outResults)
str(zoom.in.dt)
DEG.overall.gender <- outResults$ID
str(DEG.overall.gender)

temp <- aa[which(aa$adj.P.Val < 0.01),]
temp.2 <- aa[which(abs(temp$logFC) > 2),]
DEG.overall.gender <- temp.2$ID

##============================
#	Generation difference
##============================

cont.matrix = makeContrasts(O.vs.P = (-F.B.O +  F.B.P - F.C.O + F.C.P - M.B.O + M.B.P - M.C.O + M.C.P), levels =design)
fit2 = contrasts.fit(fit,cont.matrix);  
fit2 = eBayes(fit2) ;
aa <- topTable(fit2,number=Inf,adjust="BH")
str(aa[which(aa$adj.P.Val < 0.01),])


outFile <- paste (outDir, "overall_generation_diff_BH01.txt", sep ="")
outResults <- aa[which(aa$adj.P.Val < 0.01),]
write.table (outResults, outFile,  sep = "\t",  row.names= FALSE)
outResults
zoom.in.dt<-processed.dt[(rownames(processed.dt)%in%outResults$ID),]
ASE.Heatmap(zoom.in.dt, "Between generation difference")


##==========================================
#	Gender difference in parental
##==========================================

cont.matrix = makeContrasts(F.P.vs.M.P = (F.B.P + F.C.P - M.B.P - M.C.P), levels =design)
fit2 = contrasts.fit(fit,cont.matrix);  
fit2 = eBayes(fit2) ;
aa <- topTable(fit2,number=Inf,adjust="BH")
str(aa[which(aa$adj.P.Val < 0.01),])


outFile <- paste (outDir, "gender_diff_in_parental_BH01.txt", sep ="")
outResults <- aa[which(aa$adj.P.Val < 0.01),]
write.table (outResults, outFile,  sep = "\t",  row.names= FALSE)
outResults
zoom.in.dt<-processed.dt[(rownames(processed.dt)%in%outResults$ID),]
ASE.Heatmap(zoom.in.dt, "Gender difference in parentals")
DEG.gender.in.parental  <- outResults$ID

temp <- aa[which(aa$adj.P.Val < 0.01),]
temp.2 <- aa[which(abs(temp$logFC) > 2),]
DEG.gender.in.parental  <- temp.2$ID

##==========================================
#	Gender difference in off springs
##==========================================

cont.matrix = makeContrasts(F.O.vs.M.O = (F.B.O + F.C.O ) - (M.B.O +  M.C.O ), levels =design)

fit2 = contrasts.fit(fit,cont.matrix);  
fit2 = eBayes(fit2) ;
aa <- topTable(fit2,number=Inf,adjust="BH")
str(aa[which(aa$adj.P.Val < 0.01),])


outFile <- paste (outDir, "gender_diff_in_offsprings_BH01.txt", sep ="")
outResults <- aa[which(aa$adj.P.Val < 0.01),]
write.table (outResults, outFile,  sep = "\t",  row.names= FALSE)
outResults
zoom.in.dt<-processed.dt[(rownames(processed.dt)%in%outResults$ID),]
ASE.Heatmap(zoom.in.dt, "Gender difference in offsprings")
DEG.gender.in.offsprings <- outResults$ID

temp <- aa[which(aa$adj.P.Val < 0.01),]
temp.2 <- aa[which(abs(temp$logFC) > 2),]
DEG.gender.in.offsprings <- temp.2$ID

##========================================
#	Strain difference in parental
##========================================

cont.matrix = makeContrasts(C.P.vs.B.P = (F.B.P - F.C.P + M.B.P - M.C.P), levels =design)
fit2 = contrasts.fit(fit,cont.matrix);  
fit2 = eBayes(fit2) ;
aa <- topTable(fit2,number=Inf,adjust="BH")
str(aa[which(aa$adj.P.Val < 0.01),])


outFile <- paste (outDir, "strain_diff_in_parental_BH01.txt", sep ="")
outResults <- aa[which(aa$adj.P.Val < 0.01),]
write.table (outResults, outFile,  sep = "\t",  row.names= FALSE)
outResults
zoom.in.dt<-processed.dt[(rownames(processed.dt)%in%outResults$ID),]
ASE.Heatmap(zoom.in.dt, "Strain B6.vs.C3 diff in parentals")

##==============================================
##	Draw a venn diagram of gender difference
##==============================================
library(VennDiagram)

venn.plot <- draw.triple.venn (

	area1 = length(DEG.overall.gender) ,
	area2 = length(DEG.gender.in.parental ),
	area3 = length(DEG.gender.in.offsprings),

	n12 = length(get_common_list(DEG.overall.gender, DEG.gender.in.parental )),
	n13 = length(get_common_list(DEG.overall.gender,  DEG.gender.in.offsprings)) ,
	n23 = length(get_common_list(DEG.gender.in.parental ,  DEG.gender.in.offsprings)),

	n123 =  length(get_common_list(get_common_list(DEG.overall.gender, DEG.gender.in.parental ), DEG.gender.in.offsprings)) ,
	
	category = c("Overall", "Parental", "OffSprings"),
	fill = c("blue", "red", "green"), 
	lty = "blank", 
	cex = 2,
	cat.cex = 2,
	cat.col = c("blue", "red", "green"),
	
);

title ="Comparative expression gender difference"
getwd()
#tiff(filename = "DEG_genderDiff_2fc_BH01.tiff", compression = "lzw");
grid.draw(venn.plot);
#dev.off()



##===========OVERLAP list=================================================================================
DEG.gender.overlap.list <- get_common_list(get_common_list(DEG.overall.gender, DEG.gender.in.parental ), DEG.gender.in.offsprings)
length( DEG.gender.overlap.list)
DEG.gender.overlap.dt <- processed.dt[(rownames(processed.dt)%in%DEG.gender.overlap.list),]
colnames(DEG.gender.overlap.dt) <- sub (".mm9.genes.fpkm_tracking.fpkm.txt", "", colnames(DEG.gender.overlap.dt))
ASE.Heatmap (DEG.gender.overlap.dt, "Overlap list")



##==========================================
##	Consult ASE list
##==========================================
gene.names <- rownames(DEG.gender.overlap.dt)
gene.names
gene.only <- strsplit(gene.names, "\\|")
length(gene.only)
gene.ids <- c()
for (i in 1:length(gene.only))
{
	gene.ids[i] <- gene.only[[i]][1]
}

length(gene.ids)

rownames(DEG.gender.overlap.dt) <- gene.ids
ASE.genes <- rownames(ASE.dt)
rownames(DEG.gender.overlap.dt)[(rownames(DEG.gender.overlap.dt) %in% ASE.genes)]

#ONLY one gene was found ASE!!!
#"Cyp2c38"

##  No need for heatmap!!!
#ASE.Heatmap (DEG.gender.overlap.dt[(rownames(DEG.gender.overlap.dt) %in% ASE.genes),], "Overlap list")


##===========UNION list=================================================================================
DEG.gender.union.list <-  unique(c(DEG.overall.gender, DEG.gender.in.parental,DEG.gender.in.offsprings))
length(DEG.gender.union.list)
DEG.gender.union.dt <- processed.dt[(rownames(processed.dt)%in%DEG.gender.union.list),]
colnames(DEG.gender.union.dt) <- sub (".mm9.genes.fpkm_tracking.fpkm.txt", "", colnames(DEG.gender.union.dt))
ASE.Heatmap (DEG.gender.union.dt, "Union list")


gene.names <- rownames(DEG.gender.union.dt)
gene.names
gene.only <- strsplit(gene.names, "\\|")
length(gene.only)
gene.ids <- c()
for (i in 1:length(gene.only))
{
	gene.ids[i] <- gene.only[[i]][1]
}

length(gene.ids)

rownames(DEG.gender.union.dt) <- gene.ids
ASE.genes <- rownames(ASE.dt)
rownames(DEG.gender.union.dt)[(rownames(DEG.gender.union.dt) %in% ASE.genes)]
ASE.Heatmap (DEG.gender.union.dt[(rownames(DEG.gender.union.dt) %in% ASE.genes),], "Overlap between union DEGs list and ASE")

##=================================================
##	Do some assessment on DEG
##=================================================
##	Heatmap

ASE.Heatmap (DEG.gender.overlap.dt, "Overlap list")
ASE.Heatmap (DEG.gender.union.dt, "Union list")

##	PCA

library(bio3d)
data.pca <- pca.xyz(t(DEG.gender.overlap.dt))

cols <- rep("red", 24)
cols[which(colnames(DEG.gender.overlap.dt)=="_F_")] <- "blue"
cols[which(colnames(DEG.gender.overlap.dt)=="_M_")] <- "yellow"


plot.pca.w.label(data.pca,  pch = as.vector(dimnames(data.pca$z)[[1]]), col=cols,  cex=0.8)

##====================================================
#	QuasiSeq model fitting
##=====================================================
library(QuasiSeq)
library(imputation)

str(raw.FPKM)
stats(raw.FPKM)

temp <- raw.FPKM
temp[which(temp ==0)] <- NA
temp.2 <- row.filtering(temp)
temp.3 <- kNNImpute(as.matrix(temp.2),9)
temp.4 <- temp.3$x
temp.5 <- flt.any.0(temp.4)
temp.5 <- ceiling(temp.5)


threshold <- apply(temp.5,2,sum)
threshold



dt.no.NA <- temp.5 # I need data
str(dt.no.NA)
#pairs(log2(dt.no.NA))
str(p.dt)

grp = paste(p.dt$Gender, p.dt$Specie, p.dt$Generation, sep = "_")

design.list<-vector("list",2)
design.list[[1]]<-model.matrix(~as.factor(grp)) #This also could have just been ‘‘trt’’.
design.list[[2]]<-rep(1,length(rownames(p.dt)))

### Analyze using QL, QLShrink and QLSpline methods applied to quasi-Poisson model
fit<-QL.fit(dt.no.NA, design.list,log.offset=log2(threshold), Model="Poisson")
results<-QL.results(fit)
results$P.values[[3]]


### How many significant genes at FDR=.05 from QLSpline method?
## No need as there was NO significnat p-value!!

apply(results$Q.values[[3]]<.05,2,sum)

### Indexes for Top 100 most significant genes from QLSpline method
order(results$P.values[[3]])[1:100]


##===================================================================
### applied to quasi-negative binomial model

fit2<-QL.fit(dt.no.NA, design.list,log.offset=log2(threshold), Model="NegBin")
results2<-QL.results(fit2)
results$P.values[[3]]

### How many significant genes at FDR=.05 from QLSpline method?
## No need as there was NO significnat p-value!!

apply(results$Q.values[[3]]<.05,2,sum)

str(results$Q.values[[3]]<.05)
### Indexes for Top 100 most significant genes from QLSpline method
order(results$P.values[[3]])[1:100]




