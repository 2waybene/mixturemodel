library(RGtk2)
require(gWidgets)
require(tcltk)
options(guiToolkit="RGtk2")

main<- gwindow("MetAP")   
win <- ggroup(horizontal=FALSE,cont=main)
nb <- gnotebook(cont=win)

###LIST
#######"System Variables"##############


list <- glayout(cont=nb, label="System Variables")

list[1,1] <- "Directory of the code"
list[2,1] <- glabel("Code source directory:", cont=list)

list[2,2] <- gfilebrowse (text = "Select directory...", type = "selectdir", quote = FALSE,
              cont = list)

 
##radio button handler function
f<-function(h,...) { control=svalue(h$obj)
		   			 .GlobalEnv$control <- control
					if (svalue(h$obj)=="TRUE" || svalue(h$obj)=="rms" ){
			   			return (control)
                     } else { 
			   				return (control)
                         }
			 					return (control)
							}

##combining the lists and running the script
run <- function(h,...){
c(mylist,mylist2,mylist3,mylist4,mylist5)
source(file.path(mylist[1],"run_analysis.R"))

}

list[3,1] <- "Files present"
list[4,1] <- "Group 1:"
controlrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list,
		   handler=f)

svalue(controlrb,index=TRUE) <- 1
list[4,2] <- controlrb

list[5,1] <- "Group 2:"
subjectrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list,
		   handler=f)

svalue(subjectrb,index=TRUE) <- 1
list[5,2] <- subjectrb 

#list[6,1] <- "Control background :"
#input1rb <- gradio(c("TRUE","FALSE"),
#		   selected=2,index=FALSE, 
#	           horizontal=TRUE, cont=list,
#		   handler=f)

#svalue(input1rb,index=TRUE) <- 2
#list[6,2] <- input1rb 

#list[7,1] <- "Subject background :"
#input2rb <- gradio(c("TRUE","FALSE"),
#		   selected=2,index=FALSE, 
#	           horizontal=TRUE, cont=list,
#		   handler=f)

#svalue(input2rb,index=TRUE) <- 2
#list[7,2] <- input2rb

list[8,1] <- "Select Which Analysis To Perform"
list[9,1] <- "Quality Control:"
qcontrolrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list,
		   handler=f)

svalue(qcontrolrb,index=TRUE) <- 1
list[9,2] <- qcontrolrb

list[10,1] <- "Differential Methylation Analysis:"
analysisrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list,
		   handler=f)

svalue(analysisrb,index=TRUE) <- 1
list[10,2] <- analysisrb


list[11,1] <- "Cluster Analysis:"
clstrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list,
		   handler=f)

svalue(clstrb,index=TRUE) <- 1
list[11,2] <- clstrb


list[12,1] <- "Calibration:"
calibrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list,
		   handler=f)

svalue(calibrb,index=TRUE) <- 1
list[12,2] <- calibrb


#list[13,1] <- "Correlation plots:"
#corrrb <- gradio(c("TRUE","FALSE"),
#		   selected=2,index=FALSE, 
#	           horizontal=TRUE, cont=list,
#		   handler=f)

#svalue(corrrb,index=TRUE) <- 1
#list[13,2] <-l corrrb


list[14,1] <- "Data folders"
list[15,1] <- "Group 1 :"
list[15,2] <- gfilebrowse (text = "Select folder...", type = "selectdir", quote = FALSE,
              cont = list)

list[16,1] <- "Group 2 :"
list[16,2] <- gfilebrowse (text = "Select folder...", type = "selectdir", quote = FALSE,
              cont = list)

#list[17,1] <- "Background control folder"
#list[17,2] <- gfilebrowse (text = "Select folder...", type = "selectdir", quote = FALSE,
#              cont = list)
#list[18,1] <- "Background subject folder"
#list[18,2] <- gfilebrowse (text = "Select folder...", type = "selectdir", quote = FALSE,
#              cont = list)

list[19,1] <- "Analysis"
list[20,1] <- "Output Directory" ##Met_test
list[20,2] <- gfilebrowse (text = "Select folder...", type = "selectdir", quote = FALSE,
              cont = list)
list[21,1] <- "Number of Cores to use"
list[21,2] <- gedit("1", cont=list, coerce.with = as.numeric)

finalsave <- function(h,...){

	x <- svalue(list[2,2])
    	x2 <- svalue(list[15,2])
    	x3 <- svalue(list[16,2])
		#x4 <- svalue(list[17,2])
		#x5 <- svalue(list[18,2])
	x6 <- svalue(list[20,2])
	x8 <- svalue(list[21,2])
	x4 <- svalue(list2[2,2])
   	x5 <- svalue(list2[3,2])
    	x7 <- svalue(list2[4,2])
	

	y1 <- svalue(list3[2,2])
   	y2 <- svalue(list3[3,2])
    	y3 <- svalue(list3[4,2])
	y4 <- svalue(list3[5,2])
	y5 <- svalue(list3[6,2])
	y6 <- svalue(list3[7,3])
	y7 <- svalue(list3[7,5])
	y8 <- svalue(list4[4,2])
	y9 <- svalue(list4[6,2])
	y10 <- svalue(list4[7,2])
	z1 <- svalue(list5[7,2])
    	z2 <- svalue(list5[6,3])
	z3 <- svalue(list5[6,5])
    if(x!="" && x2!="" && x3!="" && x6!="" && x!="Select directory..." && x2!="Select folder..." && x3!="Select folder..." && x6!="Select folder..." && !is.na(x4)==TRUE && !is.na(x5)==TRUE && !is.na(x7)==TRUE && x5!=0 && x7!=0 && !is.na(y1)==TRUE && !is.na(y2)==TRUE && !is.na(y3)==TRUE && !is.na(y4)==TRUE && !is.na(y5)==TRUE && !is.na(y6)==TRUE && !is.na(y7)==TRUE && !is.na(y8)==TRUE && !is.na(y9)==TRUE && !is.na(y10)==TRUE && y10>=1 && !is.na(z1)==TRUE && !is.na(z2)==TRUE && !is.na(z3)==TRUE && z2!=0 && z3!=0 && x8 >=1){
    	gmessage("Data saved!", parent=win)

mylist <- list(source.dir=svalue(list[2,2]),control=svalue(controlrb),subject=svalue(subjectrb),#c.input=svalue(input1rb),s.input=svalue(input2rb),
Quality.Control=svalue(qcontrolrb),Analysis=svalue(analysisrb),calibration=svalue(calibrb),to.cluster=svalue(clstrb),#correlation.plots=svalue(corrrb), 
folder.control=svalue(list[15,2]),folder.subject=svalue(list[16,2]), #folder.control.background=svalue(list[17,2]),folder.subject.background=svalue(list[18,2]), 
base.directory=svalue(list[20,2]), num_cores=svalue(list[21,2]))
.GlobalEnv$mylist <- mylist

mylist2 <- list(shift=svalue(list2[2,2]), extend=svalue(list2[3,2]), ws=svalue(list2[4,2]),uniq=svalue(uniqrb),BSgenome=svalue(list2[6,2]),chr.select=svalue(list2[7,2]))
.GlobalEnv$mylist2 <- mylist2

	x <- svalue(list3[7,3])
	y <- svalue(list3[7,5])
	r <-c(x:y)
	
mylist3 <- list(cpg.screen=svalue(list3[2,2]), saturation.screen=svalue(list3[3,2]), seq.coverage=svalue(list3[4,2]),nit=svalue(list3[5,2]),nrit=svalue(list3[6,2]),cov.level=r)
.GlobalEnv$mylist3 <- mylist3

mylist4 <- list(p.adj=svalue(list4[2,2]), type=svalue(typerb), minRowSum=svalue(list4[4,2]),p.value=svalue(list4[6,2]), ratio=svalue(list4[7,2]))
.GlobalEnv$mylist4 <- mylist4

r <-seq(z2,z3, by=0.1)
mylist5 <- list(clust.method=svalue(list5[2,2]),dist.method=svalue(list5[3,2]),resample.depth=r, bootstrap=svalue(list5[7,2]))
.GlobalEnv$mylist5 <- mylist5
}

	else if(x=="" || x2=="" || x3=="" || x6=="" || x=="Select directory..." || x2=="Select folder..." || x3=="Select folder..." || x6=="Select folder...") {
    	gmessage("All field(s) under 'System Variables' must be filled out", parent=win)
	}

	else if(is.na(x4)==TRUE || is.na(x5)==TRUE || is.na(x7)==TRUE || x5==0 || x7==0) {
    	gmessage("Please check entered value(s) under 'Global Variables'", parent=win)
	}
	
	else if(is.na(y1)==TRUE || is.na(y2)==TRUE || is.na(y3)==TRUE || is.na(y4)==TRUE || is.na(y5)==TRUE || is.na(y6)==TRUE || is.na(y7)==TRUE) {
    	gmessage("Please check entered value(s) under 'Quality Control'", parent=win)
	}

	else if(is.na(y8)==TRUE || is.na(y9)==TRUE || is.na(y10)==TRUE || y10<1) {
    	gmessage("Please check entered value(s) under 'Profile Variables'", parent=win)
	}

	else if(is.na(z1)==TRUE || is.na(z2)==TRUE || is.na(z3)==TRUE || z2==0 || z3==0) {
    	gmessage("Please check entered value(s) under 'Clustering Analysis options'", parent=win)
	}
	
	else if(x8 == 0)
	{
	gmessage("0 Cores Selected, Please set number of cores to perform Analysis on", parent = win, icon = "error")
	}

	
}
			
#list[21,3]<-gbutton("Save", cont=list, handler=save1)
#list[21,4]<-gbutton("Cancel", cont=list, handler = function(h,...) dispose(win))


##LIST2
#######"Global Variables to be used for analysis"##########
list2 <- glayout(cont=nb, label="Global Variables")

list2[2,1] <- "Shift Reads?"
list2[2,2] <- gedit("0", cont=list2, coerce.with =as.numeric)
list2[3,1] <- "Extend Reads?"
list2[3,2] <- gedit("300", cont=list2, coerce.with =as.numeric)
list2[4,1] <- "Window Size:"
list2[4,2] <- gedit("500", cont=list2, coerce.with =as.numeric)


list2[5,1] <- "Select only Unique Reads?"
uniqrb <- gradio(c("TRUE","FALSE"),
		   selected=2,index=FALSE, 
	           horizontal=TRUE, cont=list2,
		   handler=f)

svalue(uniqrb,index=TRUE) <- 1
list2[5,2] <- uniqrb

list2[6,1] <- "Reference BSgenome:"
list2[6,2] <- gedit("BSgenome.Hsapiens.UCSC.hg19", cont=list2)

list2[7,1] <- "Select a chromosome (separate by ,):"
list2[7,2] <- gedit("all", cont=list2)

#save2 <- function(){
#   x <- svalue(list2[2,2])
#    x2 <- svalue(list2[3,2])
#    x3 <- svalue(list2[4,2])
#    if(!is.na(x)==TRUE && !is.na(x2)==TRUE && !is.na(x3)==TRUE && x2!=0 && x3!=0){

    	#gmessage("Data saved!", parent=win)
#mylist2 <- list(shift=svalue(list2[2,2]), extend=svalue(list2[3,2]), ws=svalue(list2[4,2]),uniq=svalue(uniqrb),BSgenome=svalue(list2[6,2]),chr.select=svalue(list2[7,2]))
#.GlobalEnv$mylist2 <- mylist2}
	
#	else {
#    	gmessage("Invalid input2", parent=win)
#	}

#}

#list2[43,16]<-gbutton("Save", cont=list2, handler=save2)
#list2[43,17]<-gbutton("Cancel", cont=list2, handler=function(h,...) dispose(win))


##LIST3
#####Quality Control Variables and Options #######
list3 <- glayout(cont=nb, label="Quality Control")
list3[1,1] <- "Filters"

list3[2,1] <- "Minimum CpG Enrichment Threshold :"
list3[2,2] <- gedit("1.7", cont=list3, coerce.with =as.numeric,width=3)

list3[3,1] <- "Minimum Saturation Threshold :"
list3[3,2] <- gedit("0.5", cont=list3, coerce.with =as.numeric,width=3)

list3[4,1] <- "Minimum Sequence Coverage Threshold :"
list3[4,2] <- gedit("0.05", cont=list3, coerce.with =as.numeric,width=3)

list3[5,1] <- "Number of Iterations :"
list3[5,2] <- gedit("5", cont=list3, coerce.with =as.numeric,width=3)

list3[6,1] <- "Number of Repeat Iterations :"
list3[6,2] <- gedit("1", cont=list3, coerce.with =as.numeric,width=3)

list3[7,1] <- "coverage level :"
list3[7,2] <- "start value"
list3[7,3] <- gedit("0", cont=list3, coerce.with =as.numeric, width=3)
list3[7,4] <- "end value" 
list3[7,5] <- gedit("5", cont=list3, coerce.with =as.numeric, width=3)

#save3 <- function(){
#	x1 <- svalue(list3[2,2])
#    x2 <- svalue(list3[3,2])
#    x3 <- svalue(list3[4,2])
#	x4 <- svalue(list3[5,2])
#	x5 <- svalue(list3[6,2])
#	x6 <- svalue(list3[7,3])
#	x7 <- svalue(list3[7,5])
#    if(!is.na(x1)==TRUE && !is.na(x2)==TRUE && !is.na(x3)==TRUE && !is.na(x4)==TRUE && !is.na(x5)==TRUE && !is.na(x6)==TRUE && !is.na(x7)==TRUE) {
	
    	#gmessage("Data saved!", parent=win)

#	x <- svalue(list3[7,3])
#	y <- svalue(list3[7,5])
#	r <-c(x:y)
	
#mylist3 <- list(cpg.screen=svalue(list3[2,2]), saturation.screen=svalue(list3[3,2]), seq.coverage=svalue(list3[4,2]),nit=svalue(list3[5,2]),nrit=svalue(list3[6,2]),cov.level=r)
#.GlobalEnv$mylist3 <- mylist3}

#	else{
#    	gmessage("Invalid input3", parent=win)
#	}
#}

#list3[42,19]<-gbutton("Save", cont=list3, handler=save3)
#list3[42,20]<-gbutton("Cancel", cont=list3, handler=function(h,...) dispose(win))


##LIST4
#######Methylation Profile Variables #########
list4 <- glayout(cont=nb, label="Profile Variables")

list4[2,1] <- "P-Value Adjustment Method :"
vals<-c("fdr","holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none")
list4[2,2] <- gdroplist(vals, cont=list4)

list4[3,1] <- "Data to Use:"
typerb <- gradio(c("rms","rpkm"),
		   selected=2,index=FALSE, 
	       horizontal=TRUE, cont=list2,
		   handler=f)

svalue(typerb,index=TRUE) <- 2
list4[3,2] <- typerb

list4[4,1] <- "Minimum Reads per Window :"
list4[4,2] <- gedit("1", cont=list4, coerce.with =as.numeric)

list4[5,1] <- "Differential Methylation Selection Filters"

list4[6,1] <- "P value Threshold:"
list4[6,2] <- gedit("0.01", cont=list4, coerce.with =as.numeric)

list4[7,1] <- "Minimum Expression Ratio :"
list4[7,2] <- gedit("1", cont=list4, coerce.with =as.numeric)
list4[8,2] <- "(Not less than '1')"
 
#save4 <- function(){
#	x1 <- svalue(list4[4,2])
#	x2 <- svalue(list4[6,2])
#	x3 <- svalue(list4[7,2])
#    if(!is.na(x1)==TRUE && !is.na(x2)==TRUE && !is.na(x3)==TRUE && x3>=1)  {
	
    	#gmessage("Data saved!", parent=win)
	
#mylist4 <- list(p.adj=svalue(list4[2,2]), type=svalue(typerb), minRowSum=svalue(list4[4,2]),p.value=svalue(list4[6,2]), ratio=svalue(list4[7,2]))
#.GlobalEnv$mylist4 <- mylist4
#}
#	else{
#    	gmessage("Invalid input4", parent=win)
#	}
#}


#list4[43,11]<-gbutton("Save", cont=list4, handler=save4)
#list4[43,12]<-gbutton("Cancel", cont=list4, handler=function(h,...) dispose(win))


##LIST5
#######Clustering Analysis options #########
list5 <- glayout(cont=nb, label="Clustering Analysis options")

list5[1,1] <- "General Variables"
list5[2,1] <- "Cluster Method :"
clstvals <- c("mcquitty", "ward", "single", "complete", "average", "median", "centroid")
list5[2,2] <- gdroplist(clstvals, cont=list5)

list5[3,1] <- "distance matrix method :"
distvals <- c("correlation" , " ")
list5[3,2] <- gdroplist(distvals, cont=list5)

list5[5,1] <- "Pvclust options"
list5[6,1] <- "Resample depth :"
list5[6,2] <- "start value"
list5[6,3] <- gedit("0.5", cont=list5, coerce.with =as.numeric, width=3)
list5[6,4] <- "end value"
list5[6,5] <- gedit("1.4", cont=list5, coerce.with =as.numeric, width=3)

list5[7,1] <- "Bootstrap Depth :"
list5[7,2] <- gedit("100", cont=list5, coerce.with =as.numeric, width=3)


#save5 <- function(){
	
#	z <- svalue(list5[7,2])
#	x <- svalue(list5[6,3])
#	y <- svalue(list5[6,5])
#	if(!is.na(x)==TRUE && !is.na(y)==TRUE && !is.na(z)==TRUE && y!=0 && z!=0){

		#gmessage("Data saved!", parent=win)
#		r <-seq(x,y, by=0.1)

#mylist5 <- list(clust.method=svalue(list5[2,2]),dist.method=svalue(list5[3,2]),resample.depth=r, bootstrap=svalue(list5[7,2]))
#.GlobalEnv$mylist5 <- mylist5}

#	else{
#		gmessage("Invalid input5", parent=win)
#	}
#}

#finalsave <- function(h,...){
#save1()
#save2()
#save3()
#save4()
#save5()
#}
#list5[43,5]<-gbutton("Save", cont=win, handler=save5)
#list5[43,6]<-gbutton("Run script", cont=win, handler=run)
#list5[43,7]<-gbutton("Cancel", cont=win, handler=function(h,...) dispose(win))
lyt<-ggroup( Horizontal=TRUE, cont=win)

#lyt <- glayout(cont=win)
lyt.save <-gbutton("Save", cont=lyt, handler=finalsave)
lyt.run<-gbutton("Run script", cont=lyt, handler=run)
lyt.cancel <-gbutton("Cancel", cont=lyt, handler=function(h,...) dispose(win))

gtkMain()
#gtkMainQuit()
