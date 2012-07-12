rm(list=ls())
source("K_calcs_Johnson_OS.R")
library(Hmisc)
library(MASS)
library(lhs)


##############################################################
# create pdfs of fluxes based on data frames of temperature, #
# salinity, windspeed and air and water phase concentrations #
# of a gas                                                   #
##############################################################


#list to specify ditribution type (can be normal or lognormal)
defaultdistlist<-list(
			"Salinity"="normal",
			"Temperature"="normal",
			"C_air"="lognormal",
			"C_water"="lognormal",
			"WindSpeed"="lognormal"
)

##############################################################
# helper functions


getndist_within_bounds<-function(mean,sd,conf=1.0,n=1000){
	#function to remove values above given confidence limit when generating random samples from normal distribution. Default is no substitution
	y<-rnorm(n,mean=mean,sd=sd)
	ulim<-qnorm(conf,mean,sd)
	llim<-qnorm(1-conf,mean,sd)
	count=0
	#substitute max or min value for a new one sampled from the distribution as required
	while(max(y)>ulim|min(y)<llim){
		
		while(min(y)<llim){				
			y[which(y==min(y))]<-rnorm(1,mean=mean,sd=sd)
			count<-count+1}
		while(max(y)>ulim){
			y[which(y==max(y))]<-rnorm(1,mean=mean,sd=sd)
			count<-count+1}	
	}
	print(paste(count,"substitutions"))
	y
}


getlndist_within_bounds<-function(meanlog,sdlog,conf=1.0,n=1000){
	#function to remove values above given confidence limit when generating random samples from lognormal distribution
	y<-rlnorm(n,mean=meanlog,sd=sdlog)
	ulim<-qlnorm(conf,meanlog,sdlog)
	llim<-qlnorm(1-conf,meanlog,sdlog)
	count=0
	while(max(y)>ulim|min(y)<llim){
		#substitute max value for a new one sampled from the distribution
		while(min(y)<llim){				
			y[which(y==min(y))]<-rlnorm(1,mean=meanlog,sd=sdlog)
			count<-count+1}
		while(max(y)>ulim){
			y[which(y==max(y))]<-rlnorm(1,mean=meanlog,sd=sdlog)
			count<-count+1}	
	}

	print(paste(count,"substitutions"))
	y
}

randomsamplelhc<-function(compound, dist_n, lhc_n, rdwind, rdsal, rdtemp, rdcair, rdcwater){
	#from randomly sampled distributions of variables subset a latn hypercube and calculate air-sea fluxes for each combination of variables

	#create a latin hypercube with numbers ranging from 1 to dist_n to sub-sample the distributions
	lh<-round((dist_n-1)*randomLHS(lhc_n,k=5))+1
	
	#sort each random sample distribution and then subsample by indexing against lhc
	x<-data.frame(sort(rdwind)[lh[,1]])
	colnames(x)<-c("u")
	x$S<-sort(rdsal)[lh[,2]]
	x$T<-sort(rdtemp)[lh[,3]]
	x$gconc<-sort(rdcair)[lh[,4]]
	x$swconc<-sort(rdcwater)[lh[,5]]
	# Henry's law
	x$KH<-KH(compound,x$T,x$S)
	# delta C
	x$dC<-x$swconc-(x$gconc/x$KH)
	#transfer velocity	
	x$Kw<-Kw(compound,x$T,x$u,x$S)
	#flux
	x$F<-x$Kw*x$dC
	#return x	
	x
	
}	

repeatlhc<-function(compound, dist_n, lhc_n, repeatlhc_n, rdwind, rdsal, rdtemp, rdcair, rdcwater){
	y<-randomsamplelhc(compound, dist_n, lhc_n, rdwind, rdsal, rdtemp, rdcair, rdcwater)	
	for(i in c(1:repeatlhc_n-1)){
		y<-rbind(y,randomsamplelhc(compound, dist_n, lhc_n, rdwind, rdsal, rdtemp, rdcair, rdcwater))
		print(i)
	}
	y
}


calculate_pdf<-function(filename,compound,dist_n=10000,lhc_n=1000,repeatlhc_n=100,distlist=defaultdistlist){
	#load data in 
	dataset<-read.csv(filename, header=TRUE)
	attach(dataset)
	

	#fit distributions to each data type and randomly sample them dist_n times 
	winds<-na.omit(WindSpeed)
	attributes(winds)<-NULL
	dwind<-fitdistr(winds,distlist$WindSpeed)
	if(distlist$WindSpeed=="normal")
		{
		rdwind<-getndist_within_bounds(mean=dwind[[1]][[1]],sd=dwind[[1]][[2]],n=dist_n)
		}
	else
		{
		rdwind<-getlndist_within_bounds(meanlog=dwind[[1]][[1]],sdlog=dwind[[1]][[2]],n=dist_n)
		}

	sal<-na.omit(Salinity)
	attributes(sal)<-NULL
	dsal<-fitdistr(sal,distlist$Salinity)
	if(distlist$Salinity=="normal")
			{
			rdsal<-getndist_within_bounds(mean=dsal[[1]][[1]],sd=dsal[[1]][[2]],n=dist_n)
			}
	else
			{			
			rdsal<-getlndist_within_bounds(meanlog=dsal[[1]][[1]],sdlog=dsal[[1]][[2]],n=dist_n)
			}

	temp<-na.omit(Temperature)
	attributes(temp)<-NULL
	dtemp<-fitdistr(temp,distlist$Temperature)
	if(distlist$Temperature=="normal")	
		{
		rdtemp<-getndist_within_bounds(mean=dtemp[[1]][[1]],sd=dtemp[[1]][[2]],n=dist_n)
		}
	else
		{
		rdtemp<-getlndist_within_bounds(meanlog=dtemp[[1]][[1]],sdlog=dtemp[[1]][[2]],n=dist_n)
		}


	cair<-na.omit(C_air)
	attributes(cair)<-NULL
	dcair<-fitdistr(cair,distlist$C_air)	
	if(distlist$C_air=="normal")
		{
		rdcair<-getndist_within_bounds(mean=dcair[[1]][[1]],sd=dcair[[1]][[2]],n=dist_n)
		}
	else
		{
		rdcair<-getlndist_within_bounds(meanlog=dcair[[1]][[1]],sdlog=dcair[[1]][[2]],n=dist_n)
		}

	cwater<-na.omit(C_water)
	attributes(cwater)<-NULL
	dcwater<-fitdistr(cwater,distlist$C_water)
	if(distlist$C_water=="normal")
		{
		rdcwater<-getndist_within_bounds(mean=dcwater[[1]][[1]],sd=dcwater[[1]][[2]],n=dist_n)
		}
	else
		{
		rdcwater<-getlndist_within_bounds(meanlog=dcwater[[1]][[1]],sdlog=dcwater[[1]][[2]],n=dist_n)
		}


	
	pdfdata<-repeatlhc(compound, dist_n, lhc_n, repeatlhc_n, rdwind, rdsal, rdtemp, rdcair, rdcwater)

	pdfdata
}
