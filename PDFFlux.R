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


getndist_within_bounds<-function(mean,sd,conf=1.0,n=1000){
	#function to remove values above given confidence limit when generating random samples from normal distribution
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



