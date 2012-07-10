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


randomsamplelhc<-function(compound){
	#randomly sample the distributions, then create a latin hypercube sample
	rdwind<-rlnorm(n=100000,meanlog=dwind[[1]][[1]],sdlog=dwind[[1]][[2]])
	rdsal<-rnorm(100000,mean=dsal[[1]][[1]],sd=dsal[[1]][[2]])
	rdtemp<-rnorm(100000,mean=dtemp[[1]][[1]],sd=dtemp[[1]][[2]])
	ifelse(compound=="CHBr3",gasd<-dgas,gasd<-dgas2)
	rdgas<-rlnorm(n=100000,meanlog=gasd[[1]][[1]],sdlog=gasd[[1]][[2]])
	ifelse(compound=="CHBr3",dsw_high<-dCHBr3_high,dsw_high<-dCH2Br2_high)
	ifelse(compound=="CHBr3",dsw_low<-dCHBr3_low,dsw_low<-dCH2Br2_low)	
	ifelse(highlow=="high",
		rdsw<-rlnorm(n=100000,meanlog=dsw_high[[1]][[1]],sdlog=dsw_high[[1]][[2]]),
		rdsw<-rlnorm(n=100000,meanlog=dsw_low[[1]][[1]],sdlog=dsw_low[[1]][[2]]))	
	lh<-round(999*randomLHS(1000,k=5))+1
	x<-data.frame(rdwind[lh[,1]])
	colnames(x)<-c("u")
	x$S<-rdsal[lh[,2]]
	x$T<-rdtemp[lh[,3]]
	x$gconc<-rdgas[lh[,4]]
	x$swconc<-rdsw[lh[,5]]
	# Henry's law
	x$KH<-KH(compound,x$T,x$S)
	# delta C
	x$dC<-x$swconc-(x$gconc/x$KH)
	#transfer velocity	
	x$Kw<-Kw(compound,x$T,x$u,x$S)
	#flux
	x$F<-x$Kw*x$dC
	x
	
}	



calculate_pdf<-function(filename,compound,dist_n=10000,lhc_n=1000,repeatlhc_n=100,distlist=defaultdistlist){
	#load data in 
	dataset<-read.csv(filename, sep="\t", header=TRUE)
	attach(dataset)
	attributes
	

}
