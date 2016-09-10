#Define Server for PipPrototypeApp2
#Remove all current objects in the workspace
rm(list=ls(all=TRUE))

#Load package libraries
library(shiny)
library(shinydashboard)	#Fancy GUI
library(ggplot2)  #Plotting
library(grid)  #Plotting
library(plyr)  #ddply function
library(compiler)  #Compile repeatedly called functions
library(rmarkdown)  #Generate report (Word document)

#Directories on Windows
#dir <- "//psf/Home/Desktop/PipPrototypeApp2/"	#Directory where application files are saved
#pandocdir <- "C:/Program Files/RStudio/bin/pandoc"	#Directory for pancdoc (writing to word document)

#Directories on Mac
dir <- "/Volumes/Prosecutor/PipPrototypeApp/"	#Directory where application files are saved
pandocdir <- "/Applications/RStudio.app/Contents/MacOS/pandoc"	#Directory for pancdoc (writing to word document)

#-------------------------------------------------------------------------------------
#Read in previous saved data
saved.data <- read.csv(paste(dir,"saved_data.csv",sep=""), stringsAsFactors=F, na.strings=".")

#Elements that only need to run once on application start-up
#Use custom ggplot2 theme
theme_bw2 <- theme_set(theme_bw(base_size = 16))
theme_bw2 <- theme_update(plot.margin = unit(c(1,1,1,1), "lines"),
							axis.title.x=element_text(size = 16, vjust = 0),
							axis.title.y=element_text(size = 16, vjust = 0, angle = 90),
							strip.text.x=element_text(size = 14),
							strip.text.y=element_text(size = 14, angle = 90),
							legend.direction = "vertical",
              legend.position = "bottom",
              legend.box = "horizontal")

#Describe population model
#1-compartment IV infusion
#Define parameter values
#Thetas
POPCL <- 9.14  #Clearance, L/hour
POPV <- 12.2  #Volume, L
COVCRCL <- 4.6
COVWT <- 9.49
#Omegas (as SD)
PPVCL <- 31.1/100
PPVV <- 38.0/100
#Epsilons (as SD)
ERRPRO <- 9.33/100  #Proportional residual error (log-normally distributed)

#Set up a data frame of individuals in order to calculate and plot concentration confidence intervals
n <- 50	#This is hard coded - no reason for a user to change this number
ID <- seq(from = 1, to = n, by = 1)

#Time interval
TINT <- 0.5

#Calculate concentrations at each time-point for the individual
#Function for calculating concentrations in a loop
conc.function <- function(df){
  #Set initial values in the central compartment
	for(i in 2:nrow(df)) {
		df$A1[df$TIME == 0] <- 0  #Drug amount in the central compartment at time zero
		df$IPRE[df$TIME == 0] <- 0  #Drug concentration in the central compartment at time zero
		df$CLgrp <- df$CL+COVCRCL^(df$CRCL/68.7) #Effect of CrCL on CL
		df$Vgrp <- df$V+COVWT^(df$WT/61.1)  #Effect of WT on V
		KE <- df$CLgrp[i]/df$Vgrp[i]  #Elimination rate-constant
		time <- df$TIME[i]-df$TIME[i-1]
		A1.prev <- df$A1[i-1]
		dose.rate <- df$RATE[i]
		df$A1[i] <- dose.rate/KE*(1-exp(-time*KE))+A1.prev*exp(-time*KE) #Amount
		df$IPRE[i] <- df$A1[i]/df$Vgrp[i] #Concentration
	}
	df
}
conc.function.cmpf <- cmpfun(conc.function)

#Use the patient's posterior information and the prior population model to estimate values of CL and V for the individual
bayes.function <- function(OBS.data) {
	#Bayes Objective Function Value
	bayes.ofv <- function(par) {
		CLfit <- par[1]
		Vfit <- par[2]
		OBS.data$CL <- CLfit
		OBS.data$V <- Vfit
		iCOVCRCL <- COVCRCL^(OBS.data$CRCL[1]/68.7)
		iCOVWT <- COVWT^(OBS.data$WT[1]/61.1)
			#Population model information for "bayes.function"
			CLprior <- POPCL+iCOVCRCL
			CLbsvprior <- PPVCL
			Vprior <- POPV+iCOVWT
			Vbsvprior <- PPVV
			Cprior <- 0  #Correlation between CL and V
			sigma <- ERRPRO
		#Calculate concentrations using the "fit" parameters
		CONC.data <- conc.function.cmpf(OBS.data)
		#If DV was NA, IPRE needs to be NA too (for calculating log-likelihood)
		CONC.data$IPRE[is.na(CONC.data$DV) == T] <- NA
		DV2 <- CONC.data$DV[is.na(CONC.data$DV)==F]  #DV as a vector
		IPRE2 <- CONC.data$IPRE[is.na(CONC.data$IPRE)==F] #IPRE as a vector
		#Posterior component (from the data)
		#Log densities of residuals
		#Residual error model, Y = IPRE*exp(ERR), log(Y) = log(IPRE) + ERR
		loglikpost <- dnorm(log(na.omit(DV2)), mean=log(na.omit(IPRE2)), sd=sigma, log=T)
		#Prior component (from the model)
		ETA <- log(c(CLfit+iCOVCRCL,Vfit+iCOVWT)/c(CLprior,Vprior))
		ETABSV <- c(CLbsvprior,Vbsvprior)
		loglikprior <- dnorm(ETA, mean=0, sd=ETABSV, log=T)
		#Calculate the combined likelihood
		OFVbayes <- -1*sum(loglikpost,loglikprior)
		OFVbayes  #Optimise the Bayes Objective Function Value
	}
	#Parameter estimates
	par <- c(2,2)
	#Optimise
	BAYES.data <- optim(par, bayes.ofv, method="L-BFGS-B", lower=c(0.001,0.001), upper=c(Inf,Inf), control = list(parscale=par, factr=1e3), hessian=TRUE)
	BAYES.data
}
bayes.function.cmpf <- cmpfun(bayes.function)

#Calculate concentrations at each time-point for the individual (for the next dose)
#Function for calculating concentrations in a loop
sim.conc.function <- function(SIM.data,TIME.prev,AMT.prev,CONC.prev){
	for(i in 2:nrow(SIM.data)) {
		SIM.data$A1[SIM.data$TIME == TIME.prev] <- AMT.prev  #Drug amount in the central compartment at time of next dose
		SIM.data$IPRE[SIM.data$TIME == TIME.prev] <- CONC.prev  #Concentration in the central compartment at time of next dose
		SIM.data$CLgrp <- SIM.data$CL+COVCRCL^(SIM.data$CRCL/68.7) #Effect of CrCL on CL
		SIM.data$Vgrp <- SIM.data$V+COVWT^(SIM.data$WT/61.1)  #Effect of WT on V
		KE <- SIM.data$CLgrp[i]/SIM.data$Vgrp[i]  #Elimination rate-constant
		time <- SIM.data$TIME[i]-SIM.data$TIME[i-1]
		A1.prev <- SIM.data$A1[i-1]
		dose.rate <- SIM.data$RATE[i]
		SIM.data$A1[i] <- dose.rate/KE*(1-exp(-time*KE))+A1.prev*exp(-time*KE) #Amount
		SIM.data$IPRE[i] <- SIM.data$A1[i]/SIM.data$Vgrp[i] #Concentration
	}
	SIM.data
}
sim.conc.function.cmpf <- cmpfun(sim.conc.function)

#Extract results from rBAYES.data and rSIM.data1, and make a new data frame of concentrations calculated using the new fitted parameters to create confidence intervals!
ci.function <- function(BAYES.data,OBS.data,SIM.data) {
	#Last time-point in OBS.data
	TIME.prev <- tail(OBS.data$TIME,1)
	#Repeat OBS.data by the number of pre-specified ID's to calculate confidence intervals
	REP.OBS.data <- lapply(OBS.data, rep.int, times=length(ID))
	REP.OBS.data <- as.data.frame(REP.OBS.data)
	#Remove existing CL and V columns
	REP.OBS.data <- REP.OBS.data[c(-5,-6)]
	ID.OBS <- rep(ID, times=length(OBS.data$TIME))
	ID.OBS <- sort(ID.OBS)
	REP.OBS.data$ID <- ID.OBS
		#Repeat SIM.data by the number of pre-specified ID's to calculate confidence intervals
		REP.SIM.data <- lapply(SIM.data, rep.int, times=length(ID))
		REP.SIM.data <- as.data.frame(REP.SIM.data)
		#Remove existing CL, V, A1, CLgrp and Vgrp columns
		ID.SIM <- rep(ID, times=length(SIM.data$TIME))
		ID.SIM <- sort(ID.SIM)
		REP.SIM.data$ID <- ID.SIM
		REP.SIM.data <- REP.SIM.data[c(-5,-6,-11,-12,-13)]
		#Remove the first row from each individual (should be in REP.OBS.data as the last row)
		REP.SIM.data <- REP.SIM.data[REP.SIM.data$TIME != TIME.prev,]
			#rbind REP.OBS.data and REP.SIM.data
			REP.data <- rbind(REP.OBS.data,REP.SIM.data)
			REP.data <- REP.data[with(REP.data, order(REP.data$ID,REP.data$TIME)), ]
	#Confidence intervals
	#Calculate standard errors for individual parameter estimates
	#Assign the R matrix to the Hessian returned by optim
	Rmatrix <- BAYES.data$hessian
	#Calculate the variance-covariance matrix
	VCmatrix <- solve(Rmatrix)
	#Calculate the parameter standard errors
	se.par <- sqrt(diag(VCmatrix))
	#Calculate parameter standard errors as a percentage (CL, V)
	ERR.CL <- rnorm(n, mean=BAYES.data$par[1], sd=se.par[1])
	ERR.V <- rnorm(n, mean=BAYES.data$par[2], sd=se.par[2])
	#Make a data frame
	IND.data <- data.frame(ID=ID, CL=ERR.CL, V=ERR.V)
	MULTI.data <- merge(IND.data,REP.data,by=c("ID"),all=T)
	#Calculate concentrations at each time-point for each individual
	CI.data <- conc.function.cmpf(MULTI.data)
}
ci.function.cmpf <- cmpfun(ci.function)

#Function for calculating time above MIC
time.function <- function(data) {
	time.EMIC <- length(data$MIC_FLAG[data$MIC_FLAG == 1]) #Count how many time-points are above the MIC
	percent.EMIC <- round(time.EMIC/length(data$TIME)*100, digits = 2)	#Percentage of time above the MIC
}
time.function.cmpf <- cmpfun(time.function)

#Functions for calculating confidence intervals
CI95lo <- function(x) quantile(x, probs = 0.025)
CI95hi <- function(x) quantile(x, probs = 0.975)

#Create a function that calculates the maximum concentration for the scale
max.function <- function(FIT.data,CI.data,SIM.data) max(c(FIT.data$IPRE,CI.data$IPRE,SIM.data$IPRE))

#Returns a single value for each ID of dataframe
oneperID <- function(x) tail(x, n = 1)

#-------------------------------------------------------------------------------------
#Define the "server" part of the Shiny application
shinyServer(function(input,output,session) {

###########
##_INPUT_##
###########
#Make a reactive expression that calls in the prior concentrations and their corresponding time-points
rOBS.data <- reactive({
	#Define TIME sequence dependent on the dosing frequency and the number of previous doses
		#Dosing frequency
		if (input$PFREQ == 1) {PFREQ <- 4}
		if (input$PFREQ == 2) {PFREQ <- 6}
		if (input$PFREQ == 3) {PFREQ <- 8}
		#Number of previous doses
		if (input$NPDOSE == 1) {NPDOSE <- 0}
		if (input$NPDOSE == 2) {NPDOSE <- 1}
		if (input$NPDOSE == 3) {NPDOSE <- 2}
		if (input$NPDOSE == 4) {NPDOSE <- 3}
		#TIME
		TIME <- seq(from=0, to = NPDOSE*PFREQ+PFREQ, by=TINT) #Time elapsed by previous doses + Time for sampled dose
	#Set up a data frame of individuals in order to calculate and plot concentration confidence intervals
	ID.1 <- rep(1, times=length(TIME))
	#Call in patient specific data
	AGE <- round(as.double(difftime(input$DDATE,input$BDATE))/365.25, digits = 1) #Patient's age (years)
	WT <- input$WT  #Patient's ideal body weight (kg)
	SEX <- input$SEX  #Patient's gender (1 = Male, 2 = Female)
	SECR <- input$SECR  #Patient's serum creatinine (µmol/L)
		#Calculate creatinine clearance for the individual
		if (input$SEX == 1) {CRCL <- ((140-AGE)*WT*1.23)/SECR}
		if (input$SEX == 2) {CRCL <- ((140-AGE)*WT*1.04)/SECR}
		#Calculate individual's effect of CRCL and WT to CL and V, respectively
		iCOVCRCL <- COVCRCL^(CRCL/68.7)  #Effect of CrCL on CL
		iCOVWT <- COVWT^(WT/61.1)  #Effect of WT on V
	#Call in previous concentrations and their corresponding times
	DV <- rep(NA, times = length(TIME))
	PCONC1 <- input$PCONC1	#There is always a PCONC1!
	PTIME1 <- input$PTIME1
	DV[TIME == NPDOSE*PFREQ+PTIME1] <- PCONC1
	PCONC2 <- NA
	PTIME2 <- NA
	PCONC3 <- NA
	PTIME3 <- NA
	PCONC4 <- NA
	PTIME4 <- NA
	if (input$NCONC >= 2) {
		PCONC2 <- input$PCONC2
		PTIME2 <- input$PTIME2
		DV[TIME == NPDOSE*PFREQ+PTIME2] <- PCONC2
	}
	if (input$NCONC >= 3) {
		PCONC3 <- input$PCONC3
		PTIME3 <- input$PTIME3
		DV[TIME == NPDOSE*PFREQ+PTIME3] <- PCONC3
	}
	if (input$NCONC == 4) {
		PCONC4 <- input$PCONC4
		PTIME4 <- input$PTIME4
		DV[TIME == NPDOSE*PFREQ+PTIME4] <- PCONC4
	}
  #Setting up the infusion and dosing event information
  PDOSE <- input$PDOSE	#Dose (Amount)
  PINFD <- input$PINFD	#Infusion Duration (hours)
	RATET <- PDOSE/PINFD	#"Temporary Rate"; Dose (Amount)/Infusion Duration (hours)
	INFEND <- PINFD+TINT	#Time interval when infusion ended
	END <- max(TIME)+PFREQ	#A time-point long after the specified dosing time
		#Vector marking infusion's time events (dependent on the number of previous doses)
		#PFREQ+TINT; time interval when next infusion is started
		#PFREQ+INFEND; time when infusion stopped
		#Vector marking infusion's rates (dependent on the number of previous doses)
		if (NPDOSE == 0) {	#1st dose is sampled dose
		TIMEinf <- c(0,INFEND,END)
		RATEinf <- c(RATET,0,0)
		}
		if (NPDOSE == 1) {	#2nd dose is sampled dose
		TIMEinf <- c(0,INFEND,PFREQ+TINT,PFREQ+INFEND,END)
		RATEinf <- c(RATET,0,RATET,0,0)
		}
		if (NPDOSE == 2) {	#3rd dose is sampled dose
		TIMEinf <- c(0,INFEND,PFREQ+TINT,PFREQ+INFEND,2*PFREQ+TINT,2*PFREQ+INFEND,END)
		RATEinf <- c(RATET,0,RATET,0,RATET,0,0)
		}
		if (NPDOSE == 3) {	#4th dose is sampled dose
		TIMEinf <- c(0,INFEND,PFREQ+TINT,PFREQ+INFEND,2*PFREQ+TINT,2*PFREQ+INFEND,3*PFREQ+TINT,3*PFREQ+INFEND,END)
		RATEinf <- c(RATET,0,RATET,0,RATET,0,RATET,0,0)
		}
		#Define an interpolation function that returns rate when given time - "const"
		step.doseinf <- approxfun(TIMEinf, RATEinf, method = "const")
		RATE <- step.doseinf(TIME)
	#Assign levels to the different dosing intervals
	NDOSE <- seq(from=1, to=NPDOSE+1, by=1) #Number of previous dosing intervals + sampled dosing interval
	if (NPDOSE == 0) {NDOSE <- rep(NDOSE, times=length(TIME)/length(NDOSE))}
	if (NPDOSE != 0) {NDOSE <- c(1,rep(NDOSE, times=length(TIME)/length(NDOSE))) }
	NDOSE <- sort(NDOSE)
	#Make a data frame of previous concentrations and patient's characteristics
	IPRE <- 0
	OBS.data <- data.frame("ID"=ID.1,TIME,NDOSE,RATE,"CL"=iCOVCRCL,"V"=iCOVWT,IPRE,DV,"CRCL"=CRCL,"WT"=WT)
	as.data.frame(OBS.data)
})

#Make a reactive expression that calls in the prior concentrations and their corresponding time-points and parts of the population model to estimate new parameters for the individual
rBAYES.data <- reactive({
	BAYES.data <- bayes.function.cmpf(rOBS.data())
})

#Extract results from rBAYES.data and make a new data frame of concentrations calculated using the new fitted parameters
rFIT.data <- reactive({
  BAYES.data <- rBAYES.data()
  OBS.data <- rOBS.data()

  OBS.data$CL <- BAYES.data$par[1]
  OBS.data$V <- BAYES.data$par[2]
  #Calculate concentrations at each time-point
  FIT.data <- conc.function.cmpf(OBS.data)
  as.data.frame(FIT.data)

	# fit.function <- function(BAYES.data,OBS.data) {
	# 	OBS.data$CL <- BAYES.data$par[1]
	# 	OBS.data$V <- BAYES.data$par[2]
	# 	#Calculate concentrations at each time-point
	# 	FIT.data <- conc.function.cmpf(OBS.data)
	# }
	# FIT.data <- fit.function(rBAYES.data(),rOBS.data())
	# as.data.frame(FIT.data)
})

################
##_rSIM.data1_##
################
#Extract results from rBAYES.data and simulate what the next dose would be
#Requires rBAYES.data for estimates of CL and V
#Requires rFIT.data for the last concentrations and amounts
#Requires rOBS.data for the population CL and V values
rSIM.data1 <- reactive({
	sim.function <- function(BAYES.data,FIT.data,OBS.data) {
		#Setting up the new infusion and new dosing event information
		SDOSE1 <- input$SDOSE1	#Dose (Amount)
		SINFD1 <- input$SINFD1	#Infusion Duration (hours)
		RATES <- SDOSE1/SINFD1	#"Temporary Simulated Rate"; Dose (Amount)/Infusion Duration (hours)
		INFEND <- SINFD1+TINT	#Time interval when infusion ended
		#New TIME sequence
		TIME.prev <- tail(FIT.data$TIME,1) #Last time point in FIT.data from which simulated dose starts from
		if (input$SFREQ1 == 1) {SFREQ1 <- 4}
		if (input$SFREQ1 == 2) {SFREQ1 <- 6}
		if (input$SFREQ1 == 3) {SFREQ1 <- 8}
		TIMES <- seq(from=TIME.prev, to=TIME.prev+3*SFREQ1, by=TINT)
		END <- max(TIMES)+SFREQ1	#A time-point long after the specified dosing time
			#Vector marking infusion's time events (dependent on the number of previous doses)
			#SFREQ+TINT; time interval when next infusion is started
			#SFREQ+INFEND; time when infusion stopped
			#Vector marking infusion's rates (dependent on the number of previous doses)
			TIMEinf <- TIME.prev + c(0,INFEND,SFREQ1+TINT,SFREQ1+INFEND,2*SFREQ1+TINT,2*SFREQ1+INFEND,END)
			RATEinf <- c(RATES,0,RATES,0,RATES,0,0)
			#Define an interpolation function that returns rate when given time - "const"
			step.doseinf <- approxfun(TIMEinf, RATEinf, method="const")
			RATE <- step.doseinf(TIMES)
		#Assign levels to the different dosing intervals
		NDOSE.prev <- tail(OBS.data$NDOSE,1)  #Last NDOSE value in OBS.data from which simulated doses starts from
		SNDOSE <- seq(from=NDOSE.prev+1, to=NDOSE.prev+3, by=1) #Number of simulated dosing intervals (3)
		SNDOSE <- c(NDOSE.prev+1,rep(SNDOSE, times=length(TIMES)/length(SNDOSE)))
		SNDOSE <- sort(SNDOSE)
		#Patient parameters
		CL <- BAYES.data$par[1]
		V <- BAYES.data$par[2]
		#Call in the drug amount from the previous dose
		AMT.prev <- FIT.data$A1[FIT.data$TIME == TIME.prev]
		CONC.prev <- FIT.data$IPRE[FIT.data$TIME == TIME.prev]
		#Make a data frame of previous parameters and patient's characteristics
		SIM.data1 <- data.frame(ID=1,TIME=TIMES,NDOSE=SNDOSE,RATE,CL,V,IPRE=NA,DV=NA,CRCL=OBS.data$CRCL[1],WT=OBS.data$WT[1])
		SIM.data1 <- sim.conc.function.cmpf(SIM.data1,TIME.prev,AMT.prev,CONC.prev)
	}
	SIM.data1 <- sim.function(rBAYES.data(),rFIT.data(),rOBS.data())
	as.data.frame(SIM.data1)
})

################
##_rSIM.data2_##
################
#Extract results from rBAYES.data and simulate what the next dose would be
#Requires rBAYES.data for estimates of CL and V
#Requires rFIT.data for the last concentrations and amounts
#Requires rOBS.data for the population CL and V values
rSIM.data2 <- reactive({
	sim.function <- function(BAYES.data,FIT.data,OBS.data) {
		#Setting up the new infusion and new dosing event information
		SDOSE2 <- input$SDOSE2	#Dose (Amount)
		SINFD2 <- input$SINFD2	#Infusion Duration (hours)
		RATES <- SDOSE2/SINFD2	#"Temporary Simulated Rate"; Dose (Amount)/Infusion Duration (hours)
		INFEND <- SINFD2+TINT	#Time interval when infusion ended
		#New TIME sequence
		TIME.prev <- tail(FIT.data$TIME,1) #Last time point in OBS.data from which simulated dose starts from
		if (input$SFREQ2 == 1) {SFREQ2 <- 4}
		if (input$SFREQ2 == 2) {SFREQ2 <- 6}
		if (input$SFREQ2 == 3) {SFREQ2 <- 8}
		TIMES <- seq(from=TIME.prev, to=TIME.prev+3*SFREQ2, by=TINT)
		END <- max(TIMES)+SFREQ2	#A time-point long after the specified dosing time
			#Vector marking infusion's time events (dependent on the number of previous doses)
			#SFREQ+TINT; time interval when next infusion is started
			#SFREQ+INFEND; time when infusion stopped
			#Vector marking infusion's rates (dependent on the number of previous doses)
			TIMEinf <- TIME.prev + c(0,INFEND,SFREQ2+TINT,SFREQ2+INFEND,2*SFREQ2+TINT,2*SFREQ2+INFEND,END)
			RATEinf <- c(RATES,0,RATES,0,RATES,0,0)
			#Define an interpolation function that returns rate when given time - "const"
			step.doseinf <- approxfun(TIMEinf, RATEinf, method="const")
			RATE <- step.doseinf(TIMES)
		#Assign levels to the different dosing intervals
		NDOSE.prev <- tail(OBS.data$NDOSE,1)  #Last NDOSE value in OBS.data from which simulated doses starts from
		SNDOSE <- seq(from=NDOSE.prev+1, to=NDOSE.prev+3, by=1) #Number of simulated dosing intervals (3)
		SNDOSE <- c(NDOSE.prev+1,rep(SNDOSE, times=length(TIMES)/length(SNDOSE)))
		SNDOSE <- sort(SNDOSE)
		#Patient parameters
		CL <- BAYES.data$par[1]
		V <- BAYES.data$par[2]
		#Call in the drug amount from the previous dose
		AMT.prev <- FIT.data$A1[FIT.data$TIME == TIME.prev]
		CONC.prev <- FIT.data$IPRE[FIT.data$TIME == TIME.prev]
		#Make a data frame of previous parameters and patient's characteristics
		SIM.data2 <- data.frame(ID=1,TIME=TIMES,NDOSE=SNDOSE,RATE,CL,V,IPRE=NA,DV=NA,CRCL=OBS.data$CRCL[1],WT=OBS.data$WT[1])
		SIM.data2 <- sim.conc.function.cmpf(SIM.data2,TIME.prev,AMT.prev,CONC.prev)
	}
	SIM.data2 <- sim.function(rBAYES.data(),rFIT.data(),rOBS.data())
	as.data.frame(SIM.data2)
})

###############
##_rCI.data1_##
###############
#Extract results from rBAYES.data and rSIM.data1, and make a new data frame of concentrations calculated using the new fitted parameters to create confidence intervals!
rCI.data1 <- reactive({
	CI.data1 <- ci.function(rBAYES.data(),rOBS.data(),rSIM.data1())
	as.data.frame(CI.data1)
})

###############
##_rCI.data2_##
###############
#Extract results from rBAYES.data and rSIM.data2, and make a new data frame of concentrations calculated using the new fitted parameters to create confidence intervals!
rCI.data2 <- reactive({
	CI.data2 <- ci.function(rBAYES.data(),rOBS.data(),rSIM.data2())
	as.data.frame(CI.data2)
})

################
##_rMIC.data1_##
################
#Calculate the time above MIC for the sampled dose
#Uses rFIT.data for calculating time above MIC for fitted concentrations
#Uses rCI.data1 for calculating the confidence interval for time above MIC
#Uses rSIM.data1 for calculating the time above MIC for simulated concentrations
rMIC.data1 <- reactive({
	mic.function1 <- function(FIT.data,CI.data1,SIM.data1) {
		#Assign values to the MIC input
		if (input$MIC == 1) {MIC <- 0.25}
		if (input$MIC == 2) {MIC <- 0.5}
		if (input$MIC == 3) {MIC <- 1}
		if (input$MIC == 4) {MIC <- 2}
		if (input$MIC == 5) {MIC <- 4}
		if (input$MIC == 6) {MIC <- 8}
		if (input$MIC == 7) {MIC <- 16}
		if (input$MIC == 8) {MIC <- 32}
		if (input$MIC == 9) {MIC <- 64}
		#Identify the sampled dosing interval, then subset FIT.data for this dosing interval
		#Number of previous doses
		if (input$NPDOSE == 1) {NPDOSE <- 0}
		if (input$NPDOSE == 2) {NPDOSE <- 1}
		if (input$NPDOSE == 3) {NPDOSE <- 2}
		if (input$NPDOSE == 4) {NPDOSE <- 3}
		#Time above MIC for sampled data
		SAMP_DOSE <- NPDOSE + 1 #Number of previous doses + 1 = sampled dose
		SAMP.data <- subset(FIT.data, NDOSE == SAMP_DOSE)
		SAMP.data$MIC_FLAG <- 0  #Below the MIC
		SAMP.data$MIC_FLAG[SAMP.data$IPRE > MIC] <- 1	#Above the MIC
		percent.EMIC.SAMP <- ddply(SAMP.data, .(), time.function.cmpf)
		#Confidence intervals for time above MIC for sampled data
		SAMP.CI.data <- subset(CI.data1, NDOSE == SAMP_DOSE)
		SAMP.CI.data$MIC_FLAG <- 0  #Below the MIC
		SAMP.CI.data$MIC_FLAG[SAMP.CI.data$IPRE > MIC] <- 1	#Above the MIC
		percent.EMIC.SAMP.CI <- ddply(SAMP.CI.data, .(ID), time.function.cmpf)
		percent.EMIC.SAMP.CIlo <- signif(CI95lo(percent.EMIC.SAMP.CI$V1), digits = 3)
		percent.EMIC.SAMP.CIhi <- signif(CI95hi(percent.EMIC.SAMP.CI$V1), digits = 3)
		#Time above MIC for the next simulated dose
		SIM_DOSE <- NPDOSE + 2 #Number of previous doses + 1 sampled dose + 1 = first simulated dose
		SIM.data <- subset(SIM.data1, NDOSE == SIM_DOSE)
		SIM.data <- SIM.data[-1,]
		SIM.data$MIC_FLAG <- 0  #Below the MIC
		SIM.data$MIC_FLAG[SIM.data$IPRE > MIC] <- 1	#Above the MIC
		percent.SMIC.SIM <- ddply(SIM.data, .(), time.function.cmpf)
		#Confidence intervals for time above MIC for the next simulated dose
		SIM.CI.data <- subset(CI.data1, NDOSE > SAMP_DOSE)
		SIM.CI.data$MIC_FLAG <- 0  #Below the MIC
		SIM.CI.data$MIC_FLAG[SIM.CI.data$IPRE > MIC] <- 1	#Above the MIC
		percent.SMIC.SIM.CI <- ddply(SIM.CI.data, .(ID), time.function.cmpf)
		percent.SMIC.SIM.CIlo <- CI95lo(percent.SMIC.SIM.CI$V1)
		percent.SMIC.SIM.CIlo <- signif(percent.SMIC.SIM.CIlo, digits = 3)
		percent.SMIC.SIM.CIhi <- CI95hi(percent.SMIC.SIM.CI$V1)
		percent.SMIC.SIM.CIhi <- signif(percent.SMIC.SIM.CIhi, digits = 3)
		#Make a data frame of results (MIC.data)
		EMIC <- paste(signif(percent.EMIC.SAMP$V1, digits = 3)," (",percent.EMIC.SAMP.CIlo," - ",percent.EMIC.SAMP.CIhi,")",sep="")
		SMIC <- paste(signif(percent.SMIC.SIM$V1, digits = 3)," (",percent.SMIC.SIM.CIlo," - ",percent.SMIC.SIM.CIhi,")",sep="")
		MIC.data <- data.frame(Value = c(EMIC,SMIC),
													row.names = c("Bayes Estimated (i.e., Sampled Dose)","Predicted (i.e., Next 3 Doses)"))
	}
	MIC.data1 <- mic.function1(rFIT.data(),rCI.data1(),rSIM.data1())
})

################
##_rMIC.data2_##
################
#Calculate the time above MIC for the sampled dose
#Uses rFIT.data for calculating time above MIC for fitted concentrations
#Uses rCI.data2 for calculating the confidence interval for time above MIC
#Uses rSIM.data2 for calculating the time above MIC for simulated concentrations
rMIC.data2 <- reactive({
	mic.function2 <- function(FIT.data,CI.data2,SIM.data2) {
		#Assign values to the MIC input
		if (input$MIC == 1) {MIC <- 0.25}
		if (input$MIC == 2) {MIC <- 0.5}
		if (input$MIC == 3) {MIC <- 1}
		if (input$MIC == 4) {MIC <- 2}
		if (input$MIC == 5) {MIC <- 4}
		if (input$MIC == 6) {MIC <- 8}
		if (input$MIC == 7) {MIC <- 16}
		if (input$MIC == 8) {MIC <- 32}
		if (input$MIC == 9) {MIC <- 64}
		#Identify the sampled dosing interval, then subset FIT.data for this dosing interval
		#Number of previous doses
		if (input$NPDOSE == 1) {NPDOSE <- 0}
		if (input$NPDOSE == 2) {NPDOSE <- 1}
		if (input$NPDOSE == 3) {NPDOSE <- 2}
		if (input$NPDOSE == 4) {NPDOSE <- 3}
		#Time above MIC for sampled data
		SAMP_DOSE <- NPDOSE + 1 #Number of previous doses + 1 = sampled dose
		SAMP.data <- subset(FIT.data, NDOSE == SAMP_DOSE)
		SAMP.data$MIC_FLAG <- 0  #Below the MIC
		SAMP.data$MIC_FLAG[SAMP.data$IPRE > MIC] <- 1	#Above the MIC
		percent.EMIC.SAMP <- ddply(SAMP.data, .(), time.function.cmpf)
		#Confidence intervals for time above MIC for sampled data
		SAMP.CI.data <- subset(CI.data2, NDOSE == SAMP_DOSE)
		SAMP.CI.data$MIC_FLAG <- 0  #Below the MIC
		SAMP.CI.data$MIC_FLAG[SAMP.CI.data$IPRE > MIC] <- 1	#Above the MIC
		percent.EMIC.SAMP.CI <- ddply(SAMP.CI.data, .(ID), time.function.cmpf)
		percent.EMIC.SAMP.CIlo <- signif(CI95lo(percent.EMIC.SAMP.CI$V1), digits = 3)
		percent.EMIC.SAMP.CIhi <- signif(CI95hi(percent.EMIC.SAMP.CI$V1), digits = 3)
		#Time above MIC for the next simulated dose
		SIM_DOSE <- NPDOSE + 2 #Number of previous doses + 1 sampled dose + 1 = first simulated dose
		SIM.data <- subset(SIM.data2, NDOSE == SIM_DOSE)
		SIM.data <- SIM.data[-1,]
		SIM.data$MIC_FLAG <- 0  #Below the MIC
		SIM.data$MIC_FLAG[SIM.data$IPRE > MIC] <- 1	#Above the MIC
		percent.SMIC.SIM <- ddply(SIM.data, .(), time.function.cmpf)
		#Confidence intervals for time above MIC for the next simulated dose
		SIM.CI.data <- subset(CI.data2, NDOSE > SAMP_DOSE)
		SIM.CI.data$MIC_FLAG <- 0  #Below the MIC
		SIM.CI.data$MIC_FLAG[SIM.CI.data$IPRE > MIC] <- 1	#Above the MIC
		percent.SMIC.SIM.CI <- ddply(SIM.CI.data, .(ID), time.function.cmpf)
		percent.SMIC.SIM.CIlo <- CI95lo(percent.SMIC.SIM.CI$V1)
		percent.SMIC.SIM.CIlo <- signif(percent.SMIC.SIM.CIlo, digits = 3)
		percent.SMIC.SIM.CIhi <- CI95hi(percent.SMIC.SIM.CI$V1)
		percent.SMIC.SIM.CIhi <- signif(percent.SMIC.SIM.CIhi, digits = 3)
		#Make a data frame of results (MIC.data)
		EMIC <- paste(signif(percent.EMIC.SAMP$V1, digits = 3)," (",percent.EMIC.SAMP.CIlo," - ",percent.EMIC.SAMP.CIhi,")",sep="")
		SMIC <- paste(signif(percent.SMIC.SIM$V1, digits = 3)," (",percent.SMIC.SIM.CIlo," - ",percent.SMIC.SIM.CIhi,")",sep="")
		MIC.data <- data.frame(Value = c(EMIC,SMIC),
													row.names = c("Bayes Estimated (i.e., Sampled Dose)","Predicted (i.e., Next 3 Doses)"))
	}
	MIC.data2 <- mic.function2(rFIT.data(),rCI.data2(),rSIM.data2())
})

#Extract input information AND individual bayes estimates and add to "saved_data.csv"
rNEW.data <- reactive({
	#Calculate AGE (years)
	AGE <- round(as.double(difftime(input$DDATE,input$BDATE))/365.25, digits = 1)
	#Assign labels to SEX
	if (input$SEX == 1) {SEXf <- "Male"}
	if (input$SEX == 2) {SEXf <- "Female"}
	#Assign NA's to CONCs and TIMEs
	CONC2 <- NA
	TIME2 <- NA
	CONC3 <- NA
	TIME3 <- NA
	CONC4 <- NA
	TIME4 <- NA
	#Call in previous concentrations in their corresponding times
	if (input$NCONC >= 2) {
		CONC2 <- input$PCONC2
		TIME2 <- input$PTIME2
	}
	if (input$NCONC >= 3) {
		CONC3 <- input$PCONC3
		TIME3 <- input$PTIME3
	}
	if (input$NCONC == 4) {
		CONC4 <- input$PCONC4
		TIME4 <- input$PTIME4
	}
	#Assign labels to MIC
	if (input$MIC == 1) {MIC <- 0.25}
	if (input$MIC == 2) {MIC <- 0.5}
	if (input$MIC == 3) {MIC <- 1}
	if (input$MIC == 4) {MIC <- 2}
	if (input$MIC == 5) {MIC <- 4}
	if (input$MIC == 6) {MIC <- 8}
	if (input$MIC == 7) {MIC <- 16}
	if (input$MIC == 8) {MIC <- 32}
	if (input$MIC == 9) {MIC <- 64}
	save.ind.function <- function(OBS.data,BAYES.data,MIC.data1) {
	#Make a data frame for the individual to be added to saved.data (via rbind)
	ind.data <- data.frame(Date = format(input$DDATE, "%d/%m/%Y"),
													URN = input$URN,
													First.Name = input$FNAME,
													Last.Name = input$LNAME,
													D.O.B. = format(input$BDATE, "%d/%m/%Y"),
													Age = AGE,
													SeCreat = input$SECR,
													Weight = input$WT,
													Gender = SEXf,
													CreatCL = round(OBS.data$CRCL[1], digits = 2),
													Dose.Amount = input$PDOSE,
													Infusion.Duration = input$PINFD,
													CONC1 = input$PCONC1,
													TIME1 = input$PTIME1,
													CONC2,
													TIME2,
													CONC3,
													TIME3,
													CONC4,
													TIME4,
													CL = round(BAYES.data$par[1], digits = 2),
													V = round(BAYES.data$par[2], digits = 2),
													Recommended.Dose.Amount = input$SDOSE1,
													Recommended.Infusion.Duration = input$SINFD1,
													Target.MIC = MIC,
													Predicted.Time.Above.MIC = MIC.data1$Value[2],
													Pharmacist = input$PHARMI)
	#new.data <- ind.data
	new.data <- rbind(saved.data,ind.data)
	}
	NEW.data <- save.ind.function(rOBS.data(),rBAYES.data(),rMIC.data1())
})
#Don't add individual's data to existing file until button is clicked
observe({
	if (input$ADD == T) {	#Linked to action button for adding ind.data to "saved_data.csv"
		isolate({	#Needs to be isolated so that user can select when to add data to the saved.data data frame
			write.csv(rNEW.data(), file=paste(dir,"saved_data.csv",sep=""), na=".", quote=F, row.names=F)
		})	#Brackets closing "isolate"
	}
})

#############
###_OUTPUT_##
#############
#Show calculation for creatinine clearance underneath input for serum creatinine
output$CrCLText <- renderText({
	CrCL.function <- function(OBS.data) {
		CrCL <- paste("Creatinine Clearance (mL/min): ",round(OBS.data$CRCL[1], digits = 1),sep="")
	}
	CrCL <- CrCL.function(rOBS.data())
})

#Plot the population predicted and individual Bayes predicted concentrations
output$concPlot1 <- renderPlot({
	#Plot
	plotobj1 <- NULL
	plotobj1 <- ggplot(rFIT.data())
	plotobj1 <- plotobj1 + stat_summary(aes(x = TIME, y = IPRE), data = rCI.data1(), geom = "ribbon", fun.ymin = CI95lo, fun.ymax = CI95hi, alpha = 0.2)
	plotobj1 <- plotobj1 + geom_line(aes(x = TIME, y = IPRE), size = 1, colour = "steelblue3")  #Blue line = bayes estimated
	plotobj1 <- plotobj1 + geom_line(aes(x = TIME, y = IPRE), data = rSIM.data1(), size = 1, colour = "firebrick3", linetype = "longdash")  #Red line = simulated
	#Set up a labelled dashed-line representing the MIC depending on input
	#Had to set it up like this because having the SMIC1 object in geom_hline wasn't working
	if (input$MIC == 1) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 0.25), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 2) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 0.5), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 3) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 4) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 2), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 5) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 4), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 6) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 8), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 7) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 16), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 8) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 32), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 9) {
	plotobj1 <- plotobj1 + geom_hline(aes(yintercept = 64), linetype = "dashed", colour = "black")
	}
	plotobj1 <- plotobj1 + geom_point(aes(x = TIME, y = DV), size = 4)  #Black circles = observed data used to estimate
	plotobj1 <- plotobj1 + scale_x_continuous("\nTime (hours)")
	#Find the maximum concentration for the y-axis
	max.CONC <- max.function(rFIT.data(),rCI.data1(),rSIM.data1())
	plotobj1 <- plotobj1 + scale_y_continuous("Piperacillin Concentration (mg/L)\n", lim = c(0,max.CONC))
	print(plotobj1)
})

#Plot the population predicted and individual Bayes predicted concentrations
output$concPlot2 <- renderPlot({
	#Plot
	plotobj2 <- NULL
	plotobj2 <- ggplot(rFIT.data())
	plotobj2 <- plotobj2 + stat_summary(aes(x = TIME, y = IPRE), data = rCI.data2(), geom = "ribbon", fun.ymin = CI95lo, fun.ymax = CI95hi, alpha = 0.2)
	plotobj2 <- plotobj2 + geom_line(aes(x = TIME, y = IPRE), size = 1, colour = "steelblue3")  #Blue line = bayes estimated
	plotobj2 <- plotobj2 + geom_line(aes(x = TIME, y = IPRE), data = rSIM.data2(), size = 1, colour = "firebrick3", linetype = "longdash")  #Red line = simulated
	#Set up a labelled dashed-line representing the MIC depending on input
	#Had to set it up like this because having the SMIC1 object in geom_hline wasn't working
	if (input$MIC == 1) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.25), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 2) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 0.5), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 3) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 1), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 4) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 2), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 5) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 4), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 6) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 8), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 7) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 16), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 8) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 32), linetype = "dashed", colour = "black")
	}
	if (input$MIC == 9) {
	plotobj2 <- plotobj2 + geom_hline(aes(yintercept = 64), linetype = "dashed", colour = "black")
	}
	plotobj2 <- plotobj2 + geom_point(aes(x = TIME, y = DV), size = 4)  #Black circles = observed data used to estimate
	plotobj2 <- plotobj2 + scale_x_continuous("\nTime (hours)")
	#Find the maximum concentration for the y-axis
	max.CONC <- max.function(rFIT.data(),rCI.data2(),rSIM.data2())
	plotobj2 <- plotobj2 + scale_y_continuous("Piperacillin Concentration (mg/L)\n", lim = c(0,max.CONC))
	print(plotobj2)
})

#Create table output for time above MIC
output$micTable1 <- renderTable({
	rMIC.data1()
})

#Create table output for time above MIC
output$micTable2 <- renderTable({
	rMIC.data2()
})

#Subset "saved.data.csv" by URN
output$prevTable <- renderTable({
	prev.data <- subset(saved.data, URN == input$URN & Last.Name == input$LNAME)
	prev.data <- prev.data[with(prev.data, order(prev.data$Date, decreasing = TRUE)), ]
	#prev.data <- prev.data[-c(5,7,9,11:22)]
})

#Generate a Word document of patient summary results
output$downloadReport <- downloadHandler(
	filename = function() {
		paste(format(input$DDATE,"%Y-%m-%d"),input$LNAME,input$URN,"Piperacillin_Report.docx",sep="_")
	},
	content = function(file) {
		src <- normalizePath(paste(dir,"report.Rmd",sep=""))
		inputEnv <- new.env()
		inputEnv$dir <- dir
		inputEnv$n <- n
		inputEnv$rOBS.data <- rOBS.data()
		inputEnv$rBAYES.data <- rBAYES.data()
		inputEnv$rFIT.data <- rFIT.data()
		inputEnv$rSIM.data1 <- rSIM.data1()
		inputEnv$rCI.data1 <- rCI.data1()
		inputEnv$rMIC.data1 <- rMIC.data1()
		#Temporarily switch to the temp dir, in case you don't have the permission to write to the current working directory
		owd <- setwd(tempdir())
		on.exit(setwd(owd))
		file.copy(src, "report.Rmd")
		Sys.setenv(RSTUDIO_PANDOC=pandocdir)
		out <- render("report.Rmd", word_document(fig_width = 8, fig_height = 6, reference_docx = paste(dir,"mystyles.docx",sep="")), envir = inputEnv)
		#out <- render("report.Rmd", pdf_document(fig_width = 8, fig_height = 6), envir = inputEnv)
		file.rename(out, file)
	}
)	#Brackets closing "downloadHandler"

#Make a test output table
output$testTable <- renderTable({
	head(rSIM.data1(),50)
})

#############
##_SESSION_##
#############
#Close the R session when Chrome closes
session$onSessionEnded(function() {
	stopApp()
})

#-------------------------------------------------------------------------------------
})	#Brackets closing "server"
