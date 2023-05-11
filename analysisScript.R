#####################
#Analysis of potential cloud changes from 2020 IMO regulations
#####################
rm(list = ls())
set.seed(103) #pseudo-random number generator
library("ncdf4")
library("fields")
library("viridis")
library("RColorBrewer")
library("leaps")
source("/Users/michaeldiamond/Code/IMO2020_prelim/analysisFns.R")

#Set variogram covariance model
covModel <- "exp"

#Set region
regName = "C" #Update of SE region from D+20 that better centers the cloud anomalies

#Set colors
nColors <- 10

#Load data and region information
all <- nc_open("/Users/michaeldiamond/Documents/Data/CERES/SSF1deg/SSF1deg_shipkrige_Terra.nc")
lat <- rev(ncvar_get(all, "lat"))
lon <- ncvar_get(all, "lon")
load("/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/processed/regInfo.rda")

#
###Analysis scripts
#

#
###Terra: Climatology (2002-2019)
#
covars <- c("Terra_Tskin_clim","Terra_EIS_clim","Terra_WS_clim")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_clim"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .01, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_clim"), regName, covars, 
                covModel, iniPhi = 1.25, iniSigma2 = .00022, 
                months=months, lon, lat, nColors,savenc=TRUE)


#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_clim"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0125, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_clim"), regName, covars, 
                covModel, iniPhi = 2, iniSigma2 = .000275, 
                months, lon, lat, nColors,savenc=TRUE)

#
###Terra: 2002-2004
#

covars <- c("Terra_Tskin_2002","Terra_EIS_2002","Terra_WS_2002")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2002"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0275, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2002"), regName, covars, 
                covModel, iniPhi = 1, iniSigma2 = .00035, 
                months, lon, lat, nColors,savenc=TRUE)

#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2002"), regName , covars,
                covModel, iniPhi = 1, iniSigma2 = .02, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2002"), regName, covars, 
                covModel, iniPhi = 1.5, iniSigma2 = .0004, 
                months, lon, lat, nColors,savenc=TRUE)

#
###Terra: 2005-2007
#
covars <- c("Terra_Tskin_2005","Terra_EIS_2005","Terra_WS_2005")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2005"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .03, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2005"), regName, covars, 
                covModel, iniPhi = 1, iniSigma2 = .00035, 
                months, lon, lat, nColors,savenc=TRUE)

#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2005"), regName, covars,
                covModel, iniPhi = .25, iniSigma2 = .0125, 
                months, lon, lat, nColors,savenc=TRUE)
tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2005"), regName, covars, 
                covModel, iniPhi = 1.75, iniSigma2 = .00035, 
                months, lon, lat, nColors,savenc=TRUE)


#
###Terra: 2008-2010
#
covars <- c("Terra_Tskin_2008","Terra_EIS_2008","Terra_WS_2008")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2008"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .03, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2008"), regName, covars, 
                covModel, iniPhi = 1, iniSigma2 = .0006, 
                months, lon, lat, nColors,savenc=TRUE)


#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2008"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0225, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2008"), regName, covars, 
                covModel, iniPhi = 1.5, iniSigma2 = .000375, 
                months, lon, lat, nColors,savenc=TRUE)


#
###Terra: 2011-2013
#
covars <- c("Terra_Tskin_2011","Terra_EIS_2011","Terra_WS_2011")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2011"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0175, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2011"), regName, covars, 
                covModel, iniPhi = 1.5, iniSigma2 = .00045, 
                months, lon, lat, nColors,savenc=TRUE)

#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2011"), regName, covars,
                covModel, iniPhi = .5, iniSigma2 = .0175, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2011"), regName, covars, 
                covModel, iniPhi = 1.5, iniSigma2 = .00025, 
                months, lon, lat, nColors,savenc=TRUE)

#
###Terra: 2014-2016
#
covars <- c("Terra_Tskin_2014","Terra_EIS_2014","Terra_WS_2014")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2014"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0275, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2014"), regName, covars, 
                covModel, iniPhi = 1, iniSigma2 = .0004, 
                months, lon, lat, nColors,savenc=TRUE)

#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2014"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .02, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2014"), regName, covars, 
                covModel, iniPhi = 1.5, iniSigma2 = .0004, 
                months, lon, lat, nColors,savenc=TRUE)


#
###Terra: 2017-2019
#
covars <- c("Terra_Tskin_2017","Terra_EIS_2017","Terra_WS_2017")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2017"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0175, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2017"), regName, covars, 
                covModel, iniPhi = 1, iniSigma2 = .0005, 
                months, lon, lat, nColors,savenc=TRUE)

#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2017"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .015, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2017"), regName, covars, 
                covModel, iniPhi = 1.25, iniSigma2 = .0003, 
                months, lon, lat, nColors,savenc=TRUE)


#
###Terra: 2020-2022
#
covars <- c("Terra_Tskin_2020","Terra_EIS_2020","Terra_WS_2020")

#SON values
months = c(9:11)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2020"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .02, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2020"), regName, covars, 
                covModel, iniPhi = 1.25, iniSigma2 = .00035, 
                months, lon, lat, nColors,savenc=TRUE)

#ANN values
months = c(1:12)

tran = "none"
temp <- analyze(respName = sprintf("Terra_cer_2020"), regName, covars,
                covModel, iniPhi = 1, iniSigma2 = .0125, 
                months, lon, lat, nColors,savenc=TRUE)

tran = "logit"
temp <- analyze(respName = sprintf("Terra_Acld_2020"), regName, covars, 
                covModel, iniPhi = 1.25, iniSigma2 = .00022, 
                months, lon, lat, nColors,savenc=TRUE)

