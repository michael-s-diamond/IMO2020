ilogit <- function(x){
  exp(x)/(1 + exp(x))
}

logit <- function(p) {
  log(p/(1-p))
}


analyze <- function(respName, regName, covars,
                    covModel, iniPhi, iniSigma2, 
                    months, lon, lat, nColors, savenc=FALSE) {
  pdf(sprintf("/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/Results/Figures_%s_%s_M%sto%s.pdf", respName, regName, months[1], months[length(months)]))
  #Information about the region 
  if (regName == "NE") {
    notShip <- regInfo$notShipNE; ship <- regInfo$shipNE
    notShipInd <- regInfo$notShipNEInd; shipInd <- regInfo$shipNEInd
  } else if (regName == 'NW') {
    notShip <- regInfo$notShipNW; ship <- regInfo$shipNW
    notShipInd <- regInfo$notShipNWInd; shipInd <- regInfo$shipNWInd
  } else if (regName == 'C') {
    notShip <- regInfo$notShipC; ship <- regInfo$shipC
    notShipInd <- regInfo$notShipCInd; shipInd <- regInfo$shipCInd
  } else if (regName == 'Cs') {
    notShip <- regInfo$notShipCs; ship <- regInfo$shipCs
    notShipInd <- regInfo$notShipCsInd; shipInd <- regInfo$shipCsInd
  } else if (regName == "SE") {
    notShip <- regInfo$notShipSE; ship <- regInfo$shipSE
    notShipInd <- regInfo$notShipSEInd; shipInd <- regInfo$shipSEInd
  } else {
    notShip <- regInfo$notShipSW; ship <- regInfo$shipSW
    notShipInd <- regInfo$notShipSWInd; shipInd <- regInfo$shipSWInd
  }
  
  #response variables
  resp <- ncvar_get(all, respName)
  resp <- resp[,180:1, months]
  respMean <- apply(resp, 1:2, mean)
  if (tran == "log") {
    respMean <- log(respMean)
  } else if (tran == "logit") {
    respMean <- suppressWarnings(logit(respMean)) #ok to suppress warning bc issue is values outside (0,1) on land
  } 
  
  #means of all covariates
  nCovars <- length(covars)
  covarMeans <- array(dim = c(nCovars, 360, 180))
  for (i in 1:nCovars) {
    covar <- ncvar_get(all, covars[i])
    covar <- covar[,180:1, months]
    covarMean <- apply(covar, 1:2, mean)
    covarMeans[i,,] <- covarMean
  }

  ####################################
  #Regressions to identify trend
  ####################################
  respY <- respMean[notShipInd]
  covarD <- matrix(nrow = nrow(notShipInd), ncol = nCovars)
  for (i in 1:nCovars) {
    covarD[,i] <- as.vector(covarMeans[i,,][notShipInd])
  }
  covarD <- cbind(covarD, notShip[,1], notShip[,1]^2, notShip[,2], notShip[,2]^2,
                  notShip[,1]*notShip[,2])
  covarNames <- c(covars, "lon",  "lon2", "lat", "lat2", "lon*lat")
  colnames(covarD) <- covarNames
  
  #model selection
  modSelec <- regsubsets(x = covarD, y = respY, nbest = 1, method = "exhaustive")
  sumModSelec <- summary(modSelec)
  parSel <- sumModSelec$which[which.min(sumModSelec$bic),]
  parSel <- parSel[-1] #don't want intercept recorded in which parameters we keep in the model
  parSel <- which(parSel)
  lm1 <- lm(respY ~ covarD[,parSel])
  par(mfrow = c(1, 1))
  plot(lm1, main = sprintf("Linear Regression Check: %s_%s", respName, regName)) 
  
  xBd <- max(abs(min(lm1$residuals)), abs(max(lm1$residuals))) + .01
  hist(lm1$residuals, xlim = c(-xBd, xBd), 
       main = sprintf("Linear Regression Check: %s_%s", respName, regName))
  
  #Look at observed vs. fitted values in regressions
  bbX <- min(notShipInd[,1] - 1):max(notShipInd[,1] + 1); nX <- length(bbX)
  bbY <- min(notShipInd[,2] - 1):max(notShipInd[,2] + 1); nY <- length(bbY)

  regFit <- matrix(nrow = 360, ncol = 180, data = NA)
  regFit[notShipInd] <-  lm1$fitted.values
  respMeanNoShip <- respMean
  respMeanNoShip[is.na(regFit)] <- NA
  
  zMin <- min(min(regFit[bbX, bbY], na.rm = T), 
              min(respMeanNoShip[bbX, bbY], na.rm = T)) 
  zMax <- max(max(regFit[bbX, bbY], na.rm = T), 
              max(respMeanNoShip[bbX, bbY], na.rm = T))

  par(mfrow = c(1, 2))
  image.plot(lon[bbX], lat[bbY], regFit[bbX, bbY], 
             zlim = c(zMin, zMax), col = viridis(nColors),
             xlab = "", ylab = "", main = sprintf("Estimated %s_%s", respName, regName))
  image.plot(lon[bbX], lat[bbY], respMeanNoShip[bbX, bbY], zlim = c(zMin, zMax), 
             col = viridis(nColors),
             xlab = "", ylab = "", main = sprintf("Observed %s_%s", respName, regName))
  resid <-  respMeanNoShip - regFit
  par(mfrow = c(1, 1))
  residBd <- max(abs(min(resid, na.rm = T)), abs(max(resid, na.rm = T)))
  image.plot(lon[bbX], lat[bbY], resid[bbX, bbY],  col = brewer.pal(nColors,'RdYlBu'),
             main = sprintf("Observed - Fitted  %s_%s", respName, regName), zlim = c(-residBd, residBd))
  
  ####################################
  #Variogram Fitting
  ###################################
  library("geoR")
  
  #Make geoDat object
  NdGeoDat <- as.geodata(cbind(notShip, respY,  covarD[,parSel]),
                         coords.col = 1:2, data.col = 3,
                         covar.col = 4:(4 + length(parSel) - 1))
  
  #Check that geoDat object looks like data
  par(mfrow = c(1, 2))
  points(NdGeoDat, main = sprintf("Check GeoDat:\n %s %s", respName, regName))             
  image.plot(lon[bbX], lat[bbY], respMeanNoShip[bbX, bbY], col = viridis(nColors),
             xlab = "", ylab = "",
             main = sprintf("Check GeoDat:\n %s %s", respName, regName))
  
  #Look at variogram as an exploratory tool to determine form of covar structure
  vCloud <- variog(NdGeoDat, trend=~covarD[,parSel], option = "cloud")
  vBin <- variog(NdGeoDat, trend=~covarD[,parSel], option = "bin")

  
  #model parameter estimates using WLS
  wls <- variofit(vBin, ini = c(iniSigma2, iniPhi), cov.model = covModel,
                  fix.nugget = TRUE)
  
  #Check krige fits against data
  par(mfrow = c(1, 2))
  plot(vBin, main = sprintf("Binned Variogram  \n%s_%s", respName, regName),
       pch = 20)
  lines(wls, col = 'blue')
  
  plot(vCloud, pch = 20,  cex = .1,
       main = sprintf("Cloud Variogram \n%s_%s", respName, regName))
  lines(wls, col = 'blue')

  
  ##############################
  #Perform universal kriging 
  #############################
  covarL <- matrix(nrow = nrow(shipInd), ncol = nCovars)
  for (i in 1:nCovars) {
    covarL[,i] <- as.vector(covarMeans[i,,][shipInd])
  }
  covarL <- cbind(covarL, ship[,1], ship[,1]^2, ship[,2], ship[,2]^2,
                  ship[,1]*ship[,2])
  colnames(covarL) <- covarNames
  KC <- krige.control(trend.d=~covarD[,parSel], trend.l=~covarL[,parSel],obj=wls)
  outCont <- output.control(simulations = TRUE, n.pred = 5000)
  kr <- krige.conv(NdGeoDat, loc=ship,  krige=KC, output = outCont)   
  
  #################
  #Significance testing
  respEstLB <- respEstUB <- respEst  <- matrix(nrow = 360, ncol = 180, 
                                               data = respMean)
  respEst[shipInd] <- kr$predict
  respEstLBVec <- kr$predict - qnorm(.975)*sqrt(kr$krige.var)
  respEstUBVec <- kr$predict + qnorm(.975)*sqrt(kr$krige.var)
  respEstLB[shipInd] <- respEstLBVec
  respEstUB[shipInd] <- respEstUBVec
  
  #Outside of intervals
  outIndB <- which(respMean < respEstLB, arr.ind = T)
  nOutB <- nrow(outIndB)
  outIndA <- which(respMean > respEstUB, arr.ind = T)
  nOutA <- nrow(outIndA)
  if (nOutA > nOutB) {
    dir <- "above"
    nOut <- nOutA
  } else if (nOutB > nOutA) {
    dir <- "below"
    nOut <- nOutB
  } else {
    dir <- "tied"
    nOut <- NA
  }
  nTotal <- nrow(shipInd)
  perOut <- nOut/nTotal
  
  #Simulations under the null; also store mean difference
  nSim <- outCont$n.predictive
  if (dir != "tied") {
    nOutSim <- rep(NA, nSim)
  } else {
    nOutSimA <- rep(NA, nSim)
    nOutSimB <- rep(NA, nSim)
  }
  meanDiff <- rep(NA, nSim)
  for (i in 1:nSim) {
    respSim <- kr$simulations[,i]
    if (tran == "logit") {
      meanDiff[i] <- mean(ilogit(respMean[shipInd]) - ilogit(respSim))
    } else if (tran == "log") {
      meanDiff[i] <- mean(exp(respMean[shipInd]) - exp(respSim))
    } else {
      meanDiff[i] <- mean(respMean[shipInd] - respSim)
    }
    if (dir == "above") {
      outSimCurr <- which(respSim > respEstUBVec)
      nOutSim[i] <- length(outSimCurr)
    } else if (dir == "below") {
      outSimCurr <- which(respSim < respEstLBVec)
      nOutSim[i] <- length(outSimCurr)
    } else if (dir == "tied") {
      outSimCurrA <-  which(respSim > respEstUBVec)
      outSimCurrB <-  which(respSim < respEstLBVec)
      nOutSimA[i] <- length(outSimCurrA)
      nOutSimB[i] <- length(outSimCurrB)
    }
  }
  
  if (dir == "tied") {
    test <- sum(nOutSimA >= nOutA)/nSim >= sum(nOutSimB >= nOutB)/nSim
    if (test) {
      dir <- "above"
      nOutSim <- nOutSimA
      nOut <- nOutA
    } else {
      dir <- "below"
      nOutSim <- nOutSimB
      nOut <- nOutB
    }
  }
  
  pVal <- sum(nOutSim >= nOut)/nSim
  par(mfrow = c(1, 1))
  xlim = c(0, max(c(nOutSim, nOut)) + 3)
  hist(nOutSim, xlab = "Number of Significant Grid Boxes", xlim = xlim,
       main = sprintf("Generated Null Distribution \n%s_%s", respName, regName))
  abline(v = nOut, col = "red")
    
  #plot
  zMin <- min(min(respEst[bbX, bbY]), min(respMean[bbX, bbY]), 
              min(respEstLB[bbX, bbY]), min(respEstUB[bbX, bbY])) 
  zMax <- max(max(respEst[bbX, bbY]), max(respMean[bbX, bbY]), 
              max(respEstLB[bbX, bbY]), max(respEstUB[bbX, bbY])) 
    
  par(mfrow = c(2, 2))
  image.plot(lon[bbX], lat[bbY], respEst[bbX, bbY], 
              main = sprintf("Predicted \n%s_%s", respName, regName),
              zlim = c(zMin, zMax), col = viridis(nColors), xlab = "", ylab = "")
  points(ship, pch = 0, col = 'white', cex = .25)
  points(lon[outIndA[,1]], lat[outIndA[,2]], col = 'red', pch = 20, cex = .5)
  points(lon[outIndB[,1]], lat[outIndB[,2]], col = 'blue', pch = 20, cex = .5)
  image.plot(lon[bbX], lat[bbY], respMean[bbX, bbY], 
             main = sprintf("Observed \n%s_%s", respName, regName), 
             zlim = c(zMin, zMax), col = viridis(nColors), xlab = "", ylab = "")
  points(ship, pch = 0, col = 'white', cex = .25)
  points(lon[outIndA[,1]], lat[outIndA[,2]], col = 'red', pch = 20, cex = .5)
  points(lon[outIndB[,1]], lat[outIndB[,2]], col = 'blue', pch = 20, cex = .5)
  image.plot(lon[bbX], lat[bbY], respEstLB[bbX, bbY], 
             main = sprintf("Interval Lower Bound \n%s_%s", respName, regName),
             zlim = c(zMin, zMax), col = viridis(nColors), xlab = "", ylab = "")
  points(ship, pch = 0, col = 'white', cex = .25)
  points(lon[outIndA[,1]], lat[outIndA[,2]], col = 'red', pch = 20, cex = .5)
  points(lon[outIndB[,1]], lat[outIndB[,2]], col = 'blue', pch = 20, cex = .5)
  image.plot(lon[bbX], lat[bbY], respEstUB[bbX, bbY],
             main = sprintf("Interval Upper Bound \n%s_%s", respName, regName), 
             zlim = c(zMin, zMax), col = viridis(nColors), xlab = "", ylab = "")
  points(ship, pch = 0, col = 'white', cex = .25)
  points(lon[outIndA[,1]], lat[outIndA[,2]], col = 'red', pch = 20, cex = .5)
  points(lon[outIndB[,1]], lat[outIndB[,2]], col = 'blue', pch = 20, cex = .5)
  mtext("Universal Kriging", outer = TRUE, font =2, cex = 1.2)
    
  #Check that residuals outside ship track are smaller than 
  #residuals inside ship track (if significant)
  zBd <- max(max(abs(resid[bbX, bbY]), na.rm= T), 
              max(abs(respMean[bbX, bbY] - respEst[bbX, bbY])))
  par(mfrow = c(1, 2))
  image.plot(lon[bbX], lat[bbY], resid[bbX, bbY],  col = brewer.pal(nColors,'RdYlBu'),
              main = sprintf("Observed - Fitted (Outside Track) \n%s_%s",
                            respName, regName),
              zlim = c(-zBd, zBd))
  image.plot(lon[bbX], lat[bbY], respMean[bbX, bbY] - respEst[bbX, bbY], 
              main = sprintf("Observed - Fitted (In Track) \n%s_%s",
                            respName, regName),
              col = brewer.pal(nColors,'RdYlBu'), xlab = "", ylab = "",
              zlim = c(-zBd, zBd))
  points(lon[outIndA[,1]], lat[outIndA[,2]], col = 'red', pch = 20, cex = .5)
  points(lon[outIndB[,1]], lat[outIndB[,2]], col = 'blue', pch = 20, cex = .5)
  
  par(mfrow = c(1, 1))
  hist(meanDiff, main = sprintf("Delta(%s_%s) = %.4f +/- %.4f", respName, regName, mean(meanDiff),2*sd(meanDiff))) 
  dev.off()
  
  print(pVal)
  
  if (savenc) {
    respEst[respEst == -Inf] <- -999
    respEstLB[respEstLB == -Inf] <- -999
    respEstUB[respEstUB == -Inf] <- -999
    respMean[respMean == -Inf] <- -999
    
    londim <- ncdim_def("lon","degrees_east",lon) 
    latdim <- ncdim_def("lat","degrees_north",lat)
    simdim <- ncdim_def("nsim","NA",c(1:nSim))
    shipdim <- ncdim_def("ship_dim","NA",c(1:70))
    bindim <- ncdim_def("bins","NA",vBin$u)
    fillvalue <- -999
    Est <- ncvar_def(name="Est",units="",dim=list(londim,latdim),missval=fillvalue)
    lowEst <- ncvar_def(name="lowEst",units="",dim=list(londim,latdim),missval=fillvalue)
    highEst <- ncvar_def(name="highEst",units="",dim=list(londim,latdim),missval=fillvalue)
    Obs <- ncvar_def(name="Obs",units="",dim=list(londim,latdim),missval=fillvalue)
    krSims <- ncvar_def(name="krSims",units="",dim=list(shipdim,simdim),missval=fillvalue)
    shipObs <- ncvar_def(name="shipObs",units="",dim=list(shipdim),missval=fillvalue)
    Res <- ncvar_def(name="Residual",units="",dim=list(londim,latdim),missval=fillvalue)
    semivar <- ncvar_def(name="Semivariance",units="",dim=list(bindim),missval=fillvalue)
    ncfname <- sprintf("/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/Results/Data_%s_%s_M%sto%s.nc", respName, regName, months[1], months[length(months)])
    ncout <- nc_create(filename = ncfname, vars=list(Est,lowEst,highEst,Obs,krSims,shipObs,Res,semivar), force_v4=TRUE, verbose=FALSE)
    ncvar_put(nc=ncout, varid=Est, vals=respEst)
    ncvar_put(nc=ncout, varid=lowEst, vals=respEstLB)
    ncvar_put(nc=ncout, varid=highEst, vals=respEstUB)
    ncvar_put(nc=ncout, varid=Obs, vals=respMean)
    ncvar_put(nc=ncout, varid=krSims, vals=kr$simulations)
    ncvar_put(nc=ncout, varid=shipObs, vals=respMean[shipInd])
    ncvar_put(nc=ncout, varid=Res, vals=resid)
    ncvar_put(nc=ncout, varid=semivar, vals=vBin$v)
    ncatt_put(nc=ncout, varid=0, attname="pVal", pVal)
    ncatt_put(nc=ncout, varid=0, attname="nOut", nOut)
    ncatt_put(nc=ncout, varid=0, attname="tran", tran)
    ncatt_put(nc=ncout, varid=0, attname="iniPhi", iniPhi)
    ncatt_put(nc=ncout, varid=0, attname="iniSigma2", iniSigma2)
    ncatt_put(nc=ncout, varid=0, attname="Phi", wls$cov.pars[2])
    ncatt_put(nc=ncout, varid=0, attname="Sigma2", wls$cov.pars[1])
    ncatt_put(nc=ncout, varid=0, attname="respName", respName)
    ncatt_put(nc=ncout, varid=0, attname="regName", regName)
    ncatt_put(nc=ncout, varid=0, attname="months", months)
    ncatt_put(nc=ncout, varid=0, attname="parSel", paste(covarNames[parSel],collapse=','))
    nc_close(ncout)
  }
  
  
  return(list("parSel" = covarNames[parSel], "dir" = dir, "nOut" = nOut, 
              "nTotal" = nTotal,  "perOut" = perOut, "pVal" = pVal,
              "covModel" = covModel, "meanDiff" = meanDiff))
}

