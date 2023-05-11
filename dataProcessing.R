#####################
#Prepare regions and ship masks
#####################
rm(list = ls())
library("viridis")
library("ncdf4")
library("fields")

getInd <- function(x) {
  c(which(lon == x[1]), which(lat == x[2]))
}

#load all data
all <- nc_open("/Users/michaeldiamond/Documents/Data/CERES/FluxByCldType/FBCT_shipkrige.nc")
lat <- rev(ncvar_get(all, "lat"))
lon <- ncvar_get(all, "lon")

#Read in and compute annual mean of SO2
EDGAR_SO2 <- ncvar_get(all, "EDGAR_SO2")
EDGAR_SO2 <- EDGAR_SO2[,180:1,]
EDGAR_SO2_AnnAv <- apply(EDGAR_SO2, 1:2, mean)

#regions
xNE <- which(lon >= -12.5 & lon <= 2.5); nXNE <- length(xNE)
yNE <- which(lat >= -10 & lat <= -2); nYNE <- length(yNE)
xSE <- which(lon >= -7.5 & lon <= 7.5); nXSE <- length(xSE)
ySE <- which(lat >= -18 & lat <= -10); nYSE <- length(ySE)
stopifnot(nXNE == nXSE); stopifnot(nYNE == nYSE)
nX <- nXNE; nY <- nYNE

#central region
xC <- which(lon >= -12.5 & lon <= 7.5); nXC <- length(xC)
yC <- which(lat >= -18 & lat <= -8); nYC <- length(yC)

#central region shifted
xCs <- which(lon >= -12.5 & lon <= 7.5); nXCs <- length(xCs)
yCs <- which(lat >= -18 & lat <= -8); nYCs <- length(yCs)

#define region and ship track based on where the maximum value S02 value is
possBoxNE <- EDGAR_SO2_AnnAv[xNE[1]:xNE[nX], yNE[1]:yNE[nY]]
rowMaxNE <- apply(possBoxNE, 2, which.max)
shipNE <- shipNEInd <-  matrix(ncol = 2, nrow = 0)
bdShip <- 2
for (i in 1:nY) {
  l1 <- rowMaxNE[i] - bdShip + xNE[1] - 1
  l2 <- rowMaxNE[i] + bdShip + xNE[1] - 1
  shipNE <- rbind(shipNE, cbind(lon[l1:l2], rep(lat[yNE[i]], 5)))
  shipNEInd <- rbind(shipNEInd, cbind(l1:l2, rep(yNE[i], 5)))
}

possBoxSE <- EDGAR_SO2_AnnAv[xSE[1]:xSE[nX], ySE[1]:ySE[nY]]
rowMaxSE <- apply(possBoxSE, 2, which.max)
shipSE <- shipSEInd <-  matrix(ncol = 2, nrow = 0)
for (i in 1:nY) {
  l1 <- rowMaxSE[i] - bdShip + xSE[1] - 1
  l2 <- rowMaxSE[i] + bdShip + xSE[1] - 1
  shipSE <- rbind(shipSE, cbind(lon[l1:l2], rep(lat[ySE[i]], 5)))
  shipSEInd <- rbind(shipSEInd, cbind(l1:l2, rep(ySE[i], 5)))
}

possBoxC <- EDGAR_SO2_AnnAv[xC[1]:xC[nXC], yC[1]:yC[nYC]]
rowMaxC <- apply(possBoxC, 2, which.max)
shipC <- shipCInd <-  matrix(ncol = 2, nrow = 0)
for (i in 1:nYC) {
  l1 <- rowMaxC[i] - bdShip + xC[1] - 1 -2
  l2 <- rowMaxC[i] + bdShip + xC[1] - 1
  shipC <- rbind(shipC, cbind(lon[l1:l2], rep(lat[yC[i]], 5)))
  shipCInd <- rbind(shipCInd, cbind(l1:l2, rep(yC[i], 5)))
}

possBoxCs <- EDGAR_SO2_AnnAv[xCs[1]:xCs[nXCs], yCs[1]:yCs[nYCs]]
rowMaxCs <- apply(possBoxCs, 2, which.max)
shipCs <- shipCsInd <-  matrix(ncol = 2, nrow = 0)
for (i in 1:nYCs) {
  l1 <- rowMaxCs[i] - bdShip + xCs[1] - 1 -2 -2
  l2 <- rowMaxCs[i] + bdShip + xCs[1] - 1 -2
  shipCs <- rbind(shipCs, cbind(lon[l1:l2], rep(lat[yCs[i]], 5)))
  shipCsInd <- rbind(shipCsInd, cbind(l1:l2, rep(yCs[i], 5)))
}

#Make NE and SE boxes
NEBox <- SEBox <- matrix(nrow = 360, ncol = 180, data = FALSE)
NEBoxInd <- expand.grid(xNE, yNE)
NEBox[NEBoxInd[,1], NEBoxInd[,2]] <- TRUE
SEBoxInd <- expand.grid(xSE, ySE)
SEBox[SEBoxInd[,1], SEBoxInd[,2]] <- TRUE

CBox <- matrix(nrow = 360, ncol = 180, data = FALSE)
CBoxInd <- expand.grid(xC, yC)
CBox[CBoxInd[,1], CBoxInd[,2]] <- TRUE

CsBox <- matrix(nrow = 360, ncol = 180, data = FALSE)
CsBoxInd <- expand.grid(xCs, yCs)
CsBox[CsBoxInd[,1], CsBoxInd[,2]] <- TRUE

#Remove ship track from NE and SE box
NEBox[shipNEInd] <- FALSE
SEBox[shipSEInd] <- FALSE
CBox[shipCInd] <- FALSE
CsBox[shipCsInd] <- FALSE

#Find indices of NE and SE grid boxes not in ship track
notShipNEInd <- which(NEBox, arr.ind = T)
notShipSEInd <- which(SEBox, arr.ind = T)
notShipNE <- cbind(lon[notShipNEInd[,1]], lat[notShipNEInd[,2]])
notShipSE <- cbind(lon[notShipSEInd[,1]], lat[notShipSEInd[,2]])

notShipCInd <- which(CBox, arr.ind = T)
notShipC <- cbind(lon[notShipCInd[,1]], lat[notShipCInd[,2]])

notShipCsInd <- which(CsBox, arr.ind = T)
notShipCs <- cbind(lon[notShipCsInd[,1]], lat[notShipCsInd[,2]])

#shift west by 90 and south by 5 to create second set of regions 
notShipNW <- cbind(notShipNE[,1] - 90, notShipNE[,2] - 5)
notShipSW <- cbind(notShipSE[,1] - 90, notShipSE[,2] - 5)
shipNW <- cbind(shipNE[,1] - 90, shipNE[,2] - 5)
shipSW <- cbind(shipSE[,1] - 90, shipSE[,2] - 5)

#create indices for west regions
notShipNWInd <- t(apply(notShipNW, 1, getInd))
shipNWInd <- t(apply(shipNW, 1, getInd))
notShipSWInd <- t(apply(notShipSW, 1, getInd))
shipSWInd <- t(apply(shipSW, 1, getInd))

#Plot regions to visualize
pdf("/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/regionMap.pdf", height = 4, width = 6)
zMin <- min(log(EDGAR_SO2_AnnAv + .01)); zMax <- max(log(EDGAR_SO2_AnnAv + .01))
image.plot(lon, lat, log(EDGAR_SO2_AnnAv + .01), col = viridis(10),
           zlim = c(zMin, zMax), main = "Log SO2 Emissions \n Annual Average")
points(notShipNE, pch = 20, col = 'red', cex = .01)
points(notShipSE, pch = 20, col = 'blue', cex = .01)
points(notShipNW, pch = 20, col = 'yellow', cex = .01)
points(notShipSW, pch = 20, col = 'green', cex = .01)
#dev.off()
points(notShipC, pch = 20, col = 'purple', cex = .01)


points(shipNE, pch = 20, cex = .01)
points(shipSE, pch = 20, cex = .01)
points(shipNW, pch = 20, cex = .01)
points(shipSW, pch = 20, cex = .01)


#repeat plotting with (x, y) coords to check ind values
image.plot(1:360, 1:180, log(EDGAR_SO2_AnnAv + .01), col = viridis(10),
           zlim = c(zMin, zMax))
points(notShipNEInd, pch = 20, col = 'red', cex = .01)
points(notShipSEInd, pch = 20, col = 'blue', cex = .01)
points(notShipNWInd, pch = 20, col = 'yellow', cex = .01)
points(notShipSWInd, pch = 20, col = 'green', cex = .01)

points(shipNEInd, pch = 20, cex = .01)
points(shipSEInd, pch = 20, cex = .01)
points(shipNWInd, pch = 20, cex = .01)
points(shipSWInd, pch = 20, cex = .01)


regInfo <- list("notShipNE" = notShipNE, "notShipSE" = notShipSE, 
                "notShipNW" = notShipNW, "notShipSW" = notShipSW, 
                "notShipC" = notShipC, "notShipCs" = notShipCs,
                "shipNE" = shipNE, "shipSE" = shipSE,
                "shipNW" = shipNW, "shipSW" = shipSW, 
                "shipC" = shipC, "shipCs" = shipCs,
                "notShipNEInd" = notShipNEInd, "notShipSEInd" = notShipSEInd, 
                "notShipNWInd" = notShipNWInd, "notShipSWInd" = notShipSWInd, 
                "notShipCInd" = notShipCInd, "notShipCsInd" = notShipCsInd,
                "shipNEInd" = shipNEInd, "shipSEInd" = shipSEInd,
                "shipNWInd" = shipNWInd, "shipSWInd" = shipSWInd, 
                "shipCInd" = shipCInd, "shipCsInd" = shipCsInd)
save(regInfo, file = "/Users/michaeldiamond/Documents/Projects/IMO2020_prelim/processed/regInfo.rda")

