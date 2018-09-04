### NOT PART OF THESIS ANALYSIS ###

library(GWmodel)
library(RColorBrewer)
library(rgdal)
library(leaflet)
library(psych)
library(classInt)

## Calling dataset
setwd("E\\GW_Thesis\\ThesisData")
dsn <-getwd()
ea.border<-readOGR(dsn = dsn, layer = "ea_Env_DHS")
# converting aridity index to original values
ea.border@data$Aridity <- ea.border@data$Aridity * 0.0001

## Local multidimensional outliers ############
colnames(ea.scaled.spdf@data)[22] <- "age_lessth_6"
colnames(ea.scaled.spdf@data)[23] <- "age6_11"
colnames(ea.scaled.spdf@data)[24] <- "age12_23"
colnames(ea.scaled.spdf@data)[25] <- "age24_35"
colnames(ea.scaled.spdf@data)[26] <- "age36_47"
colnames(ea.scaled.spdf@data)[27] <- "age48_59"

ea.mat <- as.matrix(ea.scaled.spdf@data[,1:27])
ea.outlier <- gwpca.cv.contrib(ea.mat,
                               Coords,
                               bw=bw.gwpca.k8,
                               k=8,
                               adaptive=TRUE)
boxplot(ea.outlier,horizontal=TRUE,xlab='Local PC Discrepancy')
which(ea.outlier > 60)

#GW PC plot 
## Editing gw.pcplot function
gw.pcplot.EA <- function (data, vars, focus, bw, adaptive = FALSE, ylim = NULL, 
                          ylab = "", fixtrans = FALSE, p = 2, theta = 0, longlat = F, 
                          dMat, ...) 
{
  if (is(data, "Spatial")) {
    p4s <- proj4string(data)
    loc <- coordinates(data)
  }
  else stop("Given data must be a Spatial*DataFrame")
  data <- as(data, "data.frame")
  dp.n <- nrow(data)
  i <- focus
  col.nm <- colnames(data)
  var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
  if (length(var.idx) == 0) 
    stop("Variables input doesn't match with data")
  x <- data[, var.idx]
  x <- as.matrix(x)
  m <- ncol(x)
  if (missing(dMat)) {
    DM.given <- F
    if (dp.n <= 5000) {
      dMat <- gw.dist(dp.locat = loc, p = p, theta = theta, 
                      longlat = longlat)
      DM.given <- T
    }
  }
  else {
    DM.given <- T
    dim.dMat <- dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
      stop("Dimensions of dMat are not correct")
  }
  if (DM.given) 
    dists <- dMat[, i]
  else dists <- gw.dist(dp.locat = loc, focus = i, p = p, theta = theta, 
                        longlat = longlat)
  if (adaptive) {
    rnk <- rank(dists, ties.method = "first")
    bw <- dists[rnk == bw]
  }
  nbrlist <- which(dists < bw)
  dists <- dists^2
  wts <- (1 - dists/(bw * bw))^12
  xss <- scale(x)
  span <- 1:m
  tsc <- 25/length(nbrlist)
  if (is.null(ylim)) 
    ylim <- c(min(xss[nbrlist, ]), max(xss[nbrlist, ]))
  plot(span, xss[i, ], type = "l", ylim = ylim, xlim = c(0.5, m + 0.5), 
       col = "red", lwd = 3, axes = F, xlab = "", 
       ylab = ylab, ...)
  axis(1, at = 1:m, labels = colnames(x), las = 2, cex.axis = 0.8)
  axis(2, at = seq(floor(ylim[1]), ceiling(ylim[2]),by = 1), 
       cex.axis = 0.8)
  abline(v = 1:m, col = grey(0.6))
  lines(c(1, m), c(0, 0), col = grey(0.2))
  if (fixtrans) {
    for (nbr in nbrlist) lines(span, xss[nbr, ], col = rgb(0.1, 
                                                           0.1, 0.1, 0.3), lwd = 1)
  }
  else {
    for (nbr in nbrlist) lines(span, xss[nbr, ], col = rgb(0.1, 
                                                           0.1, 0.1, tsc * wts[nbr]), lwd = 3)
  }
}

gw.pcplot.EA(ea.scaled.spdf,  
             vars = colnames(ea.scaled.spdf@data[,1:27]),
             focus=c(13,585),  
             bw = sqrt(bw.gwpca.k8))
