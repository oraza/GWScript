## GWPCA on East Africa data set for Childhood morbidity 

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

#### Global PCA #####
# Bartlett's test of sphericity
bart<-function(dat){ #dat is your raw data
  R<-cor(dat)
  p<-ncol(dat)
  n<-nrow(dat)
  chi2<- -((n-1)-((2*p)+5)/6 ) * log(det(R)) #this is the formula
  df<-(p*(p-1)/2)
  crit<-qchisq(.95,df) #critical value
  p<-pchisq(chi2,df,lower.tail=F) #pvalue
  cat("Bartlett's test of sphericity: X2(",
      df,")=",chi2,", p=",
      round(p,3),sep="" )   
}
# applying Bartlett's test for all potential confounders
bart(ea.border@data[,16:41]) 

# scaling all explanatory variables including main
ea.df.sc = scale(as.matrix(ea.border@data[,15:41])) 

# global PCA using psych package
ea.pca2 <- principal(as.matrix(ea.border@data[,16:41]), 
                     nfactors=6, rotate="varimax")

# inserting Global PCs into main SpatialPolygonsDataFrame
ea.border$global.PC1 <- ea.pca2$scores[,1]
ea.border$global.PC2 <- ea.pca2$scores[,2]
ea.border$global.PC5 <- ea.pca2$scores[,5]
ea.border$global.PC3 <- ea.pca2$scores[,3]
ea.border$global.PC4 <- ea.pca2$scores[,4]
ea.border$global.PC6 <- ea.pca2$scores[,6]

# standardizing th PC Score (dividing by their SDs)
ea.border@data$st.global.PC1 <- ea.border@data$global.PC1/(sd(ea.border@data$global.PC1))
ea.border@data$st.global.PC2 <- ea.border@data$global.PC2/(sd(ea.border@data$global.PC2))
ea.border@data$st.global.PC3 <- ea.border@data$global.PC3/(sd(ea.border@data$global.PC3))
ea.border@data$st.global.PC4 <- ea.border@data$global.PC4/(sd(ea.border@data$global.PC4))
ea.border@data$st.global.PC5 <- ea.border@data$global.PC5/(sd(ea.border@data$global.PC5))
ea.border@data$st.global.PC6 <- ea.border@data$global.PC6/(sd(ea.border@data$global.PC6))

# Scree plot
scree(ea.pca2$r.scores, 
      factors = F, 
      pc = T, 
      hline = -9,
      main = "Scree Plot for global PCA")

# Renaming variable
names(ea.border@data)[names(ea.border@data) == "PoorMid"] <- "Poor"
names(ea.border@data)[names(ea.border@data) == "childmorb"] <- "ChildMorbid"


## Summary of Geographical distribution of variables
round(summary(ea.border$ChildMorbid),2)
round(summary(ea.border$Poor),2)
round(summary(ea.border$Rural),2)
round(summary(ea.border$Solidfuel),2)
round(summary(ea.border$HealthAcc),2)
round(summary(ea.border$singlemom),2)
round(summary(ea.border$FEedu_NilP),2)
round(summary(ea.border$FEedu_SecH),2)
round(summary(ea.border$MAedu_NilP),2)
round(summary(ea.border$MAedu_SecH),2)
round(summary(ea.border$Moth_now),2)
round(summary(ea.border$Moth_prof),2)
round(summary(ea.border$Moth_agri),2)
round(summary(ea.border$GilrChild),2)
round(summary(ea.border$NonImprWat),2)
round(summary(ea.border$WateratHom),2)
round(summary(ea.border$Water30mAw),2)
round(summary(ea.border$CompImmu),2)
round(summary(ea.border$Aridity),2)
round(summary(ea.border$MaxTemp),2)
round(summary(ea.border$PM2_5),2)
round(summary(ea.border$DroughtHz),2)
round(summary(ea.border$agecat1),2)
round(summary(ea.border$agecat2),2)
round(summary(ea.border$agecat3),2)
round(summary(ea.border$agecat4),2)
round(summary(ea.border$agecat5),2)
round(summary(ea.border$agecat6),2)

round(sd(ea.border$ChildMorbid),2)
round(sd(ea.border$Poor),2)
round(sd(ea.border$Rural),2)
round(sd(ea.border$Solidfuel),2)
round(sd(ea.border$HealthAcc),2)
round(sd(ea.border$singlemom),2)
round(sd(ea.border$FEedu_NilP),2)
round(sd(ea.border$FEedu_SecH),2)
round(sd(ea.border$MAedu_NilP),2)
round(sd(ea.border$MAedu_SecH),2)
round(sd(ea.border$Moth_now),2)
round(sd(ea.border$Moth_prof),2)
round(sd(ea.border$Moth_agri),2)
round(sd(ea.border$GilrChild),2)
round(sd(ea.border$NonImprWat),2)
round(sd(ea.border$WateratHom),2)
round(sd(ea.border$Water30mAw),2)
round(sd(ea.border$CompImmu),2)
round(sd(ea.border$Aridity),2)
round(sd(ea.border$MaxTemp),2)
round(sd(ea.border$PM2_5),2)
round(sd(ea.border$DroughtHz),2)
round(sd(ea.border$agecat1),2)
round(sd(ea.border$agecat2),2)
round(sd(ea.border$agecat3),2)
round(sd(ea.border$agecat4),2)
round(sd(ea.border$agecat5),2)
round(sd(ea.border$agecat6),2)

#### GWPCA #####
# assigning coordinates
Coords=as.matrix(cbind(ea.border$xlong,
                       ea.border$ylat))

# merging coordinates with data set and converting it into shapefile
ea.scaled.spdf=SpatialPointsDataFrame(Coords,as.data.frame(ea.df.sc))
proj4string(ea.scaled.spdf) <- CRS("+init=epsg:4326")

##################### Base map ##########################
m <- leaflet(ea.border) %>%
  addTiles() %>%
  setView(37.701546, -6.765599, 4) %>%
  addProviderTiles("MapBox", options = providerTileOptions(
    id = "mapbox.light",
    accessToken = Sys.getenv('MAPBOX_ACCESS_TOKEN')))

#### Descriptive Results
# prevalence of diarrhea

ea.border@data$DiaProp <- as.numeric(ea.border@data$Diarrhea*100)
pal <- colorNumeric(palette = "Reds",
                    domain = ea.border@data$DiaProp)

m %>% addPolygons(
  stroke = F, smoothFactor = 0.2,
  fillOpacity = 0.7,
  color = ~pal(ea.border@data$DiaProp))  %>%
  addLegend("bottomright", 
            title = "Prevalence of Diarrhea <br> in East Africa (%)",
            pal = pal, 
            values = ea.border@data$DiaProp)

# prevalence of ARI
ea.border@data$ARIProp <- as.numeric(ea.border@data$ARI*100)
pal <- colorNumeric(palette = "Reds",
                    domain = ea.border@data$ARIProp)

m %>% addPolygons(
  stroke = F, smoothFactor = 0.2,
  fillOpacity = 0.7,
  color = ~pal(ea.border@data$ARIProp))  %>%
  addLegend("bottomright", 
            title = "Prevalence of ARI <br> in East Africa (%)",
            pal = pal, 
            values = ea.border@data$ARIProp)

# prevalence of Fever
ea.border@data$FeverProp <- as.numeric(ea.border@data$fever*100)
pal <- colorNumeric(palette = "Reds",
                    domain = ea.border@data$FeverProp)

m %>% addPolygons(
  stroke = F, smoothFactor = 0.2,
  fillOpacity = 0.7,
  color = ~pal(ea.border@data$FeverProp))  %>%
  addLegend("bottomright", 
            title = "Prevalence of fever <br> in East Africa (%)",
            pal = pal, 
            values = ea.border@data$FeverProp)

# prevalence of ChildMorb
ea.border@data$ChilMorbProp <- as.numeric(ea.border@data$ChildMorbid *100)
pal <- colorNumeric(palette = "Reds",
                    domain = ea.border@data$ChilMorbProp)

m %>% addPolygons(
  stroke = F, smoothFactor = 0.2,
  fillOpacity = 0.7,
  color = ~pal(ea.border@data$ChilMorbProp))  %>%
  addLegend("topright", 
            title = "Prevalence of childhood <br> morbidity in East Africa (%)",
            pal = pal, 
            values = ea.border@data$ChilMorbProp)

# prevalence of Poor
ea.border@data$PoorProp <- as.numeric(ea.border@data$Poor *100)
pal <- colorNumeric(palette = "Blues",
                    domain = ea.border@data$PoorProp)

m %>% addPolygons(
  stroke = F, smoothFactor = 0.2,
  fillOpacity = 0.7,
  color = ~pal(ea.border@data$PoorProp))  %>%
  addLegend("topright", 
            title = "Prevalence of poor <br> in East Africa (%)",
            pal = pal, 
            values = ea.border@data$PoorProp)

## calculting GW bandwidth to find an optimal adaptive bandwidth
## using a bi-square kernel 
bw.gwpca.k8=bw.gwpca(ea.scaled.spdf, 
                     vars=colnames(ea.scaled.spdf@data[,1:27]), 
                     k = 8,
                     robust = FALSE,
                     adaptive = TRUE,
                     kernel = "bisquare")

prop.var=function(gwpca.obj, n.components){
  return((rowSums(gwpca.obj$var[,1:n.components])/
            rowSums(gwpca.obj$var))*100)
}

# PTV within GWPC1
ea.ptv1<-(gwpca.k8$var[,1:1])/rowSums(gwpca.k8$var)*100
ea.border$ea.ptv1.k8=ea.ptv1
summary(ea.border$ea.ptv1.k8)
sd(ea.border$ea.ptv1.k8)

# mapping PTV within GWPC1
mypal.6 <- c('#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026')
bins <- c(quantile(ea.border$ea.ptv1.k8, probs = seq(0,1, 0.2),
                   type = 8))
pal <- colorBin(mypal.6, domain = ea.border$ea.ptv1.k8, bins = bins)

m %>% addPolygons(
  fillColor = ~pal(ea.border$ea.ptv1.k8),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.5)  %>%
  addLegend("bottomright", 
            title = "Percentage of total variance <br>within first component",
            pal = pal, 
            values = ea.border$ea.ptv1.k8)

# PTV within first GWPC2
ea.ptv2=prop.var(gwpca.k8, 2)
ea.border$ea.ptv2.k8=ea.ptv2
summary(ea.border$ea.ptv2.k8)
sd(ea.border$ea.ptv2.k8)

# mapping PTV within GWPC2
mypal.6 <- c('#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026')
bins <- c(quantile(ea.border$ea.ptv2.k8, probs = seq(0,1, 0.2),
                   type = 8))
pal <- colorBin(mypal.6, domain = ea.border$ea.ptv2.k8, bins = bins)

m %>% addPolygons(
  fillColor = ~pal(ea.border$ea.ptv2.k8),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.5)  %>%
  addLegend("bottomright", 
            title = "Percentage of total variance <br>within first two components",
            pal = pal, 
            values = ea.border$ea.ptv2.k8)

# PTV within first 3 comp
ea.ptv3=prop.var(gwpca.k8, 3)
ea.border$ea.ptv3.k8=ea.ptv3
summary(ea.border$ea.ptv3.k8)
sd(ea.border$ea.ptv3.k8)

# mapping PTV within GWPC3
mypal.6 <- c('#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026')
bins <- c(quantile(ea.border$ea.ptv3.k8, probs = seq(0,1, 0.2),
                   type = 8))
pal <- colorBin(mypal.6, domain = ea.border$ea.ptv3.k8, bins = bins)

m %>% addPolygons(
  fillColor = ~pal(ea.border$ea.ptv3.k8),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.5)  %>%
  addLegend("bottomright", 
            title = "Percentage of total variance <br>within first three components",
            pal = pal, 
            values = ea.border$ea.ptv3.k8)

# PTV within first 4 comp
ea.ptv4=prop.var(gwpca.k8, 4)
ea.border$ea.ptv4.k8=ea.ptv4
summary(ea.border$ea.ptv4.k8)
sd(ea.border$ea.ptv4.k8)

# mapping PTV within GWPC4
mypal.6 <- c('#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026')
bins <- c(quantile(ea.border$ea.ptv4.k8, probs = seq(0,1, 0.2),
                   type = 8))
pal <- colorBin(mypal.6, domain = ea.border$ea.ptv4.k8, bins = bins)

m %>% addPolygons(
  fillColor = ~pal(ea.border$ea.ptv4.k8),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.5)  %>%
  addLegend("bottomright", 
            title = "Percentage of total variance <br>within first four components",
            pal = pal, 
            values = ea.border$ea.ptv4.k8)

# PTV within first 5 comp
ea.ptv5=prop.var(gwpca.k8, 5)
ea.border$ea.ptv5.k8=ea.ptv5
summary(ea.border$ea.ptv5.k8)
sd(ea.border$ea.ptv5.k8)

# mapping PTV within GWPC5
mypal.6 <- c('#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026')
bins <- c(quantile(ea.border$ea.ptv5.k8, probs = seq(0,1, 0.2),
                   type = 8))
pal <- colorBin(mypal.6, domain = ea.border$ea.ptv5.k8, bins = bins)

m %>% addPolygons(
  fillColor = ~pal(ea.border$ea.ptv5.k8),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.5)  %>%
  addLegend("bottomright", 
            title = "Percentage of total variance <br>within first five components",
            pal = pal, 
            values = ea.border$ea.ptv5.k8)

# PTV within first 6 comp
ea.ptv6=prop.var(gwpca.k8, 6)
ea.border$ea.ptv6.k8=ea.ptv6
summary(ea.border$ea.ptv6.k8)
sd(ea.border$ea.ptv6.k8)

# mapping PTV within GWPC6
mypal.6 <- c('#ffffb2', '#fecc5c', '#fd8d3c', '#f03b20', '#bd0026')
bins <- c(quantile(ea.border$ea.ptv6.k8, probs = seq(0,1, 0.2),
                   type = 8))
pal <- colorBin(mypal.6, domain = ea.border$ea.ptv6.k8, bins = bins)

m %>% addPolygons(
  fillColor = ~pal(ea.border$ea.ptv6.k8),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.5)  %>%
  addLegend("bottomright", 
            title = "Percentage of total variance <br>within first six components",
            pal = pal, 
            values = ea.border$ea.ptv6.k8)

# extracting most influential variables in GWPC1
loadings.pc1.k8=gwpca.k8$loadings[,,1] # View(loadings.pc1.k8)
PC1.WinVar.k8=max.col(abs(loadings.pc1.k8)) # summary(PC1.WinVar.k8) ?Is it a position of the kth variable? 
ea.border$PC1.WinVar.k8=PC1.WinVar.k8
table(PC1.WinVar.k8)

# mapping
ea.border$PC1.WinVar.k8.cat <- as.factor(ea.border$PC1.WinVar.k8)
levels(ea.border$PC1.WinVar.k8.cat) <-c("Poor", "Rural", "Solid Fuel",
                                        "Healthcare Access",
                                        "Mother Edu < Secondary",
                                        "Mother Edu >= Secondary", 
                                        "Father Edu < Secondary", 
                                        "Father Edu >= Secondary",
                                        "Mother Occ - Prof",
                                        "Unimproved Drinking Water",
                                        "Water at Home")

factpal <- colorFactor(mypal.11, ea.border$PC1.WinVar.k8.cat)

m %>% 
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.5,
    color = ~factpal(ea.border$PC1.WinVar.k8.cat)
  ) %>%
  addLegend("bottomright", 
            title = "Variables with max loadings <br> within 1st component",
            pal = factpal, 
            values = ea.border$PC1.WinVar.k8.cat)

# extracting most influential variables in GWPC2
loadings.pc2.k8=gwpca.k8$loadings[,,2]
PC2.WinVar.k8=max.col(abs(loadings.pc2.k8))
ea.border$PC2.WinVar.k8=PC2.WinVar.k8
table(PC2.WinVar.k8)

# mapping
ea.border$PC2.WinVar.k8.cat <- as.factor(ea.border$PC2.WinVar.k8)
levels(ea.border$PC2.WinVar.k8.cat) <-c("Poor", "Rural", "Solid Fuel",
                                        "Single Mother",
                                        "Mother Occ - None",
                                        "Girl Child",
                                        "Unimproved Drinking Water", 
                                        "Complete Immunization", "Aridity",
                                        "Max Temperature","PM2.5","Drought Hazard",
                                        "6-11 months", "12-23 months",
                                        "24-35 months", "36-47 months",
                                        "48-59 months")

factpal <- colorFactor(mypal.17, ea.border$PC2.WinVar.k8.cat)
m %>% 
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.7,
    color = ~factpal(ea.border$PC2.WinVar.k8.cat)
  ) %>%
  addLegend("bottomright", 
            title = "Variables with max loadings <br> within 2nd component",
            pal = factpal, 
            values = ea.border$PC2.WinVar.k8.cat)

# extracting most influential variables in GWPC3
loadings.pc3.k8=gwpca.k8$loadings[,,3]
PC3.WinVar.k8=max.col(abs(loadings.pc3.k8))
ea.border$PC3.WinVar.k8=PC3.WinVar.k8
table(PC3.WinVar.k8)

# mapping
ea.border$PC3.WinVar.k8.cat <- as.factor(ea.border$PC3.WinVar.k8)
levels(ea.border$PC3.WinVar.k8.cat) <-c("Solid Fuel","Father Edu < Secondary",
                                        "Mother Occ - None", "Girl Child",
                                        "Complete Immunization", 
                                        "Maximum Temperature", "PM2.5",
                                        "< 6 months",
                                        "6-11 months",
                                        "12-23 months",
                                        "24-35 months",
                                        "36-47 months",
                                        "48-59 months")

factpal <- colorFactor(mypal.17, ea.border$PC3.WinVar.k8.cat)

m %>% 
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.7,
    color = ~factpal(ea.border$PC3.WinVar.k8.cat)
  ) %>%
  addLegend("bottomright", 
            title = "Variables with max loadings <br> within 3rd component",
            pal = factpal, 
            values = ea.border$PC3.WinVar.k8.cat)

# extracting most influential variables in GWPC4
loadings.pc4.k8=gwpca.k8$loadings[,,4]
PC4.WinVar.k8=max.col(abs(loadings.pc4.k8))
ea.border$PC4.WinVar.k8=PC4.WinVar.k8
table(PC4.WinVar.k8)

# mapping
ea.border$PC4.WinVar.k8.cat <- as.factor(ea.border$PC4.WinVar.k8)
levels(ea.border$PC4.WinVar.k8.cat) <-c("Solid Fuel","Healthcare Access",
                                        "Single Mother",
                                        "Mother Occ - None",
                                        "Mother Occ - Prof",
                                        "Girl Child",
                                        "Unimproved Drinking Water",
                                        "Complete Immunization", 
                                        "Maximum Temperature",
                                        "< 6 months",
                                        "6-11 months",
                                        "12-23 months",
                                        "24-35 months",
                                        "36-47 months",
                                        "48-59 months")

factpal <- colorFactor(mypal.15, ea.border$PC4.WinVar.k8.cat)

m %>% 
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.7,
    color = ~factpal(ea.border$PC4.WinVar.k8.cat)
  ) %>%
  addLegend("bottomright", 
            title = "Variables with max loadings <br> within 4th component",
            pal = factpal, 
            values = ea.border$PC4.WinVar.k8.cat)

# extracting most influential variables in GWPC5
loadings.pc5.k8=gwpca.k8$loadings[,,5]
PC5.WinVar.k8=max.col(abs(loadings.pc5.k8))
ea.border$PC5.WinVar.k8=PC5.WinVar.k8
table(PC5.WinVar.k8)

# mapping
ea.border$PC5.WinVar.k8.cat <- as.factor(ea.border$PC5.WinVar.k8)
levels(ea.border$PC5.WinVar.k8.cat) <-c("Poor","Solid Fuel", 
                                        "Healthcare Access",
                                        "Single Mother",
                                        "Mother Edu < Secondary",
                                        "Mother Occ - None",
                                        "Mother Occ - Prof",
                                        "Girl Child",
                                        "Unimproved Drinking Water",
                                        "Complete Immunization",
                                        "Maximum Temperature",
                                        "< 6 months",
                                        "6-11 months",
                                        "12-23 months",
                                        "24-35 months",
                                        "36-47 months",
                                        "48-59 months")

factpal <- colorFactor(mypal.17, ea.border$PC5.WinVar.k8.cat)

m %>% 
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.7,
    color = ~factpal(ea.border$PC5.WinVar.k8.cat)
  ) %>%
  addLegend("bottomright", 
            title = "Variables with max loadings <br> within 5th component",
            pal = factpal, 
            values = ea.border$PC5.WinVar.k8.cat)

# extracting most influential variables in GWPC6
loadings.pc6.k8=gwpca.k8$loadings[,,6]
PC6.WinVar.k8=max.col(abs(loadings.pc6.k8))
ea.border$PC6.WinVar.k8=PC6.WinVar.k8
table(PC6.WinVar.k8)

# mapping
ea.border$PC6.WinVar.k8.cat <- as.factor(ea.border$PC6.WinVar.k8)
levels(ea.border$PC6.WinVar.k8.cat) <-c("Solid Fuel", 
                                        "Healthcare Access",
                                        "Single Mother",
                                        "Mother Occ - None",
                                        "Mother Occ - Prof",
                                        "Girl Child",
                                        "Unimproved Drinking Water",
                                        "Complete Immunization",
                                        "Maximum Temperature",
                                        "< 6 months",
                                        "6-11 months",
                                        "12-23 months",
                                        "24-35 months",
                                        "36-47 months",
                                        "48-59 months")

factpal <- colorFactor(mypal.17, ea.border$PC6.WinVar.k8.cat)

m %>% 
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.7,
    color = ~factpal(ea.border$PC6.WinVar.k8.cat)
  ) %>%
  addLegend("bottomright", 
            title = "Variables with max loadings <br> within 6th component",
            pal = factpal, 
            values = ea.border$PC6.WinVar.k8.cat)

## Monte Carlo Test for first component (prespecified BW)
## Simulation (n=1000) with Monte Carlo Test for 
## first component (auto-calibirated BW )
## estimated time ~ 15 hours

gwpca.mc27.1000 <- montecarlo.gwpca.2(ea.scaled.spdf, 
                                      vars = colnames(ea.scaled.spdf@data[,1:27]),
                                      k = 27, nsims=1000,
                                      adaptive = TRUE)
with(ea.scaled.spdf@data[,1:27],{
  plot(density(gwpca.mc10.1000$sims),
       main = "Test statistics for eigenvalue nonstationarity",
       xlab = "SD of local eigenvalues from randomisations \n N = 1000")
  abline(v=2.160, col = "red", lwd = 2)
  text(2.25, 2.5, paste("Observed SD of local eigenvalues \n p value = 0.0160"), col = "red",
       srt = 90)
  abline(v = 0.8164263, col = "darkblue", lwd = 2, lty = 5)
  text(0.9, 2.5, paste("Mean of SD of local \n eigenvalues from randomisations"),
       col = "darkblue", srt = 90)
  rug(gwpca.mc10.1000$sims, col = "orange")
})
abline(h=-0.1, v=2.15900, col = "red", lwd = 2)


