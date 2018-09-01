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

