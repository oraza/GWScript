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
bart(ea.border@data[,16:41])

# scale
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
