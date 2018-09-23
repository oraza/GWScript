## GWR on Tanzania DHS data set for Childhood morbidity 
## 

library(GWmodel)
library(sp)
library(RColorBrewer) # http://colorbrewer2.org
library(rgdal)
library(leaflet)
library(psych)

# Do not run following code, set project directory first!
setwd("C:\\Users\\Owais\\Dropbox (Personal)\\GW_Thesis\\ThesisData")
dsn <-getwd()
ea.border<-readOGR(dsn = dsn, layer = "ea_Env_DHS")
TZ.data <- ea.border[ea.border$name == "Tanzania", ]
#Calculating Bandwidth
TZ.gwr <- bw.gwr(ARI ~ PoorMid + Rural+
                   FEedu_SecH+MAedu_SecH+Moth_now+
                   GilrChild+agecat6,
                 data = TZ.data, approach = 'AICc',
                 kernel = 'bisquare', adaptive = TRUE)
## GWR 
TZ.gwr.res.ari <- gwr.basic(ARI ~ PoorMid + Rural+
                          FEedu_SecH+MAedu_SecH+Moth_now+
                          GilrChild+agecat6,
                        data = TZ.data, bw = 155, 
                        kernel = 'bisquare', adaptive = TRUE, 
                        F123.test = TRUE)
## GWR output
TZ.gwr.res.ari

# Local collinearity diagnostics for basic GWR
TZ.gwr.collin<- gwr.collin.diagno(ARI ~ PoorMid + Rural+
                                    FEedu_SecH+MAedu_SecH+Moth_now+
                                    GilrChild+agecat6,
                                  data = TZ.data, bw = 155,
                                  kernel="bisquare", adaptive=TRUE)

TZ.adjust.gwt <- gwr.t.adjust(TZ.gwr.res.ari)
summary(TZ.adjust.gwt$results$bh)

n <- leaflet(TZ.data) %>%
  addTiles() %>%
  setView(37.701546, -6.765599, 6) %>%
  addProviderTiles("MapBox", options = providerTileOptions(
    id = "mapbox.light",
    accessToken = Sys.getenv('MAPBOX_ACCESS_TOKEN')))
TZ.adjust.gwt <- gwr.t.adjust(TZ.gwr.res.ari) # produces a list
names(TZ.adjust.gwt$SDF) # check the names of items within it
### poor
TZ.adjust.gwtPoorTtable <- TZ.adjust.gwt$SDF@data 
names(TZ.adjust.gwtPoorTtable)
View(TZ.adjust.gwtPoorTtable)
##### Mapping benf p values for Poor in full range
pal <-  colorNumeric(
  palette = "Reds",
  domain = TZ.adjust.gwtPoorTtable$PoorMid_p
)
n %>%
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.5,
    color = ~pal(TZ.adjust.gwtPoorTtable$PoorMid_p)
  ) %>%
  addLegend("bottomright", pal = pal, values = ~TZ.adjust.gwtPoorTtable$PoorMid_p,
            title = "p values", 
            opacity = 0.5)

# color pallate
mypal.6 = c("#4d4d4d", "#999999", "#e0e0e0", "#fddbc7", "#ef8a62","#b2182b")

## plotting GW Coeff for Poor
bins <- c(quantile(TZ.gwr.res.ari$SDF$PoorMid, probs = seq(0, 1, 0.20), type =  8))
pal <- colorBin(mypal.6, domain = TZ.gwr.res.ari$SDF$PoorMid, bins = bins)
n %>% addPolygons(
  fillColor = ~pal(TZ.gwr.res.ari$SDF$PoorMid),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.6)  %>%
  addLegend("bottomright", 
            title = "GW Coefficients for<br> Percentage of Poor",
            pal = pal, 
            #labels= c("Low", "","","","High"),
            values = TZ.gwr.res.ari$SDF$PoorMid)

#### Rural
TZ.adjust.gwtRuralTtable <- TZ.adjust.gwt$SDF@data 
names(TZ.adjust.gwtRuralTtable)
View(TZ.adjust.gwtRuralTtable)
##### Mapping benf p values for Poor in full range
pal <-  colorNumeric(
  palette = "Reds",
  domain = TZ.adjust.gwtRuralTtable$Rural_p
)
n %>%
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.5,
    color = ~pal(TZ.adjust.gwtRuralTtable$Rural_p)
  ) %>%
  addLegend("bottomright", pal = pal, values = ~TZ.adjust.gwtRuralTtable$Rural_p,
            title = "p values", 
            opacity = 0.5)

# color pallate
mypal.6 = c("#4d4d4d", "#999999", "#e0e0e0", "#fddbc7", "#ef8a62","#b2182b")

## plotting GW Coeff for Rural
bins <- c(quantile(TZ.gwr.res.ari$SDF$Rural, probs = seq(0, 1, 0.20), type =  8))
pal <- colorBin(mypal.6, domain = TZ.gwr.res.ari$SDF$Rural, bins = bins)
n %>% addPolygons(
  fillColor = ~pal(TZ.gwr.res.ari$SDF$Rural),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.6)  %>%
  addLegend("bottomright", 
            title = "GW Coefficients for<br> Percentage of Rural",
            pal = pal, 
            #labels= c("Low", "","","","High"),
            values = TZ.gwr.res.ari$SDF$Rural)
#### FEedu_SecH
TZ.adjust.gwtPoorTtable <- TZ.adjust.gwt$SDF@data 
names(TZ.adjust.gwtPoorTtable)
View(TZ.adjust.gwtPoorTtable)
##### Mapping benf p values for Poor in full range
pal <-  colorNumeric(
  palette = "Reds",
  domain = TZ.adjust.gwtPoorTtable$FEedu_SecH_p
)
n %>%
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.5,
    color = ~pal(TZ.adjust.gwtPoorTtable$FEedu_SecH_p)
  ) %>%
  addLegend("bottomright", pal = pal, values = ~TZ.adjust.gwtPoorTtable$FEedu_SecH_p,
            title = "p values", 
            opacity = 0.5)

# color pallate
mypal.6 = c("#4d4d4d", "#999999", "#e0e0e0", "#fddbc7", "#ef8a62","#b2182b")

## plotting GW Coeff for FEedu_SecH
bins <- c(quantile(TZ.gwr.res.ari$SDF$FEedu_SecH, probs = seq(0, 1, 0.20), type =  8))
pal <- colorBin(mypal.6, domain = TZ.gwr.res.ari$SDF$FEedu_SecH, bins = bins)
n %>% addPolygons(
  fillColor = ~pal(TZ.gwr.res.ari$SDF$FEedu_SecH),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.6)  %>%
  addLegend("bottomright", 
            title = "GW Coefficients for<br> Percentage of Female <br>with Higher Education",
            pal = pal, 
            #labels= c("Low", "","","","High"),
            values = TZ.gwr.res.ari$SDF$FEedu_SecH)
##### Mapping GW Coeff for FEedu_SecH in full range
pal <-  colorNumeric(
  palette = "Reds",
  domain = TZ.gwr.res.ari$SDF$FEedu_SecH
)
n %>%
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.5,
    color = ~pal(TZ.gwr.res.ari$SDF$FEedu_SecH)
  ) %>%
  addLegend("bottomright", pal = pal, values = ~TZ.gwr.res.ari$SDF$FEedu_SecH,
            title = "p values", 
            opacity = 0.5)

#### GilrChild
TZ.adjust.gwtPoorTtable <- TZ.adjust.gwt$SDF@data 
names(TZ.adjust.gwtPoorTtable)
View(TZ.adjust.gwtPoorTtable)
##### Mapping benf p values for Poor in full range
pal <-  colorNumeric(
  palette = "Reds",
  domain = TZ.adjust.gwtPoorTtable$GilrChild_p
)
n %>%
  addPolygons(
    stroke = FALSE, smoothFactor = 0.2, fillOpacity = 0.5,
    color = ~pal(TZ.adjust.gwtPoorTtable$GilrChild_p)
  ) %>%
  addLegend("bottomright", pal = pal, values = ~TZ.adjust.gwtPoorTtable$GilrChild_p,
            title = "p values", 
            opacity = 0.5)

# color pallate
mypal.6 = c("#4d4d4d", "#999999", "#e0e0e0", "#fddbc7", "#ef8a62","#b2182b")

## plotting GW Coeff for GilrChild
bins <- c(quantile(TZ.gwr.res.ari$SDF$GilrChild, probs = seq(0, 1, 0.20), type =  8))
pal <- colorBin(mypal.6, domain = TZ.gwr.res.ari$SDF$GilrChild, bins = bins)
n %>% addPolygons(
  fillColor = ~pal(TZ.gwr.res.ari$SDF$GilrChild),
  weight = 1,
  opacity = 0.3,
  color = 'grey75',
  fillOpacity = 0.6)  %>%
  addLegend("bottomright", 
            title = "GW Coefficients for<br> Percentage of Girl Child",
            pal = pal, 
            #labels= c("Low", "","","","High"),
            values = TZ.gwr.res.ari$SDF$GilrChild)
