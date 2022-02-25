prodir=getwd()

url="https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/Lab04.R"
download.file(url,"Lab04.R")
file.edit("Lab04.R")
options(repos ="http://cran.us.r-project.org")  # required to get latest libs

if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,ggplot2,patchwork)
pacman::p_load(EcoHydRology,rnoaa,curl,httr)


source("https://raw.githubusercontent.com/Miladvahidi44/BSE5304/main/functions.R")


## Question 1 .........................

url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/bky5q4i2wkdjmvm3t0d1gniq/wss_aoi_2022-02-11_09-31-24.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")



myflowgage_id="01421618"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
mysoil=readOGR("wss_aoi_2022-02-24_20-45-37/spatial/soilmu_a_aoi.shp")    
mybbox=c(mysoil@bbox)

mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)
plot(mydem)

summary(terrain(mydem, opt='slope',unit = "degrees"))
# What is this 'slope'? Use the man page for the terrain() function to answer
plot(terrain(mydem, opt='TPI',unit = "degrees"))
# What is this 'TPI'? 
summary(terrain(mydem, opt='TRI',unit = "degrees"))
plot(terrain(mydem, opt='TRI',unit = "degrees"))
# What is this 'TRI'? 

stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)


WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
  modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)]=0 # A Quick BUT sloppy removal of NAs
TMWB=modeldata


TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata

#############################################
 summary(mu2chmax)

 awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
 TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
 MidslopeZ=mean(mu2chmax$hzdepb_r)*10
 BotslopeZ=max(mu2chmax$hzdepb_r)*10

AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent

summary(terrain(mydem, opt='slope',unit = "degrees"))
 SlopeTop=1.837328 #degree
 SlopeBot=1.837328 #degree
 SlopeMid=40.600004 #degree


 TopSlope = TMWBmodel(TMWB = TopSlope,SFTmp = 1, 
                       AWCval = AWCvalTop,
                       Tlag = .5,fcres=.3,Slope = atan(SlopeTop/100))
 MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
 MidSlope = TMWBmodel(TMWB = MidSlope,SFTmp = 1, 
                       AWCval = AWCvalMid,
                       Tlag = .5,fcres=0.5,Slope = atan(SlopeMid/100))
# Low Slope and lowest ksat, $fcres=0.2
 BotSlope$P=MidSlope$Excess+BotSlope$P
 BotSlope = TMWBmodel(TMWB = BotSlope,SFTmp = 1, 
                       AWCval = AWCvalBot,
                       Tlag = .5,fcres=0.2,Slope = atan(SlopeBot/100))
##############AW Plots
 p1=ggplot() +
  geom_line(data=BotSlope,aes(x=date, y = AW,colour="BotSlope")) +
  geom_line(data=MidSlope,aes(x=date, y = AW,colour="MidSlope")) +
  geom_line(data=TopSlope,aes(x=date, y = AW,colour="TopSlope")) +
  labs(x = 'Date', y = 'AW (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlope", "MidSlope", "TopSlope"),
                      values = c("black", "blue","red"))+
   theme(legend.position = "bottom",
         legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
   ggtitle("(a), TMWB")
##############Excess
 p2=ggplot() +
  geom_line(data=BotSlope,aes(x=date, y = Excess,colour="BotSlope")) +
  geom_line(data=MidSlope,aes(x=date, y = Excess,colour="MidSlope")) +
  geom_line(data=TopSlope,aes(x=date, y = Excess,colour="TopSlope")) +
  labs(x = 'Date', y = 'Excess (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlope", "MidSlope", "TopSlope"),
                      values = c("black", "blue","red"))+
   theme(legend.position = "bottom",
         legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+  ggtitle("(b) TMWB")
 p1 + p2 + plot_layout(ncol = 1, widths = c(1, 1))


 ggplot() +
   geom_line(data=BotSlope,aes(x=date, y = Qpred )) +
   labs(x = 'Date', y = 'Flow Q (mm)')+
   theme(legend.position = "bottom",
         legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
   ggtitle("(TMWB)")
 




############################
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata


TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)

MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)
#AW Plots HW1

p1=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = AW,colour="BotSlopeCN")) +
  geom_line(data=MidSlopeCN,aes(x=date, y = AW,colour="MidSlopeCN")) +
  geom_line(data=TopSlopeCN,aes(x=date, y = AW,colour="TopSlopeCN")) +
  labs(x = 'Date', y = 'AW (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlopeCN", "MidSlopeCN", "TopSlopeCN"),
                      values = c("black", "blue","red"))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(a)")


# Excess Plots HW1

p2=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = Excess,colour="BotSlopeCN")) +
  geom_line(data=MidSlopeCN,aes(x=date, y = Excess,colour="MidSlopeCN")) +
  geom_line(data=TopSlopeCN,aes(x=date, y = Excess,colour="TopSlopeCN")) +
  labs(x = 'Date', y = 'Excess (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlopeCN", "MidSlopeCN", "TopSlopeCN"),
                      values = c("black", "blue","red"))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(b)")

p1 + p2 + plot_layout(ncol = 1, widths = c(1, 1))


# PET and ET HW2

p3=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = PET,colour="BotSlopeCN PET")) +
  geom_line(data=BotSlopeCN,aes(x=date, y = ET,colour="BotSlopeCN")) +
  geom_line(data=MidSlopeCN,aes(x=date, y = ET,colour="MidSlopeCN")) +
  geom_line(data=TopSlopeCN,aes(x=date, y = ET,colour="TopSlopeCN")) +
  
  labs(x = 'Date', y = '(P)ET (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlopeCN PET","BotSlopeCN", "MidSlopeCN", "TopSlopeCN"),
                      values = c("green","black", "blue","red"))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(c)")


p4=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = cumsum(PET),colour="BotSlopeCN PET")) +
  geom_line(data=BotSlopeCN,aes(x=date, y = cumsum(ET),colour="BotSlopeCN")) +
  geom_line(data=MidSlopeCN,aes(x=date, y = cumsum(ET),colour="MidSlopeCN")) +
  geom_line(data=TopSlopeCN,aes(x=date, y = cumsum(ET),colour="TopSlopeCN")) +
  
  labs(x = 'Date', y = '(P)ET (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlopeCN PET","BotSlopeCN", "MidSlopeCN", "TopSlopeCN"),
                      values = c("green","black", "blue","red"))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(d)")

p3 + p4 + plot_layout(ncol = 1, widths = c(1, 1))



# Cumulative Summary of QPred is very informative
p5=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = cumsum(Qpred) ,colour="BotSlopeCN")) +
  geom_line(data=MidSlopeCN,aes(x=date, y = cumsum(Qpred),colour="MidSlopeCN")) +
  geom_line(data=TopSlopeCN,aes(x=date, y = cumsum(Qpred),colour="TopSlopeCN")) +
  labs(x = 'Date', y = 'Flow Q Cumulative Summary (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlopeCN", "MidSlopeCN", "TopSlopeCN"),
                      values = c("black", "blue","red"))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(E)")


# Model Performance 
p5=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = Qpred )) +
  labs(x = 'Date', y = 'Flow Q (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("CN")


NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)
NSeff(BotSlope$Qmm,BotSlope$Qpred)




## Required packages ............................
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster)
pacman::p_load(EcoHydRology,rnoaa,curl,httr)


## Downloading streamflow dataset from USGS website ................
myflowgage_id="02038850"

myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2017-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3


## Downloading the topography and soils for the area ..................
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_02080207_HU8_Shape.zip"
download.file(url,"NHD_H_02080207_Shape.zip")
unzip("NHD_H_02080207_Shape.zip",exdir="02080207")
# Take a quick look at what is included in the NHD dataset

## next 
streams=readOGR("02080207/Shape/NHDFlowline.dbf")

plot(streams)
mystream=subset(streams,gnis_name=="Appomattox River")
lines(mystream,col="red")


mystream=subset(streams,gnis_id=="01478950")
plot(mystream,col="red")

c(mystream@bbox)
# What is this returning? Why do we care?
mybbox=c(mystream@bbox)
# This needs to be completed based on your download
mysoil=readOGR("wss_aoi_2022-02-14_12-11-49/spatial/soilmu_a_aoi.shp")

## .........
# First associate mukey with cokey from component
unique(mysoil$MUKEY)
mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
head(mu2co)


cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
View(mu2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)
summary(mu2chmax)


# downloading elevation data .........
install.packages("rasterVis")
library(rasterVis)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.1)

gplot(mydem) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(terrain.colors(225))) +
  coord_equal() + theme_void() + theme(legend.position = "none")+
  theme_set(theme_bw())




mydem@crs
mydem@ncols
res(mydem)   # Can you figure out the resolution in meters? 
plot(terrain(mydem, opt='TPI',unit = "degrees"))
# What is this 'TPI'? 
plot(terrain(mydem, opt='TRI',unit = "degrees"))
# What is this 'TRI'? 
lines(mysoil,col="black")
lines(mystream,col="red")


## ............
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)



WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
Data = data_frame("Date"=WXData$date ,
                  "PRCP"=WXData$prcp,
                  "TMAX"=WXData$tmax,
                  "TMIN"=WXData$tmin)


modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
mean(modeldata$Qmm)
mean(modeldata$P)
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MaxTemp)/2.0

modeldata[is.na(modeldata)]=0 # A Quick BUT sloppy removal of NAs



summary(modeldata)
TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata

#############################################
summary(mu2chmax)
summary(terrain(mydem, opt='slope',unit = "degrees"))

awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10

AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent

summary(terrain(mydem, opt='slope',unit = "degrees"))
SlopeTop=2.877927 #degree
SlopeBot=2.877927 #degree
SlopeMid=28.856836 #degree

TMWB=modeldata

TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata


TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN,IaFrac = 0.02, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=0.1,
                     func_z=500,fnc_fcres=0.15)

MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN,IaFrac = 0.02, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=0.1,
                     func_z=750,fnc_fcres=0.15)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN,IaFrac = 0.02, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=0.1,
                     func_z=1000,fnc_fcres=0.15)

NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)



## q4



url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/x2qg1plc5ik5fu44pft4w1pu/wss_aoi_2022-02-24_20-45-37.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")


TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata

#############################################
summary(mu2chmax)

awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10

AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent

summary(terrain(mydem, opt='slope',unit = "degrees"))
SlopeTop= 2.067430 #degree
SlopeBot= 2.067430 #degree
SlopeMid= 5.994341 #degree


TopSlope = TMWBmodel(TMWB = TopSlope,SFTmp = 1, 
                     AWCval = AWCvalTop,
                     Tlag = .5,fcres=.3,Slope = atan(SlopeTop/100))
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlope = TMWBmodel(TMWB = MidSlope,SFTmp = 1, 
                     AWCval = AWCvalMid,
                     Tlag = .5,fcres=0.5,Slope = atan(SlopeMid/100))
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlope = TMWBmodel(TMWB = BotSlope,SFTmp = 1, 
                     AWCval = AWCvalBot,
                     Tlag = .5,fcres=0.2,Slope = atan(SlopeBot/100))
NSeff(BotSlope$Qmm,BotSlope$Qpred)





TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata


TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)

MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)



NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)
NSeff(BotSlope$Qmm,BotSlope$Qpred)


p1=ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y = Qpred )) +
  labs(x = 'Date', y = 'Flow Q (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(CN) LICK RUN")




p2=ggplot() +
  geom_line(data=BotSlope,aes(x=date, y = Qpred )) +
  labs(x = 'Date', y = 'Flow Q (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TMWB) LICK RUN")


p1 + p2 + plot_layout(ncol = 1, widths = c(1, 1))


