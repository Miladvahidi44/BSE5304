### Lab 4 ...............
rm(list = ls())


url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/cybakmopm0fhsg44mzzwzkgr/wss_aoi_2022-02-14_12-11-49.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
list.files()

list.files(pattern = "wss")

 
list.files("wss_aoi_2022-02-14_12-11-49.zip/",pattern="shp")
list.files(paste0(list.files(pattern="wss"),"/spatial/"))
#Note you will need to replace the ?? with the values from your soils directory
list.files("wss_aoi_2022-02-14_12-11-49/tabular/")

  
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
list.files("02080207/Shape/",pattern = "dbf")
    
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
# Explore the mysoil dataset which is returned
 head(mysoil@data)
 class(mysoil)
mysoil@data
plot(mysoil)
  
bg = "grey92"
AG <- fortify(mysoil)
ggplot() + 
  geom_polygon(data = AG, aes(long, lat, group = group, fill = hole), 
colour = alpha("green3", 1/2), size = 0.7) + 
  scale_fill_manual(values = c("white", bg)) + 
  theme(panel.background = element_rect(fill = bg),
   legend.position = "none")


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
 
 
 
 
 plot(mydem)
 lines(mysoil,col="black")
 lines(mystream,col="red")
 class(mydem)
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

ggplot(Data) +
  geom_line(aes(x = Date , y= PRCP))+
  ylab("")+
  ggtitle('HOLIDAY CREEK NEAR ANDERSONVILLE, VA')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0.5,size=10, angle = 45 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank()) 
        




## Question 1 ...............................


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
TMWB=modeldata
soildrying<-function(AWprev,dP,AWC){
  AW<-AWprev*exp(dP/AWC)
  excess<-0.0
  c(AW,excess)
}

soil_wetting_above_capacity<-function(AWprev,dP,AWC){
  AW<-AWC
  excess<-AWprev+dP-AWC
  c(AW,excess)
}

soilwetting<-function(AWprev,dP,AWC){
  AW<-AWprev+dP
  excess<-0.0
  c(AW,excess)
}         
         

TMWBmodel=function(TMWB=TMWB,fcres=.25,SFTmp=0,bmlt6=2.5,bmlt12=1,Tlag=.5,AWCval=200,Slope=0){
  attach(TMWB)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = Slope, aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
  # Note that the -3 in the above 
  detach(TMWB)
  TMWB$SNO=SNO_Energy$SnowWaterEq_mm
  TMWB$SNOmlt=SNO_Energy$SnowMelt_mm
  attach(TMWB)
  TMWB$Albedo=.23
  TMWB$Albedo[TMWB$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,lat_radians = myflowgage$declat*pi/180) * 1000
  TMWB$PET=PET
  detach(TMWB)
  # add in rm
  rm(list=c("PET"))
  
  
  #TMWB$AWC=(0.45-0.15)*1000 #Fld Cap = .45, Wilt Pt = .15, z=1000mm
  TMWB$AWC=AWCval
  # Oh, this we want to vary some of these around our watershed!
  TMWB$dP = 0 # Initializing Net Precipitation
  TMWB$ET = 0 # Initializing ET
  TMWB$AW = 0 # Initializing AW
  TMWB$Excess = 0 # Initializing Excess
  # Functions for the Thornthwaite-Mather
  #
  
  # Loop to calculate AW and Excess
  attach(TMWB)
  for (t in 2:length(AW)){
    # This is where ET and Net Precipitation is now calculated
    ET[t] = min (AW[t-1],PET[t])
    ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model
    if(AvgTemp[t] >= SFTmp){
      dP[t] = P[t] - ET[t] + SNOmlt[t]
    }  else {
      dP[t] = ET[t]
    }
    # From here onward, everything is the same as Week2’s lab
    if (dP[t]<=0) {
      values<-soildrying(AW[t-1],dP[t],AWC[t])
    } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
      values<-soilwetting(AW[t-1],dP[t],AWC[t])
    } else {
      values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
    }
    AW[t]<-values[1]
    Excess[t]<-values[2]
  }
  TMWB$AW=AW
  TMWB$Excess=Excess
  # TMWB$dP=Excess  # This was in error originally
  TMWB$dP=dP
  TMWB$ET=ET
  detach(TMWB) # IMPORTANT TO DETACH
  rm(list=c("AW","Excess","dP","ET"))
  
  TMWB$Qpred=NA
  TMWB$Qpred[1]=0
  TMWB$S=NA
  TMWB$S[1]=0
  attach(TMWB)
  #fcres=.3        # Oh, this we want to vary in different areas
  for (t in 2:length(Qpred)){
    S[t]=S[t-1]+Excess[t]     
    Qpred[t]=fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  TMWB$S=S
  TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  detach(TMWB) # IMPORTANT TO DETACH
  rm(list=c("S","Qpred"))
  return(TMWB)
}

TopSlope=TMWB
MidSlope=TMWB
BotSlope=TMWB  

top1=TMWBmodel(TMWB=TopSlope,Slope=0.0759)
MidSlope$P=top1$Excess+TMWB$P
mid1=TMWBmodel(TMWB=MidSlope,Slope = 0.5036)
BotSlope$P=mid1$Excess+TMWB$P
bot1=TMWBmodel(TMWB=BotSlope,Slope=0)



colors <- c("Top" = "green3", "Mid" = "blue","Bottom"="red")

ggplot()+
  
  geom_line(data=top1, aes(x=date,y=Excess,color='Top'), size=.5) +
  geom_line(data=mid1, aes(x=date,y=Excess,color='Mid'), size=.5)+
  geom_line(data=bot1, aes(x=date,y=Excess,color='Bottom'), size=.5)+
  
  theme(
    axis.title.y = element_text(color = 'black', size=12),
    axis.title.y.right = element_text(color = 'black', size=12),
    plot.title = element_text(hjust = 0.5,size=14),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
    axis.ticks.length=unit(-0.1, "cm")
  ) +
  xlab('Year')+
  ylab('Excess[mm/day]')+
  ggtitle('Excess at HOLIDAY CREEK NEAR ANDERSONVILLE, VA')+
  scale_color_manual('Parameters:',values = colors)

colors <- c("Top" = "red", "Mid" = "blue","Bottom"="darkgreen")

ggplot()+
  
  geom_line(data=top1, aes(x=date,y=AW,color='Top'), size=.5) +
  geom_line(data=mid1, aes(x=date,y=AW,color='Mid'), size=.5)+
  geom_line(data=bot1, aes(x=date,y=AW,color='Bottom'), size=.5)+
  
  theme(
    axis.title.y = element_text(color = 'black', size=12),
    axis.title.y.right = element_text(color = 'black', size=12),
    plot.title = element_text(hjust = 0.5,size=14),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
    axis.ticks.length=unit(-0.1, "cm")
  ) +
  xlab('Year')+
  ylab('AW[mm/day]')+
  ggtitle('Excess at HOLIDAY CREEK NEAR ANDERSONVILLE, VA')+
  scale_color_manual('Parameters:',values = colors)
         
         
         ## Question2  ,..............
## .........
# First associate mukey with cokey from component

install.packages("httr")
install.packages("aqp")
install.packages("soilDB")
install.packages("rgdal")
install.packages("raster")
install.packages("rgeos")
install.packages("knitr")

library(aqp)
library(soilDB)
library(sp)
library(rgdal)
library(raster)
library(rgeos)


s <- fetchKSSL(bbox=c(-78.84301, 37.08077 , -77.27340  , 37.54968))
coordinates(s) <- ~ x + y
       
proj4string(s) <- '+proj=longlat +datum=WGS84'

# extract "site data"
s.sp <- as(s, 'SpatialPointsDataFrame')

# perform SDA query on collection of points
mu.data <- SDA_spatialQuery(s.sp, what = 'geom')

# use local spatial intersection to link source + SDA data
s.sp$mukey <- over(s.sp, mu.data)$mukey

# join results to original SoilProfileCollection using 'pedlabsampnum'
site(s) <- as.data.frame(s.sp)[, c('pedlabsampnum', 'mukey')]



par(mar=c(0,2,4,0))
groupedProfilePlot(s, groups='mukey', group.name.cex=1, color='clay', name='hzn_desgn', id.style='side', label='pedon_id', max.depth=100)
# describe IDs
mtext('user pedon ID', side=2, line=-1.5)
mtext('mukey', side=3, line=-1, at = c(0,0), adj = 0)



par(mar=c(0,2,4,0))
groupedProfilePlot(s, groups='mukey', group.name.cex=1, color='sand', name='hzn_desgn', id.style='side', label='pedon_id', max.depth=100)
# describe IDs
mtext('user pedon ID', side=0.25, line=-1)
mtext('mukey', side=1, line=-10, at = c(0,0), adj = 0)



par(mar=c(0,2,4,0))
groupedProfilePlot(s, groups='mukey', group.name.cex=1, color='silt', name='hzn_desgn', id.style='side', label='pedon_id', max.depth=100)
# describe IDs
mtext('user pedon ID', side=0.5, line=-1)
mtext('mukey', side=1, line=-10, at = c(0,0), adj = 0)
