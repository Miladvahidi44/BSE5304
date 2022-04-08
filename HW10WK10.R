if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)
setwd("~/Week10/")



make_usgs_gage_list=function(siteNo,parameterCd ,start.date,end.date){
  
  USGSlist=list()   # Organize the data in a nice list as in previous labs
  USGSlist[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  head(USGSlist$flowdata)  # Note that we have 00060 and 00065...
  # And of course we want to work in SI units so:
  USGSlist$flowdata$depth_m=USGSlist$flowdata$X_00065_00000*0.3048-.8
  # m/ft depth
  USGSlist$flowdata$cms=USGSlist$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  USGSlist[["site"]]=readNWISsite(siteNo)
  head(USGSlist$site)
  class(USGSlist$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  USGSlist$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(USGSlist$site)=~dec_long_va+dec_lat_va
  #
  return(USGSlist)
}

#USGS08374550=make_usgs_gage_list(siteNo = "08374550", parameterCd = c("00060","00065"),start.date = "2017-05-01",  end.date = "2017-11-01")
USGS07010000 =make_usgs_gage_list(siteNo = "07010000",parameterCd = c("00060","00065"),start.date = "2017-05-01",  end.date = "2017-11-01" )
USGS07010035  =make_usgs_gage_list(siteNo = "07010035", parameterCd = c("00060","00065"),start.date = "2017-05-01",  end.date = "2017-11-01")
 


ab_ll=rbind(USGS07010000$site,
            USGS07010035 $site)

class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
       trunc((180+coordinates(USGS07010035 $site)[1])/6+1), 
       " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
          z = 12, prj = proj4_utm,expand=1)
#
# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# From Lab02, I know I can get an overview of streams with the 
# USGS H
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)



USGS07010000$flowdata=USGS07010000$flowdata[,c(1,2,3,4,5,8,10)]


USGS07010000[["rating"]]=readNWISrating(USGS07010000$site$site_no) 
plot(USGS07010000$rating$DEP,USGS07010000$rating$INDEP,xlab="DEP",ylab="INDEP")


USGS07010000$flowdata$X_00065_00000=approx(USGS07010000$rating$DEP,
                                           USGS07010000$rating$INDEP, xout = USGS07010000$flowdata$X_00060_00000, ties = min)$y
points(USGS07010000$flowdata$X_00060_00000,USGS07010000$flowdata$X_00065_00000,
       col="red")
#
USGS07010000$flowdata$depth_m=USGS07010000$flowdata$X_00065_00000*0.3048
# m/ft depth
#

A=SpatialPoints(USGS07010035 $site)# Up gradient site Lick Run
B=SpatialPoints(USGS07010000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS07010035 $site$L=SpatialLinesLengths(AtoB) # km to m
USGS07010035 $site$L # reach length in m
#
#
# Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
# estimate of slope than taking the point slopes at the gage site
#
USGS07010035 $site$slope=(raster::extract(mydem,A_utm)-
                           raster::extract(mydem,B_utm))/USGS07010035 $site$L
USGS07010035 $site$slope



#B=(n*Q)/(y^(5/3)*sqrt(So))
USGS07010035 $flowdata$B=(USGS07010035 $site$man_n*
                             USGS07010035 $flowdata$cms)/(USGS07010035 $flowdata$depth_m^(5/3)*
                                                             sqrt(USGS07010035 $site$slope))
head(USGS07010035 $flowdata)
#  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
#1  	USGS 05267000 2017-05-01 04:00:00      	6.38            	A
#2  	USGS 05267000 2017-05-01 04:05:00      	6.38            	A
#  X_00065_00000 X_00065_00000_cd tz_cd   	cms  depth_m    	B
#1      	2.74            	A   UTC 0.1806816 0.835152 0.103032
#2      	2.74            	A   UTC 0.1806816 0.835152 0.103032
#
# Lets look at how B changes with flow.
plot(USGS07010035 $flowdata$dateTime,USGS07010035 $flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")


plot(USGS07010035 $flowdata$dateTime,USGS07010035 $flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.


plot(USGS07010035 $flowdata$cms,USGS07010035 $flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")




#USGS07010035 $flowdata$ck = ???
# ANS
USGS07010035 $flowdata$ck =
  5/3*sqrt(USGS07010035 $site$slope)/USGS07010035 $site$man_n*
  (USGS07010035 $flowdata$depth_m^(2/3))

#USGS07010035 $flowdata$dt = ???
USGS07010035 $flowdata$dt =
  USGS07010035 $site$L/USGS07010035 $flowdata$ck

plot(USGS07010035 $flowdata$dateTime,USGS07010035 $flowdata$dt)
USGS07010035 $flowdata$outTime=USGS07010035 $flowdata$dateTime+
  USGS07010035 $flowdata$dt

# Find beginning of  Waves
USGS07010035 $flowdata$newwave=
  USGS07010035 $flowdata$cms *1.1 <
  data.table::shift(USGS07010035 $flowdata$cms)
summary(USGS07010035 $flowdata$newwave)
# Add plot of the point found
len=length(USGS07010035 $flowdata$newwave)
USGS07010035 $flowdata$newwave[is.na(USGS07010035 $flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS07010035 $flowdata$newwave[i]==T &
     USGS07010035 $flowdata$newwave[i-1]==T){
    USGS07010035 $flowdata$newwave[i]=F
  }
}
plot(USGS07010035 $flowdata$dateTime,USGS07010035 $flowdata$cms,type="l")
points(USGS07010035 $flowdata$dateTime[USGS07010035 $flowdata$newwave],
       USGS07010035 $flowdata$cms[USGS07010035 $flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS07010035 $flowdata$newwave == TRUE)
plot(USGS07010035 $flowdata$dateTime,USGS07010035 $flowdata$cms,
     type="l",xlim=c(USGS07010035 $flowdata$dateTime[1109],
                     USGS07010035 $flowdata$dateTime[1109+200]), main="Engelholm Creek near Wellston, MO, TO Mississippi River at St. Louis, MO")
lines(USGS07010035 $flowdata$outTime,USGS07010035 $flowdata$cms,col=2)

Data = data.frame("Date" =USGS07010035 $flowdata$dateTime, "B" = USGS07010035 $flowdata$B, "cms"=USGS07010035 $flowdata$cms,
                  "depth_m"=USGS07010035 $flowdata$depth_m, "outTime"=USGS07010035 $flowdata$outTime, "dt"=USGS07010035 $flowdata$dt)

png("cms.png", width = 6, height = 4, units = 'in', res = 600) 
ggplot(Data) +
  geom_line(aes(x = as.Date(Date)  , y= cms))+
  ylab("cms")+
  xlab("Date")+
  scale_x_date(date_breaks = "30 days", expand = c(0.005,0.005))+
  ggtitle('Engelholm Creek near Wellston, MO, TO Mississippi River at St. Louis, MO')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0.5,size=10, angle = 45 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman", size =15 ),
        axis.ticks = element_blank()) 
dev.off()

Data_subset= Data[1109:1309,]
png("D.png", width = 6, height = 4, units = 'in', res = 600) 

ggplot(Data_subset )+
  geom_line(aes(x = Date  , y= cms,color = 'black'))+
  geom_line(aes(x = outTime  , y= cms, color = 'red'))+
  ylab("cms")+
  xlab("Date")+
  #scale_x_date(date_breaks = "30 days", expand = c(0.005,0.005))+
  ggtitle('River Des Peres near University City, MO, TO Mississippi River at St. Louis, MO')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0.5,size=10, angle = 45 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman", size =15 ),
        axis.ticks = element_blank()) 
dev.off()

