Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
options(repos ="http://cran.r-project.org")  # required to get latest libs
if (!require("pacman")) install.packages("pacman")
pacman::p_load(aqp,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB,circlize,topmodel,DEoptim)
install.packages("EcoHydRology", repos="http://R-Forge.R-project.org")
pacman::p_load(EcoHydRology)

myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2022-03-24")

myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

#
declat=myflowgage$declat
declon=myflowgage$declon
WXData=FillMissWX(declat, declon,30,
                  date_min="2010-01-01",
                  date_max="2022-03-16")
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")

proj4_utm = paste0("+proj=utm +zone=", trunc((180+declon)/6+1), " +datum=WGS84 +units=m +no_defs")

proj4_ll = "+proj=longlat"

crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)

area=myflowgage$area   
latlon <- cbind(declon,declat)
gagepoint_ll <- SpatialPoints(latlon)
proj4string(gagepoint_ll)=proj4_ll
gagepoint_utm=spTransform(gagepoint_ll,crs_utm)
pourpoint=SpatialPoints(gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=gagepoint_utm@coords
searchlength=sqrt(area*4)*1000 
bboxpts=rbind(bboxpts,bboxpts)
bboxpts[2,1]=bboxpts[2,1]-searchlength
bboxpts[2,2]=bboxpts[2,2]+searchlength
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
mydem=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
Sys.getenv("PATH")
if(!grepl("TauDEM",Sys.getenv("PATH"))){
  old_path <- Sys.getenv("PATH")
  old_path
  Sys.setenv(PATH = paste(old_path,
                          paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
                          sep = ":"))
}
# Test to make sure TauDEM runs
system("mpirun aread8")

writeRaster(mydem,filename = "mydem.tif",overwrite=T)
# remember our intro to terminal
# ls; cd ~; pwd;  #Linux/Mac
# dir; cd ; 

z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)

plot(z-fel)
# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
p=raster("mydemp.tif")
plot(p)
sd8=raster("mydemsd8.tif")
plot(sd8)

# Contributing area
system("mpiexec -n 8 aread8 -p mydemp.tif -ad8 mydemad8.tif")
ad8=raster("mydemad8.tif")
plot(log(ad8))

# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
plot(ang)
slp=raster("mydemslp.tif")
plot(slp)

# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
plot(log(sca))
#
# This area could be useful in 
threshold=area*10^6/(res(mydem)^2)[1]/10
threshold
# Threshold
system(paste0("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh ",threshold))
src=raster("mydemsrc.tif")
plot(src,xlim=c(pourpoint@coords[1]-5*res(src)[1],pourpoint@coords[1]+5*res(src)[1]),
     ylim=c(pourpoint@coords[2]-5*res(src)[2],pourpoint@coords[2]+5*res(src)[2]))
plot(pourpoint,add=T)

pourpointSPDF=SpatialPointsDataFrame(pourpoint,data.frame(Outlet=row(pourpoint@coords)[1]))
writeOGR(pourpointSPDF, dsn=".", "ApproxOutlets", driver="ESRI Shapefile", overwrite_layer = TRUE)


# Move Outlets
system("mpiexec -n 8  moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=readOGR("Outlet.shp")
points(outpt,col=2,add=T)

# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 

# Threshold
system(paste0("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh ",threshold))
src1=raster("mydemsrc1.tif")
plot(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
mydemw=raster("mydemw.tif")
plot(mydemw)
# Plot streams using stream order as width
plot(readOGR("mydemnet.dbf"),add=T)
streams=readOGR("mydemnet.dbf")
# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi,add=T)

mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
plot(TIC)
# Make a poly with command line gdal (fast)
# 
system("gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp")
mydemw_poly=readOGR("mydemw_poly_gdal.shp")
plot(mydemw_poly,add=T,border="black")
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)

#source CNmodel function
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
# Use the estimated S for our watershed (Lab06)
Sest = 157
nTIclass=5
VSAsol=data.table(TIClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
#
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
# Calculate TI Class localized sigma and Curve Number
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol

TIC01=modeldata
TIC02=modeldata
TIC03=modeldata
TIC04=modeldata
TIC05=modeldata
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAsol$CN[1], 
                declat=declat,declon=declon)
TIC02$P=TIC01$Qpred+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[2], 
                declat=declat,declon=declon)
TIC03$P=TIC02$Qpred+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3], 
                declat=declat,declon=declon)
TIC04$P=TIC03$Qpred+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[4], 
                declat=declat,declon=declon)
TIC05$P=TIC04$Qpred+TIC05$P
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[5], 
                declat=declat,declon=declon)
DPTI05=data.frame(date=TIC05$date,
                  Rt=TIC05$Qpred/1000, 
                  Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)

PLossFunc=function(DPTIdf=DPTI01
                   ,tau=9.3  # days
                   ,dt=1     # days time step
                   ,kF=.015  # Table 2
                   ,TICn=1){

  DPTIdf$MF=0
  DPTIdf$DF=0
  DPTIdf$MF[(format(DPTIdf$date,"%j") %in% c(121,213,274))]=5.4*10^-4
  attach(DPTIdf)
  #
  # Loop to solve MF and DF
  for (i in 2:length(date)){
    if(MF[i]<=MF[i-1]){
      MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
    }
    DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
  }
  DPTIdf$MF=MF
  DPTIdf$DF=DF
  detach(DPTIdf)
  rm(list=c("MF","DF")) # Clean up the environment 
  
  plot(DPTI05$date,DPTI05$MF)
  plot(DPTI05$date,DPTI05$DF)
  
  
  muTS_TI=(((520-0)/5)*TICn+0) 
  MS_TI=(((18.5-3.7)/5)*TICn+3.7)*2000   

  QS= 3.0 
  TR=20   
  DPTIdf$muS= muTS_TI*QS^((DPTIdf$Tavg-TR)/10) 
  DPTIdf$DS=(DPTIdf$muS*MS_TI*DPTIdf$Rt)/10^9     
  return(DPTIdf)
}

### Question 2 ............

TIC01$Rt=TIC01$Qpred/1000
TIC01$Tavg=(TIC01$MaxTemp+TIC01$MinTemp)/2
TIC01=PLossFunc(DPTIdf = TIC01,tau=2,dt=1,kF=.015,TICn=1)

TIC02$Rt=TIC02$Qpred/1000
TIC02$Tavg=(TIC02$MaxTemp+TIC02$MinTemp)/2
TIC02=PLossFunc(DPTIdf = TIC02,tau=2,dt=1,kF=.015,TICn=2)

TIC03$Rt=TIC03$Qpred/1000
TIC03$Tavg=(TIC03$MaxTemp+TIC03$MinTemp)/2
TIC03=PLossFunc(DPTIdf = TIC03,tau=2,dt=1,kF=.015,TICn=3)

TIC04$Rt=TIC04$Qpred/1000
TIC04$Tavg=(TIC04$MaxTemp+TIC04$MinTemp)/2
TIC04=PLossFunc(DPTIdf = TIC04,tau=2,dt=1,kF=.015,TICn=4)

TIC05$Rt=TIC05$Qpred/1000
TIC05$Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2
TIC05=PLossFunc(DPTIdf = TIC05,tau=2,dt=1,kF=.015,TICn=5)

DPLT=data.frame(date=TIC05$date,
                Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9

TIC01$area=myflowgage$area/5*10^6
TIC02$area=myflowgage$area/5*10^6
TIC03$area=myflowgage$area/5*10^6
TIC04$area=myflowgage$area/5*10^6
TIC05$area=myflowgage$area/5*10^6

DPLT$LT=DPLT$LB +
  TIC01$area*(TIC01$DF + TIC01$DS)+
  TIC02$area*(TIC02$DF + TIC02$DS)+
  TIC03$area*(TIC03$DF + TIC03$DS)+
  TIC04$area*(TIC04$DF + TIC04$DS)+
  TIC05$area*(TIC05$DF + TIC05$DS)
mean(DPLT$LT)
# ...............

## Question 3 changing the tau parameter to 9.3. .........
TIC01$Rt=TIC01$Qpred/1000
TIC01$Tavg=(TIC01$MaxTemp+TIC01$MinTemp)/2
TIC01=PLossFunc(DPTIdf = TIC01,tau=9.3,dt=1,kF=.015,TICn=1)

TIC02$Rt=TIC02$Qpred/1000
TIC02$Tavg=(TIC02$MaxTemp+TIC02$MinTemp)/2
TIC02=PLossFunc(DPTIdf = TIC02,tau=9.3,dt=1,kF=.015,TICn=2)

TIC03$Rt=TIC03$Qpred/1000
TIC03$Tavg=(TIC03$MaxTemp+TIC03$MinTemp)/2
TIC03=PLossFunc(DPTIdf = TIC03,tau=9.3,dt=1,kF=.015,TICn=3)

TIC04$Rt=TIC04$Qpred/1000
TIC04$Tavg=(TIC04$MaxTemp+TIC04$MinTemp)/2
TIC04=PLossFunc(DPTIdf = TIC04,tau=9.3,dt=1,kF=.015,TICn=4)

TIC05$Rt=TIC05$Qpred/1000
TIC05$Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2
TIC05=PLossFunc(DPTIdf = TIC05,tau=9.3,dt=1,kF=.015,TICn=5)

DPLT=data.frame(date=TIC05$date,
                Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
mean(DPLT$LB)


##Plotting P loss from soil (DS) and baseflow (LB)
ggplot(DPLT, aes(x=date)) +
  geom_line(aes( y=LB),size = 0.2) + 
  xlab("")+
  scale_x_date(date_breaks = "60 days", expand = c(0.005,0.005))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0,size=10, angle = 90 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank()) 
ggplot(TIC05, aes(x=date)) +
  geom_line(aes( y=DS),size = 0.2) + 
  xlab("")+
  scale_x_date(date_breaks = "60 days", expand = c(0.005,0.005))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0,size=10, angle = 90 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank())

TIC01$area=myflowgage$area/5*10^6
TIC02$area=myflowgage$area/5*10^6
TIC03$area=myflowgage$area/5*10^6
TIC04$area=myflowgage$area/5*10^6
TIC05$area=myflowgage$area/5*10^6

DPLT$LT=DPLT$LB +
  TIC01$area*(TIC01$DF + TIC01$DS)+
  TIC02$area*(TIC02$DF + TIC02$DS)+
  TIC03$area*(TIC03$DF + TIC03$DS)+
  TIC04$area*(TIC04$DF + TIC04$DS)+
  TIC05$area*(TIC05$DF + TIC05$DS)

mean(DPLT$LT)
plot(DPLT$date,DPLT$LT, type="l")
ggplot(DPLT, aes(x=date)) +
  geom_line(aes( y=LT),size = 0.2) + 
  xlab("")+
  scale_x_date(date_breaks = "60 days", expand = c(0.005,0.005))+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0,size=10, angle = 90 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank())

## In 20cm of soil depth,

amount = (0.2)*mean(c(18.5,3.7))*2000/1000
#To lose 10% of it we just need to 
Updated_amount = 0.1*amount
.1 * .2 * mean(c(18.5,3.7))*2000/1000 # 20% of 10% of 1m , .02m * kg/m^3 = [kg/m^2]

#Whole the amount that we lose from 2016 to 2022 can be calculated using sum function
Total = sum(TIC03$DS)
Years = Updated_amount/Total*(2022-2016)


rcp01=mean(TIC01$DF+TIC01$DS)*365
rcp02=mean(TIC02$DF+TIC02$DS)*365
rcp03=mean(TIC03$DF+TIC03$DS)*365
rcp04=mean(TIC04$DF+TIC04$DS)*365
rcp05=mean(TIC05$DF+TIC05$DS)*365

### 


p_loss<-function(Runoffmodel=TIC05,tau=9.3, dt=1, kF=.015, TI=1){
  DPTI05=data.frame(date=TIC05$date,
                    Rt=TIC05$Qpred/1000, 
                    Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
  DPTI05$MF=0
  DPTI05$DF=0
  DPTI05$MF[(format(DPTI05$date,"%j") %in% c(121,213,274))]=5.4*10^-4
  # Remember what we say about attaching! 
  attach(DPTI05)
  #
  # Loop to solve MF and DF
  for (i in 2:length(date)){
    if(MF[i]<=MF[i-1]){
      MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
    }
    DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
  }
  DPTI05$MF=MF
  DPTI05$DF=DF
  detach(DPTI05)
  rm(list=c("MF","DF")) # Clean up the environment 
  muTS_TI05=(((520-0)/5)*VSAsol$TIClass[TI]+0) # use VSAsol$TIClass table to 

  MS_TI05=(((18.5-3.7)/5)*VSAsol$TIClass[TI]+3.7)*2000  # range from Easton et 

  QS= 3.0 # A guess using the middle of the range 1-5
  TR=20   # reference Temperature from Table 2.
  DPTI05$muS= muTS_TI05*QS^((DPTI05$Tavg-TR)/10)  # Eq. 5
  DPTI05$DS=(DPTI05$muS*MS_TI05*DPTI05$Rt)/10^6          # Eq. 4
  
  DPTI05$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
  muTB=2.1*10^(-5) # Easton Table 2
  QB=2.2           # Easton Table 2
  TB=17            # Easton Table 2
  DPTI05$muB=muTB*QB^((DPTI05$Tavg-TB)/10)  # Easton eq. 10
  DPTI05$LB=DPTI05$muB*DPTI05$B     # Easton eq. 9
  return(DPTI05)
  
}
