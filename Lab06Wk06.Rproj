prodir=getwd()

options(repos ="http://cran.us.r-project.org")  # required to get latest libs

if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,ggplot2,patchwork)
pacman::p_load(EcoHydRology,rnoaa,curl,httr)

source("https://raw.githubusercontent.com/Miladvahidi44/BSE5304/main/functions.R")
download.file(laburl,"Lab05Sol.R")
file.edit("Lab05Sol.R")

myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

CNmodeldf=modeldata


mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)
proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)

stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
#summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
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
modeldata[is.na(modeldata)]=0 # A Quick BUT sloppy removal of NAs
TMWB=modeldata

summary(terrain(mydem, opt='slope',unit = "degrees"))

TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = 11,Tlag = .5,AWCval = 100)

CNmodeldf=modeldata

f1 <- function (x) {
  CNopt       =x[1]
  Iaopt       =x[2]
  func_DAWCopt=x[3]
  func_zopt   =x[4]
  fnc_fcresopt=x[5]
  CNmodelnew=CNmodel(CNmodeldf =CNmodeldf,CNavg=CNopt, IaFrac=Iaopt,fnc_aspect=0,func_DAWC=func_DAWCopt, func_z=func_zopt, fnc_fcres=fnc_fcresopt)
  1-NSE(CNmodelnew$Qmm,CNmodelnew$Qpred)  
}

lower <- c(35,.01,0.1,500,0.1)
upper <- c(99,.25,0.35,1000,0.3)

DEoptim(f1, lower, upper,control = DEoptim.control(itermax=40))


f2 <- function (x) {
  fcresopt=x[1]
  SFTmpopt=x[2]
  Tlagopt = x[3]
  AWCvalopt=x[4]
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcresopt,SFTmp = SFTmpopt,Tlag = Tlagopt,AWCval = AWCvalopt)
  1-NSE(TMWBnew$Qmm,TMWBnew$Qpred)  
}
lower <- c(0.1,6,0.1,100)
upper <- c(0.5,12,1,300)


DEoptim(f2, lower, upper,control = DEoptim.control(itermax=40))



## TMWB with optimized parameters
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.33,SFTmp = 9.10,Tlag =0.93,AWCval = 100)


## CN model with optimized parametrs
CNmodelnew=CNmodel(CNmodeldf =CNmodeldf,CNavg =96.53,IaFrac=0.026,fnc_slope=0, 
                   fnc_aspect=0,func_DAWC=0.33,func_z=988.446,fnc_fcres=0.29)



## Plotting Qmm and Flow 

NSE(TMWBnew$Qmm,TMWBnew$Qpred)
NSE(CNmodelnew$Qmm,CNmodelnew$Qpred)
ggplot() +
  geom_line(data=TMWBnew,aes(x=date, y = Qmm,colour="Qmm")) +
  geom_line(data=TMWBnew,aes(x=date, y = Qpred,colour="Qpred_TMWB,NSE=0.35")) +
  geom_line(data=CNmodelnew,aes(x=date, y = Qpred,colour="Qpred_CN,NSE=0.52")) +
  labs(x = 'Date', y = 'Flow (mm)')+
  scale_colour_manual("", 
                      breaks = c("Qmm", "Qpred_TMWB,NSE=0.35", "Qpred_CN,NSE=0.52"),
                      values = c("black", "blue","red"))+
  theme(legend.position = "bottom",
        plot.title=element_text(hjust=0.5),
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("Discharge Comparison between CN model 
          and TMWB against Observed data")





## 
pacman::p_load(lubridate, data.table)
BasinTMWB_JO=TMWBnew[(month(TMWBnew$date) > 5 
                      & month(TMWBnew$date) < 11),]
attach(BasinTMWB_JO)
plot(dP,Qmm)
detach(BasinTMWB_JO)

attach(BasinTMWB_JO)
plot(dP,Qmm)
points(dP,dP^2/(dP+45),col="red")  # S guestimates in bold
points(dP,dP^2/(dP+260),col="blue")# S guestimates in bold
NSE(Qmm,dP^2/(dP+260))

NSE(Qmm,dP^2/(dP+45))


f <- function (x) {
  Sest=x
  NSE(Qmm,dP^2/(dP+Sest))
}
optimize(f, c(50,500), tol = 0.0001,maximum = TRUE)$maximum
Sest=150.7943
plot(dP,Qmm)
points(dP,dP^2/(dP+Sest),col="red") 

detach(BasinTMWB_JO)



nTIclass=5
VSAsol=data.table(WetClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
VSAsol 
#
# Now fill in the missing value
#
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
VSAsol
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol
plot(VSAsol$As,VSAsol$sigma)
lines(VSAsol$As,VSAsol$sigma)
plot(VSAsol$As,VSAsol$CN)
lines(VSAsol$As,VSAsol$CN)


TIC05=modeldata
TIC04=modeldata
TIC03=modeldata
TIC02=modeldata
TIC01=modeldata

TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[1],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[2],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[4],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAsol$CN[5],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC04$P=TIC05$Excess+TIC04$P
TIC03$P=TIC04$Excess+TIC03$P
TIC02$P=TIC03$Excess+TIC02$P
TIC01$P=TIC02$Excess+TIC01$P

# Qmm Plots 
png("TIC05.png", width = 6.12, height = 3.5, units = 'in', res = 600) 
ggplot() +
  geom_line(data=TIC05,aes(x=date, y = Qpred)) +
  labs(x = 'Date', y = 'Qpred (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class5)")
dev.off()


png("TIC04.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC04,aes(x=date, y = Qpred)) +
  labs(x = 'Date', y = 'Qpred (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class4)")
dev.off()

png("TIC03.png", width = 6.12, height = 3.5, units = 'in', res = 600) 
ggplot() +
  geom_line(data=TIC03,aes(x=date, y = Qpred)) +
  labs(x = 'Date', y = 'Qpred (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class3)")
dev.off()

png("TIC02.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC02,aes(x=date, y = Qpred)) +
  labs(x = 'Date', y = 'Qpred (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class2)")
dev.off()

png("TIC01.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC01,aes(x=date, y = Qpred)) +
  labs(x = 'Date', y = 'Qpred (mm)')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class1)")
dev.off()






# AW Plots 
png("AW5.png", width = 6.12, height = 3.5, units = 'in', res = 600) 
ggplot() +
  geom_line(data=TIC05,aes(x=date, y = AW)) +
  labs(x = 'Date', y = 'AW')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class5)")
dev.off()


png("AW4.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC04,aes(x=date, y = AW)) +
  labs(x = 'Date', y = 'AW')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class4)")
dev.off()

png("AW3.png", width = 6.12, height = 3.5, units = 'in', res = 600) 
ggplot() +
  geom_line(data=TIC03,aes(x=date, y = AW)) +
  labs(x = 'Date', y = 'AW')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class3)")
dev.off()

png("AW2.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC02,aes(x=date, y = AW)) +
  labs(x = 'Date', y = 'AW')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class2)")
dev.off()

png("AW1.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC01,aes(x=date, y = AW)) +
  labs(x = 'Date', y = 'AW')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class1)")
dev.off()




# ET Plots 
png("ET5.png", width = 6.12, height = 3.5, units = 'in', res = 600) 
ggplot() +
  geom_line(data=TIC05,aes(x=date, y = ET)) +
  labs(x = 'Date', y = 'ET')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class5)")
dev.off()


png("ET4.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC04,aes(x=date, y = ET)) +
  labs(x = 'Date', y = 'ET')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class4)")
dev.off()

png("ET3.png", width = 6.12, height = 3.5, units = 'in', res = 600) 
ggplot() +
  geom_line(data=TIC03,aes(x=date, y = ET)) +
  labs(x = 'Date', y = 'ET')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class3)")
dev.off()

png("ET2.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC02,aes(x=date, y = ET)) +
  labs(x = 'Date', y = 'ET')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class2)")
dev.off()

png("ET1.png", width = 6.12, height = 3.5, units = 'in', res = 600) 

ggplot() +
  geom_line(data=TIC01,aes(x=date, y = ET)) +
  labs(x = 'Date', y = 'ET')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),,text = element_text(family="Times New Roman",size = 15))+
  ggtitle("(TI class1)")
dev.off()


