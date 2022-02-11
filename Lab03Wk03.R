source("https://goo.gl/Cb8zGn")
if (!require("pacman"))
install.packages("pacman")
pacman::p_load(rnoaa,scales,extrafont,EcoHydRology,extrafontdb,tidyverse,lattice,ggplot2,dplyr,patchwork,hrbrthemes)#load library

myflowgage_id="14216500"
myflowgage=get_usgs_gage(myflowgage_id,
                         begin_date="2017-02-01",end_date="2022-02-01")
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 20,
  limit = NULL
)
Data = data_frame("Date"=myflowgage$flowdata$mdate ,
                  "Flow"=myflowgage$flowdata$flow,
                  "Quality"=myflowgage$flowdata$quality,
                  "Site Number"=myflowgage$flowdata$flow)

WXData=meteo_pull_monitors(
  monitors=stns[13,1],   
  keep_flags = FALSE,
  date_min = "2017-02-01",
  date_max = "2022-02-01",
  var = c("TMAX","TMIN","PRCP") 
)

myflowgage$flowdata$flowmm=myflowgage$flowdata$flow/myflowgage$area/10^3

modeldata=merge(WXData, myflowgage$flowdata, by.x="date", by.y="mdate")

modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$Qmm=modeldata$flowmm
modeldata$P=modeldata$prcp/10 # Converting to mm

TMWB=modeldata

##**Ploting ET, QMM and P

SFTmp = 1  # referred to as SFTMP in SWAT input (Table 1)
bmlt6 = 4.5   # referred to as SMFMX in SWAT input (Table 1)
bmlt12 = 0.0  # referred to as SMFMN in SWAT input adjusted for season
Tmlt = SFTmp  # Assumed to be same as SnowFall Temperature
Tlag = 1  # referred to as TIMP in SWAT input (Table 1)
TMWB$AvgTemp= (TMWB$MaxTemp + TMWB$MinTemp)/2
TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))
# Initialize SNO, Tsno as well as the first values of each
TMWB$SNO = 0  # Snow Depth (mm)
TMWB$Tsno = 0  # Snow Temp (C)
TMWB$SNOmlt = 0  # Snow Melt (mm)
TMWB$AvgTemp[is.na(TMWB$AvgTemp)]=0

attach(TMWB)
for (t in 2:length(date)){
  Tsno[t]= Tsno[t-1] * (1.0-Tlag) +  AvgTemp[t] * Tlag
  if(AvgTemp[t] < SFTmp){
    SNO[t]= SNO[t-1] + P[t]
    }  else {
      SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t]+MaxTemp[t])/2 - Tmlt) 
      SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
      SNO[t]= SNO[t-1] -SNOmlt[t]
      }
}
detach(TMWB)
TMWB$Tsno=Tsno
TMWB$SNO=SNO
TMWB$SNOmlt=SNOmlt
rm(list=c("SNO", "SNOmlt", "Tsno"))


 ggplot(TMWB, aes(x=date)) +
   geom_line(aes( y=SNO),size = 0.2) + 
   xlab("")+
   scale_x_date(date_breaks = "60 days", expand = c(0.005,0.005))+
   ggtitle('MUDDY RIVER BELOW CLEAR CREEK NEAR COUGAR, WA')+
   theme(legend.position = "bottom",
         legend.justification = c("center"),
         legend.box.just = "left",
         legend.background = element_blank(),
         axis.text.x.bottom = element_text(vjust = 0,size=10, angle = 90 ),
         plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
         text=element_text(family="Times New Roman"),
         axis.ticks = element_blank()) 
 
 
 ##Describe in words what you see in the plot of Snow Water Equivalent, which slope aspect combination has the greatest SWE, which has the least, why?

TMWB[is.na(TMWB)] = 0
SS = cbind(atan(c(0,10,10,45,45)/100) , pi/180*c(0,0,180,315,225) )
sne <- matrix(0, NROW(TMWB), 5)
colnames(sne) <- paste("SNE", 1:5, sep = "") 

attach(TMWB)

for (i in 1:5) {
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope =  SS[i,1] ,
                      aspect = SS[i,2] , tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
  sne[,i] = SNO_Energy$SnowWaterEq_mm
  }

detach(TMWB)

SNE = as.data.frame(sne)
TMWB = cbind(TMWB , SNE)

colors <- c("SNE1" = "green",
            "SNE2" = "red",
            "SNE3" = "blue",
            "SNE4" = "black",
            "SNE5" = "gray")

ggplot(TMWB, aes(x= date)) +
  geom_line(aes( y= SNE1,color = "SNE1"),size = 1) +
  geom_line(aes( y= SNE2,color = "SNE2"),size = 1) + 
  geom_line(aes( y= SNE3,color = "SNE3"),size = 1) + 
  geom_line(aes( y= SNE4,color = "SNE4"),size = 1) + 
  geom_line(aes( y= SNE5,color = "SNE5"),size = 1) + 
  xlab("")+
  ylab("SNE ")+
  labs(x = "", colour = "Legend")+
  scale_color_manual("SNE", values = colors)+
  ggtitle('MUDDY RIVER BELOW CLEAR CREEK NEAR COUGAR, WA')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0,size=10, angle = 90 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank()) 

##**Plotting mean and commulative of SNE

df1 = apply(sne,2, mean  )
df2 = apply(sne,2, sum  )

D1 = data_frame("SNE" = c("flat slope (0% slope)", "10% slope North facing", "10% slope South facing",
                          "45% slope Northwest facing","45% slope Southwest facing"), "SNEaverage" = df1)
D2 = data_frame("SNE" = c("flat slope (0% slope)", "10% slope North facing", "10% slope South facing",
                   "45% slope Northwest facing","45% slope Southwest facing"), "SNESum" = df2)


ggplot(D1, aes(x = SNE , y= SNEaverage)) +
  geom_bar(stat="identity", fill="green3", width = 0.05)+
  ylab("SNE Average")+
  ggtitle('MUDDY RIVER BELOW CLEAR CREEK NEAR COUGAR, WA')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0.5,size=10, angle = 45 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank()) 
        
  ggplot(D2, aes(x = SNE , y= SNESum)) +
  geom_bar(stat="identity", fill="green3", width = 0.05)+
  ylab("SNE Sum")+
  ggtitle('MUDDY RIVER BELOW CLEAR CREEK NEAR COUGAR, WA')+
  theme(legend.position = "bottom",
        legend.justification = c("center"),
        legend.box.just = "left",
        legend.background = element_blank(),
        axis.text.x.bottom = element_text(vjust = 0.5,size=10, angle = 45 ),
        plot.title = element_text(hjust = 0.5, vjust = 0, size=12),
        text=element_text(family="Times New Roman"),
        axis.ticks = element_blank())      
        




