# functions ...............
NSE = function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}


soil_wetting_above_capacity=function(AWprev,dP,AWC){
  AW<-AWC
  excess<-AWprev+dP-AWC
  c(AW,excess)
}


soildrying =function(AWprev,dP,AWC){
  AW<-AWprev*exp(dP/AWC)
  excess<-0.0
  c(AW,excess)
}

soilwetting = function(AWprev,dP,AWC){
  AW<-AWprev+dP
  excess<-0.0
  c(AW,excess)
}

TMWBmodel = function(TMWB=TMWB,fcres=.25,SFTmp=0,bmlt6=2.5,bmlt12=1,Tlag=.5,AWCval=200,Slope=0){
  # Now complete the model… what flows from TopSlope to MidSlope, and down to 
  # BotSlope. How will these be connected?
  
  # notice that there is an Energy Balance based Snow Accumulation 
  # and Melt model in the EcoHydRology package.
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


CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0, 
                  fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  
  # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2, 
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN   
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a 
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t]    # CN Solution
    # Is the soil saturated, and thus can't take more dP? 
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction? 
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}


