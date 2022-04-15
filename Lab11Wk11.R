if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)
# KEEP OPEN AS YOU WILL BE WALKING THROUGH IT FOR LAB	
vignette("multisensi-vignette")
#

#PET function .....................................
PET_fromTemp <- function (Jday, Tmax_C, Tmin_C, lat_radians, AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, TerrestEmiss = 0.97, aspect = 0, slope = 0, forest = 0, PTconstant=1.26, AEparams=list(vp=NULL, opt="linear"))
{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_radians, Jday, Tmax_C, Tmin_C, albedo, forest, slope, aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}

J <- seq(from = 1, to = 365, by = 5)

PET_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp(lat=X$lat[i],
                              Jday=Jday, Tmax=X$Tmax[i], 
                              Tmin=(X$Tmax[i]-X$Trange[i]),
                              X$slope[i],X$aspect[i],X$albedo[i],X$TerrestEmiss[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}




n <- 10
set.seed(1234)
X <- data.frame(Tmax_C = runif(n, min = 5, max = 25),Trange = runif(n, min = 0,max = 10),
                slope = runif(n, min = 0.0, max = 0.2),aspect = runif(n, min = 0.0, max = 0.2),
                albedo = runif(n, min = 0.0, max = 1),TerrestEmiss = runif(n, min = 0.0, max = 1),
                lat_radians=runif(n, min = 0.0, max = 1.1))

Y <- PET_Looped (X,Jday = J)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Surface Short Wave Rad(W/m^2)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
} 

X <- expand.grid(Tmax = c(5,15,25), 
                 Trange = c(2,6,10), 
                 slope = c(0.1,0.15,0.2),
                 aspect = c(0.1,.15,0.2),
                 albedo = c(0.1,.5,1.0),
                 TerrestEmiss=c(0.1,.5,1.0),
                 lat=c(0.1,.66,1.1))

Y <- PET_Looped(X,Jday=J) 
PET_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

dev.off() # Clean up previous par()
plot(PET_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(PET_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")

PET_Looped.pca <- multisensi(model=PET_Looped, reduction=basis.ACP, scale=FALSE,
                             design.args = list(Tmax = c(5,15,25), 
                                                Trange = c(2,6,10), 
                                                slope = c(0.1,0.15,0.2),
                                                aspect = c(0.1,.15,0.2),
                                                albedo = c(0.1,.5,1.0),
                                                TerrestEmiss=c(0.1,.5,1.0),
                                                lat=c(0.1,.66,1.1)))

summary(PET_Looped.pca, digits = 2)
dev.off()
plot(PET_Looped.pca, graph = 1)
dev.off()
plot(PET_Looped.pca, graph = 2)
dev.off()
plot(PET_Looped.pca, graph = 3)

#.............................sobol2007..........................................................

m <- 10000
X_Sobo <- data.frame(Tmax = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 albedo = runif(n, min = 0.0, max = 1),
                 TerrestEmiss = runif(n, min = 0.0, max = 1),
                 lat=runif(m, min = 0.0, max = 1.1))

PET_Looped.seq.sobol <- multisensi(design = sobol2007, model = PET_Looped,
                                   reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                   design.args = list(X1 = X_Sobo[1:(m/2), ], X2 = X_Sobo[(1 + m/2):m, ], nboot = 100),
                                   analysis.args = list(keep.outputs = FALSE))

print(PET_Looped.seq.sobol, digits = 2)
dev.off()
plot(PET_Looped.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

#fast99
PET_Looped.seq.fast <- multisensi(design = fast99, model = PET_Looped,
                                  center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                  design.args=list( factors=c("Tmax","Trange","slope","aspect","lat",'albedo','TerrestEmiss'), 
                                                    n=1000, q = "qunif",
                                                    q.arg = list(list(min=5, max=30), 
                                                                 list(min=2, max=16),
                                                                 list(min=0, max=.2),
                                                                 list(min=0, max=.2),
                                                                 list(min = 0.0, max = 1.1),
                                                                 list(min = 0.0, max = 1),
                                                                 list(min = 0.0, max = 1))),
                                  analysis.args=list(keep.outputs=FALSE))

print(PET_Looped.seq.fast,digits=2)

plot(PET_Looped.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

############ Net

J <- seq(from = 1, to = 365, by = 5)

NetRad_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- NetRad(lat=X$lat[i],
                       Jday=Jday, Tx=X$Tx[i], 
                       Tn=(X$Tx[i]-X$Trange[i]), 
                       X$slope[i],X$aspect[i],units="Wm2")
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}
n <- 10
set.seed(1234)
X <- data.frame(Tx = runif(n, min = 5, max = 30), Trange = runif(n, min = 2,max = 16),
                slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                lat=runif(n, min = 0.0, max = 1.1))

Y <- NetRad_Looped(X,Jday = J)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Surface Short Wave Rad(W/m^2)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}  

NetRad_Looped.seq <- multisensi(model=NetRad_Looped, reduction=NULL, center=FALSE,
                                design.args = list( Tx = c(5,15,25), 
                                                    Trange = c(2,9,16), 
                                                    slope = c(0.1,0.2,0.3),
                                                    aspect = c(0.1,.5,1.0),
                                                    lat=c(0.1,.77,1.1)))

dev.off() # Clean up previous par()
plot(NetRad_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(NetRad_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")

X <- expand.grid(Tx = c(5,15,25), 
                 Trange = c(2,9,16), 
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 lat=c(0.1,.77,1.1))
Y <- NetRad_Looped(X,Jday=J) # can be performed outside R if necessary
NetRad_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

NetRad_Looped.pca <- multisensi(model=NetRad_Looped, reduction=basis.ACP, scale=FALSE,
                                design.args = list( Tx = c(5,15,25), 
                                                    Trange = c(2,9,16), 
                                                    slope = c(0.1,0.2,0.3),
                                                    aspect = c(0.1,.5,1.0),
                                                    lat=c(0.1,.77,1.1)))


dev.off()
plot(NetRad_Looped.pca, graph = 1)
plot(NetRad_Looped.pca, graph = 2)
plot(NetRad_Looped.pca, graph = 3)

#sobol2007
m <- 10000
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 lat=runif(m, min = 0.0, max = 1.1))

NetRad_Looped.seq.sobol <- multisensi(design = sobol2007, model = NetRad_Looped,
                                      reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                      design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                      analysis.args = list(keep.outputs = FALSE))

print(NetRad_Looped.seq.sobol, digits = 2)
dev.off()
#plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)
plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

#fast99
NetRad_Looped.seq.fast <- multisensi(design = fast99, model = NetRad_Looped,
                                     center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                     design.args=list( factors=c("Tx","Trange","slope","aspect","lat"), 
                                                       n=1000, q = "qunif",
                                                       q.arg = list(list(min=5, max=30), 
                                                                    list(min=2, max=16),
                                                                    list(min=0, max=.2),
                                                                    list(min=0, max=.2),
                                                                    list(min = 0.0, max = 1.1))),
                                     analysis.args=list(keep.outputs=FALSE))

print(NetRad_Looped.seq.fast,digits=2)
plot(NetRad_Looped.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

################## question 2. 

CN1<- seq(from = 20, to = 100, by = 5)

SoilStorage_Looped <- function(X, CN = CN1) {
  S_av = (1000/CN-10)*25.4
  out <- matrix(nrow = nrow(X), ncol = length(S_av), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- SoilStorage(S_avg=S_avg[i],field_capacity=X$field_capacity[i],
                              soil_water_content=X$soil_water_content[i], 
                            porosity=X$porosity[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("CN", CN, sep = "")
  return(out)
}



n <- 10
set.seed(1234)
X <- data.frame(field_capacity = runif(n, min = 0.11, max = 0.43),
                soil_water_content = runif(n, min = 0.1,max = 0.3),
                porosity=runif(n, min = 0.2, max = 0.4))

Y <- SoilStorage_Looped (X,CN = CN1)


X <- expand.grid(field_capacity = c(0.11,0.27,0.43), 
                 soil_water_content = c(0.1,0.2,0.3), 
                 porosity = c(0.2,0.3,0.4))

Y <- SoilStorage_Looped(X,CN = CN1) 
SoilStorage_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

dev.off() # Clean up previous par()
plot(SoilStorage_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Curve Number")
plot(SoilStorage_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Average S")


SoilStorage.pca <- multisensi(model=SoilStorage_Looped, reduction=basis.ACP, scale=FALSE,
                             design.args = list(field_capacity = c(0.11,0.27,0.43), 
                                                soil_water_content = c(0.1,0.2,0.3), 
                                                porosity = c(0.2,0.3,0.4)))

summary(SoilStorage.pca, digits = 2)
dev.off()
plot(SoilStorage.pca, graph = 1)
dev.off()
plot(SoilStorage.pca, graph = 2)
dev.off()
plot(SoilStorage.pca, graph = 3)





#sobol2007
m <- 1000
Xb <- data.frame(field_capacity = runif(m, min = 0.11, max = 0.43),
                 soil_water_content = runif(m, min = 0.1,max = 0.3),
                 porosity=runif(m, min = 0.2, max = 0.4))

SoilStorage_Looped.seq.sobol <- multisensi(design = sobol2007, model = SoilStorage_Looped,
                                      reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                      design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                      analysis.args = list(keep.outputs = FALSE))

print(SoilStorage_Looped.seq.sobol, digits = 2)
dev.off()
plot(SoilStorage_Looped.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)

#fast99
SoilStorage_Looped.seq.fast <- multisensi(design = fast99, model = SoilStorage_Looped,
                                     center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                     design.args=list( factors=c("field_capacity","soil_water_content","porosity"), 
                                                       n=100, q = "qunif",
                                                       q.arg = list(list(min=0.11, max=0.43), 
                                                                    list(min=0.1, max=0.3),
                                                                    list(min=0.2, max=0.4))),
                                     analysis.args=list(keep.outputs=FALSE))

print(SoilStorage_Looped.seq.fast,digits=2)
plot(SoilStorage_Looped.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
