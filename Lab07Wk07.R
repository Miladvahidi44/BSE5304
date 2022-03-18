objects()
 rm(list=objects())

 if (!require("pacman")) install.packages("pacman")
 pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr)

# Get our gold standard flow data from USGS 02038850 `
 myflowgage_id="02038850"
 myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                           end_date = "2019-01-01")

# We want Q in mm/day for the basin
 myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

 
trunc((180+myflowgage$declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
 print(proj4_utm)
 
 # Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
 
 # Now we will build our proj4strings which define our “Coordinate 
 # Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)




myflowgage$area   # area in km2

# If the watershed was square, which it is not, the size would be the 
# square root of the area. Also, the gage/pour point is not in the center
# so we will search around the gage.
# Build sp point for USGS gage location, in both _ll and _utm
latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
             myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)
# We are going to over estimate our area
sqrt(myflowgage$area)   # guestimating square watershed
# For our search we are going to multiply the area by 6 and
# to get the distance
sqrt(myflowgage$area*8)
 searchlength=sqrt(myflowgage$area*8)*1000 
 pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)
 bboxpts=myflowgage$gagepoint_utm@coords
 bboxpts=rbind(bboxpts,bboxpts+searchlength)
 bboxpts=rbind(bboxpts,bboxpts-searchlength)
 box = bboxpts
 
 bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
 bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
 bboxpts
 bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
 mydem=get_aws_terrain(locations=bboxpts@coords, 
                        z = 12, prj = proj4_utm,src ="aws",expand=1)
 res(mydem)
 plot(mydem)
 plot(bboxpts,add=T)
 plot(pourpoint,add=T,col="red")

# Write our raster to a geotiff file that can be used with
# OS level hydrological models 
writeRaster(mydem,filename = "mydem.tif",overwrite=T)

zoom(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")
   


old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
                        sep = ":"))

system("mpirun aread8")


library(raster)
library(shapefiles)

# Set working directory to your location
setwd("C:/Users/dtarb/Scratch/mydem")

z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)


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
zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
zoom(gord)

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
zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 4000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint, add = T)
zoom(src)
plot(pourpoint, add = T)
zoom(src)
plot(pourpoint, add = T)

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("qutlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

zoom(src)


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 4000")
src1=raster("mydemsrc1.tif")
plot(src1)
zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"), add = T)
zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"))

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)

#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Deduce that the optimal threshold is 300 
# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 4000")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


mydemw=raster("mydemw.tif")

mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasinmask)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasinmask)

Demplt = projectRaster(mydem , mydemw)
demmask=trim(Demplt,padding=2)
Demplt = crop(Demplt,mybasinmask)
Demplt = mask(Demplt,mybasinmask)
plot(Demplt)

## Downloading NLCD

Boarder =list(x=c(-78.590822,-78.725756 ), y=c( 37.381884,37.516446 ))
ANDERSONVILLE <- FedData::polygon_from_extent(raster::extent(Boarder), proj4string='+proj=longlat +datum=WGS84 +no_defs')
NLCD <- FedData::get_nlcd(template = ANDERSONVILLE,extraction.dir = "~",year = 2019,
label = "ANDERSONVILLE VA")
plot(NLCD)

NLCDplt = projectRaster(NLCD , mydemw)
demmask=trim(NLCDplt,padding=2)
NLCDplt = crop(NLCDplt,mybasinmask)
NLCDplt = mask(NLCDplt,mybasinmask)
plot(NLCDplt)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
zoom(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#
# A series of plots to show all of the components
#
par(mfrow = c(4, 2))
plot(TI, MAIN ="TI")
plot(TIC, MAIN ="TIC")
plot(mybasinsca)
plot(mybasinslp, MAIN ="slp")
plot(TIC,col=rainbow(5), MAIN ="TIC")
plot(Demplt, main ="DEM")
plot(NLCDplt , main = "NLCD")
plot(TIC,col=rainbow(5))
dev.off()

plot(TIC,col=rainbow(5))
# Make a poly with raster library (slow)
# or from thee command line gdal (fast)
# gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp
mydemw_poly=rasterToPolygons(mydemw,dissolve = T,na.rm = T)
plot(mydemw_poly,add=T,border="red")
mydemw_poly
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)


zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/4tv24bkq1czxdkcvlsxvfxtg/wss_aoi_2022-03-17_18-03-58.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
mysoil=readOGR("wss_aoi_2022-03-17_18-03-58/spatial/soilmu_a_aoi.shp")  

# Or, if you want to make your life easier, 
mydemw_poly_gdal=readOGR("mydemw.dbf")
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)

x <- mapunit_geom_by_ll_bbox(bbox(spTransform(mydemw_poly_gdal,crs_ll)))
x=spTransform(x,crs_utm) # convert LL to UTM
out <- gIntersection(x, mydemw_poly_gdal, byid=TRUE)
plot(out,col=rainbow(100))


ff = ISSR800.wcs(mydemw_poly, var, res = 800, quiet = FALSE)


# Get the data for a specific area
bboxpts = data.frame("X" = c(547847.47, 549011.30, 549115.76,550286.39), "Y" = c(4118538.37, 4117188.47,4118822.74,4117761.91))
pourpoint=SpatialPoints(bboxpts,proj4string = crs_utm) 

bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
res(mydem)
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")




# Write our raster to a geotiff file that can be used with
# OS level hydrological models 
writeRaster(mydem,filename = "mydem.tif",overwrite=T)

zoom(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")



old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
                        sep = ":"))

system("mpirun aread8")


library(raster)
library(shapefiles)

# Set working directory to your location
setwd("C:/Users/dtarb/Scratch/mydem")

z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)


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
zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
zoom(gord)

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
zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 100")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint, add = T)
zoom(src)
plot(pourpoint, add = T)
zoom(src)
plot(pourpoint, add = T)

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("qutlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

zoom(src)


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 100")
src1=raster("mydemsrc1.tif")
plot(src1)
zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"), add = T)
zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"))

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)

#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Deduce that the optimal threshold is 300 
# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 100")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


mydemw=raster("mydemw.tif")

mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasinmask)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
zoom(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#
# A series of plots to show all of the components
#
par(mfrow = c(2, 2))
plot(TI, main = "Stream Lab TI")
plot(TIC , main = "Stream Lab TIC")
plot(mybasinsca , main = "Stream Lab mybasinsca")
plot(mybasinslp, main = "Stream Lab mybasinslp")
dev.off()

plot(TIC,col=rainbow(5))
# Make a poly with raster library (slow)
# or from thee command line gdal (fast)
# gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp
mydemw_poly=rasterToPolygons(mydemw,dissolve = T,na.rm = T)
plot(mydemw_poly,add=T,border="red")
mydemw_poly
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)
# We will use this ESRI shape file, a zipped group of it, to download 
# our soil extent from the WebSoilSurvey Website
zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
# Download to your local machine mydemw.zip from the "Files" tab
# Open the WebSoilSurvey site to: 
browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Creat AOI from a zipped shapefile"
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/s4cjk1wdxfbo0i4n3eltrsbs/wss_aoi_2022-03-17_22-38-37.zip"
download.file(url,"wss_aoi.zip")
unzip("wss_aoi.zip")

url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_02080207_HU8_Shape.zip"
download.file(url,"NHD_H_02080207_HU8_Shape.zip")
unzip("NHD_H_02080207_HU8_Shape.zip",exdir="02080207")


streams=readOGR("02080207/Shape/NHDFlowline.dbf")
plot(streams)
mystream=subset(streams,gnis_id=="01478950")
plot(mystream)

mybbox=c(mystream@bbox)
mysoil = mapunit_geom_by_ll_bbox(mybbox)
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# THE ONLY DIFFERENCE BETWEEN LAST WEEKS LAB AND THE HW3 IS:
# We modify these following lines from our Lab4 exercise to access 
# the "resdept" from the "corestrictions" table as described in the PDF.
# So, from last week we modify this line:
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# To HW3 needs (here we are using "cr" for corestrictions where before 
# we used "ch" for "chorizon"... naming should reflect the data (though 
# many would use the full name so they can remember when read their code)
q_co2cr = paste("SELECT ksat_r,resdept_r FROM corestrictions WHERE cokey IN ", cokey_statement, sep="")
co2cr = SDA_query(q_co2ch)
View(co2cr)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2cr=merge(mu2co,co2cr)
mu2crmax=aggregate(mu2cr,list(mu2cr$mukey),max)
depth2restrict<-merge(mysoil,mu2crmax,by="mukey")
# And then a simple plot will do...
plot(mystream,col="red",lwd=4)
plot(depth2restrict,col=topo.colors(10),add=T)
plot(mystream,col="red",lwd=4,add=T)
# Though you want to explore nicer plots
pols1 <- list("sp.lines", as(mystream, 'SpatialLines'), col = "red", lwd = 4)
spplot(depth2restrict, sp.layout=list(pols1), zcol="ksat_r", xlab="Longitude", ylab="Latitude", main="Depth to Restrictive Layer", colorkey=T, col.regions=colorRampPalette(c("grey87","royalblue4"))(100), checkEmptyRC=T, add=T,xlim=mystream@bbox[1,],ylim=mystream@bbox[2,])
# or other colormap
spplot(depth2restrict, sp.layout=list(pols1), zcol="ksat_r", xlab="Longitude", ylab="Latitude", main="Depth to Restrictive Layer", colorkey=T, col.regions=topo.colors(100), checkEmptyRC=T, add=T,xlim=mystream@bbox[1,],ylim=mystream@bbox[2,])

