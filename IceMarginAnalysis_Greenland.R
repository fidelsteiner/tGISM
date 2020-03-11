################################################################################
# Spatial Analysis of the Ice Margin over Greenland
#
#
# ReadMe:
# 
#
# Created:          2018/02/05
# Latest Revision:  2017/02/05
#
# Jakob F Steiner| PhD candidate | Faculty of Geosciences | Universiteit Utrecht | Princetonlaan 8a, 3584 CB Utrecht 
# Vening Meinesz building, room 4.30 | P.O. Box 80.115, 3508 TC Utrecht | j.f.steiner@uu.nl | www.uu.nl/staff/jfsteiner | www.mountainhydrology.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# packages (if not installed yet: install.packages('examplePackage')
#install.packages('pacman')
library(pacman)
library(gstat)
p_load(rgdal,rgeos,maptools,raster,foreach,lubridate,compare,colorRamps,data.table,circular,parallel,snowfall,truncnorm,rlecuyer,forecast,rasterVis,R.utils,zoo,rlist,ggplot2)
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(ggplot2)
library(foreach)
library(ncdf4)


##################
# paths
##################

# Paths
path_outlines<-'F:\\OtherResearch\\Greenland\\PROMICE_database\\Glacier Outline'  # PROMICE
fnMarginOutline<-'PROMICdiss'
path_outlines_RGI<-'F:\\OtherResearch\\Greenland\\Outlines\\RGI\\05_rgi60_GreenlandPeriphery'                          # CCI
fnMarginOutline_RGI<-'05_rgi60_GreenlandPeriphery'
pathDEM_topo <- 'F:\\OtherResearch\\Greenland\\Bathymetry'
path_figs <- 'F:\\OtherResearch\\Greenland\\Figures'
path_mankoff <- 'F:\\OtherResearch\\Greenland\\PROMICE_database\\MankoffData\\sector_D.csv'
path_mouginot <- 'F:\\OtherResearch\\Greenland\\Mouginot\\'

path_sectorMargins <- 'F:\\OtherResearch\\Greenland\\Code\\SectorMargins'

fnOuterBuffer<-'Outer_Buffer2'

# Bedmachine data
bedmachineRAW <- 'F:\\OtherResearch\\Greenland\\Bedmachine\\5000000451838\\160281892\\BedMachineGreenland-2017-09-20.nc'

# basin data
GLbasins <- 'F:\\OtherResearch\\Greenland\\Mouginot\\drainagebasins\\GRE_Basins_IMBIE2_v1.3.shp'
# sector data
GLsectors <- 'F:\\OtherResearch\\Greenland\\Mouginot\\drainagebasins\\MouginotRignot\\doi_10.7280_D1WT11__v1\\Greenland_Basins_PS_v1.4.2.shp'

#read PROMICE margin  including all ice caps
ogrInfo(path_outlines,fnMarginOutline)
margin_pEXT<-readOGR(dsn=path_outlines,layer=fnMarginOutline)
margin_pEXT<-SpatialPolygons(margin_pEXT@polygons,proj4string=margin_pEXT@proj4string)
projection(margin_pEXT)<-CRS("+init=epsg:4326")
margin_pEXT_transf<-spTransform(margin_pEXT,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

margin_actual <- margin_pEXT_transf[1]
margin_glaciers <- margin_pEXT_transf[2]
margin_peary <- margin_pEXT_transf[3]
margin_icecaps <- margin_pEXT_transf[4]

margin_line <- as(margin_actual,'SpatialLines')
#margin_line_buffered <- gBuffer(margin_line,width=1000,byid=T)

#read Rastner margin (Rastner2012)
ogrInfo(path_outlines_RGI,fnMarginOutline_RGI)
margin_RGI<-readOGR(dsn=path_outlines_RGI,layer=fnMarginOutline_RGI)
margin_RGI<-SpatialPolygons(margin_RGI@polygons,proj4string=margin_RGI@proj4string)
projection(margin_RGI)<-CRS("+init=epsg:4326")
margin_RGI_transf<-spTransform(margin_RGI,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

##################
#################
###################
#CONTINUE HERE TO REMOVE THE RASTER DATA FROM THE PROMICE DATA
rgeos::gDifference

# read basin outlines
ogrInfo(GLbasins)
basins_pEXT<-readOGR(dsn=GLbasins)
basins_pEXT<-SpatialPolygons(basins_pEXT@polygons,proj4string=basins_pEXT@proj4string)
projection(basins_pEXT)<-CRS("+init=epsg:4326")
basins_pEXT_transf<-spTransform(basins_pEXT,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

# read sector outlines (Mouginot2019)
ogrInfo(GLsectors)
sectors_pEXT<-readOGR(dsn=GLsectors)
sectors_pEXT<-SpatialPolygons(sectors_pEXT@polygons,proj4string=sectors_pEXT@proj4string)
projection(sectors_pEXT)<-CRS("+init=epsg:4326")
#sectors_pEXT_transf<-spTransform(sectors_pEXT,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
projection(sectors_pEXT)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

# Read Bedmachine (Citation!)
bedmachine_bed <- raster(bedmachineRAW, varname ='bed')
bedmachine_bed[bedmachine_bed<=-1000]<-NA
projection(bedmachine_bed) <- '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

coarse_bed <- aggregate(bedmachine_bed,fact=1)

################### 
# read and process velocities (Mouginot2019)
################### 
ncin_x <- raster(path_mouginot&'Velocities\\2010_2017\\vel_2016-07-01_2017-06-31.nc',varname='VX')
ncin_y <- raster(path_mouginot&'Velocities\\2010_2017\\vel_2016-07-01_2017-06-31.nc',varname='VY')
velmap <- sqrt(ncin_x^2 + ncin_y^2)
projection(velmap)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
projection(ncin_x)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
projection(ncin_y)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
velmap[velmap[]<0] <- NA
velmap_bed <- raster::resample(velmap,coarse_bed)
xvel <- raster::resample(ncin_x,coarse_bed)
yvel <- raster::resample(ncin_y,coarse_bed)

dir_v <- atan2(yvel,xvel) * 180 / pi
velmap_bed[coarse_bed[]<0]<-NA
xvel[coarse_bed[]<0]<-NA
yvel[coarse_bed[]<0]<-NA
dir_v[coarse_bed[]<0]<-NA

writeRaster(xvel, path_outlines2&'\\dumper\\xvel.tif','GTiff',overwrite=T)
writeRaster(yvel, path_outlines2&'\\dumper\\yvel.tif','GTiff',overwrite=T)
writeRaster(dir_v, path_outlines2&'\\dumper\\dir_v.tif','GTiff',overwrite=T)
writeRaster(velmap_bed, path_outlines2&'\\dumper\\mag.tif','GTiff',overwrite=T)

# extract velocities along margin
velmap_coarse <- aggregate(velmap,fact= 20)
bedmachine_coarse <- aggregate(bedmachine_bed,fact=20)
velmap_bed <- raster::resample(bedmachine_coarse,velmap_coarse)
margin_line <- as(margin_actual,'SpatialLines')
velMarg <- rasterize(margin_line,velmap_bed)
velMarg_holed <- boundaries(velMarg, type='inner', classes=FALSE, directions=8,asNA=T)

velMarg_holed[which(velMarg_holed[]>0&!is.na(velMarg_holed[]))] <- 2
velMarg_holed[which(bedmachine_coarse[]<=0&!is.na(velMarg_holed[]))] <- 3



# Read Mankoff Fluxes
s.sf <- st_read(GLsectors)
s.sf$NAME
mankoffFlux <- read.csv(path_mankoff,sep=',',header=T)

fluxVal <- vector()
for(i in 1:length(s.sf$NAME)){
  
  if(c(levels(s.sf$NAME)[i]) %in% colnames(mankoffFlux)==TRUE){
    fluxVal[i] <- mean(mankoffFlux[,c(levels(s.sf$NAME)[i])])
  }
  else{
    fluxVal[i] <- NA
  }
}


# Read Mankoff Fluxes
p<-sectors_pEXT

( pid <- sapply(slot(p, "polygons"), function(x) slot(x, "ID")) )
( p.df <- data.frame( ID=1:length(p), row.names = pid) )  
p <- SpatialPolygonsDataFrame(p, p.df)
p@data$flux <- fluxVal
colVals <- round(p@data$flux*length(pal(100)))
colVals[is.na(colVals)]<-1
colRamp <- c('#ffffff',col=pal(99))

png(file=path_figs&'\\iceDischarge.png', res = 160,width=1800,height=1800)
par(mar=c(2,2,2,2),cex.lab=1.2,cex.axis=1.2)
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(p,col=colRamp[colVals])
legend_image <- as.raster(matrix(colRamp, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'discharge [Gt yr-1]')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,38,l=5))
rasterImage(rev(legend_image[0:100]), 0, 0, 1,1)
dev.off()



# Test one subbasin
ltermFrac <- vector()
for(i in 1:length(sectors_pEXT)){
sector_line <- as(sectors_pEXT[i],'SpatialLines')
xx <- gBuffer(sector_line,width=2000,byid=T)
op <- crop(margin_line,extent(sector_line))

if(!is.null(op)){
margin_Sector <- gIntersection(xx,op)

df<-SpatialLinesDataFrame(margin_Sector, data.frame(id=1:length(margin_Sector)))
writeOGR(df, dsn=path_sectorMargins&'\\marginSector'&i&'.shp' ,layer="id",driver="ESRI Shapefile")

opt <- rasterize(margin_Sector,coarse_bed)
opt[which(coarse_bed[]>0&!is.na(opt[]))] <- 2
opt[which(coarse_bed[]<=0&!is.na(opt[]))] <- 3
opt[which(opt[]==1)] <- NA
writeRaster(opt, path_sectorMargins&'\\marginSector'&i&'.tif','GTiff',overwrite=T)

#marginasPoints <- rasterToPoints(opt)

ltermFrac[i] <- length(which(opt[]==2)) / length(which(opt[]>1))
}
print(i)
}

p<-sectors_pEXT

( pid <- sapply(slot(p, "polygons"), function(x) slot(x, "ID")) )
( p.df <- data.frame( ID=1:length(p), row.names = pid) )  
p <- SpatialPolygonsDataFrame(p, p.df)
ltermFrac2<- seq(1,260,1)*0 + 1
ltermFrac2[1:length(ltermFrac)] <- ltermFrac
p@data$frac <- ltermFrac2
pal <- colorRampPalette(c("white", "red"))

colRamp <- c('#ffffff',col=pal(99))
colVals <- round(p@data$frac*length(pal(100)))


png(file=path_figs&'\\relativeTerrestrialMargin.png', res = 160,width=1800,height=1800)
par(mar=c(2,2,2,2),cex.lab=1.2,cex.axis=1.2)
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(p,col=colRamp[colVals])
legend_image <- as.raster(matrix(colRamp, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'terrestrial fraction')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0.3,1,l=5))
rasterImage(rev(legend_image[30:100]), 0, 0, 1,1)
dev.off()

p@data$flux <- fluxVal
colVals <- round(p@data$flux*length(pal(100)))
colVals[is.na(colVals)]<-1

pal <- terrain.colors(1000)

png(file=path_figs&'\\IceMargin.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(bedmachine_bed,xlab = 'Westing [m]', ylab ='Northing [m]',legend.args=list(text=expression('Elevation [m]')),col=pal)
plot(margin_line,add=TRUE)
grid(NULL,NULL)
dev.off()

# Rasterize Margin to Bedmachine
rasterized_margin_icesheet <- rasterize(margin_actual,coarse_bed)
# Remove the inside cells (ice sheet)
rasterized_margin_holed <- rasterized_margin_icesheet
rasterized_margin_holed <- boundaries(rasterized_margin_icesheet, type='inner', classes=FALSE, directions=8,asNA=T)
writeRaster(rasterized_margin, path_outlines2&'\\RasterizedMargin_icesheet_promicee150m.tif','GTiff')


rasterized_margin_holed[which(rasterized_margin_holed[]==0)] <- NA
rasterized_margin_holed[which(coarse_bed[]>0&!is.na(rasterized_margin_holed[]))] <- 2
rasterized_margin_holed[which(coarse_bed[]<=0&!is.na(rasterized_margin_holed[]))] <- 3
rasterized_margin_holed[which(rasterized_margin_holed[]==1)] <- NA
marginasPoints <- rasterToPoints(rasterized_margin_holed)

spdf <- SpatialPointsDataFrame(marginasPoints[,1:2],as.data.frame(marginasPoints))
projection(spdf) <- '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

writeOGR(obj=spdf, dsn=path_outlines2&"\\margin_icesheet_150m.shp", layer='margin', driver="ESRI Shapefile") # this is in geographical projection

percLand <- 100 * sum(rasterized_margin_holed[] == 2,na.rm=T) / (sum(rasterized_margin_holed[] == 2,na.rm=T) + sum(rasterized_margin_holed[] == 3,na.rm=T))
percLand_icecaps <- 100 * sum(rasterized_margin_holed_promice[] == 2,na.rm=T) / (sum(rasterized_margin_holed_promice[] == 2,na.rm=T) + sum(rasterized_margin_holed_promice[] == 3,na.rm=T))

rbPal <- colorRampPalette(c('red','blue'))
col <- rbPal(2)[as.numeric(cut(marginasPoints[,3],breaks = 2))]
png(file=path_figs&'\\IceMargin_landocean.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(coarse_bed,xlab = 'Westing [m]', ylab ='Northing [m]',legend.args=list(text=expression('Elevation [m]')),col=pal)
points(marginasPoints,col=col)
#plot(rasterized_margin_holed,add=TRUE,legend=F,col=rainbow(2,alpha=0.8))
#plot(margin_pEXT_transf,add=T)
grid(NULL,NULL)
dev.off()

col <- rbPal(2)[as.numeric(cut(marginasPoints_promice[,3],breaks = 2))]
png(file=path_figs&'\\IceMargin_landocean_plusicecaps.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(coarse_bed,xlab = 'Westing [m]', ylab ='Northing [m]',legend.args=list(text=expression('Elevation [m]')),col=pal)
points(marginasPoints_promice,add=T,col=col)
#plot(rasterized_margin_holed,add=TRUE,legend=F,col=rainbow(2,alpha=0.8))
plot(margin_pEXT_transf,add=T)
grid(NULL,NULL)
dev.off()

#levelplot(op2)

#bin <- matrix(0,11,2)
#startK<-0
#for(iExtent in 1:length(seq(extent(op2)[1],extent(op2)[2],150000))){
#  kount <- 150000
#e <- extent(extent(op2)[1] + startK,extent(op2)[1] + startK + kount,-3392600,-632600)
#startK <- startK + kount + 1
#allCells <- extract(op2, e)
#bin[iExtent,] <- cbind(sum(allCells == 2,na.rm=T),sum(allCells == 3,na.rm=T))
#}

#barplot(bin[,1]/rowSums(bin)*100)