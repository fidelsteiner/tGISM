################################################################################
# Extract Ice Margin Morphology for separate segments
#
# 
# IceMarginExtractor.R
#
# ReadMe: 
#
# Input:

# 
#
# Created:          2021/10/15
# Latest Revision:  2021/10/05
#
# Jakob F Steiner| ICIMOD | jakob.steiner@icimod.org | x-hydrolab.org 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')


##################
#to change
##################
GlacierName <- 'RedRock'

fnDEM_1 <- 'PGO_2017-08-15_RedRockCliff_DEM_2m.tif'
fnDEM_2 <- 'c0_reg_dem_ProjectRaster1.tif'


HiawathaDEM_1 <- '29_36_1_1_2m_v3.0_reg_dem_proj.tif'
HiawathaDEM_2 <- '29_36_1_2_2m_v3.0_reg_dem_proj.tif'
HiawathaDEM_3 <- '29_36_2_2_2m_v3.0_reg_dem_proj.tif'
HiawathaDEM_4 <- '29_35_2_1_2m_v3.0_reg_dem_proj.tif'
HiawathaDEM_5 <- '29_35_2_2_2m_v3.0_reg_dem_proj.tif'
HiawathaDEM_6 <- '29_35_1_2_2m_v3.0_reg_dem_proj.tif'

RussellDEM_1 <- '16_38_2_1_2m_v3.0_reg_dem_proj.tif'
RussellDEM_2 <- '15_38_2_2_2m_v3.0_reg_dem_proj.tif'

KronprinsDEM_1 <- '31_43_2_2_2m_v3.0_reg_dem_proj.tif'
KronprinsDEM_2 <- '31_43_2_1_2m_v3.0_reg_dem_proj.tif'

fnMarginOutline <- 'RedRock_ManualBuffer'
fnMarginOutlineGrIS <- 'RedRockGrIS_ManualBuffer'

fnMarginHiawatha <- 'marginSector176_buffer'
fnMarginRussell <- 'marginSector71_buffer'
fnMarginKronprins <- 'marginSector57_buffer'
##################
#to change
##################

# Paths
path_outlines<-'E:\\Research\\OtherResearch\\Greenland\\Code\\SectorMargins\\ManualRedRock'
path_figs<-'E:\\Research\\OtherResearch\\Greenland\\Figures'
pathDEM_Pleiades<-'E:\\Research\\OtherResearch\\Greenland\\Pleiades\\DEM\\2017-08-15_RedRockCliff\\PGO\\2017-08-15_RedRockCliff'
pathDEM_Arctic<-'D:\\Work\\Research\\tGISM\\DEMData\\ArcticDEM'

projec<-'+proj=utm +datum=WGS84'

projec_utm<-'+proj=utm +zone=19N +datum=WGS84'
projec_27N<-'+proj=utm +zone=27N +datum=WGS84'

# packages (if not installed yet: install.packages('examplePackage')
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(ggplot2)
library(dplyr)
library(ggridges)
library(SpaDES)
library(RColorBrewer)

# read Margin Outlines
ogrInfo(path_outlines,fnMarginOutline)
margin_pEXT<-readOGR(dsn=path_outlines,layer=fnMarginOutline)
#projection(margin_pEXT)<-projec
margin_pEXT<-SpatialPolygons(margin_pEXT@polygons,proj4string=margin_pEXT@proj4string)

ogrInfo(path_outlines,fnMarginOutlineGrIS)
margin_pEXTGrIS<-readOGR(dsn=path_outlines,layer=fnMarginOutlineGrIS)
#projection(margin_pEXT)<-projec
margin_pEXTGrIS<-SpatialPolygons(margin_pEXTGrIS@polygons,proj4string=margin_pEXTGrIS@proj4string)

#######################
# Read DEM Data
#######################
DEM_temp <- raster(pathDEM_Pleiades&'/'&fnDEM_1)
projection(DEM_temp) <- projec_utm

DEM_Arctic <- raster(pathDEM_Arctic&'\\27_35_2_2_2m_v3.0'&'/'&fnDEM_2)
projection(DEM_Arctic) <- projec_utm


DEM_masked <- mask(crop(DEM_temp,margin_pEXT),margin_pEXT)
DEM_maskedGrIS <- mask(crop(DEM_temp,margin_pEXTGrIS),margin_pEXTGrIS)

DEMArctic_masked <- mask(crop(DEM_Arctic,margin_pEXT),margin_pEXT)
DEMArctic_maskedGrIS <- mask(crop(DEM_Arctic,margin_pEXTGrIS),margin_pEXTGrIS)

RR_slopeGrIS <- terrain(DEM_maskedGrIS, opt="slope", unit="degrees", neighbors=8)
RR_slope <- terrain(DEM_masked, opt="slope", unit="degrees", neighbors=8)

RR_slopeArcticGrIS <- terrain(DEMArctic_maskedGrIS, opt="slope", unit="degrees", neighbors=8)
RR_slopeArctic <- terrain(DEMArctic_masked, opt="slope", unit="degrees", neighbors=8)


set.seed(1234)
wdata_IceCap = data.frame(
  DEM = factor(rep(c("ArcticDEM", "Pleiades"), each=length(which(!is.na(RR_slope[]))))),
  slope = c(RR_slopeArctic[which(!is.na(RR_slope[]))], RR_slope[which(!is.na(RR_slope[]))])
)
head(wdata_IceCap, 4)

mu <- wdata_IceCap %>% 
  group_by(DEM) %>%
  summarise(grp.mean = median(slope,na.rm=T))

png(file=path_figs&'\\RedRock_ICeCap_DEMmatch.png', res = 300,width=3600,height=900)
par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
ggplot(wdata_IceCap, aes(x = slope, y = DEM)) +
  geom_density_ridges(aes(fill = DEM)) +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) +
  geom_vline(aes(xintercept = grp.mean, color = DEM),
             data = mu, linetype = "dashed") +
  scale_x_continuous(name=expression('Slope [?]'), limits=c(0, 80)) +
  theme(axis.text.x = element_text(size=14, angle=90),
        axis.text.y = element_text(size=14, angle=0)) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank()) +
  theme(legend.text=element_text(size=rel(1)))
theme_bw()
dev.off()


set.seed(1234)
wdata_IceSheet = data.frame(
  DEM = factor(rep(c("ArcticDEM", "Pleiades"), each=length(which(!is.na(RR_slopeGrIS[]))))),
  slope = c(RR_slopeArcticGrIS[which(!is.na(RR_slopeGrIS[]))], RR_slopeGrIS[which(!is.na(RR_slopeGrIS[]))])
)
head(wdata_IceSheet, 4)

mu <- wdata_IceSheet %>% 
  group_by(DEM) %>%
  summarise(grp.mean = median(slope,na.rm=T))

png(file=path_figs&'\\RedRock_IceSheet_DEMmatch.png', res = 300,width=3600,height=900)
par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
ggplot(wdata_IceSheet, aes(x = slope, y = DEM)) +
  geom_density_ridges(aes(fill = DEM)) +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) +
  geom_vline(aes(xintercept = grp.mean, color = DEM),
             data = mu, linetype = "dashed") +
  scale_x_continuous(name=expression('Slope [?]'), limits=c(0, 80)) +
  theme(axis.text.x = element_text(size=14, angle=90),
        axis.text.y = element_text(size=14, angle=0)) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank()) +
  theme(legend.text=element_text(size=rel(1)))
theme_bw()
dev.off()

#######################
# Read Hiawatha Site
#######################

ogrInfo('D:\\Work\\Research\\tGISM\\Code','marginSector176_buffer')
margin_HIA<-readOGR(dsn='D:\\Work\\Research\\tGISM\\Code',layer='marginSector176_buffer')

margin_HIA<-SpatialPolygons(margin_HIA@polygons,proj4string=margin_HIA@proj4string)
margin_HIA <- spTransform(margin_HIA, CRS('+proj=utm +zone=19N +datum=WGS84'))

DEM_Hiawatha1 <- raster(pathDEM_Arctic&'\\29_36_1_1_2m_v3.0'&'/'&HiawathaDEM_1)
projection(DEM_Hiawatha1) <- projec_utm

DEM_Hiawatha2 <- raster(pathDEM_Arctic&'\\29_36_1_2_2m_v3.0'&'/'&HiawathaDEM_2)
projection(DEM_Hiawatha2) <- projec_utm

DEM_Hiawatha3 <- raster(pathDEM_Arctic&'\\29_36_2_2_2m_v3.0'&'/'&HiawathaDEM_3)
projection(DEM_Hiawatha3) <- projec_utm

DEM_Hiawatha4 <- raster(pathDEM_Arctic&'\\29_35_2_1_2m_v3.0'&'/'&HiawathaDEM_4)
projection(DEM_Hiawatha4) <- projec_utm

DEM_Hiawatha5 <- raster(pathDEM_Arctic&'\\29_35_2_2_2m_v3.0'&'/'&HiawathaDEM_5)
projection(DEM_Hiawatha5) <- projec_utm

DEM_Hiawatha6 <- raster(pathDEM_Arctic&'\\29_35_1_2_2m_v3.0'&'/'&HiawathaDEM_6)
projection(DEM_Hiawatha6) <- projec_utm


DEM_HIA1 <- mask(crop(DEM_Hiawatha1,crop(margin_HIA,extent(DEM_Hiawatha1))),crop(margin_HIA,extent(DEM_Hiawatha1)))
HIA_slope1 <- terrain(DEM_HIA1, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_HIA1)
file.remove(c(f, extension(f, '.gri')))
DEM_HIA2 <- mask(crop(DEM_Hiawatha2,crop(margin_HIA,extent(DEM_Hiawatha2))),crop(margin_HIA,extent(DEM_Hiawatha2)))
HIA_slope2 <- terrain(DEM_HIA2, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_HIA2)
file.remove(c(f, extension(f, '.gri')))
DEM_HIA3 <- mask(crop(DEM_Hiawatha3,crop(margin_HIA,extent(DEM_Hiawatha3))),crop(margin_HIA,extent(DEM_Hiawatha3)))
HIA_slope3 <- terrain(DEM_HIA3, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_HIA3)
file.remove(c(f, extension(f, '.gri')))
DEM_HIA4 <- mask(crop(DEM_Hiawatha4,crop(margin_HIA,extent(DEM_Hiawatha4))),crop(margin_HIA,extent(DEM_Hiawatha4)))
HIA_slope4 <- terrain(DEM_HIA4, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_HIA4)
file.remove(c(f, extension(f, '.gri')))
DEM_HIA5 <- mask(crop(DEM_Hiawatha5,crop(margin_HIA,extent(DEM_Hiawatha5))),crop(margin_HIA,extent(DEM_Hiawatha5)))
HIA_slope5 <- terrain(DEM_HIA5, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_HIA5)
file.remove(c(f, extension(f, '.gri')))
DEM_HIA6 <- mask(crop(DEM_Hiawatha6,crop(margin_HIA,extent(DEM_Hiawatha6))),crop(margin_HIA,extent(DEM_Hiawatha6)))
HIA_slope6 <- terrain(DEM_HIA6, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_HIA6)
file.remove(c(f, extension(f, '.gri')))
gc()
rm(DEM_HIA1, DEM_HIA2, DEM_HIA3, DEM_HIA4, DEM_HIA5, DEM_HIA6)
rm(DEM_Hiawatha1, DEM_Hiawatha2, DEM_Hiawatha3, DEM_Hiawatha4,DEM_Hiawatha5,DEM_Hiawatha6)

#######################
# Split up margin
#######################

grid <- raster(extent(margin_HIA), resolution = c(10000,10000), crs = proj4string(margin_HIA))
grid <- raster::extend(grid, c(1,1))
gridPolygon <- rasterToPolygons(grid)
plot(margin_HIA)
plot(gridPolygon, add = T)

gridPolygon$id <- 1:nrow(gridPolygon)

intersectGridClipped <- raster::intersect(gridPolygon, margin_HIA)
intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]
rgdal::writeOGR(intersectGridClipped,'D:\\Work\\Research\\tGISM\\Code\\marginSector176_clipped',layer = 'intersectGridClipped',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(intersectGrid,'D:\\Work\\Research\\tGISM\\Code\\marginSector176_Grid',layer = 'intersectGrid',driver = 'ESRI Shapefile',overwrite_layer = T)

slopeVals_segment <- list()
steepSlopes <- matrix(NA, nrow = 1,ncol = 2)

for(k in 1:dim(intersectGridClipped)[1]){
  save.image(file="temp.RData")
  rm(list=ls())
  load(file="temp.RData")
  gc()
  slopeVals_segment[[k]] <- NA
if(tryCatch(!is.null(crop(HIA_slope1,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
  segmentMargin <- mask(crop(HIA_slope1,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
  segmentMargin[segmentMargin[]<5]<-NA
  
  slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
  steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])
} 
if(tryCatch(!is.null(crop(HIA_slope2,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(HIA_slope2,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
    
    slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
    steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])} 
if(tryCatch(!is.null(crop(HIA_slope3,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(HIA_slope3,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
    
  slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
  steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])} 
if(tryCatch(!is.null(crop(HIA_slope5,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(HIA_slope5,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
    
  slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
  steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])
}
if(tryCatch(!is.null(crop(HIA_slope4,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(HIA_slope4,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
    
    slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
    steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])
  } 
if(tryCatch(!is.null(crop(HIA_slope6,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(HIA_slope6,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
    
  slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
  steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])
  } 
  
  # Save slope points >33deg

}

df_steepcoord <- data.frame(steepSlopes
)

write.csv(df_steepcoord,'D:\\Work\\Research\\tGISM\\Code\\steeppixels_hiawatha.csv', row.names = FALSE)

do.call(file.remove, list(list.files(tmpDir(), full.names = TRUE)))

slopeclass <- matrix(NA, nrow = dim(intersectGridClipped)[1],ncol = 4)
for(l in 1:dim(intersectGridClipped)[1]){
slopeclass[l,1] <- length(which(slopeVals_segment[[l]][]<10))/length(slopeVals_segment[[l]][])
slopeclass[l,2] <- length(which(slopeVals_segment[[l]][]>=10 & slopeVals_segment[[l]][]<25))/length(slopeVals_segment[[l]][])
slopeclass[l,3] <- length(which(slopeVals_segment[[l]][]>=25 & slopeVals_segment[[l]][]<40))/length(slopeVals_segment[[l]][])
slopeclass[l,4] <- length(which(slopeVals_segment[[l]][]>=40))/length(slopeVals_segment[[l]][])
}

png(file=path_figs&'\\Hiawatha_slopeclasses.png', res = 300,width=2000,height=900) 
barplot(t(slopeclass),col=brewer.pal(n = 4, name = "Reds"),cex.names=1,
        names.arg=seq(1,29,1),ylim=c(0,1))
legend("top", fill = brewer.pal(n = 4, name = "Reds"), legend = c('<10','10-25','25-40','>40'), 
       horiz = TRUE, inset = c(0,-0.5), xpd = TRUE)
dev.off()

#######################
# Make Figures
#######################

  library(RColorBrewer)
  set.seed(1234)

    wdata_Hiawatha = data.frame(
    location = factor(rep(c("S01","S09","S11","S13","S18","S21","S19"),c(length(slopeVals_segment[[1]]),length(slopeVals_segment[[9]]),length(slopeVals_segment[[11]]),length(slopeVals_segment[[13]]),length(slopeVals_segment[[18]]),length(slopeVals_segment[[21]]),length(slopeVals_segment[[19]])))),
    slope = c(slopeVals_segment[[1]],slopeVals_segment[[9]],slopeVals_segment[[11]],slopeVals_segment[[13]],slopeVals_segment[[18]],slopeVals_segment[[21]],slopeVals_segment[[19]])
  )
  head(wdata_Hiawatha, 4)

  mu <- wdata_Hiawatha%>% 
   # arrange(order)  %>% 
    group_by(location) %>%
    summarise(grp.mean = median(slope,na.rm=T))

  png(file=path_figs&'\\HiawathaSlope distributions.png', res = 300,width=2000,height=900)  
  par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
  ggplot(wdata_Hiawatha, aes(x = slope, y = location)) +
    geom_density_ridges(aes(fill = location)) +
    scale_fill_manual(values = brewer.pal(n = 7, name = "Dark2")) +
    geom_vline(aes(xintercept = grp.mean, color = location),
               data = mu, linetype = "dashed") +
    scale_x_continuous(name=expression('Slope [?]'), limits=c(0, 70)) +
    theme(axis.text.x = element_text(size=14, angle=90),
          axis.text.y = element_text(size=14, angle=0)) +
    theme(axis.title.y =element_blank(),
          axis.text.y=element_blank()) +
    theme(legend.text=element_text(size=rel(1)))
  theme_bw()
  dev.off()
  
  
#######################
# Read Russell Site
#######################
  
ogrInfo('D:\\Work\\Research\\tGISM\\Code','marginSector71_buffer')
margin_RUS<-readOGR(dsn='D:\\Work\\Research\\tGISM\\Code',layer='marginSector71_buffer')
  
margin_RUS<-SpatialPolygons(margin_RUS@polygons,proj4string=margin_RUS@proj4string)
margin_RUS <- spTransform(margin_RUS, CRS('+proj=utm +zone=19N +datum=WGS84'))
  
DEM_Russell1 <- raster(pathDEM_Arctic&'\\16_38_2_1_2m_v3.0'&'/'&RussellDEM_1)
projection(DEM_Russell1) <- projec_utm
  
DEM_Russell2 <- raster(pathDEM_Arctic&'\\15_38_2_2_2m_v3.0'&'/'&RussellDEM_2)
projection(DEM_Russell1) <- projec_utm
  
DEM_RUS1 <- mask(crop(DEM_Russell1,crop(margin_RUS,extent(DEM_Russell1))),crop(margin_RUS,extent(DEM_Russell1)))
RUS_slope1 <- terrain(DEM_RUS1, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_RUS1)
file.remove(c(f, extension(f, '.gri')))
DEM_RUS2 <- mask(crop(DEM_Russell2,crop(margin_RUS,extent(DEM_Russell2))),crop(margin_RUS,extent(DEM_Russell2)))
RUS_slope2 <- terrain(DEM_RUS2, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_RUS2)
file.remove(c(f, extension(f, '.gri')))

gc()
rm(DEM_Russell1, DEM_Russell2)
rm(DEM_RUS1, DEM_RUS2)
  
#######################
# Split up margin
#######################
  
grid <- raster(extent(margin_RUS), resolution = c(10000,10000), crs = proj4string(margin_RUS))
grid <- raster::extend(grid, c(1,1))
gridPolygon <- rasterToPolygons(grid)
plot(margin_RUS)
plot(gridPolygon, add = T)
  
gridPolygon$id <- 1:nrow(gridPolygon)
  
intersectGridClipped <- raster::intersect(gridPolygon, margin_RUS)
intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]
rgdal::writeOGR(intersectGridClipped,'D:\\Work\\Research\\tGISM\\Code\\marginSector71_clipped',layer = 'intersectGridClipped',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(intersectGrid,'D:\\Work\\Research\\tGISM\\Code\\marginSector71_Grid',layer = 'intersectGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
  
slopeVals_segment <- list()
steepSlopes <- matrix(NA, nrow = 1,ncol = 2)
  
for(k in 1:dim(intersectGridClipped)[1]){
  save.image(file="temp.RData")
  rm(list=ls())
  load(file="temp.RData")
  gc()
  slopeVals_segment[[k]] <- NA
  if(tryCatch(!is.null(crop(RUS_slope1,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(RUS_slope1,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
      
    slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
    steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])
  } 
  if(tryCatch(!is.null(crop(RUS_slope2,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(RUS_slope2,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
      
    slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
    steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])} 

    
}
  
df_steepcoord <- data.frame(steepSlopes
)
  
write.csv(df_steepcoord,'D:\\Work\\Research\\tGISM\\Code\\steeppixels_russell.csv', row.names = FALSE)
  
do.call(file.remove, list(list.files(tmpDir(), full.names = TRUE)))
  
slopeclass <- matrix(NA, nrow = dim(intersectGridClipped)[1],ncol = 4)
for(l in 1:dim(intersectGridClipped)[1]){
  slopeclass[l,1] <- length(which(slopeVals_segment[[l]][]<10))/length(slopeVals_segment[[l]][])
  slopeclass[l,2] <- length(which(slopeVals_segment[[l]][]>=10 & slopeVals_segment[[l]][]<25))/length(slopeVals_segment[[l]][])
  slopeclass[l,3] <- length(which(slopeVals_segment[[l]][]>=25 & slopeVals_segment[[l]][]<40))/length(slopeVals_segment[[l]][])
  slopeclass[l,4] <- length(which(slopeVals_segment[[l]][]>=40))/length(slopeVals_segment[[l]][])
}
  
png(file=path_figs&'\\Russell_slopeclasses.png', res = 300,width=2000,height=900) 
barplot(t(slopeclass),col=brewer.pal(n = 4, name = "Reds"),cex.names=1,
        names.arg=seq(1,dim(intersectGridClipped)[1],1),ylim=c(0,1))
legend("top", fill = brewer.pal(n = 4, name = "Reds"), legend = c('<10','10-25','25-40','>40'), 
       horiz = TRUE, inset = c(0,-0.5), xpd = TRUE)
dev.off()
  
#######################
# Make Figures
#######################
  
set.seed(1234)
  
wdata_Russell = data.frame(
  location = factor(rep(c("S01","S02","S04","S05","S08"),c(length(slopeVals_segment[[1]]),length(slopeVals_segment[[2]]),length(slopeVals_segment[[4]]),length(slopeVals_segment[[5]]),length(slopeVals_segment[[8]])))),
  slope = c(slopeVals_segment[[1]],slopeVals_segment[[2]],slopeVals_segment[[4]],slopeVals_segment[[5]],slopeVals_segment[[8]])
)
head(wdata_Russell, 4)
  
mu <- wdata_Russell%>% 
  # arrange(order)  %>% 
  group_by(location) %>%
  summarise(grp.mean = median(slope,na.rm=T))
  
png(file=path_figs&'\\RussellSlope distributions.png', res = 300,width=2000,height=900)  
par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
ggplot(wdata_Russell, aes(x = slope, y = location)) +
  geom_density_ridges(aes(fill = location)) +
  scale_fill_manual(values = brewer.pal(n = 5, name = "Dark2")) +
  geom_vline(aes(xintercept = grp.mean, color = location),
             data = mu, linetype = "dashed") +
  scale_x_continuous(name=expression('Slope [?]'), limits=c(0, 70)) +
  theme(axis.text.x = element_text(size=14, angle=90),
        axis.text.y = element_text(size=14, angle=0)) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank()) +
  theme(legend.text=element_text(size=rel(1)))
theme_bw()
dev.off()
  
#######################
# Read Kronprins Site
#######################
  
ogrInfo('D:\\Work\\Research\\tGISM\\Code','marginSector57_buffer')
margin_KRO<-readOGR(dsn='D:\\Work\\Research\\tGISM\\Code',layer='marginSector57_buffer')
  
margin_KRO <-SpatialPolygons(margin_KRO@polygons,proj4string=margin_KRO@proj4string)
margin_KRO <- spTransform(margin_KRO, CRS('+proj=utm +zone=27N +datum=WGS84'))
  
DEM_Kron1 <- raster(pathDEM_Arctic&'\\31_43_2_2_2m_v3.0'&'/'&KronprinsDEM_1)
projection(DEM_Kron1) <- projec_27N
  
DEM_Kron2 <- raster(pathDEM_Arctic&'\\31_43_2_1_2m_v3.0'&'/'&KronprinsDEM_2)
projection(DEM_Kron2) <- projec_27N
  
DEM_KRO1 <- mask(crop(DEM_Kron1,crop(margin_KRO,extent(DEM_Kron1))),crop(margin_KRO,extent(DEM_Kron1)))
KRO_slope1 <- terrain(DEM_KRO1, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_KRO1)
file.remove(c(f, extension(f, '.gri')))
DEM_KRO2 <- mask(crop(DEM_Kron2,crop(margin_KRO,extent(DEM_Kron2))),crop(margin_KRO,extent(DEM_Kron2)))
KRO_slope2 <- terrain(DEM_KRO2, opt="slope", unit="degrees", neighbors=8)
f <- filename(DEM_KRO2)
file.remove(c(f, extension(f, '.gri')))
  
gc()
rm(DEM_Kron1, DEM_Kron2)
rm(DEM_KRO1, DEM_KRO2)
  
#######################
# Split up margin
#######################
  
grid <- raster(extent(margin_KRO), resolution = c(10000,10000), crs = proj4string(margin_KRO))
grid <- raster::extend(grid, c(1,1))
gridPolygon <- rasterToPolygons(grid)
plot(margin_KRO)
plot(gridPolygon, add = T)
  
gridPolygon$id <- 1:nrow(gridPolygon)
  
intersectGridClipped <- raster::intersect(gridPolygon, margin_KRO)
intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]
rgdal::writeOGR(intersectGridClipped,'D:\\Work\\Research\\tGISM\\Code\\marginSector57_clipped',layer = 'intersectGridClipped',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(intersectGrid,'D:\\Work\\Research\\tGISM\\Code\\marginSector57_Grid',layer = 'intersectGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
  
slopeVals_segment <- list()
steepSlopes <- matrix(NA, nrow = 1,ncol = 2)
  
for(k in 1:dim(intersectGridClipped)[1]){
  save.image(file="temp.RData")
  rm(list=ls())
  load(file="temp.RData")
  gc()
  slopeVals_segment[[k]] <- NA
  if(tryCatch(!is.null(crop(KRO_slope1,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(KRO_slope1,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
      
    slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
    steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])
  } 
  if(tryCatch(!is.null(crop(KRO_slope2,extent(intersectGridClipped[k,]))), error=function(e) return(FALSE))){
    segmentMargin <- mask(crop(KRO_slope2,extent(intersectGridClipped[k,])),intersectGridClipped[k,])
    segmentMargin[segmentMargin[]<5]<-NA
      
    slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
    steepSlopes <- rbind(steepSlopes,coordinates(segmentMargin)[which(segmentMargin[]>33),])} 
    
}
  
df_steepcoord <- data.frame(steepSlopes
)
  
write.csv(df_steepcoord,'D:\\Work\\Research\\tGISM\\Code\\steeppixels_kronprins.csv', row.names = FALSE)
  
do.call(file.remove, list(list.files(tmpDir(), full.names = TRUE)))
  
  slopeclass <- matrix(NA, nrow = dim(intersectGridClipped)[1],ncol = 4)
  for(l in 1:dim(intersectGridClipped)[1]){
    slopeclass[l,1] <- length(which(slopeVals_segment[[l]][]<10))/length(slopeVals_segment[[l]][])
    slopeclass[l,2] <- length(which(slopeVals_segment[[l]][]>=10 & slopeVals_segment[[l]][]<25))/length(slopeVals_segment[[l]][])
    slopeclass[l,3] <- length(which(slopeVals_segment[[l]][]>=25 & slopeVals_segment[[l]][]<40))/length(slopeVals_segment[[l]][])
    slopeclass[l,4] <- length(which(slopeVals_segment[[l]][]>=40))/length(slopeVals_segment[[l]][])
  }
  
  png(file=path_figs&'\\Kronprins_slopeclasses.png', res = 300,width=2000,height=900) 
  barplot(t(slopeclass),col=brewer.pal(n = 4, name = "Reds"),cex.names=1,
          names.arg=seq(1,dim(intersectGridClipped)[1],1),ylim=c(0,1))
  legend("top", fill = brewer.pal(n = 4, name = "Reds"), legend = c('<10','10-25','25-40','>40'), 
         horiz = TRUE, inset = c(0,-0.5), xpd = TRUE)
  dev.off()
  
  #######################
  # Make Figures
  #######################
  
  
  
  
  set.seed(1234)
  
  wdata_Kronprins = data.frame(
    location = factor(rep(c("S01","S05","S06","S08","S10"),c(length(slopeVals_segment[[1]]),length(slopeVals_segment[[5]]),length(slopeVals_segment[[6]]),length(slopeVals_segment[[8]]),length(slopeVals_segment[[10]])))),
    slope = c(slopeVals_segment[[1]],slopeVals_segment[[5]],slopeVals_segment[[6]],slopeVals_segment[[8]],slopeVals_segment[[10]])
  )
  head(wdata_Kronprins, 4)
  
  mu <- wdata_Kronprins%>% 
    # arrange(order)  %>% 
    group_by(location) %>%
    summarise(grp.mean = median(slope,na.rm=T))
  
  png(file=path_figs&'\\KronprinsSlope distributions.png', res = 300,width=2000,height=900)  
  par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
  ggplot(wdata_Kronprins, aes(x = slope, y = location)) +
    geom_density_ridges(aes(fill = location)) +
    scale_fill_manual(values = brewer.pal(n = 5, name = "Dark2")) +
    geom_vline(aes(xintercept = grp.mean, color = location),
               data = mu, linetype = "dashed") +
    scale_x_continuous(name=expression('Slope [?]'), limits=c(0, 70)) +
    theme(axis.text.x = element_text(size=14, angle=90),
          axis.text.y = element_text(size=14, angle=0)) +
    theme(axis.title.y =element_blank(),
          axis.text.y=element_blank()) +
    theme(legend.text=element_text(size=rel(1)))
  theme_bw()
  dev.off()