################################################################################
# Produce figures 
#
# IceMargin_SubbasinPlots.R
#
# ReadMe: 
#
# https://github.com/fidelsteiner/tGISM
#
# Input:
# 
#
# Created:          2021/10/15
# Latest Revision:  2024/02/16
#
# Jakob F Steiner | jakob.steiner@uni-graz.at | fidelsteiner.github.io 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# packages (if not installed yet: install.packages('examplePackage')
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(ggplot2)
library(dplyr)
library(ggridges)
library(RColorBrewer)
library(terra)
library(sf)
library(ncdf4)
library(plotly)
##################
# Multiple subsites for analysis of the margin
##################

data_basepath <- 'E:\\GeospatialData'                     # path for heavy input data
data_ArcticDEM <- 'E:\\ArcticDEM'                         # path for ArcticDEM data
output_basepath <- 'C:\\Work\\Research\\tGISM'            # path for code/output   

##################
# Paths
##################

path_outlines <- 'E:\\Research\\OtherResearch\\Greenland\\Code\\SectorMargins\\ManualRedRock' # Not used, remove
path_figs <- paste(output_basepath,'\\Figures',sep = "")
pathDEM_Pleiades <- paste(data_basepath,'DEMData\\ValidationDEMs\\RedRock_Pleiades\\Pleiades\\DEM\\2017-08-15_RedRockCliff\\PGO\\2017-08-15_RedRockCliff',sep = "")
pathDEM_Arctic <- data_ArcticDEM
data_Bedmachine <- 'E:\\GeospatialData\\Bedmachine\\5000000451838\\160281892' # path Bedmachine
data_lakes <- 'E:\\GeospatialData\\Lakes\\Greenland' # path Lakes
folderSectors <- paste(output_basepath,'\\Sectors',sep = "")


projec <-'+proj=utm +datum=WGS84'

projec_utm <-'+proj=utm +zone=19 +north +datum=WGS84'
projec_27N <-'+proj=utm +zone=27N +datum=WGS84'

regionSW <- c(137,136,135,134,133,132,131,130,71,70,69,40,22,21,20,19,18,17,16,15,14,13,12,11)
regionSE <- c(63,64,65,67,68,73,79,80,83,88,89,90,92,93,94,106,107,108,110,111,113,114,115,117,118,119,120,121,125,126,128,186,187,188,189,190,191,192,193,199,200,201,202,209,210,211,212,213,231,232,233,237,243,258)
regionCE <- c(60,61,62,74,75,76,95,96,100,102,122,124,203,204,248,249,250,251,252,255,256)
regionCW <- c(3,4,5,6,7,8,9,10,72,81,82,85,218,219,220,221,222,230)
regionNE <- c(58,59,87,129,138,139,140,141,142,144,145,217,244,245,246,247)
regionNO <- c(39,50,51,52,53,54,55,56,57,77,84,109,146,176,215,216,223,224)
regionNW <- c(1,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,78,86,103,104,105,127,147,148,159,160,161,162,163,164,165,166,167,168,169,170,172,173,174,175,180,185,195,196,197,207,208,214,225,238,239,240)
regionPGICNO <- c(302,303,304)
regionPGICNE <- c(305,306,307)
regionPGICCE <- c(308,309,310)
regionPGICSE <- c(311,312,313)
regionPGICSW <- c(314,315)

PGICCW <- 300
PGICNW <- 301
PGICNOW <- c(302)
PGICNOC <- c(303)
PGICNOE <- c(304)
PGICNEN <- c(305)
PGICNEC <- c(306)
PGICNES <- c(307)
PGICCEN <- c(308)
PGICCEC <- c(309)
PGICCES <- c(310)
PGICSEN <- c(311)
PGICSEC <- c(312)
PGICSES <- c(313)
PGICSWS <- c(314)
PGICSWN <- c(315)

# Read general margin statistics
marginStats <- 'C:\\Work\\Research\\tGISM\\Sectors\\marginstatistics.csv'

margStats <- read.csv(marginStats)


colorMargins <- c('#ffa8fcff','#edf0ffff','#a8fff9ff','#abffa7ff','#fffda6ff','#ffbbadff','#cfa8bbff',
                  '#ffa8fcff','#edf0ffff','#a8fff9ff','#a8fff9ff','#a8fff9ff','#abffa7ff','#abffa7ff','#abffa7ff',
                  '#fffda6ff','#fffda6ff','#fffda6ff','#ffbbadff','#ffbbadff','#ffbbadff','#cfa8bbff','#cfa8bbff',
                  '#abffa7ff','#fffda6ff','#ffbbadff','#cfa8bbff')
nameMargins <- c('CW','NW','NO','NE', 'CE', 'SE', 'SW', 
                 'pCW','pNW','pNOW','pNOC','pNOE',
                 'pNEN','pNEC','pNES','pCEN','pCEC','pCES','pSEN','pSEC',
                 'pSES','pSWS','pSWN',
                 'pNO','pNE','pCE','pSE','pSW')

# Plot fraction plots for each subbasin
MarginDonutsGreenland <- function(regionDat, countMargin){
  
  marginDat <- margStats[match(regionDat,margStats[,1]),]
  
  dataML <- data.frame(
    category = c('terrestrial','marine','lake'),
    count=c(as.vector(colSums(marginDat[,3:5]))) 
  )
  print(sqrt(sum(dataML$count)/1000000)*10)
  dataML$fraction <- dataML$count / sum(dataML$count)
  dataML$ymax <- cumsum(dataML$fraction)
  dataML$ymin <- c(0, head(dataML$ymax, n=-1))
  dataML$labelPosition <- (dataML$ymax + dataML$ymin) / 2
  dataML$label <- paste0(round(dataML$fraction*100, digits = 1),"%")
  print(dataML$label)
  print(sum(dataML$count)/1000)
  pGM <- ggplot(dataML, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect(color =colorMargins[countMargin],lwd=8) +
    #geom_text( x=c(4.3,4.3,4.3), aes(y=labelPosition, label=label), size=30) +
    scale_fill_manual(values=c( "#008000", "#0000FF","#554400")) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA),
          #panel.background = element_blank(),
          panel.grid = element_blank(),
          title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0,"null"),
          panel.border = element_blank()) +
        
    geom_text(aes(label = nameMargins[countMargin]), x=2,y=500,size=30,color ='black')
  
  
  
  png(file=paste(output_basepath,'\\Figures\\Donuts\\',nameMargins[countMargin],'_marginlengths.png',sep = ""), res = 300,width=3600,height=3600,bg = "transparent")
  par(bty = 'n') 
  print(pGM,axes=FALSE)
  dev.off() 
}
MarginDonutsGreenland(regionCW, 1)
MarginDonutsGreenland(regionNW, 2)
MarginDonutsGreenland(regionNO, 3)
MarginDonutsGreenland(regionNE, 4)
MarginDonutsGreenland(regionCE, 5)
MarginDonutsGreenland(regionSE, 6)
MarginDonutsGreenland(regionSW, 7)
MarginDonutsGreenland(PGICCW, 8)
MarginDonutsGreenland(PGICNW, 9)
MarginDonutsGreenland(PGICNOW, 10)
MarginDonutsGreenland(PGICNOC, 11)
MarginDonutsGreenland(PGICNOE, 12)
MarginDonutsGreenland(PGICNEN, 13)
MarginDonutsGreenland(PGICNEC, 14)
MarginDonutsGreenland(PGICNES, 15)
MarginDonutsGreenland(PGICCEN, 16)
MarginDonutsGreenland(PGICCEC, 17)
MarginDonutsGreenland(PGICCES, 18)
MarginDonutsGreenland(PGICSEN, 19)
MarginDonutsGreenland(PGICSEC, 20)
MarginDonutsGreenland(PGICSES, 21)
MarginDonutsGreenland(PGICSWS, 22)
MarginDonutsGreenland(PGICSWN, 23)


# Plot slope statistics

MarginSlopesGreenland <- function(regionDat, countMargin){
  
med <- vector()
hig <- vector()
low <- vector()
cnt <- vector()

png(file=paste(output_basepath,'\\Figures\\SlopeGraphs\\',nameMargins[countMargin],'_marginslopes.png',sep = ""), res = 300,width=3600,height=3600)
par(bty = 'n') 
plot(c(0,0), c(0,0), type = "n",xlim = c(0 , 1), ylim = c( 0, 1),xaxs='i', xlab ='',ylab='',xaxt = "n",yaxt = "n") 
axis(side = 1, lwd = 4,labels =F)
axis(side = 2, lwd = 4,labels =F)
axis(side = 3, lwd = 4,labels =F)
axis(side = 4, lwd = 4,labels =F)
abline(a=0,b=1,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")
abline(a=0,b=33/45,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")

if(countMargin<10000){
for(i in 1:length(regionDat)){
    slopeStats <- read.csv(paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'_slopePercentiles.csv',sep=''))

  med[i] <- slopeStats[dim(slopeStats)[1],4]
  hig[i] <- slopeStats[dim(slopeStats)[1],6]
  low[i] <- slopeStats[dim(slopeStats)[1],2]
  cnt[i] <- slopeStats[dim(slopeStats)[1],8]
  abline(a=0,b=med[i]/45,col='red',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=hig[i]/45,col='black',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=low[i]/45,col='blue',lwd=1,xaxt = "n",yaxt = "n")
}
abline(a=0,b=sum((cnt/sum(cnt))*med)/45,col='red',lwd=8,lty=2,xaxt = "n",yaxt = "n")
text(0.3,0.9,paste0(round(sum((cnt/sum(cnt))*med),1),"°"),cex=10,col='red')
dev.off() 
}

png(file=paste(output_basepath,'\\Figures\\SlopeGraphs\\dummy_marginslopes.png',sep = ""), res = 300,width=3600,height=3600)
par(bty = 'n') 
plot(c(0,0), c(0,0), type = "n",xlim = c(0 , 1), ylim = c( 0, 1),xaxs='i', xlab ='',ylab='',xaxt = "n",yaxt = "n") 
axis(side = 1, lwd = 4,labels =F)
axis(side = 2, lwd = 4,labels =F)
axis(side = 3, lwd = 4,labels =F)
axis(side = 4, lwd = 4,labels =F)
abline(a=0,b=1,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")
abline(a=0,b=33/45,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")


  abline(a=0,b=23/45,col='red',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=65/45,col='black',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=8/45,col='blue',lwd=1,xaxt = "n",yaxt = "n")

abline(a=0,b=25/45,col='red',lwd=8,lty=2,xaxt = "n",yaxt = "n")
text(0.3,0.9,paste0('median',"°"),cex=6,col='black')
dev.off() 
}

MarginSlopesGreenland(regionCW, 1)
MarginSlopesGreenland(regionNW, 2)
MarginSlopesGreenland(regionNO, 3)
MarginSlopesGreenland(regionNE, 4)
MarginSlopesGreenland(regionCE, 5)
MarginSlopesGreenland(regionSE, 6)
MarginSlopesGreenland(regionSW, 7)
MarginSlopesGreenland(PGICCW, 8)
MarginSlopesGreenland(PGICNW, 9)
MarginSlopesGreenland(PGICNOW, 10)
MarginSlopesGreenland(PGICNOC, 11)
MarginSlopesGreenland(PGICNOE, 12)
MarginSlopesGreenland(PGICNEN, 13)
MarginSlopesGreenland(PGICNEC, 14)
MarginSlopesGreenland(PGICNES, 15)
MarginSlopesGreenland(PGICCEN, 16)
MarginSlopesGreenland(PGICCEC, 17)
MarginSlopesGreenland(PGICCES, 18)
MarginSlopesGreenland(PGICSEN, 19)
MarginSlopesGreenland(PGICSEC, 20)
MarginSlopesGreenland(PGICSES, 21)
MarginSlopesGreenland(PGICSWS, 22)
MarginSlopesGreenland(PGICSWN, 23)
MarginSlopesGreenland(regionPGICNO, 24)
MarginSlopesGreenland(regionPGICNE, 25)
MarginSlopesGreenland(regionPGICCE, 26)
MarginSlopesGreenland(regionPGICSE, 27)
MarginSlopesGreenland(regionPGICSW, 28)



# Slope violin plots per grid cell

MarginGridSections <- function(regionDat,LOC){
  steepLength <- vector()
  shallowLength <- vector()
  cliffLength <- vector()
  truecliffLength <- vector()
  
  for(subk in 1:length(regionDat)){ 

    SubBasinID <- regionDat[subk]
    sPerc <- read.csv(paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_slopePercentiles.csv',sep = ""))
    numGrid <- dim(sPerc)[1]-1
    
    sGrid <- ogrInfo(paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""))
    sGrid <- readOGR(dsn=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""),layer = sGrid$layer)
    
    shallowRamp <- vector()
    steepRamp <- vector()
    cliffRamp <- vector()
    truecliffRamp<- vector()
    slopeclass <- matrix(NA, nrow = numGrid,ncol = 4)
    
    # Save output for this subbasin
    tabSB <- matrix(NA, ncol = 8, nrow = numGrid+1)
    tabSB[,1] <- c(1:numGrid,'TOTAL')
    
    for(segk in 1:numGrid){
      gridslVal <- read.table(paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',sPerc[segk,1],'_slopeValues.txt',sep = ""),fill = TRUE)
      if(mean(unlist(gridslVal),na.rm=T)!='NaN'){
        # store grid cells according to their steepness class
        if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<10))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.75){
          shallowRamp <- c(shallowRamp,sPerc[segk,1])
          shallowLength <- c(shallowLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
          
        }
        if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<45&unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.1){
          steepRamp <- c(steepRamp,sPerc[segk,1])
          steepLength <- c(steepLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<45&unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
        }
        if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.1){
          cliffRamp <- c(cliffRamp,sPerc[segk,1])
          cliffLength <- c(cliffLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
        }
        
        if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>65))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.05){
          truecliffRamp <- c(truecliffRamp,sPerc[segk,1])
          truecliffLength <- c(truecliffLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>65))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
          
        }
        
        
        tabSB[segk,2:6] <- quantile(unlist(gridslVal), probs = c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T)
        tabSB[segk,7] <- sPerc[segk,7]
        
        tabSB[segk,8] <- sPerc[segk,8]
        
        
        slopeclass[segk,1] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<10))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        slopeclass[segk,2] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>=10 & unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        slopeclass[segk,3] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>=20 & unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        slopeclass[segk,4] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>=45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        
        wdata_violin = data.frame(
          location = factor(rep(c(SubBasinID),c(length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])))),
          slope = c(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        )
        
        mu <- wdata_violin%>% 
          # arrange(order)  %>% 
          group_by(location) %>%
          summarise(grp.mean = median(slope,na.rm=T))
        
        tryCatch({
          png(file=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',sPerc[segk,1],'_violin.png',sep=''), res = 300,width=2000,height=900)  
          par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
          print(ggplot(wdata_violin, aes(x = slope, y = location)) +
                  geom_density_ridges(aes(fill = location)) +
                  scale_fill_manual(values = brewer.pal(n = 5, name = "Dark2")) +
                  geom_vline(aes(xintercept = grp.mean, color = location),
                             data = mu, linetype = "dashed") +
                  scale_x_continuous(name=expression('Slope [°]'), limits=c(0, 95)) + 
                  #scale_y_continuous(name='', limits=c(1)) +
                  theme_bw() +
                  theme(axis.text.x = element_text(size=14, angle=270),
                        axis.text.y = element_text(size=14, angle=0)) +
                  theme(axis.title.y =element_blank(),
                        axis.text.y=element_blank()) +
                  theme(legend.position="none"))
          dev.off()
        }, error=function(e){})
      }
    }
    
    png(file=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_slopeclasses.png',sep=''), res = 300,width=2000,height=900) 
    barplot(t(slopeclass),col=brewer.pal(n = 4, name = "Reds"),cex.names=1,
            names.arg=sPerc[1:numGrid,1],ylim=c(0,1))
    legend("top", fill = brewer.pal(n = 4, name = "Reds"), legend = c('<10','10-20','20-45','>45'), 
           horiz = TRUE, inset = c(0,-0.5), xpd = TRUE)
    dev.off()
    
    shallowGrid <- sGrid[which(!is.na(match(sGrid$id,as.numeric(shallowRamp)))),]
    rgdal::writeOGR(shallowGrid,paste(output_basepath,'\\Sectors\\shallowGrids\\',SubBasinID,'ShallowGrid.shp',sep = ""),layer = 'shallowGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
    steepGrid <- sGrid[which(!is.na(match(sGrid$id,as.numeric(steepRamp)))),]
    rgdal::writeOGR(steepGrid,paste(output_basepath,'\\Sectors\\steepGrids\\',SubBasinID,'SteepGrid.shp',sep = ""),layer = 'steepGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
    cliffGrid <- sGrid[which(!is.na(match(sGrid$id,as.numeric(cliffRamp)))),]
    rgdal::writeOGR(cliffGrid,paste(output_basepath,'\\Sectors\\cliffGrids\\',SubBasinID,'CliffGrid.shp',sep = ""),layer = 'cliffGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
    truecliffGrid <- sGrid[which(!is.na(match(sGrid$id,as.numeric(truecliffRamp)))),]
    rgdal::writeOGR(truecliffGrid,paste(output_basepath,'\\Sectors\\cliffGrids\\',SubBasinID,'trueCliffGrid.shp',sep = ""),layer = 'truecliffGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
    
    
    
    # Hiawatha example
    if(SubBasinID==176){
      gridslVal <- list()
      for(segk in c(14,15,23,25)){
        gridslVal[[segk]] <- read.table(paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',sPerc[segk,1],'_slopeValues.txt',sep = ""),fill = TRUE)
      }
      
      wdata_violin = data.frame(
        location = factor(rep(c(14,15,23,25),c(length(unlist(gridslVal[[14]])),length(unlist(gridslVal[[23]])),length(unlist(gridslVal[[15]])),length(unlist(gridslVal[[25]]))))),
        slope = c(unlist(gridslVal[[14]]),unlist(gridslVal[[23]]),unlist(gridslVal[[15]]),unlist(gridslVal[[25]]))
      )
      
      mu <- wdata_violin%>% 
        # arrange(order)  %>% 
        group_by(location) %>%
        summarise(grp.mean = median(slope,na.rm=T))
      
      png(file=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_Hiawatha_violin.png',sep=''), res = 300,width=2000,height=900)  
      par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
      print(ggplot(wdata_violin, aes(x = slope, y = location)) +
              geom_density_ridges(aes(fill = location)) +
              scale_fill_manual(values = brewer.pal(n = 5, name = "Dark2")) +
              geom_vline(aes(xintercept = 10, color = 'black'),
                         data = mu, linetype = "solid") +
              geom_vline(aes(xintercept = 20, color = 'black'),
                         data = mu, linetype = 'solid') +
              geom_vline(aes(xintercept = 45, color = 'black'),
                         data = mu, linetype = 'solid') +
              scale_x_continuous(name=expression('Slope [°]'), limits=c(0, 80)) + 
              #scale_y_continuous(name='', limits=c(1)) +
              theme_bw() +
              theme(axis.text.x = element_text(size=14, angle=270),
                    axis.text.y = element_text(size=14, angle=0)) +
              theme(axis.title.y =element_blank(),
                    axis.text.y=element_blank()) +
              theme(legend.position="none"))
      dev.off()
      
    }
    
    
  }
 
  write.csv(steepLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_steeplengths.csv',sep = ""))
  write.csv(cliffLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_clifflengths.csv',sep = ""))
  write.csv(shallowLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_shallowlengths.csv',sep = ""))
  write.csv(truecliffLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_trueclifflengths.csv',sep = ""))
  
}


MarginGridSections(regionNO,'NO')

MarginGridSections(regionNW,'NW')

MarginGridSections(regionCW,'CW')

MarginGridSections(regionNE,'NE')

MarginGridSections(regionCE,'CE')

MarginGridSections(regionSE,'SE')

MarginGridSections(regionSW,'SW')

MarginGridSections(PGICCW,'pCW')

MarginGridSections(PGICNW,'pNW')