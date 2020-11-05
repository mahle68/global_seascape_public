#Codes for reproducing the results of GAMM for the spatio-temporal map of energy seascape.
#Elham Nourani, PhD. Nov. 4. 2020. enourani@ab.mpg.de
#use a parallel computing package (e.g. Parallel) to speed up the modeling and prediciton processes. Codes for parallelization are commented out with ##.

library(tidyverse)
library(mgcv)
library(itsadug)
library(mgcv)
library(raster)
library(dplyr)
library(sf)
library(fields)
library(sp)
library(TeachingDemos) #for subplot
#library(parallel) #for parallelization


#LOAD FILES and DEFINE VARIABLES ####
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#load the delta-t data. Make sure to set your working directory correctly.
load("data_ls_regional_gam.RData") #a list of 4 data frame, each representing one region: names(data_ls) will give you the names of the regions

#load spatial layer that will be used for masking the prediction raster
load("masking_layer.RData")
load("eur_sea.RData")

#load sample species tracks for plotting
load("tracks_for_global_map.RData")

#load the spatial layer with the extent of the study
load("region.RData")

#MODELING ####
#suggestion: use a parallel computing package (e.g. Parallel)


models_ls <- lapply(data_ls, function(x){
  
  gamm(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
         s(yday, by = sun_elev_f, bs = "cc") +
         s(year, bs = "re") +
         sun_elev_f , method = "REML", data = x, 
       weights = varPower(form = ~lat))
})



#SUPPLEMENTSRY TABLE 2: GAMM results ####
#create output tables in Latex
latex_ls <- lapply(names(models_ls), function(x){
  m <- models_ls[[x]]
  gamtabs(m, caption = x)
})



### FIGURE 2: regional seascapes and delta_t trends ####

#--- STEP 1: make predictions ----
#migration timing OVER THE SEA for each species (based on the empirical data and the literature)
timing <- list(OHB = c(260:294),
               GFB = c(277:299),
               AF = c(319:323), 
               EF = c(288:301),
               PF_EU = c(261:289),
               O_EU = c(222:277),
               PF_A = c(279:305),
               O_A = c(244:299))

timing_areas <- list(East_asia = c(min(c(timing$GFB,timing$OHB)):max(c(timing$GFB,timing$OHB))),
                     Americas = c(min(c(timing$O_A,timing$PF_A)):max(c(timing$O_A,timing$PF_A))),
                     Indian_ocean = timing$AF,
                     Europe = c(min(c(timing$O_EU,timing$PF_EU,timing$EF)):max(c(timing$O_EU,timing$PF_EU,timing$EF)))) #consider separating this into regions


#make predictions
## mycl <- makeCluster(4) 
## 
## clusterExport(mycl, list("timing_areas","models_ls", "data_ls", "wgs"))
## 
## clusterEvalQ(mycl, {
##   library(mgcv)
##   library(raster)
##   library(dplyr)
##   library(sf)
##   library(fields)
## })

preds <- parLapply(cl = mycl, c(names(models_ls)),function(x){ #use lapply if not parallelizing
  d <- data_ls[[x]] %>% 
    filter(yday %in% timing_areas[[x]])
  m <- models_ls[[x]]
  
  pred <- data.frame(pred = as.numeric(predict(m,d)), lon = d$lon, lat = d$lat)
  
  coordinates(pred) <- ~lon+lat
  gridded(pred) <- T
  r <- raster(pred)
  proj4string(r) <- wgs
  
  #interpolate. for visualization purposes
  surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])
  
  grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 0.1),
                     y = seq(from = extent(r)[3],to = extent(r)[4],by = 0.1))
  
  grd$coords <- matrix(c(grd$x,grd$y),ncol=2)
  
  surf.1.pred <- predict.Krig(surf.1,grd$coords)
  interpdf <- data.frame(grd$coords,surf.1.pred)
  
  colnames(interpdf)<-c("lon","lat","delta_t")
  
  coordinates(interpdf) <- ~lon+lat
  gridded(interpdf) <- TRUE
  interpr <- raster(interpdf)
  proj4string(interpr) <- wgs
  
  return(interpr)
  
})

## stopCluster(mycl)

#assign names to the list elements
names(preds) <- names(models_ls)

#--- STEP 2: mask the prediction rasters to keep regions of interest ----

preds_filt <- lapply(preds, function(x){
  x_f <- raster::mask(x,as(masking_layer,"Spatial"), inverse = T) 
  x_f
})

#only keep the Mediterranean and the Baltic seas for Europe
preds_filt$Europe <- mask(preds_filt$Europe, as(eur_sea,"Spatial"), inverse = F)

#--- STEP 3:  create a plotting function for subplots in Fig 2 ----
names <- c("South-East Asia", "The Americas", "Indian Ocean", "Europe")

reg_gam_plot <- function(x){
  m <- models_ls[[x]]
  t <- timing_areas[[x]]
  
  plot(0, type = "n", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "", ylab = "")
  abline(h = 0, col = "gray60",lty = 1, lwd = 0.5)
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view="yday", plot_all="sun_elev_f", rm.ranef=F, lwd = 1.5, #ylim=c(-2,5),
              col = "grey50", hide.label = TRUE, 
              legend_plot_all =  F, 
              h0 = NULL, add = T, lty = c(1,5,3))
  
  axis(side = 1, at = c(100, 200, 300), line = 0, labels = c(100,200,300), 
       tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015)
  axis(side = 2, at = c(-2,0,2,4), line = 0.15, labels = c(-2,0,2,4),
       tick = T , col.ticks = 1,col = NA, lty = NULL, tck = -.015, 
       las = 2) # text perpendicular to axis label 
  mtext(expression(italic(paste(Delta,"T"))), 2, line = 1.2 ,las = 2.5, cex = 0.5)
  mtext("Day of year", 1, line = 1.2 , cex = 0.5)
  mtext(names[[x]], 3, line = 0 , cex = 0.7)
  
}

#--- STEP 4:  create color palette ----
#create a raster for legend with values from -1 to 5 (range of values in the prediction rasters)
imaginary_r<-preds_filt[[1]]
imaginary_r@data@values[3:9]<- rep(5,7)
imaginary_r@data@values[10:16]<- rep(-1,7)

#create a color palette
cuts<-seq(-1,5,0.01) #set breaks
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)

#--- STEP 5:  plot everything ----
X11(width = 11, height = 5.2) #make the window proportional to region

par(mfrow=c(1,1),
    fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    #cex = 0.5,
    oma = c(0,0,0,0),
    mar = c(0, 0, 0, 0),
    lend = 1  #rectangular line endings (trick for adding the rectangle to the legend)
)


plot(region, col="#e5e5e5",border="#e5e5e5")
plot(preds_filt[[1]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[2]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[3]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[4]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 

#add tracks
lwd <- 1.2
col <- "grey25"
plot(st_geometry(sp_samples$EF_sample), add= T, lty = 2, lwd = lwd, col = col)
plot(st_geometry(sp_samples$PF_sample), add= T, lty = 3, lwd = lwd, col = col)
plot(st_geometry(sp_samples$O_sample), add= T, lty = 5, lwd = lwd, col = col)
plot(st_geometry(sp_samples$OHB_sample), add= T, lty = 4, lwd = lwd, col = col)
plot(st_geometry(sp_samples$GFB_sample), add= T, lty = 1, lwd = lwd, col = col)
plot(st_geometry(sp_samples$AF_sample), add= T, lty = 6, lwd = lwd, col = col)

#add latitude lines
clip(-130, 157, -50, 73)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")

#add a frame for the sub-plots
rect(xleft = -85,
     xright = 153,
     ybottom =  -52,
     ytop = -4,
     col="white",
     border = NA)

#add subplots...
centers_x <- c(124,-56,64,4) #distance between centers = 60

for(i in 1:length(centers_x)){
  rect(xleft = centers_x[i] - 27,
       xright = centers_x[i] + 27,
       ybottom =  -42, #-38
       ytop = -6, #-2
       col="white")
  
  subplot(reg_gam_plot(i), x = centers_x[i] + 4,y = -23, size = c(1.6,0.7),  #-23 was -19
          pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
}

#add legend
plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal,
     smallplot = c(0.059,0.169, 0.145,0.16), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                      labels = c(-1,0,2,4), 
                      col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                      col.ticks = NA,
                      line = -1.3),
     horizontal = T,
     legend.args = list(text = expression(italic(paste(Delta,"T", "(°C)"))), side = 3, font = 2, line = 0.1, cex = 0.7)
)

text(x = -105,y = -8, "Map legend", cex = 0.8)
legend(-126,-10, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                            "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)

legend(-30,-43, legend = c("sea-crossing period","High sun elevation", "Low sun elevation", "Night"),
       lty = c(1,1,2,3), lwd = c(9,1,1,1), col = c("#99CC0060", rep("black",3)),cex = 0.55, bty = "n", seg.len = 3, horiz = T)



