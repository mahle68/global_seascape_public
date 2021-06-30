# Scripts for estimating the energy seascape; i.e. GAMs for predicting uplift as a function of location and time
# This is script 3 of 4 for reproducing the results of Nourani et al 2021, ProcB.
# session info is provided at the end of script 4 (all_figures.R)
# Elham Nourani, PhD. June. 2021; enourani@ab.mpg.de
#-----------------------------------------------------------------

library(mgcv)
library(parallel)
library(fields) #for Tps

# ---------- STEP 1: load data #####
# data was downloaded from the ECMWF Era-interim data set (see manuscript for details). 
#In addition to sea-surface temperature, temperature at 2m above ground, sun elevation categories were assigned

load("input_regional_gams.RData") #data_ls_sun. contains one list per region


# ---------- STEP 2: Modeling #####

mycl <- makeCluster(5) #adjust the number of cores based on your machine

clusterExport(mycl, "data_ls_sun") 

clusterEvalQ(mycl, {
  library(mgcv)
})

(b <- Sys.time())
models_ls <- lapply(data_ls_sun, function(x){ #one model per region
  
  x$sun_elev_f <- as.factor(x$sun_elev)
  
  gamm(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
         s(yday, by = sun_elev_f, bs = "cc") +
         s(year, bs = "re") +
         sun_elev_f , method = "REML", data = x, 
       weights = varPower(form = ~lat))
  
})

Sys.time() - b 

stopCluster(mycl)


# ---------- STEP 3: Predictions #####

#load additonal data from the Dryad repository
load("timing_for_gam_preds.RData") #timing_areas; contains information for migration timing (in julian dates) OVER THE SEA for each species


(b <- Sys.time())

preds <- lapply(c(names(models_ls)),function(x){
  d <- data_ls_sun[[x]] %>% 
    filter(yday %in% timing_areas[[x]]) %>% 
    mutate(sun_elev_f = as.factor(sun_elev))
  
  m <- models_ls[[x]]
  
  pred <- data.frame(pred = as.numeric(predict(m$gam, d)), lon = d$lon, lat = d$lat)
  
  coordinates(pred) <- ~ lon + lat
  proj4string(pred) <- wgs
  
  pred_sp <- SpatialPixelsDataFrame(pred, tolerance = 0.916421, pred@data)
  r <- raster(pred_sp)
  
  #interpolate. for visualization purposes
  surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])
  
  grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 0.1),
                     y = seq(from = extent(r)[3],to = extent(r)[4],by = 0.1))
  
  grd$coords <- matrix(c(grd$x,grd$y),ncol=2)
  
  surf.1.pred <- predict.Krig(surf.1,grd$coords)
  interpdf <- data.frame(grd$coords,surf.1.pred)
  
  colnames(interpdf)<-c("lon","lat","delta_t")
  
  coordinates(interpdf) <- ~ lon + lat
  gridded(interpdf) <- TRUE
  interpr <- raster(interpdf)
  proj4string(interpr) <- wgs
  
  return(interpr)
  
})

Sys.time() -b



names(preds) <- names(models_ls)

#mask the rasters with the relevant land layer (not shown)
#the result can be found in the Dryad repository under name: "predictions_regional_gams.RData"
#see all_figures.R for plotting the seascape map (Fig. 2)
