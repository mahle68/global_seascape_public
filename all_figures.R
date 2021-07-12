# Scripts for producing the figures (main text and supplementary)
# This is script 4 of 4 for reproducing the results of Nourani et al 2021, ProcB.
# session info is provided at the end of script 4 (all_figures.R)
# Elham Nourani, PhD. June. 2021; enourani@ab.mpg.de
#-----------------------------------------------------------------

library(tidyverse)
library(TeachingDemos) #for subplot
library(itsadug) #for gam plots
library(mgcv)
library(sf)
library(move)
library(scales)
library(maptools)

source("functions.R") #as a part of the repository on Github
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")


# ---------- Fig 1: w_star #####

#Here we build the GAM model and generate the resulting plot

#load data (annotated sea-crossing points. available in the Dryad repository). This file contains annotations using ECMWF ERA5 (columns with _5 suffix) and ERA-interim (columns with _i suffix). 
#For the w_star estimations, we use the ERA_interim data. The variables of interest wer not available from ERA5 via the Movebank annotation service at the time this study was conducted. 
#We will use the ERA5 annotations later on for creating maps for Fig S3.

load("annotated_points.RData") #ann_pts

#only keep rows with positive delta_t
data <- ann_pts %>% 
  filter(delta_t>= 0)

data$w_star <- w_star(blh = data$blh, T2m = data$t2m, 
                        s_flux = data$s_flux, m_flux = data$m_flux)


#model and estimate confidence intervals
M <- gam(w_star ~ s(delta_t, bs = "cr", k = 10), data = data)
fit <- predict(M, se = T)$fit
se <- predict(M, se = T)$se.fit

lcl <- fit - 1.96* se
ucl <- fit + 1.96* se

i.for <- order(data$delta_t)
i.back <- order(data$delta_t, decreasing = TRUE)

x.polygon <- c(data$delta_t[i.for] , data$delta_t[i.back])
y.polygon <- c(ucl[i.for] , lcl[i.back])


#create a color palette
Pal <- colorRampPalette(c("darkgoldenrod1","lightpink1", "mediumblue")) 
Cols <- paste0(Pal(3), "80")

data$color <- as.factor(data$sun_elev)
levels(data$color) <- Cols

data$shape <- as.factor(data$species)
levels(data$shape) <- c(0,1,2,3,4,5)


#Plot
X11(width = 5.2, height = 3) #inches

par(mfrow=c(1,1), 
    bty = "l",
    cex.axis = 0.7,
    font.lab = 3,
    cex.lab = 0.9,
    mgp=c(2,0.5,0), #margin line for the axis titles
    mar = c(3.5,3.5,0.5,6.35),
    #oma = c(0,0,0,4),
    xpd=  TRUE) # allows drawing legend in outer margins

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(0,8.3), ylim = c(0.5,5), xlab = expression(italic(paste(Delta,"T", "(°C)"))), ylab = "w* (m/s)")

with(data,points(delta_t, w_star, col= as.character(data$color), pch = as.numeric(data$shape), cex = 0.6))

polygon(x.polygon , y.polygon , col = alpha("grey", 0.5) , border = NA) #confidence intervals

lines(data$delta_t[i.for] , fit[i.for], col = "black" , lwd = 1.1)

axis(side = 1, at = c(0,2,4,6,8), labels = c(0,2,4,6,8), 
     tick = T , col.ticks = 1, col = NA, tck = -.015,lwd = 0, lwd.ticks = 1)
axis(side= 2, at= seq(1,5,1), labels= seq(1,5,1),
     tick=T , col.ticks = 1, col = NA, tck=-.015,
     las=2) 

text(9.8,4.6, "Time of day", cex = 0.7)
legend(x = 8.7, y = 4.5, legend = c("daytime: high sun", "daytime: low sun", "night"), col = Cols, #coords indicate top-left
       cex = 0.7, pt.cex = 0.9, bg = "white", bty = "n", pch = 20)


text(9.5,3, "Species", cex = 0.7)
legend(x = 8.7, y = 2.9, legend = c("Pernis ptilorhynchus", "Pandion haliaetus","Butastur indicus", "Falco peregrinus", "Falco eleonorae"), #species : unique(data$species), #coords indicate top-left
       cex = 0.7, pt.cex = 0.5, bg = "white", bty = "n", pch = as.numeric(unique(data$shape)),
       text.font = 3)



# ---------- Fig 2: global seascape map #####

#Bird figures are not included in this code for copyright reasons.

#load the data  
load("predictions_regional_gams.RData") #preds_filt; rasters of predictions using GAMs produced in seascape_GAM_analysis.R
load("models_ls_reg_GAMs.RData") #models_ls; list of regional gams produced in seascape_GAM_analysis.R
load("timing_for_gam_preds.RData") #timing_areas; timing of sea-crossing
load("sample_tracks.RData") #sp_samples; sample tracks for each species to add to the map

#load continents shapefile as the map background
region <- st_read("continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 158, ymin = -74, ymax = 71) %>%
  st_union()

#names of regions
names <- c("South-East Asia", "The Americas", "Indian Ocean", "Europe", "Mozambique Channel")


#imaginary raster for legend. The legend should go from -1 to 5 (range of values in the prediction rasters)
imaginary_r <- preds_filt[[1]]
imaginary_r@data@values[3:9] <- rep(5,7)
imaginary_r@data@values[10:16] <- rep(-1,7)

#create a color palette
cuts <- seq(-1,5,0.01) #set breaks
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)

X11(width = 11, height = 6) #in inches

par(mfrow=c(1,1),
    fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    oma = c(0,0,0,0),
    mar = c(0, 0, 0, 0),
    lend = 1  #rectangular line endings (trick for adding the rectangle to the legend)
)

#plot the background map
plot(region, col="#e5e5e5",border="#e5e5e5")
plot(preds_filt[[1]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[2]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[3]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[4]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[5]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 

#add tracks
lwd <- 1.2
col <- "grey25"
plot(st_geometry(sp_samples$EF_sample), add= T, lty = 2, lwd = lwd, col = col)
plot(st_geometry(sp_samples$PF_sample), add= T, lty = 3, lwd = lwd, col = col)
plot(st_geometry(sp_samples$O_sample), add= T, lty = 5, lwd = lwd, col = col)
plot(st_geometry(sp_samples$OHB_sample), add= T, lty = 4, lwd = lwd, col = col)
plot(st_geometry(sp_samples$GFB_sample), add= T, lty = 1, lwd = lwd, col = col)
plot(st_geometry(sp_samples$AF_sample), add= T, lty = 6, lwd = lwd, col = col)
plot(st_geometry(sp_samples$EF_S_sample), add = T, lty = 2, lwd = lwd, col = col)

#add latitude lines
clip(-130, 157, -50, 73)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
text(x = -125, y = c(32,62), labels = c("30° N", "60° N"), cex = 0.6, col = "grey65")

#add a frame for the sub-plots and legend
rect(xleft = -130,
     xright = 159,
     ybottom =  -74.5,
     ytop = -28,
     col="white",
     border = NA)

rect(xleft = -130,
     xright = -87,
     ybottom =  -15,
     ytop = 7,
     col="white",
     border = NA)

#add subplots
centers_x <- c(130,-102,72,-44,14) #distance between centers = 58

for(i in 1:length(centers_x)){
  rect(xleft = centers_x[i] - 27.5,
       xright = centers_x[i] + 28,
       ybottom =  -66,#-42
       ytop = -30, #-4
       col="white")
  
  subplot(reg_gam_plot(i), x = centers_x[i] + 4,y = -47, size = c(1.6,0.7),  #-23 was -19
          pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
}

#add legend
plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal, 
     smallplot = c(0.04,0.15,0.37,.385), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                      labels = c(-1,0,2,4), 
                      col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                      col.ticks = NA,
                      line = -1.3),
     horizontal = T,
     legend.args = list(text = expression(italic(paste(Delta,"T", "(°C)"))), side = 3, font = 2, line = 0.1, cex = 0.7)
)


text(x = -118,y = 10, "Map legend", cex = 0.8)
legend(-130,8, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                                  "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)

 legend(-45,-67, legend = c("sea-crossing period","High sun elevation", "Low sun elevation", "Night"),
        lty = c(1,1,2,3), lwd = c(9,1,1,1), col = c("#99CC0060", rep("black",3)),cex = 0.55, bty = "n", seg.len = 3, horiz = T)



 
 
# ---------- Fig 3: INLA results: coefficients #####

#load the model and data used to build it (produced in step_selection_analysis.R; also available on the Dryad repository)
load("INLA_model.RData") #M; final INLA model
load("annotated_steps.RData") #ann_cmpl; data used for INLA modeling. This dataframe includes used and alternative steps and can be reproduced using step_generation.R

#calculate z_scores for predictor variables
all_data <- ann_cmpl %>% 
  mutate_at(c("delta_t", "wind_speed", "wind_support", "wind_support_var", "abs_cross_wind", "delta_t_var"),
            list(z = ~(scale(.)))) 

# posterior means of coefficients
graph <- as.data.frame(summary(M)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)

VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

min <- min(graph$Lower,na.rm = T)
max <- max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

#plot
X11(width = 4.1, height = 2.7)
par(cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 4.15, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-2,3), ylim = c(0.7,4.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph$Estimate, graph$Factor_n, col = "cornflowerblue", pch = 20, cex = 2)
arrows(graph$Lower, graph$Factor_n,
       graph$Upper, graph$Factor_n,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line

#add axes
axis(side= 1, at = c(-2, -1, 0, 1, 2, 3), labels = c(-2, -1, 0, 1, 2, 3), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4),
     labels = c(expression(paste(italic(paste(Delta,"T"))," : Wind support")),
                "Wind support var","Wind support", expression(italic(paste(Delta,"T")))),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

# ---------- Fig 4: INLA results: interaction effect plot #####
#plot the interaction between wind support and delta-t

#load data (prepared in step_selection_analysis.R)
load("INLA_model_preds.RData") #M_pred; INLA model used for predictions
load("new_data_for_modeling.RData") #new_data; data generated to plot predictions (interaction between delta_t and wind support)

#extract predicted values
used_na <- which(is.na(new_data$used))
 
 
#extract information for rows that had NAs as response variables
preds <- data.frame(delta_t = new_data[is.na(new_data$used) ,"delta_t_z"],
                    wind_support = new_data[is.na(new_data$used) ,"wind_support_z"],
                    wind_support_var = new_data[is.na(new_data$used) ,"wind_support_var_z"],
                    preds = M_pred$summary.fitted.values[used_na,"mean"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds))) #this should be between 0-1

#create a raster of predictions
avg_preds <- preds %>% 
  group_by(delta_t_z, wind_support_z) %>%  
  summarise(avg_pres = mean(prob_pres)) %>% 
  ungroup() %>% 
  mutate(wspt_backtr = wind_support_z * attr(all_data$wind_support_z, 'scaled:scale') + attr(all_data$wind_support_z, 'scaled:center'),
         dt_backtr = delta_t_z * attr(all_data$delta_t_z, 'scaled:scale') + attr(all_data$delta_t_z, 'scaled:center')) %>% 
  dplyr::select(-c("delta_t_z","wind_support_z")) %>% 
  as.data.frame()
 
 
coordinates(avg_preds) <-~ wspt_backtr + dt_backtr 
gridded(avg_preds) <- TRUE
r <- raster(avg_preds)
 

#interpolate. for visualization purposes
surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])

grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 2),
                   y = seq(from = extent(r)[3],to = extent(r)[4],by = 2))

grd$coords <- matrix(c(grd$x,grd$y),ncol=2)

surf.1.pred <- predict.Krig(surf.1,grd$coords)
interpdf <- data.frame(grd$coords,surf.1.pred)

colnames(interpdf) <- c("wind_support","delta_t","prob_pres")

coordinates(interpdf) <- ~ wind_support + delta_t
gridded(interpdf) <- TRUE
interpr <- raster(interpdf)


#create a color palette
cuts <- c(0, 0.25,0.5,0.75,1) #set breaks
pal <- colorRampPalette(c("aliceblue", "lightskyblue1", "khaki2", "sandybrown", "salmon2","tomato"))
colpal <- pal(200)
 


#plot
X11(width = 5, height = 4)

par(cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(0, 0, 0, 1.5),
    bty = "n",
    mgp = c(1,0.5,0)
)


plot(interpr, col = colpal, axes = F, box = F, legend = F, ext = extent(c(-22, 28.9, -9.7, 14))) #crop to the extent of observed data

#add axes
axis(side = 1, at = seq(-20,30,10),
     labels = seq(-20,30,10),
     tick = T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck = -.015, line = -5.75, cex.axis = 0.7) #tick marks smaller than default by this proportion

axis(side = 2, at = c(-5, 0,5, 10), labels = c(-5, 0,5, 10), 
     tick = T ,col = NA, col.ticks = 1, tck = -.015, las = 2, cex.axis = 0.7)

lines(x = c(-21.9, -21.9), y = c(-9.9,13.9))
abline(h =-10)

#axis titles
mtext("Wind support (m/s)", 1, line = -4, cex = 0.9, font = 3)
mtext(expression(italic(paste(Delta,"T", "(°C)"))), 2, line = 1.2, cex = 0.9)

#add legend
plot(interpr, legend.only = T, horizontal = F, col = colpal, legend.args = list("Probability of use", side = 4, font = 1, line = 1.5, cex = 0.7),
     legend.shrink = 0.4,
     #smallplot= c(0.12,0.7, 0.06,0.09),
     axis.args = list(at = seq(0,1,0.25), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                      labels = seq(0,1,0.25), 
                      col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                      col.ticks = NA,
                      line = -0.8, cex.axis = 0.7))


# ---------- Fig S1: species-specific coefficients #####

#original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)

#load the INLA model
load("INLA_model.RData") #M

#species
species_names <- c("O", "PF", "EF", "OHB")
species_names <- c("Osprey", "Peregrine \n falcon", "Eleonora's \n falcom", "Oriental \n honey buzzard")

tab_dt <- data.frame(ID = as.factor(M$summary.random$species1$ID),
                     mean = M$summary.random$species1$mean,
                     IClower = M$summary.random$species1[, 4],
                     ICupper = M$summary.random$species1[, 6])


tab_wspt <- data.frame(ID = as.factor(M$summary.random$species2$ID),
                       mean = M$summary.random$species2$mean,
                       IClower = M$summary.random$species2[, 4],
                       ICupper = M$summary.random$species2[, 6])


tab_wspt_var <- data.frame(ID = as.factor(M$summary.random$species3$ID),
                       mean = M$summary.random$species3$mean,
                       IClower = M$summary.random$species3[, 4],
                       ICupper = M$summary.random$species3[, 6])

#plot
X11(width = 4, height = 4)

par(mfrow = c(1,1), bty="n",
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 2, 0.5, 1)
)


plot(0, bty = "l", labels = FALSE, tck = 0, xlim = c(-2.5,3), ylim = c(0.5,4.5), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_dt$mean, as.numeric(tab_dt$ID) - 0.2, col = "darkgoldenrod2", pch = 19, cex = 1.3)
arrows(tab_dt$IClower, as.numeric(tab_dt$ID) - 0.2,
       tab_dt$ICupper, as.numeric(tab_dt$ID) - 0.2,
       col = "darkgoldenrod2", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line

points(tab_wspt$mean, as.numeric(tab_wspt$ID) , col = "cornflowerblue", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) ,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) ,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 2) 

points(tab_wspt_var$mean, as.numeric(tab_wspt_var$ID) + 0.2, col = "pink1", pch = 19, cex = 1.3)
arrows(tab_wspt_var$IClower, as.numeric(tab_wspt_var$ID) + 0.2,
       tab_wspt_var$ICupper, as.numeric(tab_wspt_var$ID) + 0.2,
       col = "pink1", code = 3, length = 0.03, angle = 90, lwd = 2) 

axis(side= 1, at = c(-2,-1,0,1,2), labels = c(-2,-1,0,1,2), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4), 
     labels =  tab_dt$ID, 
     tick = T ,col = NA, col.ticks = 1, 
     tck = -.015 , 
     las = 2) 

#add legend
legend(x = 1.25, y = 4.5, legend = c( "Wind support var", "Wind support", expression(italic(paste(Delta,"T")))), 
       col = c("pink1", "cornflowerblue","darkgoldenrod1"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.9)


# ---------- Fig S2: boxplot comparing used and available steps #####

load("annotated_steps.RData") #ann_cmpl; This dataframe includes used and alternative steps and can be reproduced using step_generation.R

X11(width = 6, height = 7)

par(mfrow= c(2,2), 
    oma = c(0,0,3,0), 
    las = 1,
    mgp = c(0,1,0))

labels <- c(expression(italic(paste(Delta,"T"))), "Wind support")
variables <- c("delta_t", "wind_support")
v_variables <- c("delta_t_var", "wind_support_var")

for(i in 1:length(variables)){
  
  boxplot(ann_cmpl[,variables[i]] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c(alpha("darkgoldenrod1", 0.9),"gray"), bty = "n", cex = 0.8)
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1, variables[i]] ~ ann_cmpl[ann_cmpl$used == 1,"species"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = alpha("darkgoldenrod1", 0.9),  lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1, "species"])) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0, variables[i]] ~ ann_cmpl[ann_cmpl$used == 0, "species"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1 , "species"])) + 0.15)
  
}
mtext("Instantaneous values at each step", side = 3, outer = T, cex = 1.2)

for(i in 1:length(v_variables)){
  
  boxplot(ann_cmpl[,v_variables[i]] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c(alpha("darkgoldenrod1", 0.9),"gray"), bty = "n", cex = 0.8)
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1,v_variables[i]] ~ ann_cmpl[ann_cmpl$used == 1,"species"], outcol = alpha("black", 0.2),
          yaxt = "n",xaxt = "n", add = T, boxfill = alpha("darkgoldenrod1", 0.9), lwd = 0.7,  outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,v_variables[i]] ~ ann_cmpl[ann_cmpl$used == 0,"species"], outcol = alpha("black", 0.2),
          yaxt = "n",xaxt = "n", add = T, boxfill = "grey", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
} 

mtext("40-year variances at each step", side = 3, outer = T, cex = 1.2, line = -20)


# ---------- Fig S3: maps with annotated tracking points #####

#load data: annotated sea-crossing points as a move object
load("raw_sea_points_for_maps_mv.RData") #mv

#add shapefile as the map background
region <- st_read("continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -99, xmax = 144, ymin = -30, ymax = 71) %>%
  st_union()

#add a categorical variable for wind levels
breaks_w <- c(-20,-10,-5,0,5,10,15,35)
tags_w <- c("< -10","-10 to -5","-5 to 0","0 to 5","5 to 10","10 to 15", "> 15")
#add a categorical variable for delta_t levels
breaks_dt <- c(-5,-2,0,2,5,10)
tags_dt <- c("< -5","-5 to -2","0 to 2","2 to 5", "> 5")

#create color palettes and select colors for positive and negative values.
Pal_p <- colorRampPalette(c("darkgoldenrod2", "indianred1")) #colors for positive values
Pal_n <- colorRampPalette(c("mediumblue", "cornflowerblue")) #colors for negative values
Cols_w <- paste0(c(Pal_n(3),Pal_p(4)), "80") #add transparency
Cols_dt <- paste0(c(Pal_n(2),Pal_p(3)), "80")



df <- mv %>% 
  as.data.frame() %>% 
  drop_na(c("heading","delta_t")) %>% 
  mutate(wind_support= wind_support(u = u925, v = v925, heading = heading),
         cross_wind= cross_wind(u = u925, v = v925, heading = heading)) %>% 
  mutate(binned_w = cut(wind_support, breaks = breaks_w, include.lowest = T, right = F, labels = tags_w),
         binned_dt = cut(delta_t, breaks = breaks_dt, include.lowest = T, right = F, labels = tags_dt)) %>% 
  mutate(cols_w = as.factor(binned_w),
         cols_dt = as.factor(binned_dt)) 

levels(df$cols_w) <- Cols_w
levels(df$cols_dt) <- Cols_dt

df_sp <- SpatialPointsDataFrame(coords = df[,c("coords.x1", "coords.x2")], proj4string = wgs, data = df)

#plot
X11(width = 12, height = 11.5) #in inches
par(mfrow=c(2,1),
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font.axis = 3,
    cex.lab = 0.6,
    cex = 0.9,
    oma = c(0,0,1.5,0),
    mar = c(0, 0, 0.3, 0),
    lend = 1  #rectangular line endings (trick for adding the rectangle to the legend)
)


plot(region, col="#e5e5e5",border="#e5e5e5")
points(df_sp, pch = 1, col = as.character(df_sp$cols_w), cex = 0.2)

#add latitudes
clip(-99, 144, -30, 71)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
#text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")
text(x = -95, y = c(32,62), labels = c("30° N", "60° N"), cex = 0.6, col = "grey65")


#add a frame for the sub-plots and legend
rect(xleft = -100,
     xright = -71,
     ybottom =  -30,
     ytop = 2,
     col="white",
     border = NA)

text(x = -87,y = 0, "Wind support (m/s)", cex = 0.8, font = 3)
legend(x = -100, y = 0, legend = levels(df_sp$binned_w), col = Cols_w, pch = 20, 
       bty = "n", cex = 0.8, text.font = 3)
mtext("Sea-crossing tracks annotated with wind support", 3, outer = F, cex = 1.3, line = -0.5)

plot(region, col = "#e5e5e5",border = "#e5e5e5")
points(df_sp, pch = 1, col = as.character(df_sp$cols_dt), cex = 0.2)

#add latitudes
clip(-99, 144, -30, 71)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
text(x = -95, y = c(32,62), labels = c("30° N", "60° N"), cex = 0.6, col = "grey65")

#add a frame for the sub-plots and legend
rect(xleft = -100,
     xright = -80,
     ybottom =  -30,
     ytop = 2,
     col="white",
     border = NA)

text(x = -93,y = 0,  expression(italic(paste(Delta,"T", "(°C)"))), cex = 0.8)
legend(x = -100, y = -1, legend = levels(df_sp$binned_dt), col =Cols_dt, pch = 20, 
       bty = "n", cex = 0.8, text.font = 3)
mtext(bquote(paste('Sea-crossing tracks annotated with', italic(~ Delta *"T"))), 3, outer = F, cex = 1.3, line = -0.5)



# ---------- Fig S4: boxplots for conditions at sea-crossing initiation #####

#load data: annotated sea-crossing points. This is the same file as used in Fig. S3.
load("raw_sea_points_for_maps_mv.RData") #mv

#extract first point of each track
starts_all <- mv %>% 
  as.data.frame() %>% 
  drop_na(c("heading","delta_t")) %>% 
  mutate(wind_support = wind_support(u = u925, v = v925, heading = heading),
         cross_wind = cross_wind(u = u925, v = v925, heading = heading)) %>% 
  mutate(group = ifelse(species == "OHB", "OHB",
                        ifelse(species == "GFB", "GFB",
                               ifelse(species == "O" & coords.x1 < -30, "O_A",
                                      ifelse(species == "O" & coords.x1 > -30, "O_E",
                                             ifelse(species == "EF" & coords.x2 > 0, "EF_med",
                                                    ifelse(species == "EF" & coords.x2 < 0, "EF_moz",
                                                           ifelse(species == "PF" & coords.x1 < -30, "PF_A",
                                                                  "PF_E"))))))))  %>% 
  group_by(track) %>% 
  arrange(timestamp) %>% 
  slice(1) %>%
  ungroup() %>% 
  st_as_sf(coords = c("coords.x1", "coords.x2"), crs = wgs) %>% 
  mutate(s_elev_angle = solarpos(st_coordinates(.), timestamp, proj4string = CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                           ifelse(s_elev_angle > 40, "high", "low"))) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  mutate(group = as.factor(group)) %>% 
  as.data.frame()

#relevel the group variable, for nice plots; longitudinal order
starts_all$group <- fct_relevel(starts_all$group, c("O_A", "PF_A", "O_E", "PF_E", "EF_med", "EF_moz", "GFB", "OHB"))

#color palette
Pal <- colorRampPalette(c("darkgoldenrod1","lightpink1", "mediumblue")) 
Cols <- paste0(Pal(3), "80") #add transparency.

labels <- c("Pandion haliaetus \n (America)", "Falco peregrinus \n (America)", "Pandion haliaetus \n (Europe)", "Falco peregrinus \n (Europe)",
            "Falco eleonorae \n (Mediterranean)",  "Falco eleonorae \n (Mozambique)", "Butastur indicus \n" ,"Pernis ptilorhynchus \n")

variables <- c("wind_support", "delta_t")

#plot
X11(width = 12, height = 5.5) #inches

par(mfrow= c(2,1), 
    oma = c(3,0,2.5,0), 
    mar = c(0.5,4,0.2,0.2),
    las = 1,
    bty = "l",
    cex.axis = 0.8,
    font.axis = 3,
    tck = -0.015,
    mgp = c(1,1,0))



for(i in 1:length(variables)){
  
  boxplot(starts_all[,variables[i]] ~ starts_all[,"group"], data = starts_all, boxfill = NA, border = NA, xlab = "", ylab = "", xaxt = "n")
  abline(h = 0, lty = 2, col = alpha("black", 0.6), lwd = 0.3)

  points(jitter((as.numeric(starts_all[starts_all$sun_elev == "high","group"]) -0.28),0.3), starts_all[starts_all$sun_elev == "high", variables[i]], 
         yaxt = "n", xaxt = "n", pch = 20, cex = 0.8, col = alpha("black", 0.6))
  
  boxplot(starts_all[starts_all$sun_elev == "high", variables[i]] ~ starts_all[starts_all$sun_elev == "high","group"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = Cols[1], outline = FALSE, lwd = 0.5, 
          boxwex = 0.25, at = 1:length(unique(starts_all$group)) -0.28)
  
  points(jitter(as.numeric(starts_all[starts_all$sun_elev == "low","group"]),0.3), starts_all[starts_all$sun_elev == "low", variables[i]], 
         yaxt = "n", xaxt = "n", pch = 20, cex = 0.8, col = alpha("black", 0.6))
  
  boxplot(starts_all[starts_all$sun_elev == "low", variables[i]] ~ starts_all[starts_all$sun_elev == "low", "group"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = Cols[2], outline = FALSE,  lwd = 0.5, 
          boxwex = 0.25, at = 1:length(unique(starts_all$group))+ 0)
  
  points(jitter((as.numeric(starts_all[starts_all$sun_elev == "night","group"]) +0.28),0.3), starts_all[starts_all$sun_elev == "night", variables[i]], 
         yaxt = "n", xaxt = "n", pch = 20, cex = 0.8, col = alpha("black", 0.6))
  
  boxplot(starts_all[starts_all$sun_elev == "night", variables[i]] ~ starts_all[starts_all$sun_elev == "night", "group"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = Cols[3], outline=FALSE,  lwd = 0.5, 
          boxwex = 0.25, at = 1:length(unique(starts_all$group))+ 0.28)
  
  
  if(i == length(variables)){
    axis(side = 1, at = 1:length(levels(starts_all$group)), labels = labels, 
         tick = T , col.ticks = 1, col = NA, tck = -.015, lwd = 0, lwd.ticks = 1, cex = 0.9)
  }
  
  if(variables[i] == "wind_support"){
    mtext("Wind support (m/s)", side = 2, las = 3, line = 2, font = 3)
  }
  
  if(variables[i] == "delta_t"){
    mtext(expression(italic(paste(Delta,"T", "(°C)"))), side = 2, las = 3, line = 2)
  }
  
}
mtext("Atmospheric conditions at the start of sea-crossing tracks", side = 3, outer = T, cex = 1.2, font = 1, line = 0.7)

legend(x = 7.5, y = 21.6, legend = c("daytime: high sun", "daytime: low sun", "night"), col = Cols,
       cex = 0.78, pt.cex = 1.1, bg = "white", bty = "n", pch = 15, horiz = F, xpd = NA)





# SESSION INFO
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 20.10
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scales_1.1.1       itsadug_2.4        plotfunctions_1.4  TeachingDemos_2.12 fields_11.6        spam_2.6-0         dotCall64_1.0-1    mgcv_1.8-34        nlme_3.1-152       corrr_0.4.3       
# [11] INLA_21.02.23      foreach_1.5.1      Matrix_1.3-2       maptools_1.1-1     sf_0.9-8           move_4.0.6         rgdal_1.5-23       raster_3.4-5       sp_1.4-5           geosphere_1.5-10  
# [21] lubridate_1.7.10   forcats_0.5.1      stringr_1.4.0      dplyr_1.0.5        purrr_0.3.4        readr_1.4.0        tidyr_1.1.3        tibble_3.1.0       ggplot2_3.3.3      tidyverse_1.3.0   
# 
# loaded via a namespace (and not attached):
#   [1] httr_1.4.2         maps_3.3.0         jsonlite_1.7.2     splines_4.0.3      modelr_0.1.8       assertthat_0.2.1   cellranger_1.1.0   pillar_1.5.1       backports_1.2.1    lattice_0.20-41   
# [11] glue_1.4.2         rvest_1.0.0        colorspace_2.0-0   pkgconfig_2.0.3    broom_0.7.6        haven_2.3.1        proxy_0.4-25       farver_2.1.0       generics_0.1.0     ellipsis_0.3.1    
# [21] cachem_1.0.4       withr_2.4.1        cli_2.4.0          magrittr_2.0.1     crayon_1.4.1       readxl_1.3.1       memoise_2.0.0      fs_1.5.0           fansi_0.4.2        xml2_1.3.2        
# [31] foreign_0.8-81     class_7.3-18       tools_4.0.3        hms_1.0.0          lifecycle_1.0.0    munsell_0.5.0      reprex_2.0.0       packrat_0.6.0      compiler_4.0.3     e1071_1.7-6       
# [41] rlang_0.4.10       classInt_0.4-3     units_0.7-1        iterators_1.0.13   rstudioapi_0.13    gtable_0.3.0       codetools_0.2-18   DBI_1.1.1          R6_2.5.0           fastmap_1.1.0     
# [51] utf8_1.2.1         KernSmooth_2.23-18 stringi_1.5.3      Rcpp_1.0.6         vctrs_0.3.7        dbplyr_2.1.1       tidyselect_1.1.0  