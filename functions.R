
#used in step_generation.R
#taken from NCEP pacakge, but produces NA values instead of an error when heading cannot be calculated
NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}

#used in all_figures.R
CubeRoot<-function(x){
  sign(x)*abs(x)^(1/3)
}

#used in all_figures.R
w_star <- function(g = 9.81, blh, T2m, s_flux, m_flux) {
  
  z <- blh
  T_k_2m <- T2m
  T_c_2m <- T_k_2m - 273.15
  Thetav_k_z <- (T_k_2m) + 0.006 * z
  wT <- (s_flux * -1) / 1013 / 1.2 #reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1) *1000 /1.2 #reverse the sign. ECMWF upward fluxes are negative
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}

#used in all_figures.R
reg_gam_plot <- function(x){
  m <- models_ls[[x]]$gam
  t <- timing_areas[[x]]
  
  plot(0, type = "n", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "", ylab = "")#, main = names[[x]]) #expression(paste(Delta,"T")), main = x
  abline(h = 0, col = "gray60",lty = 1, lwd = 0.5)
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view = "yday", plot_all = "sun_elev_f", rm.ranef = F, lwd = 1.5, #ylim=c(-2,5),
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
  #title(ylab = expression(italic(paste(Delta,"T"))), line = 2.5 ,las = 2, cex = 0.7)
  
}



#-------------------
#wind support and crosswind codes. from Safi et al (2013; Movement Ecology)
#used in all_figures.R and step_generation.R

check.libs <- function()
{
  if (!require(move))
  {install.packages(move)}
  if (!require(tcltk))
  {install.packages(tcltk)}
  if (!require(mapdata))
  {install.packages(mapdata)}
  if (!require(mgcv))
  {install.packages(mgcv)}
}
check.track <- function(x)
{
  if(sum(length(grep("U.Component", names(x))), length(grep("V.Component", names(x))))!=2)
  {
    stop("No U and V components found!")
  }
  if(sum(as.numeric(names(x) %in% c("ground.speed", "heading")))!=2)
  {
    message("No heading and ground speed were provided and were calculated from next locations!")
    x$heading <- unlist(lapply(lapply(split(x), angle), "c", NA))
    x$ground.speed <- unlist(lapply(lapply(split(x), speed), "c", NA))
    return(x)
  }else{
    message("Analysis will be based on heading and ground speed provided.")
    return(x)
  }
}
wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
  14
}
cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}
airspeed <- function(x)
{
  va <- sqrt((x$ground.speed - x$ws)^2 + (x$cw)^2)
  return(va)
}
make.plot <- function(f, x.df, modgs, modas)
{
  # Get the predicted values depending on one single variable while
  # keeping the other wind component at 0 and the longitude and latitude
  #at the mean of the entire data set
  # For ground speed
  pcw.gs <- predict(modgs$gam,
                    newdata=data.frame(ws=rep(0,100),
                                       cw=seq(min(x.df$cw, na.rm=T),
                                              max(x.df$cw, na.rm=T),
                                              length.out=100),
                                       location.long=rep(mean(x.df$location.long, na.rm=T), 100),
                                       location.lat=mean(x.df$location.lat, na.rm=T)))
  pws.gs <- predict(modgs$gam,
                    newdata=data.frame(ws=seq(min(x.df$ws, na.rm=T),
                                              max(x.df$ws, na.rm=T),
                                              length.out=100),
                                       cw=rep(0,100),
                                       location.long=rep(mean(x.df$location.long, na.rm=T), 100),
                                       location.lat=mean(x.df$location.lat, na.rm=T)))
  15
  # For airspeed
  pcw.as <- predict(modas$gam,
                    newdata=data.frame(ws=rep(0,100),
                                       cw=seq(min(x.df$cw, na.rm=T),
                                              max(x.df$cw, na.rm=T),
                                              length.out=100),
                                       location.long=rep(mean(x.df$location.long, na.rm=T), 100),
                                       location.lat=mean(x.df$location.lat, na.rm=T)))
  pws.as <- predict(modas$gam,
                    newdata=data.frame(ws=seq(min(x.df$ws, na.rm=T),
                                              max(x.df$ws, na.rm=T),
                                              length.out=100), cw=rep(0,100),
                                       location.long=rep(mean(x.df$location.long, na.rm=T), 100),
                                       location.lat=mean(x.df$location.lat, na.rm=T)))
  outfile <- paste(tk_choose.dir(default = "", caption = "Select directory for output"),
                   sub(".csv", ".pdf",
                       paste("Report_", grep(".csv",
                                             unlist(strsplit(f, "/")), value=T),
                             sep="")),
                   sep="/")
  # Open a pdf file and plot all into it
  pdf(outfile)
  # First page split in 4 areas
  layout(matrix(c(1,1,1,1,5,2,3,6,5,4,4,6),
                ncol=4, byrow=T),
         height=c(0.1,0.45,0.45),
         width=c(0.1,0.4,0.4,0.1))
  par(mar=c(0,0,0,0))
  # Write title
  plot.new()
  plot.window(xlim=c(-1,1), ylim=c(-0.1,0.1))
  text(0,-0.1,paste("Report for file ",
                    grep(".csv", unlist(strsplit(f, "/")), value=T),
                    sep=""), pos=3, font=2, cex=1.5)
  16
  # Make a histogram of the airspeed raw values
  par(mar=c(4,4,2,4))
  hist((data.frame(x.df)[,"airspeed"]), xlab="Airspeed", main=NA)
  # Make a histogram of ground speed raw values
  hist((data.frame(x.df)[,"ground.speed"]), xlab="Ground speed", main=NA)
  # Write the means of these measures in the bottom section
  plot(0,0, type="n", xlim=c(-10,10), ylim=c(-10,10),
       xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
  text(0,0,paste("Mean airspeed: ", round(mean(x.df$airspeed, na.rm=T), 4),
                 " m/s.\n", "Mean ground speed: ",
                 round(mean(x.df$ground.speed, na.rm=T), 4),
                 " m/s.", sep=""), cex=1.5)
  # Next page split in 4 regions
  layout(matrix(c(1,2,3,4), ncol=2, byrow=T))
  par(mar=c(5, 4, 4, 2) + 0.1)
  # plot airspeed (square root) transformed and add the
  # predicted model estimates as a function of wind support as a line
  plot(sqrt(airspeed)~ws, data=x.df, pch=16, col="grey95",
       xlab="Wind support", ylab=expression(paste(sqrt("Airspeed"), sep="")))
  points(sqrt(airspeed)~ws, data=x.df, col="grey60")
  lines(seq(min(x.df$ws, na.rm=T), max(x.df$ws, na.rm=T), length.out=100), pws.as, col="red")
  mtext(substitute(paste(sqrt("Airspeed="), i%*%"ws", " + ", j, sep=""),
                   list(i=round(summary(modas$gam)$p.coeff["ws"], 4),
                        j=round(summary(modas$gam)$p.coeff["(Intercept)"], 4))),
        side=4, cex=0.5, , col="red", line=0)
  # plot airspeed (square root) transformed and add the
  # predicted model estimates as a function of cross wind as a line
  plot(sqrt(airspeed)~cw, data=x.df, pch=16, col="grey95",
       xlab="Cross wind", ylab=expression(paste(sqrt("Airspeed"), sep="")))
  points(sqrt(airspeed)~cw, data=x.df, col="grey60")
  lines(seq(min(x.df$ws, na.rm=T), max(x.df$ws, na.rm=T), length.out=100), pcw.as, col="red")
  mtext(substitute(paste(sqrt("Airspeed="), i%*%"ws", " + ", j, sep=""),
                   list(i=round(summary(modas$gam)$p.coeff["ws"], 4),
                        j=round(summary(modas$gam)$p.coeff["(Intercept)"], 4))),
        side=4, cex=0.5, , col="red", line=0)
  # plot ground speed (square root) transformed and add the
  # predicted model estimates as a function of wind support as a line
  plot(sqrt(ground.speed)~ws, data=x.df, pch=16, col="grey95",
       xlab="Wind support", ylab=expression(paste(sqrt("Ground speed"), sep="")))
  points(sqrt(ground.speed)~ws, data=x.df, col="grey60")
  lines(seq(min(x.df$ws, na.rm=T), max(x.df$ws, na.rm=T),
            length.out=100), pws.gs, col="red")
  mtext(substitute(paste(sqrt("Ground speed="), i%*%"ws", " + ", j, sep=""),
                   list(i=round(summary(modgs$gam)$p.coeff["ws"], 4),
                        j=round(summary(modgs$gam)$p.coeff["(Intercept)"], 4))),
        side=4, cex=0.5, , col="red", line=0)
  # plot ground speed (square root) transformed and add the
  # predicted model estimates as a function of cross wind as a line
  plot(sqrt(ground.speed)~cw, data=x.df, pch=16, col="grey95",
       xlab="Cross wind", ylab=expression(paste(sqrt("Ground speed"), sep="")))
  points(sqrt(ground.speed)~cw, data=x.df, col="grey60")
  lines(seq(min(x.df$cw, na.rm=T), max(x.df$cw, na.rm=T),
            length.out=100), pcw.gs, col="red")
  mtext(substitute(paste(sqrt("Ground speed="), i%*%"cw", " + ", j, sep=""),
                   list(i=round(summary(modgs$gam)$p.coeff["cw"], 4),
                        j=round(summary(modgs$gam)$p.coeff["(Intercept)"], 4))),
        side=4, cex=0.5, col="red", line=0)
  # next pages
  layout(c(1))
  par(mar=c(5, 4, 4, 2))
  # capture stats output and write it into a plot region
  out <- capture.output(summary(modas$gam))
  plot(0,0, type="n", xlim=c(-10,10), ylim=c(-10,10),
       xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
  title("GAMM summary for airspeed")
  text(0,0,"Unsupervised !", cex=5, col="grey70", srt=45)
  text(-10,0,paste(out, collapse="\n"), pos=4, cex=0.5, family="mono")
  # capture stats output and write it into a plot region
  18
  out <- capture.output(summary(modas$lme))
  grep("Fixed effects:", out)
  plot(0,0, type="n", xlim=c(-10,10), ylim=c(-10,10),
       xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
  title("LME summary 1 of 2 for airspeed")
  text(0,0,"Unsupervised !", cex=5, col="grey70", srt=45)
  text(-10,0, paste(out[1:19], collapse="\n"), pos=4, cex=0.5, family="mono")
  # capture stats output and write it into a plot region
  plot(0,0, type="n", xlim=c(-10,10), ylim=c(-10,10),
       xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
  title("LME summary 2 of 2 for airspeed")
  text(0,0,"Unsupervised !", cex=5, col="grey70", srt=45)
  text(-10,0, paste(out[20:length(out)], collapse="\n"), pos=4, cex=0.5, family="mono")
  dev.off()
  # report successful production
  message(paste("Output report was saved to ", outfile, ".", sep=""))
}
make.report <- function(zipfile="missing", vlim=c(4,25))
{
  check.libs
  require(move)
  require(mapdata)
  require(tcltk)
  require(mgcv)
  if(zipfile=="missing")
  {
    zipfile <- tk_choose.files(default = "",
                               caption = "Select annotated zip file")
  }
  f <- unzip(zipfile, exdir=tempdir(), overwrite=T)
  x <- move(f[grep(".csv", f)])
  x <- check.track(x)
  df <- data.frame(x)
  19
  x$ws <- wind_support(df[,grep("U.Component", names(df))],
                       df[,grep("V.Component", names(df))], df$heading)
  x$cw <- cross_wind(df[,grep("U.Component", names(df))],
                     df[,grep("V.Component", names(df))], df$heading)
  x$airspeed <- airspeed(x)
  x.df <- as(x,"data.frame")
  message(paste("The ground speed range will be limited to between ",
                vlim[1], " and ", vlim[2], " m/s.", sep=""))
  x.df <- x.df[x.df$ground.speed>vlim[1] & x.df$ground.speed<vlim[2],]
  modgs <- gamm(sqrt(ground.speed)~s(location.long, location.lat)+cw*ws,
                random=list(individual.local.identifier=~1),
                data=x.df)
  modas <- gamm(sqrt(airspeed)~s(location.long, location.lat)+cw*ws,
                random=list(individual.local.identifier=~1),
                data=x.df)
  make.plot(f, x.df, modgs, modas)
}