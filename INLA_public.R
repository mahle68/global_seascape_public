#Codes for reproducing the results of step selection function for sea-crossing of terrestrial birds using INLA.
#Elham Nourani, PhD. Nov. 4. 2020. enourani@ab.mpg.de

library(tidyverse)
library(INLA)

#DATA PREP ####
#load the annotated dataset prepared for step selection function estimation. Make sure to set your working directory correctly.
load("INLA_input_public.RData") #data frame called annotated_data 

#z-transform
z_data <- annotated_data %>% 
  mutate_at(c("delta_t","wind_speed","wind_support","wind_support_var","delta_t_var"),
            list(z = ~scale(.))) %>%
  as.data.frame() 

#repeat variabels that will be used as random slopes (requirement of INLA package. Each variable can be used only once)
z_data <- z_data %>% 
  mutate(species1 = species,
         species2 = species,
         species3 = species,
         species4 = species,
         species5 = species,
         ind1 = factor(ind),
         ind2 = factor(ind),
         ind3 = factor(ind),
         ind4 = factor(ind),
         ind5 = factor(ind),
         stratum = factor(stratum)) 


#MODELING ####
#see Muff et al 2019 (Journal of Animal Ecology) for details of choice of priors and model fitting.

#set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

formulaM1 <- used ~ -1 + delta_t_z * wind_speed_z + delta_t_var_z + wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species5, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind5, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


M1 <- inla(formulaM1, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = z_data,
            num.threads = 10, #This depends on your computer how many threads you can let INLA use
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))


summary(M1)

#remove variance of delta t
formulaM2 <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


M2 <- inla(formulaM2, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = z_data,
            num.threads = 10,
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))

summary(M2)

#remove variance of wind support
formulaM3 <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


M3 <- inla(formulaM3, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = z_data,
            num.threads = 10,
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = F))


summary(M3)


#FIGURE 2: posterior means of fixed effects ####

ModelList <- list(M1,M2,M3)
graphlist<-list()
for(i in 1:length(ModelList)){
  model<-ModelList[[i]]
  
  graph<-as.data.frame(summary(model)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  graph$Model<-i
  graph$Factor<-rownames(graph)
  
  graphlist[[i]]<-graph
}

graph <- bind_rows(graphlist)

graph$Sig <- with(graph, ifelse(Lower*Upper>0, "*", ""))

graph$Model <- as.factor(graph$Model)

position <- ifelse(length(unique(graph$Model))  ==  1, "none", "right")

VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

min<-min(graph$Lower,na.rm = T)
max<-max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

X11(width = 3.5, height = 3)

par(mfrow=c(1,1), bty="n",
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 3.5, 0.5, 1),
    bty = "l"
)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-6,8), ylim = c(0,6.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph[graph$Model == 1, "Estimate"], graph[graph$Model == 1,"Factor_n"] - 0.2, col = "salmon2", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 1, "Lower"], graph[graph$Model == 1,"Factor_n"] - 0.2,
       graph[graph$Model == 1, "Upper"], graph[graph$Model == 1,"Factor_n"] - 0.2,
       col = "salmon2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(graph[graph$Model == 2, c("Estimate","Factor")], col = "palegreen3", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 2, "Lower"], graph[graph$Model == 2,"Factor_n"],
       graph[graph$Model == 2, "Upper"], graph[graph$Model == 2,"Factor_n"],
       col = "palegreen3", code = 3, length = 0.03, angle = 90)

points(graph[graph$Model == 3, "Estimate"], graph[graph$Model == 3,"Factor_n"] + 0.2, col = "steelblue1", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 3, "Lower"], graph[graph$Model == 3,"Factor_n"] + 0.2,
       graph[graph$Model == 3, "Upper"], graph[graph$Model == 3,"Factor_n"] + 0.2,
       col = "steelblue1", code = 3, length = 0.03, angle = 90)
#add axes
axis(side= 1, at= c(-5,0,5), labels= c("-5", "0", "5"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:6), 
     labels= c( expression(paste(Delta,"t"," : wind speed")),
                "wind support var","wind support",expression(paste(Delta,"t"," var")),"wind speed", expression(paste(Delta,"t"))),
     tick=T ,col = NA, col.ticks = 1, 
     tck=-.015 ,
     las=2 ) 

#add legend
legend(x = 5.3, y = 0.8, legend=c("model 3", "model 2", "model 1"), col = c("steelblue1","palegreen3","salmon2"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.75)


#SUPPLEMENTARY FIGURE 1: species-specific coefficients ####
#for the best model (M3); original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)
species_names <- unique(z_data$species)

tab_dt <- data.frame(ID = as.factor(M3$summary.random$species1$ID),
                     mean = M3$summary.random$species1$mean,
                     IClower = M3$summary.random$species1[, 4],
                     ICupper = M3$summary.random$species1[, 6])

tab_wspd <- data.frame(ID = as.factor(M3$summary.random$species2$ID),
                       mean = M3$summary.random$species2$mean,
                       IClower = M3$summary.random$species2[, 4],
                       ICupper = M3$summary.random$species2[, 6])

tab_wspt <- data.frame(ID = as.factor(M3$summary.random$species3$ID),
                       mean = M3$summary.random$species3$mean,
                       IClower = M3$summary.random$species3[, 4],
                       ICupper = M3$summary.random$species3[, 6])

X11(width = 3.5, height = 3)

par(mfrow = c(1,1), bty="n"
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 3.5, 0.5, 1),
    bty = "l"
)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-4,4), ylim = c(0,4.3), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_dt$mean, as.numeric(tab_dt$ID) - 0.2, col = "lightcoral", pch = 19, cex = 1.3)
arrows(tab_dt$IClower, as.numeric(tab_dt$ID) - 0.2,
       tab_dt$ICupper, as.numeric(tab_dt$ID) - 0.2,
       col = "lightcoral", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspd$mean, as.numeric(tab_wspd$ID), col = "yellowgreen", pch = 19, cex = 1.3)
arrows(tab_wspd$IClower, as.numeric(tab_wspd$ID),
       tab_wspd$ICupper, as.numeric(tab_wspd$ID),
       col = "yellowgreen", code = 3, length = 0.03, angle = 90) 

points(tab_wspt$mean, as.numeric(tab_wspt$ID) + 0.2, col = "paleturquoise2", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) + 0.2,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) + 0.2,
       col = "paleturquoise2", code = 3, length = 0.03, angle = 90) 

axis(side= 1, at= c(-2,0,2), labels= c("-2", "0", "2"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4),
     labels= c(expression(italic("F. eleonorae")), expression(italic("P. haliaetus")),
               expression(italic("P. ptilorhynchus")),
               expression(italic("F. peregrinus"))),
     tick=T ,col = NA, col.ticks = 1, 
     tck=-.015 , 
     las=2 ) 

#add legend
legend(x = 1.8, y = 0.6, legend=c("wind support", "wind speed", expression(paste(Delta,"t"))), 
       col = c("paleturquoise2","yellowgreen","lightcoral"),
       pch = 19, bg="white",bty="n", cex = 0.75)

