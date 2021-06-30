# Scripts for step selection function analysis
# This is script 2 of 6 for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------

#to do: add the wind support functions to functions.R

library(tidyverse)
library(lubridate)
library(INLA)
library(corrr)
library(raster)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

# ---------- STEP 1: load data #####

#load annotated data prepared in step_generation.R (available on the Dryad repository)
#for public: load(annotated_steps.RData)
load("2021/public/ssf_input_ann_cmpl_60_15.RData") #ann_cmpl; This dataframe includes used and alternative steps and can be reproduced using step_generation.R

# ---------- STEP 2: look into autocorrelation #####

ann_cmpl %>% 
  dplyr::select(c("delta_t", "wind_speed", "wind_support", "wind_support_var", "abs_cross_wind", "delta_t_var","step_length")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6)

#corr test to include in the paper
cor.test(ann_cmpl$wind_support_var, ann_cmpl$delta_t_var)

# ---------- STEP 3: z-transformation (i.e. scale all predictor variables) #####

#z-transform
all_data <- ann_cmpl %>% 
  mutate_at(c("delta_t", "wind_speed", "wind_support", "wind_support_var", "abs_cross_wind", "delta_t_var"),
            list(z = ~(scale(.)))) 


#remove: save(all_data, file = "2021/public/ssf_input_ann_60_15_z.RData")

# ---------- STEP 4: step selection analysis in INLA #####

#repeat variabels that will be used as random slopes
all_data <- all_data %>% 
  mutate(species1 = factor(species),
         species2 = factor(species),
         species3 = factor(species),
         species4 = factor(species),
         ind1 = factor(ind),
         ind2 = factor(ind),
         ind3 = factor(ind),
         ind4 = factor(ind))


# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

#prepare new data to be used for predictions (for effect plots)
#see Gómez-Rubio 2020 for details of prediction with INLA models (i.e. imputation of missing values)


#for the new data, randomly select 4 species and one stratum per species

n <- 2448 #corresponding to 48 unique strata, from 4 individuals per species, 3 strata per individual. 

inds <- all_data %>% 
  group_by(species) %>%
  distinct(ind) %>% 
  sample_n(4, replace = F)

strata <- all_data %>% 
  filter(ind %in% c(inds$ind)) %>% 
  group_by(ind) %>% 
  distinct(stratum) %>% 
  sample_n(3, replace = F)


new <-  data.frame(used = rep(NA,n),
                   delta_t_z = sample(seq(min(all_data$delta_t_z),max(all_data$delta_t_z), length.out = 5), n, replace = T), #regular intervals for wind support and delta t, so we can make a raster later on
                   wind_support_z = sample(seq(min(all_data$wind_support_z),max(all_data$wind_support_z), length.out = 10), n, replace = T),
                   stratum = factor(rep(strata$stratum, 51)),
                   species1 = factor(rep(inds$species, 51)),
                   ind1 = factor(rep(inds$ind, 51))) %>% 
  mutate(species2 = species1,
         species4 = species1,
         ind2 = ind1,
         ind4 = ind1) 

new_data <- new %>% 
  full_join(all_data)

save(new_data, file = "2021/public/new_data_48str_lowres.RData")


#alternative new data: add one new row to unique strata instead of entire empty copies of strata
n <- 500
new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA) %>% 
  full_join(all_data)

save(new_data, file = "2021/public/new_data_500n_irregular.RData")

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA) %>% 
  full_join(all_data)

save(new_data, file = "2021/public/new_data_500n_irregular.RData")


#Model formula
formulaM <- used ~ -1 + delta_t_z * wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

#Model
      (b <- Sys.time())
      M_pred <- inla(formulaM, family = "Poisson", 
                          control.fixed = list(
                            mean = mean.beta,
                            prec = list(default = prec.beta)),
                          data = new_data, 
                          num.threads = 10,
                          control.predictor = list(compute = TRUE), #this means that NA values will be predicted. link can also be set. but will make the predictions Inf (response is binary but family is poisson.. what is the correct link!!??) # “apply the first link function to everything”.
                          control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))#, mlik = T, waic = T)) 
      Sys.time() - b #1.574914 hours; 49 min; 1.6 min; 2.2 hrs
      #with link = 1, all NaN and Inf values

save(M_pred, file = "2021/public/inla_models/M_preds_all_48_lowres.RData")

#extract predicted values
used_na <- which(is.na(new_data$used))
M_pred$summary.fitted.values[used_na,]


#create a raster of predictions
preds <- data.frame(delta_t = new_data[is.na(new_data$used) ,"delta_t_z"],
                    wind_support = new_data[is.na(new_data$used) ,"wind_support_z"],
                    preds = M_pred$summary.fitted.values[used_na,"mean"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds)))


plot(preds$delta_t, preds$wind_support, col = as.factor("preds"))


avg_preds <- preds %>% 
  group_by(delta_t, wind_support) %>% 
  summarise(avg_pres = mean(prob_pres)) %>% 
  ungroup() %>% 
  mutate(wspt_backtr = wind_support * attr(all_data$wind_support_z, 'scaled:scale') + attr(all_data$wind_support_z, 'scaled:center'),
         dt_backtr = delta_t * attr(all_data$delta_t_z, 'scaled:scale') + attr(all_data$delta_t_z, 'scaled:center')) %>% 
  dplyr::select(-c("delta_t","wind_support")) %>% 
  as.data.frame()



coordinates(avg_preds) <-~ dt_backtr + wspt_backtr 
gridded(avg_preds) <- TRUE
r <- raster(avg_preds)


plot(r, ylab = "wind support (m/s)", xlab = "delta_t (°C)")

avg <- preds %>% 
  group_by(delta_t, wind_support) %>% 
  summarise(avg_pres = mean(prob_pres)) %>% 
  ungroup() %>% 
  as.data.frame()


coordinates(avg) <-~ delta_t + wind_support 
gridded(avg) <- TRUE
rr <- raster(avg)


plot(rr, ylab = "wind support (m/s)", xlab = "delta_t (°C)")

#plot

#create a color palette
cuts <- seq(-1,5,0.01) #set breaks
#pal <- colorRampPalette(c("dodgerblue","darkturquoise","goldenrod1","coral","firebrick1"))
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)



#Model
formulaM2 <- used ~ -1 + delta_t_z * wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
M2 <- inla(formulaM2, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           control.inla = list(force.diagonal = T),
           data = all_data,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = T))
Sys.time() - b #44.05078 mins

#Eextract cpo and waic
#-sum(log(M2$cpo$cpo))
#M2$waic$waic

save(M2, file = "2021/public/inla_models/m2_60_15.RData")


#make predictions to visualize the interaction term
load("2021/public/inla_models/m2_60_15.RData")

#create a sample for testing model formulas in short time
str_to_keep <- sample(unique(all_data$stratum),200)
sample <- all_data[all_data$stratum %in% str_to_keep,]


n <- 204
#for the new data, randomly select 4 species and one stratum per species
inds <- all_data %>% 
  group_by(species) %>%
  summarize(ind = sample(ind, 1))
strata <- all_data %>% 
  filter(ind %in% c(inds$ind)) %>% 
  group_by(ind) %>% 
  summarise(stratum = sample(stratum, 1))
  
# new <-  data.frame(used = rep(NA,n),
#                    delta_t_z = as.numeric(sample(all_data$delta_t_z, n, replace = T)),
#                    wind_support_z = sample(all_data$wind_support_z, n, replace = T),
#                    stratum = factor(rep(strata$stratum, 51)),
#                    species1 = factor(rep(inds$species, 51)),
#                    ind1 = factor(rep(inds$ind, 51))) %>% #choose a handful of individuals
#   mutate(species2 = species1,
#          species4 = species1,
#          ind2 = ind1,
#          ind4 = ind1) 
  
#regular intervals for wind support and delta t
new <-  data.frame(used = rep(NA,n),
                   delta_t_z = sample(seq(min(all_data$delta_t_z),max(all_data$delta_t_z), length.out = 6), n, replace = T), #sample from 10 equally spaced values of delta_t
                   wind_support_z = sample(seq(min(all_data$wind_support_z),max(all_data$wind_support_z), length.out = 6), n, replace = T),
                   stratum = factor(rep(strata$stratum, 51)),
                   species1 = factor(rep(inds$species, 51)),
                   ind1 = factor(rep(inds$ind, 51))) %>% #choose a handful of individuals
  mutate(species2 = species1,
         species4 = species1,
         ind2 = ind1,
         ind4 = ind1) 
new_data <- new %>% 
  full_join(all_data)


#method 1: (Virgilio's book section 12.3)



#transform marginals manually
fitted.values.mean <- numeric(n)

for(i in 1:n) {
  
  if (is.na(new_data$used[i])) {
    
    if (FALSE) {
      
      ## either like this, which is slower...
      
      marg = inla.marginal.transform( 
        
        function(x) exp(x)/(1+exp(x)), 
        
        result$marginals.fitted.values[[i]] )
      
      fitted.values.mean[i] = inla.expectation(
        
        function(x) x, marg)
      
    } else {
      
      ## or like this,  which is much faster...
      
      fitted.values.mean[i] = inla.expectation(
        
        function(x) exp(x)/(1 +exp(x)), 
        
        result$marginals.fitted.values[[i]])
      
    }
    
  } else {
    
    fitted.values.mean[i] = result$summary.fitted.values[i,"mean"]
    
  }
  
}


plot(fitted.values.mean)




m <- M2_pred$summary.linear.predictor[used_na,]

round(M2_pred$summary.linear.predictor[used_na,],3)


#plot interaction effects
ws_low <- which(is.na(new_data$used) & new_data$wind_speed_z == -1)
ws_high <- which(is.na(new_data$used) & new_data$wind_speed_z == 1)
ws_zero <- which(is.na(new_data$used) & new_data$wind_speed_z == 0)

X11();par(mfrow = c(1,3))
plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_sample$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = -1" ,
      xlab="delta_t" )
points(new_data[ws_low,"delta_t_z"], m1d_sample$summary.fitted.values[ws_low,"mean"])

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_sample$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = 0" ,
      xlab="delta_t" )
points(new_data[ws_zero,"delta_t_z"], m1d_sample$summary.fitted.values[ws_zero,"mean"])

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_sample$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = 1" ,
      xlab="delta_t" )
points(new_data[ws_high,"delta_t_z"], m1d_sample$summary.fitted.values[ws_high,"mean"])





#method2: https://www.flutterbys.com.au/stats/tut/tut12.10.html



new_data <- data.frame(delta_t_z = sample(all_data$delta_t_z, 204, replace = T), 
                       wind_support_z = sample(all_data$wind_support_z, 200, replace = T),
                       stratum = sample(all_data$stratum,4, replace = F),)

newdata <- data.frame(x = seq(min(data$x, na.rm=TRUE), max(data$x, na.rm=TRUE), len=100))

Xmat <- model.matrix(~ -1 + delta_t_z * wind_support_z + wind_support_var_z +
                       f(stratum, model = "iid", 
                         hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
                       f(species1, delta_t_z, model = "iid", 
                         hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
                       f(species2, wind_support_z,  model = "iid",
                         hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
                       f(species4, wind_support_var_z, model = "iid",
                         hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
                       f(ind1, delta_t_z, model = "iid",
                         hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
                       f(ind2, wind_support_z,  model = "iid",
                         hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
                       f(ind4, wind_support_var_z, model = "iid",
                         hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))), 
                     data = newdata)


lincomb <- inla.make.lincombs(as.data.frame(Xmat))
#fit the model
data.inla1 <- inla(y~x, lincomb=lincomb,  #include the linear combinations
                   data=data,
                   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                   control.inla=list(lincomb.derived.only=FALSE))
#examine the regular summary - not there are no changes to the first fit
#summary(data.inla1)



#----------------------------------------------------
#Model 3: remove wind speed and delta t
formulaM3 <- used ~ -1 + delta_t_z * wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind4, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


(b <- Sys.time())
M3 <- inla(formulaM3, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           control.inla = list(force.diagonal = T),
           data = all_data,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = T))
Sys.time() - b #27.32572 mins


save(M3, file = "2021/public/inla_models/m3_60_15.RData")

-sum(log(M3$cpo$cpo))
M3$waic$waic

#------------------------------
#method 3:
#instead of making predictions, we can also sample from the posterior distribution. (following Virgilio's book  section 2.7, 12.4)
#make sure config = TRUE in the inla call
samp1 <- inla.posterior.sample(100, M2, selection = list(wind_support_z= 1, delta_t_z = 1)) #selection can be used to determine which variables we want
#1 in the above call means that we want to keep the first element in effect 
samp2 <- inla.posterior.sample.eval(function(...) {wind_support_z * delta_t_z}, samp1)
summary(as.vector(samp2))
