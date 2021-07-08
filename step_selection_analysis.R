# Scripts for step selection function analysis
# This is script 2 of 4 for reproducing the results of Nourani et al 2021, ProcB.
# session info is provided at the end of script 4 (all_figures.R)
# Elham Nourani, PhD. June. 2021; enourani@ab.mpg.de
#-----------------------------------------------------------------

#to do: add the wind support functions to functions.R

library(tidyverse)
library(lubridate)
library(INLA)
library(corrr)
library(raster)


# ---------- STEP 1: load data #####

#load annotated data (available on the Dryad repository)
load("annotated_steps.RData") #ann_cmpl; This dataframe includes used and alternative steps and can be reproduced using step_generation.R


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

#to be able to produce effect plots for the interaction of wind support and delta-t, we need to add rows to the dataset where the dependent variable is set to NA
#see Gómez-Rubio 2020 for details of prediction with INLA models (i.e. imputation of missing values)

#add one new row to unique strata instead of entire empty copies of strata. assign wind and delta t values on a regular grid (tried with irregular, but range of predictions was off)
set.seed(200)

n <- 500
new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA,
         delta_t_z = sample(seq(min(all_data$delta_t_z),max(all_data$delta_t_z), length.out = 10), n, replace = T), #regular intervals for wind support and delta t, so we can make a raster later on
         wind_support_z = sample(seq(min(all_data$wind_support_z),max(all_data$wind_support_z), length.out = 10), n, replace = T),
         wind_support_var_z = sample(seq(min(all_data$wind_support_var_z),max(all_data$wind_support_var_z), length.out = 10), n, replace = T)) %>% 
  full_join(all_data)


#The new_data dataframe is available on the Dryad repository under name: new_data_for_modeling.RData

#Model formula
formulaM <- used ~ -1 + delta_t_z * wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

#Model
(b <- Sys.time())
M <- inla(formulaM, family = "Poisson", 
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = all_data, 
               num.threads = 10,
               control.predictor = list(compute = TRUE, link = 1), 
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b 

#This model is available on the Dryad repository under name: INLA_model.RData

#Model for predictions
(b <- Sys.time())
M_pred <- inla(formulaM, family = "Poisson", 
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data, 
               num.threads = 10,
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b 


#This model is available on the Dryad repository under name: INLA_model_preds.RData

#to plot the predictions and coefficients, see all_figures.R (Fig. 3 and 4)




