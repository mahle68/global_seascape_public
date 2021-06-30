# global_seascape_public
This repository contains R scripts and input files to reproduce the results and figures presented in "The interplay of wind and uplift facilitates over-water flight in facultative soaring birds" by Nourani et al. 

# ABSTRACT
Flying over the open sea is energetically costly for terrestrial birds. Despite this, over-water journeys of many terrestrial birds, sometimes hundreds of kilometers long, are uncovered by bio-logging technology. To understand how these birds afford their flights over the open sea, we investigated the role of atmospheric conditions, in the form of wind and uplift, in subsidizing over-water flight at the global scale. We first established that Delta_t, the temperature difference between sea surface and air, is a meaningful proxy for uplift over the open sea. Using this proxy, we showed that the spatio-temporal patterns of sea-crossing in terrestrial migratory birds is associated with favorable uplift conditions. We then analyzed route selection over the open sea for four facultative soaring species, representing all major migratory flyways worldwide. Our results showed that favorable uplift conditions, albeit not as common and powerful as over land, are not rare over the open sea. As expected, birds maximized wind support when selecting their sea-crossing routes, but  also selected higher uplift when suitable wind support is available. Our findings suggest that uplift may play a key role in the energy seascape for bird migration that in turn determines strategies and associated costs for birds  crossing ecological barriers such as the open sea.

# This repository consists of 5 R scripts:
1) step_generation.R: privdes the script necessary for generating alternative steps for step-selection analysis
2) step_selection_analysis.R: uses the data generated in the previous script to run the step selection analysis
3) seascape_GAM_analysis.R: reproduces the GAMM analysis for creating regional energy seascapes (predicting delta-t as a function of latitude, longitude and time)
4) all_figures.R: provides the code for reproducing all the figures in the manuscript, including the main text and the supplementary material. All figures are plotted using R Base Graphics.
5) functions.R: includes the functions used in various other scripts. 

# All input data are provided in the Dryad repository accessible via: 
