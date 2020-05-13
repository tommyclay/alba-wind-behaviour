# Sex-specific effects of wind on the flight decisions of a sexually-dimorphic soaring bird

Thomas A. Clay, Rocío Joo, Henri Weimerskirch, Richard A. Phillips, Olivier den Ouden, Mathieu Basille, Susana Clusella-Trullas, Jelle D. Assink, Samantha C. Patrick


## Overview 

The following repository contains codes and data to reproduce the main results of the paper, published in the Journal of Animal Ecology.

There are 2 folders: "Codes" and "Data_inputs". The provided code should create "Data_outputs" and "Plots" folders. 

This repository can be cited through its Zenodo link (Clay et al. 2020; https://zenodo.org/record/3824065#.XrvM-mhKg2x).
 

## Codes

Here are the main codes in R to reproduce the results in the aforementioned paper. 

- **wind_experienced_analysis_plotting.r**: Code to assess how wind speeds and directions encountered vary by population and sex. It contains mixed model fitting and plots, including Figure 2 in the manuscript.

In the hmms sub-folder:

- **1_run_simple_hmm_test_initial_values.r**: Since the parameter estimation in hidden Markov models (HMMs) can be sensitive to the choice of initial values, here we load in data, run initial 3-state HMMs from a range of ecologically realistic values and determine which parameter values were inferred consistently for the 3 main states. 

- **2_hmm_run_all_models.r**: Since we had a set of possible covariates that could have an effect in the transition probabilities between states, here we ran all combinations of HMMs based on covariates specified in Section 2.5 of Methods. Codes for all model outputs, including AIC tables, and autocorrelation and goodness-of-fit plots, are included in this file.

- **3_randomizing_covariates_best_model.r**: This code is related to calculations described in the Appendix S3. Overall, we wanted to double-check that the effect of the covariates of the best model selected in part 2, were indeed significant. For that purpose, we run HMM models with the covariate configuration that was selected as "the best" from the previous code, but reshuffling the values of the covariates, and compared model performance to the "best model".

- **4_plotting_predicted_transition_probabilities.r**: Here we use the model from part 2 to plot and summarize the transition probabilities. These results correspond to Figure 4 and Table 3 in the manuscript.

- **5_plotting_predicted_stationary_probabilities.r**: Here we use the model from part 2 to plot the stationary probabilities corresponding to Figure 5 in the manuscript.

* The combination of a large dataset and a relatively long list of potential covariates mean  these scripts take several days to run. For that reason, most of the codes contain lines to select random subsets of the data and do all the calculations with the subset. We recommend though that to fully reproduce the results, that analyses be conducted on the full dataset.


# Input data
 
These are the datasets used in the analysis, in a straightforward format allowing the fitting of mixed models and HMMs.

- **BothSites_forwindexperienced.csv**: data, with each value corresponding to a tracking location, to load into the wind_experienced_analysis_plotting.r code. The columns are: 
    - ID - factor encoding individual bird identity; 
    - TripID - factor encoding trip identity; 
    - Site - factor with two levels: "Crozet" or "South Georgia"; 
    - Sex - factor with two levels: "F" (female) or "M" (male); 
    - Year - years from 2010-2016; 
    - WindSp - numeric values of wind speed for each tracking location; 
    - WindDir - numeric values of wind direction for each tracking location.

- **GPS_Crozet_2010-2016.csv**: Crozet data to load into hmm codes. Columns are: 
    - ID - factor encoding individual trip identity; 
    - sex - factor with two levels: "F" (female) or "M" (male); 
    - x - longitude; 
    - y - latitude; 
    - ws - numeric value of wind speed for tracking location; 
    - dir - numeric value of wind direction relative to bird trajectory for each location; 
    - lod - factor with two levels: "L" (daylight) or "D" (darkness).


- **GPS_SouthGeorgia_2012.csv**: South Georgia data to load into hmm codes. Columns are:
    - ID - factor encoding individual trip identity; 
    - sex - factor with two levels: "F" (female) or "M" (male); 
    - x - longitude; 
    - y - latitude; 
    - ws - numeric value of wind speed for tracking location; 
    - dir - numeric value of wind direction relative to bird trajectory for each location; 
    - lod - factor with two levels: "L" (daylight) or "D" (darkness).
