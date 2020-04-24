# Sex-specific effects of wind on the flight decisions of a sexually-dimorphic soaring bird

Thomas A. Clay, Rocío Joo, Henri Weimerskirch, Richard A. Phillips, Olivier den Ouden, Mathieu Basille, Susana Clusella-Trullas, Jelle D. Assink, Samantha C. Patrick


## Overview 

The following repository contains codes and data to reproduce main results of the paper published in the Journal of Animal Ecology.

There are 4 folders: Codes, Data_inputs, Data_outputs, and Plots. 


    [RJ]: # (Do you really need to add Data_outputs and Plots?)


For now, this repository can be cited using its github link (https://github.com/tommyclay/alba-wind-behaviour), but will soon have a DOI.
 

## Codes

Here are the main codes in R to reproduce the results in the aforementioned paper. 

- **wind_experienced_analysis_plotting.r**: Code to assess how wind speeds and directions encountered vary by population and sex. It contains mixed model fitting and plots, including Figure 2 in the manuscript.

In the hmms sub-folder:

- **1_run_simple_hmm_test_initial_values.r**: Since the parameter estimation in Hidden Markov Models (HMMs) can be sensitive to the choice of initial values, here we load in data, run initial 3-state HMMs from a range of ecologically realistic values and determine which parameter values were inferred consistently for the 3 main states. 

- **2_hmm_run_all_models.r**: Since we had a set of possible covariates that could have an effect in the transition probabilities from one state to the other, here we ran all combinations of HMMs based on covariates specified in Section 2.5 of Methods. Codes for all model outputs, including AIC tables, and autocorrelation and goodness-of-fit plots, are included in this file.

- **3_randomizing_covariates_best_model.r**: This code is rather associated to calculations described in the appendix. Overall, we wanted to double-check that the effect of the covariates of the best model in part 2 were indeed significant. For that purpose, we run HMM models with the covariate configuration that was selected as "the best" from the previous code, but reshuffling the values of the covariates, and compared there performance in relation to the "best model".

- **4_plotting_predicted_transition_probabilities.r**: Here we use the model from part 2 to plot and summarize the transition probabilities. These results correspond to Figure 4 and Table 3 in the manuscript.

- **5_plotting_predicted_stationary_probabilities.r**: Here we use the model from part 2 to plot the stationary probabilities corresponding to Figure 5 in the manuscript.

* The combination of a large dataset and a relatively long list of potential covariates would make these scripts take several days to run. For that reason, most of the codes contain lines to select random subsets of the data and do all the calculations with the subset.


# Input data
 
These are the datasets used in the analysis, in a format that makes it straightforward to fit the mixed models and HMMs.

- **BothSites_forwindexperienced.csv**: data to load into the wind_experienced_analysis_plotting.r code. The columns are: 
    - ID - factor encoding individual bird identity; 
    - TripID - factor encoding trip identity; 
    - Site - factor with two levels: "Crozet" or "South Georgia"; 
    - Sex - factor with two levels: "F" (female) or "M" (male); 
    - Year - years from 2010-2016; 
    - WindSp - it takes numeric values of wind speed for tracking location; 
    - WindDir - it takes numeric values of wind direction for tracking location.


    [RJ]: # (Tommy, this description is not clear. Are these average wind values?)



- **GPS_Crozet_2010-2016.csv**: Crozet data to load into hmm codes. Columns are: 
    - ID - factor encoding individual trip identity; 
    - sex - factor with two levels: "F" (female) or "M" (male); 
    - x - longitude; 
    - y - latitude; 
    - ws - numeric value of wind speed for tracking location; 
    - dir - numeric value of wind direction relative to bird trajectory for each location; 
    - lod - factor with two levels: "L" (daylight) or "D" (darkness).


- **GPS_SouthGeorgia_2012.csv**: South Georgia to load into hmm codes. Columns are:
    - ID - factor encoding individual trip identity; 
    - sex - factor with two levels: "F" (female) or "M" (male); 
    - x - longitude; 
    - y - latitude; 
    - ws - numeric value of wind speed for tracking location; 
    - dir - numeric value of wind direction relative to bird trajectory for each location; 
    - lod - factor with two levels: "L" (daylight) or "D" (darkness).
