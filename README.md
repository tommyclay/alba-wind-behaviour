# Sex-specific effects of wind on the flight decisions of a sexually dimorphic soaring bird

Thomas A. Clay, Rocío Joo, Henri Weimerskirch, Richard A. Phillips, Olivier den Ouden, Mathieu Basille, Susana Clusella-Trullas, Jelle D. Assink, Samantha C. Patrick


## Overview 

The following repository contains codes and data to reproduce the main results of the paper, published in the Journal of Animal Ecology (2020). DOI: 10.1111/1365-2656.13267.

There are 2 folders: "Codes" and "Data_inputs". The provided code should create "Data_outputs" and "Plots" folders. 

This repository can be cited through its Zenodo link (Clay et al. 2020; https://zenodo.org/record/3824065#.XrvM-mhKg2x).
 

## Codes

Here are the main codes in R to reproduce the results in the aforementioned paper. 

- **wind_experienced_analysis_plotting.r**: Code to assess how wind speeds and directions encountered vary by population and sex. It contains mixed model fitting and plots, including Figure 2 in the manuscript.

In the hmms sub-folder:

- **1_run_simple_hmm_test_initial_values.r**: Since the parameter estimation in hidden Markov models (HMMs) can be sensitive to the choice of initial values, here we load in data, run initial 3-state HMMs from a range of ecologically realistic values and determine which parameter values were inferred consistently for the 3 main states. 

- **2_hmm_run_all_models.r**: Since we had a set of possible covariates that could have an effect in the transition probabilities between states, here we ran all combinations of HMMs based on covariates specified in Section 2.5 of Methods. Codes for all model outputs, including AIC tables, and autocorrelation and goodness-of-fit plots, are included in this file.

- **3_randomizing_covariates_best_model.r**: This code is related to calculations described in the Appendix S3. Overall, we wanted to double-check that the effect of the covariates of the best model selected in script 2 were indeed significant. For that purpose, we run HMMs with the covariate configuration that was selected as "the best" from the previous code, but reshuffling the values of the covariates, and compared model performance to the "best model".

- **4_plotting_predicted_transition_probabilities.r**: Here we use the best model from script 2 to plot and summarize the transition probabilities. These results correspond to Figure 4 and Table 3 in the manuscript.

- **5_plotting_predicted_stationary_probabilities.r**: Here we use the best model from script 2 to plot the stationary probabilities corresponding to Figure 5 in the manuscript.

* The combination of a large dataset and a relatively long list of potential covariates mean  these scripts take several days to run. For that reason, most of the codes contain lines to select random subsets of the data and do all the calculations with the subset. We recommend though that to fully reproduce the results, that analyses be conducted on the full dataset.

# General statement (please read before using the data)

The attached archived file(s) contain data derived from long-term field projects involving monitoring and tracking of individual wandering albatrosses at Bird Island (South Georgia) and Possession Island (Crozet Islands).

This is a request to please let us know if you use them. Several people have spent many years involved in the data collection and processing.

If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact Richard Phillips (raphil@bas.ac.uk) or Henri Weimerskirch (henri.weimerskirch@cebc.cnrs.fr) before doing so:

- Occasionally we discover and correct errors in the data.
- The data are complex and workers who do not know the study system are likely to benefit from advice on interpretation.
- At any one time, other people within the project or collaborators may be analysing data from this project. Someone else may already be conducting the analysis that you have in mind and it is desirable to prevent duplication of effort.
- In order to maintain funding for the project and for further analyses, every few years we submit proposals to funding agencies. It is therefore very helpful for those running the project to know which data analyses are in progress elsewhere.

If you are interested in analysing the detailed project data in any depth, you may find it helpful to have access to a larger dataset than the files available here. If so, we are always open to further collaboration.

The individual bird identities have been recoded and should therefore not be linked with data archived from other papers on this study population.


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

- **GPS_Crozet_2010-2016.csv**: Crozet data to load into HMM codes. Columns are: 
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
