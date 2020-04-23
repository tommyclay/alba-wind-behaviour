# Sex-specific effects of wind on the flight decisions of a sexually-dimorphic soaring bird
Thomas A. Clay, Rocío Joo, Henri Weimerskirch, Richard A. Phillips, Olivier den Ouden, Mathieu Basille, Susana Clusella-Trullas, Jelle D. Assink, Samantha C. Patrick

The following repository contains codes and data to reproduce main results of the paper published in the Journal of Animal Ecology - insert url here!

Codes:

- wind_experienced_analysis_plotting.r - code to calculate how wind speeds and directions encounterd vary by population and sex, incluidng code to plot Figure 2

In hmms folder:

- 1_run_simple_hmm_test_initial_values.r - load in data, run initial 3-state hmm from range of realistic values and determine if models consistently select similar parameter values for the 3 main states 

- 2_hmm_run_all_models - to run all combinations of hmms based on predictors specified in Section 2.5 of Methods, and outputting models, AIC tables and autocorrelation and goodness of fit plots

- 3_randomizing_covariates_best_model - to run the "best" model but with covariate values reshuffled to test for over-parameterization

- 4_plotting_predicted_transition_probabilities - plot predicted transition probabilities for Figure 4 and Table 3

- 5_plotting_predicted_stationary_probabilities - plot predicted stationary probabilities for Figure 5


Input data:

BothSites_forwindexperienced.csv - data to load into wind_experienced_analysis_plotting.r code. Columns are: "ID" - factor encoding individual bird identity; "TripID" - factor encoding trip identity; "Site" - factor with two levels: "Crozet" or "South Georgia"; "Sex" - factor with two levels: "F" (female) or "M" (male); Year - years from 2010-2016; WindSp - numeric value of wind speed for tracking location; WindDir -numeric value of wind direction for tracking location.

ID	sex	x	y	ws	dir	lod

GPS_Crozet_2010-2016.csv - Crozet data to load into hmm codes. Columns are: "ID" - factor encoding individual trip identity; "sex" - factor with two levels: "F" (female) or "M" (male); x - longitude; y - latitude; ws - numeric value of wind speed for tracking location; dir - numeric value of wind direction relative to bird trajectory for each location; lod - factor with two levels: "L" (daylight) or "D" (darkness).

GPS_SouthGeorgia_2012.csv - South Georgia to load into hmm codes. Columns are: "ID" - factor encoding individual trip identity; "sex" - factor with two levels: "F" (female) or "M" (male); x - longitude; y - latitude; ws - numeric value of wind speed for tracking location; dir - numeric value of wind direction relative to bird trajectory for each location; lod - factor with two levels: "L" (daylight) or "D" (darkness).
