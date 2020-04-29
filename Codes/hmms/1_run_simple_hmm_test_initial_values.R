
library(momentuHMM)
library(ggplot2)


# AIM OF SCRIPT - to load in data, run initial 3-state HMM from a range of sensible values for each state (directed flight, search and rest) and determine if models
# consistently settle on same parameter values for each state.  

# Need to run separately for each population as differences in sampling interval (South Georgia = 25 mins, Crozet = 15 mins) mean that the suitable 
# values (particularly of step lengths) for the main states will be different. 



########  1. SOUTH GEORGIA ##########


#### 1A.  LOAD IN DATA  ####

path.sg <- "./Data_inputs/GPS_SouthGeorgia_2012.csv"
sg <- read.csv(file = path.sg, na.strings = "NA")
sg$ID <- as.factor(as.character(sg$ID))
nlevels(sg$ID) # 38



#### 1B. TO REDUCE COMPUTATION TIME - RANDOMLY SAMPLE 20 INDIVIDUALS AND RUN ON SUBSET #####

# ideally this script should be run on the full sample of individuals

set.seed(1) # fixing the seed allows for comparable results, but you can remove this to make the sampling random
ids <- sample(unique(sg$ID), 20, replace = FALSE)
# create a new subset dataframe sg2
sg2 <- subset(sg, ID %in% ids,)
sg2$ID <- as.factor(as.character(sg2$ID))
nlevels(sg2$ID) # 20

# plot
ggplot(sg2, aes(x, y)) + geom_path(aes(colour = ID))



##### 1C. TURN TO MOMENTUHMM DATAFRAME #######

# Use prepData to create input object - here using subset (sg2) rather than full (sg) dataset 
dat <-prepData(sg2, type= "LL", # longs and lats
               coordNames = c("x", "y")) ## these are our names
head(dat)

# plotting step lengths and turning angles
hist(dat$step) # looks like there's a bimodal relationship - peak at around 15-20 - presumably associated with travel
# then search asociated with lower values (0-10), while large peak for very low steps (<5) likely indicative of rest.
hist(dat$angle)
# lots of lower values of turning angle and few higher values



##### 1D. TRY OUT SIMPLE HMM WITH 3 STATES #######

# The order of these labels correspond to the one in our final results. You'll may have to alter it depending on your results.
stateNames<-c("travel","search", "rest")


# STEP LENGTH INITIAL VALUES #

shape_0 <- c(15, 5, 1) # setting shape (i.e. means for the 3 states; thinking in terms of travel, search and rest, in that order)
scale_0 <- c(5, 4, 0.5)  # setting scale (i.e. variance)
stepPar0 <- c(shape_0,scale_0)
# artifact for zero step lengths for which there are only 2/24442 = <0.01% of values, so just set zero values to very small numbers
# A proper way to deal with zeros is to use the zeromass parameters. 
# As here there are only 2 zero values, this will not change the results
ind_zero <- which(dat$step == 0)
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}
ind_zero <- which(dat$step == "NA")
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}

# TURNING ANGLE INITIAL VALUES #

# we would assume they all have mean of around 0 (i.e. no bias in a particular direction), but concentration differs = travel and rest more concentrated, search large variation in angles.
location_0 <-    c(0.001,0.001, 0.001) # mean values
concentration_0 <- c(60,2,90) # concentration (large values = more concentrated angles)
anglePar0 <- c(location_0,concentration_0)


# RUNNING BASIC 3-STATE MODEL #


# gamma distribution for step lengths and von mises for turning angles
m1 <- fitHMM(data=dat, nbStates=3,
             dist=list(step="gamma",angle="vm"),
             Par0=list(step=stepPar0, angle=anglePar0),
             estAngleMean = list(angle=TRUE), # no need for fix par list beta as all transitions are possible. 
             stateNames = stateNames)
print(m1) # model outputs 
# 
# step parameters:
#   ----------------
#   travel   search      rest
# mean 18.184348 5.295932 0.6253831
# sd    5.873511 6.122766 0.3283630
# 
# angle parameters:
#   -----------------
#   travel      search        rest
# mean           0.005895415 0.007634065  0.03301047
# concentration 43.454467024 1.264779824 52.26784760

# shows what model selects as "best" parameters for different states

plot(m1) # plot distribution of states and states assigned to tracks

# seems to make sense: 
# travel has large step lengths and concentrated angles
# search has moderate step lengths and wide/variable turning angles
# rest has low step lengths and concentrated angles
# when trips are plotted, birds conduct long periods of travel, search over a restricted area interspersed with rest periods.



##### 1E. TRYING A RANGE OF 100 INITIAL VALUES  #######


# CHOOSING A RANGE OF REALISTIC PARAMETERS #

# 100 combinations should be sufficient to know which are "best" starting parameters, plus models take a long time to run.
# For step lengths = choose a range of 10 values for the mean and standard deviation of each state based on model outputs (above) and 
# histograms of of data, and then calculate all combinations of these
hist(dat$step) 
# For angles = give all states similar mean angles as we expect them all to centre around zero. For concentration, we might give the model a bit
# more freedom to choose from a large range for each state - but feel free to change these to more state-specific values. 
hist(dat$angle)

# steps 
travel_step <- expand.grid(Step_m_t = seq(13, 22, 1), Step_sd_t = seq(1, 10, 1))
search_step <- expand.grid(Step_m_s = seq(2, 11, 1), Step_sd_s = seq(2, 11, 1))
rest_step <- expand.grid(Step_m_r = seq(0.2, 2, 0.2), Step_sd_r = seq(0.1, 1, 0.1))

# angles - trying out same for all states
travel_angle <- expand.grid(Angle_loc_t =  seq(-0.2, 0.25, 0.05), Angle_con_t =  seq(10, 100, 10))
search_angle <- expand.grid(Angle_loc_s =  seq(-0.2, 0.25, 0.05), Angle_con_s =  seq(1, 19, 2))
rest_angle <- expand.grid(Angle_loc_r =  seq(-0.2, 0.25, 0.05), Angle_con_r =  seq(10, 100, 10))

# combine into dataframe and randomize order of all step and angle parameters
# it will generate a "decent" set of combinations of randomized but realistic initial values for step and angle parameters
all_vals <- cbind(travel_step[sample(nrow(travel_step)),], 
                  search_step[sample(nrow(search_step)),], 
                  rest_step[sample(nrow(rest_step)),], 
                  travel_angle[sample(nrow(travel_angle)),], 
                  search_angle[sample(nrow(search_angle)),], 
                  rest_angle[sample(nrow(rest_angle)),])



## ITERATE THROUGH AND CHOOSE ONE ROW OF EACH OF THESE AT RANDOM 100 TIMES 

# go through and each time choose a row of randomized values, run hmm and then output parameter values and AIC into dataframe
# due to initial value combinations, sometimes the model cannot converge and so the tryCatch function is a workaround that allows you to skip
# to the next combination. 
# RJ: Is that why we get the warnings? If so, you should add a line about it.

out.mod <- list()
out.df <- list()
for (i in 1:nrow(all_vals)) {
  print(paste("iter", i, sep = "..."))
  # step
  shape_0 <- c(all_vals$Step_m_t[i], all_vals$Step_m_s[i], all_vals$Step_m_r[i])
  scale_0 <- c(all_vals$Step_sd_t[i], all_vals$Step_sd_s[i], all_vals$Step_sd_r[i])   
  stepPar0 <- c(shape_0,scale_0)
  # angle par
  location_0 <- c(all_vals$Angle_loc_t[i], all_vals$Angle_loc_s[i], all_vals$Angle_loc_r[i])
  concentration_0 <- c(all_vals$Angle_con_t[i], all_vals$Angle_con_s[i], all_vals$Angle_con_r[i])
  anglePar0 <- c(location_0,concentration_0)
  tryCatch({ # to skip if model doesn't run
    out.mod[[i]] <- fitHMM(data=dat, nbStates=3,
                           dist=list(step="gamma",angle="vm"),
                           Par0=list(step=stepPar0, angle=anglePar0),
                           estAngleMean = list(angle=TRUE), 
                           stateNames = stateNames)
    # outputting dataframe
    out.df[[i]] <- data.frame(Iter = i, AIC = AIC(out.mod[[i]]),
                              # steps
                              T_step_m_IN = all_vals$Step_m_t[i], T_step_sd_IN = all_vals$Step_sd_t[i], 
                              T_step_m_OUT = out.mod[[i]][[2]]$step[1], T_step_sd_OUT = out.mod[[i]][[2]]$step[2],
                              S_step_m_IN = all_vals$Step_m_s[i], S_step_sd_IN = all_vals$Step_sd_s[i],
                              S_step_m_OUT = out.mod[[i]][[2]]$step[3], S_step_sd_OUT = out.mod[[i]][[2]]$step[4],
                              R_step_m_IN = all_vals$Step_m_r[i], R_step_sd_IN = all_vals$Step_m_r[i],
                              R_step_m_OUT = out.mod[[i]][[2]]$step[5], R_step_sd_OUT = out.mod[[i]][[2]]$step[6],
                              # angles
                              T_ang_m_IN = all_vals$Angle_loc_t[i], T_ang_sd_IN = all_vals$Angle_con_t[i],
                              T_ang_m_OUT = out.mod[[i]][[2]]$angle[1], T_ang_sd_OUT = out.mod[[i]][[2]]$angle[2],
                              S_ang_m_IN = all_vals$Angle_loc_s[i], S_ang_sd_IN = all_vals$Angle_con_s[i],
                              S_ang_m_OUT = out.mod[[i]][[2]]$angle[3], S_ang_sd_OUT = out.mod[[i]][[2]]$angle[4],
                              R_ang_m_IN = all_vals$Angle_loc_r[i], R_ang_sd_IN = all_vals$Angle_con_r[i],
                              R_ang_m_OUT = out.mod[[i]][[2]]$angle[5], R_ang_sd_OUT = out.mod[[i]][[2]]$angle[6])
  }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n")})
}
out.all <- do.call(rbind, out.df)
out.all <- out.all[order(out.all$AIC),]

# paste out

out.path <- "./Data_outputs/"
if (dir.exists(out.path) == FALSE){
  dir.create(out.path)
}
filename <- paste0(out.path,"Testing_initial_vals_SG.csv")
write.csv(out.all, filename, row.names=T)


head(out.all)
# there are several cases within 2 units of the top AIC score
# choosing just the "top" models that are within 2 AIC units of the highest value - as a fairly arbritary cut-off
top <- subset(out.all, AIC < out.all$AIC[1]+2)
nrow(top)

# plotting histogram of step lengths based on "best" models
hist(top$T_step_m_OUT, breaks = 10)
hist(top$S_step_m_OUT, breaks = 10)
hist(top$R_step_m_OUT, breaks = 10)

# plotting step lengths together
top.df <- data.frame(Vals = c(top$T_step_m_OUT, top$S_step_m_OUT, top$R_step_m_OUT),
                     State = c(rep("T", length(top$T_step_m_OUT)), rep("S", length(top$S_step_m_OUT)), rep("R", length(top$R_step_m_OUT))))

ggplot(top.df, aes(Vals)) + geom_histogram(aes(fill = State)) + xlab('Step length')
# you can see that the "best" models consistently choose the same inital means for step length, but occassionally they are labelled in the wrong order
# this is not problem - the purpose of this code is mainly to check that the model selects 3 peaks of parameter values. We know that of the 3 states
# the model selects, from our knowledge of albatross ecology we know that travel will be the state with the greatest step lengths and rest the smallest. 

# the values below are just from the subset (of 20 individuals)
# RJ: "average values?"
# travel = 18.73
# search = 5.59
# rest = 0.36

# plotting trips from "best" model to see if they make sense visually 

head(out.all) # model number 82 was "best" for me, but will change each time as order will be randomized
plot(out.mod[[82]])
# look good although states are labelled in the incorrect order. 







########  2. CROZET ##########

# This exactly the same analysis ran with the South Georgia dataset

#### 2A.  LOAD IN DATA  ####

path.cro <- "./Data_inputs/GPS_Crozet_2010-2016.csv"
cro <- read.csv(file = path.cro, na.strings = "NA")
cro$ID <- as.factor(as.character(cro$ID))
nlevels(cro$ID) # 267


#### 2B. TO REDUCE COMPUTATION TIME - RANDOMLY SAMPLE 20 INDIVIDUALS AND RUN ON SUBSET #####

# ideally this script should be run on the full sample of individuals

set.seed(2)

ids <- sample(unique(cro$ID), 20, replace = FALSE)
# create a new subset dataframe sg2
cro2 <- subset(cro, ID %in% ids)
cro2$ID <- as.factor(as.character(cro2$ID))
nlevels(cro2$ID) # 20

# plot
ggplot(cro2, aes(x, y)) + geom_path(aes(colour = ID))



##### 2C. TURN TO MOMENTUHMM DATAFRAME #######

# Use prepData to create input object - here using subset (cro2) rather than full (cro) dataset 
dat <-prepData(cro2, type= "LL", # longs and lats
               coordNames = c("x", "y")) ## these are our names
head(dat)

# plotting step lengths and turning angles
hist(dat$step) # looks like there's a bimodal relationship - peak at around 10-15 - presumably associated with travel
# then search asociated with lower values (0-10), while large peak for very low steps (<3) likely indicative of rest.
hist(dat$angle)
# lots of lower values of turning angle and few higher values




##### 2D. TRY OUT SIMPLE HMM WITH 3 STATES #######


stateNames<-c("travel","search", "rest")


# STEP LENGTH INITIAL VALUES #

shape_0 <- c(12, 3, 1) # setting shape (i.e. means for the 3 states)
scale_0 <- c(4, 3, 0.5)  # setting scale (i.e. variance)
stepPar0 <- c(shape_0,scale_0)
# artifact for zero step lengths for which there are only 124/232713 = 0.05% values, so just set zero values to small numbers
# ideally use the zero mass parameter
ind_zero <- which(dat$step == 0)
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}
ind_zero <- which(dat$step == "NA")
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}

# TURNING ANGLE INITIAL VALUES #

# we would assume they all have mean of around 0 (i.e. no bias in a particular direction), but concentration differs = travel and rest more concentrated, search large variation in angles.
location_0 <-    c(0.001,0.001, 0.001) # mean values
concentration_0 <- c(60,2,90) # concentration (large values = more concentrated angles)
anglePar0 <- c(location_0,concentration_0)



# RUNNING BASIC 3-STATE MODEL #


# gamma distribution for step lengths and von mises for turning angles
m1 <- fitHMM(data=dat, nbStates=3,
             dist=list(step="gamma",angle="vm"),
             Par0=list(step=stepPar0, angle=anglePar0),
             estAngleMean = list(angle=TRUE), # no need for fix par list beta as all transitions are possible. 
             stateNames = stateNames)
print(m1) # model outputs 

# step parameters:
#   ----------------
#   travel   search      rest
# mean 12.808237 4.384151 0.3195857
# sd    3.497444 4.996611 0.1690926
# 
# angle parameters:
#   -----------------
#   travel     search        rest
# mean           0.002588744 0.01749804  0.02294263
# concentration 47.711634938 1.30130186 43.64997801

# shows what model selects as "best" parameters for different states

plot(m1) # plot distribution of states and states assigned to tracks

# seems to make sense: 
# travel has large step lengths and concentrated angles
# search has moderate step lengths and wide/variable turning angles
# rest has low step lengths and concentrated angles
# when trips are plotted, birds conduct long periods of travel, search over a restricted area interspersed with rest periods.



##### 2E. TRYING A RANGE OF 100 INITIAL VALUES  #######


# CHOOSING A RANGE OF REALISTIC PARAMETERS #

# 100 combinations should be sufficient to know which are "best" starting parameters, plus models take a long time to run.
# For step lengths = choose a range of 10 values for the mean and standard deviation of each state based on model outputs (above) and 
# histograms of of data, and then calculate all combinations of these
hist(dat$step) 
# For angles = give all states similar mean angles as we expect them all to centre around zero. For concentration, we might give the model a bit
# more freedom to choose from a large range for each state - but feel free to change these to more state-specific values. 
hist(dat$angle)

# steps 
travel_step <- expand.grid(Step_m_t = seq(8, 17, 1), Step_sd_t = seq(1, 10, 1))
search_step <- expand.grid(Step_m_s = seq(2, 6.5, 0.5), Step_sd_s = seq(1, 10, 1))
rest_step <- expand.grid(Step_m_r = seq(0.2, 2, 0.2), Step_sd_r = seq(0.1, 1, 0.1))

# angles - trying out same for all states
travel_angle <- expand.grid(Angle_loc_t =  seq(-0.2, 0.25, 0.05), Angle_con_t =  seq(10, 100, 10))
search_angle <- expand.grid(Angle_loc_s =  seq(-0.2, 0.25, 0.05), Angle_con_s =  seq(1, 19, 2))
rest_angle <- expand.grid(Angle_loc_r =  seq(-0.2, 0.25, 0.05), Angle_con_r =  seq(10, 100, 10))

# combine into dataframe and randomize order of all step and angle parameters
all_vals <- cbind(travel_step[sample(nrow(travel_step)),], 
                  search_step[sample(nrow(search_step)),], 
                  rest_step[sample(nrow(rest_step)),], 
                  travel_angle[sample(nrow(travel_angle)),], 
                  search_angle[sample(nrow(search_angle)),], 
                  rest_angle[sample(nrow(rest_angle)),])


## ITERATE THROUGH AND CHOOSE ONE ROW OF EACH OF THESE AT RANDOM 100 TIMES 


# go through and each time choose a row of randomized values, run hmm and then output parameter values and AIC into dataframe
# due to initial value combinations, sometimes the model cannot converge and so the tryCatch function is a workaround that allows you to skip
# to the next combination. 

out.mod <- list()
out.df <- list()
for (i in 1:nrow(all_vals)) {
  print(paste("iter", i, sep = "..."))
  # step
  shape_0 <- c(all_vals$Step_m_t[i], all_vals$Step_m_s[i], all_vals$Step_m_r[i])
  scale_0 <- c(all_vals$Step_sd_t[i], all_vals$Step_sd_s[i], all_vals$Step_sd_r[i])   
  stepPar0 <- c(shape_0,scale_0)
  # angle par
  location_0 <- c(all_vals$Angle_loc_t[i], all_vals$Angle_loc_s[i], all_vals$Angle_loc_r[i])
  concentration_0 <- c(all_vals$Angle_con_t[i], all_vals$Angle_con_s[i], all_vals$Angle_con_r[i])
  anglePar0 <- c(location_0,concentration_0)
  tryCatch({ # to skip if model doesn't run
    out.mod[[i]] <- fitHMM(data=dat, nbStates=3,
                           dist=list(step="gamma",angle="vm"),
                           Par0=list(step=stepPar0, angle=anglePar0),
                           estAngleMean = list(angle=TRUE), 
                           stateNames = stateNames)
    # outputting dataframe
    out.df[[i]] <- data.frame(Iter = i, AIC = AIC(out.mod[[i]]),
                              # steps
                              T_step_m_IN = all_vals$Step_m_t[i], T_step_sd_IN = all_vals$Step_sd_t[i], 
                              T_step_m_OUT = out.mod[[i]][[2]]$step[1], T_step_sd_OUT = out.mod[[i]][[2]]$step[2],
                              S_step_m_IN = all_vals$Step_m_s[i], S_step_sd_IN = all_vals$Step_sd_s[i],
                              S_step_m_OUT = out.mod[[i]][[2]]$step[3], S_step_sd_OUT = out.mod[[i]][[2]]$step[4],
                              R_step_m_IN = all_vals$Step_m_r[i], R_step_sd_IN = all_vals$Step_m_r[i],
                              R_step_m_OUT = out.mod[[i]][[2]]$step[5], R_step_sd_OUT = out.mod[[i]][[2]]$step[6],
                              # angles
                              T_ang_m_IN = all_vals$Angle_loc_t[i], T_ang_sd_IN = all_vals$Angle_con_t[i],
                              T_ang_m_OUT = out.mod[[i]][[2]]$angle[1], T_ang_sd_OUT = out.mod[[i]][[2]]$angle[2],
                              S_ang_m_IN = all_vals$Angle_loc_s[i], S_ang_sd_IN = all_vals$Angle_con_s[i],
                              S_ang_m_OUT = out.mod[[i]][[2]]$angle[3], S_ang_sd_OUT = out.mod[[i]][[2]]$angle[4],
                              R_ang_m_IN = all_vals$Angle_loc_r[i], R_ang_sd_IN = all_vals$Angle_con_r[i],
                              R_ang_m_OUT = out.mod[[i]][[2]]$angle[5], R_ang_sd_OUT = out.mod[[i]][[2]]$angle[6])
  }, error = function(e) { cat("ERROR :", conditionMessage(e), "\n")})
}
out.all <- do.call(rbind, out.df)
out.all <- out.all[order(out.all$AIC),]

# paste out
out.path <- "./Data_outputs/"
if (dir.exists(out.path) == FALSE){
  dir.create(out.path)
}
filename <- paste0(out.path,"Testing_initial_vals_Cro.csv")
write.csv(out.all, filename, row.names=T)


head(out.all)
# there are several cases within 2 units of the top AIC score
# choosing just the "top" models that are within 2 AIC units of the highest value - as a fairly arbritary cut-off
top <- subset(out.all, AIC < out.all$AIC[1]+2)
nrow(top)

# plotting histogram of best values 
hist(top$T_step_m_OUT, breaks = 10)
hist(top$S_step_m_OUT, breaks = 10)
hist(top$R_step_m_OUT, breaks = 10)

# plotting together
top.df <- data.frame(Vals = c(top$T_step_m_OUT, top$S_step_m_OUT, top$R_step_m_OUT),
                     State = c(rep("T", length(top$T_step_m_OUT)), rep("S", length(top$S_step_m_OUT)), rep("R", length(top$R_step_m_OUT))))

ggplot(top.df, aes(Vals)) + geom_histogram(aes(fill = State)) + xlab('Step length')
# you can see that the "best" models consistently choose the same inital means for step length, but occassionally they are labelled in the wrong order
# this is not problem - the purpose of this code is mainly to check that the model selects 3 peaks of parameter values. We know that of the 3 states
# the model selects, from our knowledge of albatross ecology we know that travel will be the state with the greatest step lengths and rest the smallest. 

# the values below are just for the subset (of 20 individuals)
# travel = 12.46
# search = 3.95
# rest = 0.34

# plotting trips from "best" model to see if they make sense visually 

head(out.all) # model number 44 was "best" for me, but will change each time as order will be randomized
plot(out.mod[[44]])
# look good although states are labelled in the incorrect order. 


