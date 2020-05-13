library(ggplot2)
library(momentuHMM)


# AIM OF SCRIPT - to plot predicted stationary probabilites for a range of specified covariates and produce Figure 5 



##### LOAD IN MODEL OUTPUTS #####


# load in best models for each site - for the manuscript it was model 30
file.in <- paste0("./Data_outputs/", paste0("SG_mod_", 30, ".RData"))
load(file = file.in)
m.sg <- model

file.in <- paste0("./Data_outputs/", paste0("Cro_mod_", 30, ".RData"))
load(file = file.in)
m.cro <- model



# initial plotting to look at stationary probabilities


plotStationary(m.sg,plotCI = TRUE)
plotStationary(m.cro,plotCI = TRUE)


# plotStationary() function only allows you to predict for a row of covariate values within the plotting architecture of momentuhmm
# to extract confidence intervals of the stationary transition probabilities, R. Joo wrote a function based on momentuhmm internal functions (in Apr. 2019)
# for a range of specified covariates and output to plot in ggplot

ci.stationary <- function(model, cov, alpha) {
  gridLength <- nrow(cov)
  model2 <- momentuHMM:::delta_bc(model) # extended format allowing for backwards compatibility
  nbStates <- length(model2$stateNames)
  Sigma <- model2$mod$Sigma
  
  lci <- matrix(NA,gridLength,nbStates)
  uci <- matrix(NA,gridLength,nbStates)
  
  formula <- model2$conditions$formula # for covariates
  newForm <- momentuHMM:::newFormulas(formula, nbStates) # a formula for each state transition
  newformula <- newForm$newformula
  
  nbCovs <- ncol(model.matrix(newformula, model2$data)) - 1 # model.matrix gives values of covariate terms of the formula for each observation
  # in my opinion, nbCovs gives exactly the same number of terms as formula (-1 is minus intercept)
  
  gamInd <- (length(model2$mod$estimate) - (nbCovs + 1) * nbStates * (nbStates - 1) 
             + 1):(length(model2$mod$estimate)) - ncol(model2$covsDelta) * (nbStates - 1) * (!model2$conditions$stationary) # if the mode is stationary there wouldn't be any covariates
  # here we just enumerate parameters leaving out the first ones which are step and angle related, and the last ones which are delta related
  
  rawCovs <- model2$rawCovs
  # changing order of variables tomatch with raw data
  rawCovs <- rawCovs[,order(names(rawCovs))]
  tempCovs <- cov
  tempCovs <- tempCovs[,sort(names(rawCovs))]
  
  # fixing format problems
  for (i in which(unlist(lapply(rawCovs, is.factor)))) {
    tempCovs[[i]] <- factor(tempCovs[[i]], levels = levels(rawCovs[, 
                                                                   i])) 
  }
  tmpSplineInputs <- momentuHMM:::getSplineFormula(newformula, model2$data, 
                                                   tempCovs) # just a format for spline use later
  desMat <- model.matrix(tmpSplineInputs$formula, data = tmpSplineInputs$covs) # tmpSplineInputs$covs is tempCovs
  # model.matrix gives the design matrix for the formula in the argument
  
  probs <- as.data.frame(stationary(model2, covs=desMat)) # the stationary probability is computed based on desMat, which has only a range of values for one of the covariates
  
  for(state in 1:nbStates) {
    dN <- t(apply(desMat, 1, function(x)
      numDeriv::grad(momentuHMM:::get_stat,model2$mod$estimate[gamInd][unique(c(model2$conditions$betaCons))],covs=matrix(x,nrow=1),nbStates=nbStates,i=state,betaRef=model2$conditions$betaRef,betaCons=model2$conditions$betaCons,workBounds=model2$conditions$workBounds$beta)))
    
    se <- t(apply(dN, 1, function(x)
      suppressWarnings(sqrt(x%*%Sigma[gamInd[unique(c(model2$conditions$betaCons))],gamInd[unique(c(model2$conditions$betaCons))]]%*%x))))
    
    lci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) -
                                  qnorm(1-(1-alpha)/2) * (1/(probs[,state]-probs[,state]^2)) * se)))
    uci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) +
                                  qnorm(1-(1-alpha)/2) * (1/(probs[,state]-probs[,state]^2)) * se)))
  }
  lci <- as.data.frame(lci)
  names(lci) <- model$stateNames
  uci <- as.data.frame(uci)
  names(uci) <- model$stateNames
  ci.list <- list(lci, uci)
  names(ci.list) <- c("lower", "upper")
  return(ci.list)
}




##### 1. PLOTTING STATIONARY PROBABILITIES FOR COVARIATE VALUES OF WIND SPEED - FOR FIGURE 5 ###############



min(m.sg$data$ws)
#0.08082135 # these values correspond to the full dataset
min(m.cro$data$ws)
# 0.04510045
max(m.sg$data$ws)
# 23.6461
max(m.cro$data$ws)
# 23.1515

# predict wind speed for average wind direction which is 75 degrees relative to bird direction, and for daylight hours when birds are most active. 

cov.sg.m <- data.frame(lod = rep("L"), dir = rep(75), ws = seq(from=0.5, 23,by=0.2), sex = "M")
cov.sg.f <- data.frame(lod = rep("L"), dir = rep(75), ws = seq(from=0.5, 23,by=0.2), sex = "F")
cov.cro.m <- data.frame(lod = rep("L"), dir = rep(75), ws = seq(from=0.5, 23,by=0.2), sex = "M")
cov.cro.f <- data.frame(lod = rep("L"), dir = rep(75), ws = seq(from=0.5, 23,by=0.2), sex = "F")

s.cro.m <- as.data.frame(stationary(model = m.cro, covs = cov.cro.m))
s.cro.f <- as.data.frame(stationary(model = m.cro, covs = cov.cro.f))
s.sg.m <- as.data.frame(stationary(model = m.sg, covs = cov.sg.m))
s.sg.f <- as.data.frame(stationary(model = m.sg, covs = cov.sg.f))

# creating confidence intervals around stationary probabilities

stat.cro.m <- ci.stationary(model = m.cro, cov = cov.cro.m, alpha = 0.95) 
stat.cro.f <- ci.stationary(model = m.cro, cov = cov.cro.f, alpha = 0.95) 
stat.sg.m <- ci.stationary(model = m.sg, cov = cov.sg.m, alpha = 0.95) 
stat.sg.f <- ci.stationary(model = m.sg, cov = cov.sg.f, alpha = 0.95) 


# combine into dataframe for plotting

stat.df <- data.frame(ws = rep(cov.sg.m$ws, 12),
                      prob = c(s.cro.m$travel, s.cro.m$search, s.cro.m$rest,
                               s.cro.f$travel, s.cro.f$search, s.cro.f$rest,
                               s.sg.m$travel, s.sg.m$search, s.sg.m$rest,
                               s.sg.f$travel, s.sg.f$search, s.sg.f$rest),
                      lwr = c(stat.cro.m[[1]]$travel, stat.cro.m[[1]]$search, stat.cro.m[[1]]$rest,
                              stat.cro.f[[1]]$travel, stat.cro.f[[1]]$search, stat.cro.f[[1]]$rest,
                              stat.sg.m[[1]]$travel, stat.sg.m[[1]]$search, stat.sg.m[[1]]$rest,
                              stat.sg.f[[1]]$travel, stat.sg.f[[1]]$search, stat.sg.f[[1]]$rest),
                      upr = c(stat.cro.m[[2]]$travel, stat.cro.m[[2]]$search, stat.cro.m[[2]]$rest,
                              stat.cro.f[[2]]$travel, stat.cro.f[[2]]$search, stat.cro.f[[2]]$rest,
                              stat.sg.m[[2]]$travel, stat.sg.m[[2]]$search, stat.sg.m[[2]]$rest,
                              stat.sg.f[[2]]$travel, stat.sg.f[[2]]$search, stat.sg.f[[2]]$rest),
                      state = rep(rep(c("travel", "search", "rest"), each = length(cov.sg.m$ws)), 4),
                      sex = rep(c("Male", "Female", "Male", "Female"), each = length(cov.sg.m$ws)*3),
                      site = rep(c("Crozet", "South Georgia"), each = length(cov.sg.m$ws)*6))

stat.df$sex_state <- as.factor(as.character(paste(stat.df$sex, stat.df$state, sep = "_")))


# create plot

ggplot(stat.df, aes(ws, prob)) +facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lwr, ymax=upr, fill = sex_state), alpha=0.15)+ 
  geom_line(aes(colour = state, linetype = sex), size = 1.3) + 
  ylim(0, 1)+theme_bw() + #ylab("Stationary probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab(expression(paste("Wind speed (", ms^-1, ")", sep="")))+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))+
  # check colours below refer to right states!
  scale_colour_manual(values = c("yellow", "red", "blue"))





##### 2. PLOTTING STATIONARY PROBABILITIES FOR COVARIATE VALUES OF WIND DIRECTION - FOR FIGURE 5 ###############


min(m.sg$data$dir)
# 0.02182379 # reference valeues for the whole dataset
min(m.cro$data$dir)
# 0.001087263
max(m.sg$data$dir)
# 179.9744
max(m.cro$data$dir)
# 179.9997

# predict wind direction for average wind speed which is 9 ms-1, and for daylight hours when birds are most active. 

cov.sg.m <- data.frame(lod = rep("L"), dir = seq(from=0, 180, by= 1), ws = rep(9), sex = "M")
cov.sg.f <- data.frame(lod = rep("L"), dir = seq(from=0, 180, by= 1), ws = rep(9), sex = "F")
cov.cro.m <- data.frame(lod = rep("L"), dir = seq(from=0, 180, by= 1), ws = rep(9), sex = "M")
cov.cro.f <- data.frame(lod = rep("L"), dir = seq(from=0, 180, by= 1), ws = rep(9), sex = "F")

s.cro.m <- as.data.frame(stationary(model = m.cro, covs = cov.cro.m))
s.cro.f <- as.data.frame(stationary(model = m.cro, covs = cov.cro.f))
s.sg.m <- as.data.frame(stationary(model = m.sg, covs = cov.sg.m))
s.sg.f <- as.data.frame(stationary(model = m.sg, covs = cov.sg.f))

# creating confidence intervals around stationary probabilities

stat.cro.m <- ci.stationary(model = m.cro, cov = cov.cro.m, alpha = 0.95) 
stat.cro.f <- ci.stationary(model = m.cro, cov = cov.cro.f, alpha = 0.95) 
stat.sg.m <- ci.stationary(model = m.sg, cov = cov.sg.m, alpha = 0.95) 
stat.sg.f <- ci.stationary(model = m.sg, cov = cov.sg.f, alpha = 0.95) 

# combine into dataframe for plotting

stat.df <- data.frame(dir = rep(cov.sg.m$dir, 12),
                      prob = c(s.cro.m$travel, s.cro.m$search, s.cro.m$rest,
                               s.cro.f$travel, s.cro.f$search, s.cro.f$rest,
                               s.sg.m$travel, s.sg.m$search, s.sg.m$rest,
                               s.sg.f$travel, s.sg.f$search, s.sg.f$rest),
                      lwr = c(stat.cro.m[[1]]$travel, stat.cro.m[[1]]$search, stat.cro.m[[1]]$rest,
                              stat.cro.f[[1]]$travel, stat.cro.f[[1]]$search, stat.cro.f[[1]]$rest,
                              stat.sg.m[[1]]$travel, stat.sg.m[[1]]$search, stat.sg.m[[1]]$rest,
                              stat.sg.f[[1]]$travel, stat.sg.f[[1]]$search, stat.sg.f[[1]]$rest),
                      upr = c(stat.cro.m[[2]]$travel, stat.cro.m[[2]]$search, stat.cro.m[[2]]$rest,
                              stat.cro.f[[2]]$travel, stat.cro.f[[2]]$search, stat.cro.f[[2]]$rest,
                              stat.sg.m[[2]]$travel, stat.sg.m[[2]]$search, stat.sg.m[[2]]$rest,
                              stat.sg.f[[2]]$travel, stat.sg.f[[2]]$search, stat.sg.f[[2]]$rest),
                      state = rep(rep(c("travel", "search", "rest"), each = length(cov.sg.m$ws)), 4),
                      sex = rep(c("Male", "Female", "Male", "Female"), each = length(cov.sg.m$ws)*3),
                      site = rep(c("Crozet", "South Georgia"), each = length(cov.sg.m$ws)*6))


stat.df$sex_state <- as.factor(as.character(paste(stat.df$sex, stat.df$state, sep = "_")))

# create plot

ggplot(stat.df, aes(dir, prob)) +facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lwr, ymax=upr, fill = sex_state), alpha=0.15)+ 
  geom_line(aes(colour = state, linetype = sex), size = 1.3) + 
  ylim(0, 1)+theme_bw() + #ylab("Stationary probability")+ xlab("Relative wind direction (degrees)")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0,180,60)) +  
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))+
  scale_colour_manual(values = c("yellow", "red", "blue"))




##### 3. PLOTTING STATIONARY PROBABILITIES BY SEX - FOR FIGURE 5 ###############


#  for average wind speed and direction

cov.sg.m <- data.frame(lod = "L", dir = 75, ws = 9, sex = "M")
cov.sg.f <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")
cov.cro.m <- data.frame(lod = "L", dir = 75, ws = 9, sex = "M")
cov.cro.f <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")

s.cro.m <- as.data.frame(stationary(model = m.cro, covs = cov.cro.m))
s.cro.f <- as.data.frame(stationary(model = m.cro, covs = cov.cro.f))
s.sg.m <- as.data.frame(stationary(model = m.sg, covs = cov.sg.m))
s.sg.f <- as.data.frame(stationary(model = m.sg, covs = cov.sg.f))

# creating confidence intervals around stationary probabilities

stat.cro.m <- ci.stationary(model = m.cro, cov = cov.cro.m, alpha = 0.95) 
stat.cro.f <- ci.stationary(model = m.cro, cov = cov.cro.f, alpha = 0.95) 
stat.sg.m <- ci.stationary(model = m.sg, cov = cov.sg.m, alpha = 0.95) 
stat.sg.f <- ci.stationary(model = m.sg, cov = cov.sg.f, alpha = 0.95) 


stat.df <- data.frame(dir = rep(cov.sg.m$dir, 4),
                      prob = c(s.cro.m$travel, s.cro.m$search, s.cro.m$rest,
                               s.cro.f$travel, s.cro.f$search, s.cro.f$rest,
                               s.sg.m$travel, s.sg.m$search, s.sg.m$rest,
                               s.sg.f$travel, s.sg.f$search, s.sg.f$rest),
                      lwr = c(stat.cro.m[[1]]$travel, stat.cro.m[[1]]$search, stat.cro.m[[1]]$rest,
                              stat.cro.f[[1]]$travel, stat.cro.f[[1]]$search, stat.cro.f[[1]]$rest,
                              stat.sg.m[[1]]$travel, stat.sg.m[[1]]$search, stat.sg.m[[1]]$rest,
                              stat.sg.f[[1]]$travel, stat.sg.f[[1]]$search, stat.sg.f[[1]]$rest),
                      upr = c(stat.cro.m[[2]]$travel, stat.cro.m[[2]]$search, stat.cro.m[[2]]$rest,
                              stat.cro.f[[2]]$travel, stat.cro.f[[2]]$search, stat.cro.f[[2]]$rest,
                              stat.sg.m[[2]]$travel, stat.sg.m[[2]]$search, stat.sg.m[[2]]$rest,
                              stat.sg.f[[2]]$travel, stat.sg.f[[2]]$search, stat.sg.f[[2]]$rest),
                      state = rep(rep(c("travel", "search", "rest"), each = length(cov.sg.m$ws)), 4),
                      sex = rep(c("M", "F", "M", "F"), each = length(cov.sg.m$ws)*3),
                      site = rep(c("Crozet", "South Georgia"), each = length(cov.sg.m$ws)*6))

pd <- position_dodge()
ggplot(stat.df, aes(sex, prob)) + facet_wrap(~site)+
  geom_errorbar(aes(ymin=lwr, ymax=upr, group = sex), width=.25, position=pd) +
  geom_point(aes(fill = state), pch = 21, size = 2.7)+
  ylim(0, 1)+theme_bw() + #ylab("Stationary probability")+
  #xlab("Sex")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 16),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("yellow", "red", "blue"))




##### 4. PLOTTING STATIONARY PROBABILITIES BY DAYLIGHT VS DARKNESS - FOR FIGURE 5 ###############


# as there was no sex interaction with lod, need to calculate separately for both sexes and average the two
# for average wind speed and direction

# males

cov.sg.l.m <- data.frame(lod = "L", dir = 75, ws = 9, sex = "M")
cov.sg.d.m <- data.frame(lod = "D", dir = 75, ws = 9, sex = "M")
cov.cro.l.m <- data.frame(lod = "L", dir = 75, ws = 9, sex = "M")
cov.cro.d.m <- data.frame(lod = "D", dir = 75, ws = 9, sex = "M")

s.cro.l.m <- as.data.frame(stationary(model = m.cro, covs = cov.cro.l.m))
s.cro.d.m <- as.data.frame(stationary(model = m.cro, covs = cov.cro.d.m))
s.sg.l.m <- as.data.frame(stationary(model = m.sg, covs = cov.sg.l.m))
s.sg.d.m <- as.data.frame(stationary(model = m.sg, covs = cov.sg.d.m))

stat.cro.l.m <- ci.stationary(model = m.cro, cov = cov.cro.l.m, alpha = 0.95) 
stat.cro.d.m <- ci.stationary(model = m.cro, cov = cov.cro.d.m, alpha = 0.95) 
stat.sg.l.m <- ci.stationary(model = m.sg, cov = cov.sg.l.m, alpha = 0.95) 
stat.sg.d.m <- ci.stationary(model = m.sg, cov = cov.sg.d.m, alpha = 0.95) 

# females

cov.sg.l.f <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")
cov.sg.d.f <- data.frame(lod = "D", dir = 75, ws = 9, sex = "F")
cov.cro.l.f <- data.frame(lod = "L", dir = 75, ws = 9, sex = "F")
cov.cro.d.f <- data.frame(lod = "D", dir = 75, ws = 9, sex = "F")

s.cro.l.f <- as.data.frame(stationary(model = m.cro, covs = cov.cro.l.f))
s.cro.d.f <- as.data.frame(stationary(model = m.cro, covs = cov.cro.d.f))
s.sg.l.f <- as.data.frame(stationary(model = m.sg, covs = cov.sg.l.f))
s.sg.d.f <- as.data.frame(stationary(model = m.sg, covs = cov.sg.d.f))

stat.cro.l.f <- ci.stationary(model = m.cro, cov = cov.cro.l.f, alpha = 0.95) 
stat.cro.d.f <- ci.stationary(model = m.cro, cov = cov.cro.d.f, alpha = 0.95) 
stat.sg.l.f <- ci.stationary(model = m.sg, cov = cov.sg.l.f, alpha = 0.95) 
stat.sg.d.f <- ci.stationary(model = m.sg, cov = cov.sg.d.f, alpha = 0.95) 

# combine both
s.cro.l <- as.data.frame.list(colMeans(rbind(s.cro.l.m, s.cro.l.f)))
s.cro.d <- as.data.frame.list(colMeans(rbind(s.cro.d.m, s.cro.d.f)))
s.sg.l <- as.data.frame.list(colMeans(rbind(s.sg.l.m, s.sg.l.f)))
s.sg.d <- as.data.frame.list(colMeans(rbind(s.sg.d.m, s.sg.d.f)))

stat.cro.l <- stat.cro.l.f 
stat.cro.l$lower <- as.data.frame.list(colMeans(rbind(stat.cro.l.m$lower, stat.cro.l.f$lower)))
stat.cro.l$upper <- as.data.frame.list(colMeans(rbind(stat.cro.l.m$upper, stat.cro.l.f$upper)))
stat.sg.l <- stat.sg.l.f 
stat.sg.l$lower <- as.data.frame.list(colMeans(rbind(stat.sg.l.m$lower, stat.sg.l.f$lower)))
stat.sg.l$upper <- as.data.frame.list(colMeans(rbind(stat.sg.l.m$upper, stat.sg.l.f$upper)))

stat.cro.d <- stat.cro.d.f 
stat.cro.d$lower <- as.data.frame.list(colMeans(rbind(stat.cro.d.m$lower, stat.cro.d.f$lower)))
stat.cro.d$upper <- as.data.frame.list(colMeans(rbind(stat.cro.d.m$upper, stat.cro.d.f$upper)))
stat.sg.d <- stat.sg.d.f 
stat.sg.d$lower <- as.data.frame.list(colMeans(rbind(stat.sg.d.m$lower, stat.sg.d.f$lower)))
stat.sg.d$upper <- as.data.frame.list(colMeans(rbind(stat.sg.d.m$upper, stat.sg.d.f$upper)))


# output dataframe

stat.df <- data.frame(dir = rep(cov.sg.l.m$dir, 4),
                      prob = c(s.cro.l$travel, s.cro.l$search, s.cro.l$rest,
                               s.cro.d$travel, s.cro.d$search, s.cro.d$rest,
                               s.sg.l$travel, s.sg.l$search, s.sg.l$rest,
                               s.sg.d$travel, s.sg.d$search, s.sg.d$rest),
                      lwr = c(stat.cro.l[[1]]$travel, stat.cro.l[[1]]$search, stat.cro.l[[1]]$rest,
                              stat.cro.d[[1]]$travel, stat.cro.d[[1]]$search, stat.cro.d[[1]]$rest,
                              stat.sg.l[[1]]$travel, stat.sg.l[[1]]$search, stat.sg.l[[1]]$rest,
                              stat.sg.d[[1]]$travel, stat.sg.d[[1]]$search, stat.sg.d[[1]]$rest),
                      upr = c(stat.cro.l[[2]]$travel, stat.cro.l[[2]]$search, stat.cro.l[[2]]$rest,
                              stat.cro.d[[2]]$travel, stat.cro.d[[2]]$search, stat.cro.d[[2]]$rest,
                              stat.sg.l[[2]]$travel, stat.sg.l[[2]]$search, stat.sg.l[[2]]$rest,
                              stat.sg.d[[2]]$travel, stat.sg.d[[2]]$search, stat.sg.d[[2]]$rest),
                      state = rep(rep(c("travel", "search", "rest"), each = length(cov.sg.l.m$ws)), 4),
                      lod = rep(c("L", "D", "L", "D"), each = length(cov.sg.l.m$ws)*3),
                      site = rep(c("Crozet", "South Georgia"), each = length(cov.sg.l.m$ws)*6))

# plot

ggplot(stat.df, aes(lod, prob)) + facet_wrap(~site)+
  geom_errorbar(aes(ymin=lwr, ymax=upr, group = lod), width=.25, position=pd) +
  geom_point(aes(fill = state), pch = 21, size = 2.7) +
  ylim(0, 1)+theme_bw() + #ylab("Stationary probability")+
  #xlab("Photoperiod")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 16),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("yellow", "red", "blue"))
