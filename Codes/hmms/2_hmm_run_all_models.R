library(momentuHMM)
library(ggplot2)
library(plyr)


# AIM OF SCRIPT - to run all combinations of hmms based on predictors specified in Section 2.5 of Methods and outputting models, AIC tables and goodness of fit and
# autocorrelation plots


############# 1 - RUNNING ALL COVARIATE COMBINATIONS OF HMMS ###############


########## 1A. SOUTH GEORGIA ######


### LOAD IN TRACKS ####

path.sg <- "./Data_inputs/GPS_SouthGeorgia_2012.csv"
sg <- read.csv(file = path.sg, na.strings = "NA")
sg$ID <- as.factor(as.character(sg$ID))
nlevels(sg$ID) # 38


#### As it takes a very long time (i.e hours) to run on the full datasets - to test code we'll subset to 2 individuals per sex  ##### 

males <- subset(sg, sex == "M",)
males$ID <- as.factor(as.character(males$ID))
females <- subset(sg, sex == "F",)
females$ID <- as.factor(as.character(females$ID))
m.ids <- sample(levels(males$ID), 2, replace = FALSE)
f.ids <- sample(levels(females$ID), 2, replace = FALSE)
ids <- c(m.ids, f.ids)
sg <- subset(sg, ID %in% ids,)
sg$sex <- as.factor(as.character(sg$sex))
sg$ID <- as.factor(as.character(sg$ID))
nlevels(sg$ID) # 4


### PREPARE DATA FOR HMMS ####


dat <-prepData(sg, type= "LL", # longs and lats
               coordNames = c("x", "y")) ## these are our names
head(dat)

# assigning step lengths based on values associated with "best" models in the first code "1_run_simple_hmm_test_initial_values.r"
shape_0 <- c(18.73,5.59, 0.62) 
scale_0 <- c(5.59, 6.47, 0.32) 
stepPar0 <- c(shape_0,scale_0)

# artifact for zero step lengths of which there are only 2 values
table(which(dat$step == 0))
ind_zero <- which(dat$step == 0)
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}
ind_zero <- which(dat$step == "NA")
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}

# assigning turning angles based on first code
location_0 <- c(0.0022,  -0.004, 0.035) 
concentration_0 <- c(52.97,  1.31,  51.27)
anglePar0 <- c(location_0,concentration_0)




######  RUNNING A NULL MODEL ########

# first run a null model with no covariates on transition probabilities. You need to use this to set up parameter values based on your covariates
# to input into models with covariates (below)

stateNames<-c("travel","search", "rest")

m1 <- fitHMM(data=dat, nbStates=3,
             dist=list(step="gamma",angle="vm"),
             Par0=list(step=stepPar0, angle=anglePar0),
             estAngleMean = list(angle=TRUE),
             stateNames = stateNames)
# store model as an .rdata object so you don't have to run from scratch each time
file.out <- paste0("./Data_outputs/", paste0("SG_mod_", 1, ".RData"))
save(m1, file = file.out)



####### SETTING UP FORMULA AND PAR0 FOR ALL COMBINATIONS OF MODELS ##########


# The following code manually specifies each of 40 covariate combinations tested, which is detailed Methods section 2.5 of paper
formula <- list()	
formula[[2]] <- ~ ws	
formula[[3]] <- ~ ws + dir
formula[[4]] <- ~ (ws + I(ws^2))
formula[[5]] <- ~ ws + I(ws^2) + dir + ws:dir + I(ws^2):dir
formula[[6]] <- ~ sex	
formula[[7]] <- ~ lod	
formula[[8]] <- ~ ws + sex
formula[[9]] <- ~ ws + dir + sex
formula[[10]] <- ~ ws + I(ws^2) + sex	
formula[[11]] <- ~ ws + I(ws^2) + dir + sex	+ ws:dir + I(ws^2):dir
formula[[12]] <- ~ ws + lod	
formula[[13]] <- ~ ws+ ws:dir + dir + lod	
formula[[14]] <- ~ ws + I(ws^2) + lod
formula[[15]] <- ~ ws + I(ws^2) + dir + lod + ws:dir + I(ws^2):dir
formula[[16]] <- ~ lod + sex	
formula[[17]] <- ~ ws:sex + ws + sex	
formula[[18]] <- ~ ws:sex + ws:dir + ws+ dir + sex	+ ws:sex:dir
formula[[19]] <- ~ ws:sex + I(ws^2):sex + ws + I(ws^2) + sex	
formula[[20]] <- ~ ws:sex + I(ws^2):sex + ws + I(ws^2) + sex	+ dir + ws:dir + I(ws^2):dir + ws:sex:dir + I(ws^2):sex:dir
formula[[21]] <- ~ ws:lod + ws + lod	
formula[[22]] <- ~ ws:lod + ws + lod + ws:dir + dir + ws:dir:lod	
formula[[23]] <- ~ ws:lod + I(ws^2):lod + ws + I(ws^2) + lod	
formula[[24]] <- ~ ws:lod + I(ws^2):lod + ws + I(ws^2) + lod	+ dir + ws:lod:dir + I(ws^2):lod:dir + ws:dir + I(ws^2):dir
formula[[25]] <- ~ ws + lod + sex	
formula[[26]] <- ~ ws + lod + sex	+ dir + ws:dir
formula[[27]] <- ~ ws + I(ws^2) + lod + sex	
formula[[28]] <- ~ ws + I(ws^2) + lod + sex + dir + ws:dir + I(ws^2):dir
formula[[29]] <- ~ ws:sex + ws + sex + lod	
formula[[30]] <- ~ ws:sex + ws + sex + lod	+ dir + ws:sex:dir + ws:dir
formula[[31]] <- ~  ws:sex + I(ws^2):sex + ws + I(ws^2) + sex + lod	
formula[[32]] <- ~  ws:sex + I(ws^2):sex + ws + I(ws^2) + sex + lod	+ dir + ws:sex:dir + I(ws^2):sex:dir + ws:dir + I(ws^2):dir
formula[[33]] <- ~  ws:lod + ws + lod + sex	
formula[[34]] <- ~ ws:lod + ws + lod + sex	+ dir + ws:lod:dir + ws:dir
formula[[35]] <- ~ ws:lod + I(ws^2):lod + ws + I(ws^2) + lod + sex	
formula[[36]] <- ~  ws:lod + I(ws^2):lod + ws + I(ws^2) + lod + sex	+ dir + ws:lod:dir + I(ws^2):lod:dir + ws:dir + I(ws^2):dir
formula[[37]] <- ~  ws:lod + ws:sex + ws + sex + lod 	
formula[[38]] <- ~  ws:lod + ws:sex + ws + sex + lod + ws:lod:dir + ws:sex:dir + ws:dir
formula[[39]] <- ~  ws:lod + I(ws^2):lod + ws:sex + I(ws^2):sex + ws + I(ws^2) + lod + sex	
formula[[40]] <- ~  ws:lod + I(ws^2):lod + ws:sex + I(ws^2):sex + ws + I(ws^2) + lod + sex + dir + ws:dir + I(ws^2):dir + I(ws^2):sex:dir + ws:lod:dir + I(ws^2):lod:dir + ws:sex:dir 

# this function gets starting values for each model from existing null model fit for a specified covariate formula
Par <- list()
for (i in 2:length(formula)){
  Par[[i]] <- getPar0(model=m1, nbStates=3, formula = formula[[i]])
}



####### MODEL FORMULATION ########


stateNames<-c("travel","search", "rest")


##### RUNNING MODELS ###### 


# the following code iterates through and runs all 40 models, pasting each out in turn

m.list <- list()
for (i in 2:length(formula)) {
  print(i)
  m.list[[i]] <- fitHMM(data=dat, nbStates=3,
                        dist=list(step="gamma",angle="vm"),
                        Par0=list(step=Par[[i]]$Par$step, angle=Par[[i]]$Par$angle,delta0 = Par[[i]]$delta),
                        estAngleMean = list(angle=TRUE), beta0 = Par[[i]]$beta,
                        stateNames = stateNames, 
                        formula = formula[[i]])
  model <- m.list[[i]]
  file.out <- paste0("./Data_outputs/", paste0("SG_mod_", i, ".RData"))
  save(model, file = file.out)
}




############# 1B. CROZET  ###############


### LOAD IN TRACKS ####

path.cro <- "./Data_inputs/GPS_Crozet_2010-2016.csv"
cro <- read.csv(file = path.cro, na.strings = "NA")
cro$ID <- as.factor(as.character(cro$ID))
nlevels(cro$ID) # 267


#### As it takes a very long time (i.e days) to run on the full datasets - to test code we'll subset to 2 individuals per sex  ##### 

males <- subset(cro, sex == "M",)
males$ID <- as.factor(as.character(males$ID))
females <- subset(cro, sex == "F",)
females$ID <- as.factor(as.character(females$ID))
m.ids <- sample(levels(males$ID), 2, replace = FALSE)
f.ids <- sample(levels(females$ID), 2, replace = FALSE)
ids <- c(m.ids, f.ids)
cro <- subset(cro, ID %in% ids,)
cro$sex <- as.factor(as.character(cro$sex))
cro$ID <- as.factor(as.character(cro$ID))
nlevels(cro$ID) # 4


### PREPARE DATA FOR HMMS ####


dat <-prepData(cro, type= "LL", # longs and lats
               coordNames = c("x", "y")) ## these are our names
head(dat)

# assigning step lengths based on values associated with "best" models in the first code "1_run_simple_hmm_test_initial_values.r"
shape_0 <- c(12.46,3.95, 0.34) 
scale_0 <- c(3.734, 4.44, 	0.19) 
stepPar0 <- c(shape_0,scale_0)

table(which(dat$step == 0))
ind_zero <- which(dat$step == 0)
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}
ind_zero <- which(dat$step == "NA")
if (length(ind_zero)>0){
  dat$step[ind_zero] <- runif(length(ind_zero))/10000
}

# assigning turning angles based on first code
location_0 <- c(0.0033,  -0.016, 0.03) 
concentration_0 <- c(47.15,  1.16,  39.00)
anglePar0 <- c(location_0,concentration_0)





######  RUNNING A NULL MODEL ########

# first run a null model with no covariates on transition probabilities. You need to use this to set up parameter values based on your covariates
# to input into models with covariates (below)

stateNames<-c("travel","search", "rest")

m1 <- fitHMM(data=dat, nbStates=3,
             dist=list(step="gamma",angle="vm"),
             Par0=list(step=stepPar0, angle=anglePar0),
             estAngleMean = list(angle=TRUE),
             stateNames = stateNames)
# store model as an .rdata object so you don't have to run from scratch each time
file.out <- paste0("./Data_outputs/", paste0("Cro_mod_", 1, ".RData"))
save(m1, file = file.out)




####### SETTING UP FORMULA AND PAR0 FOR ALL COMBINATIONS OF MODELS ##########


# The following code manually specifies each of 40 covariate combinations tested, which is detailed Methods section 2.5 of paper
formula <- list()	
formula[[2]] <- ~ ws	
formula[[3]] <- ~ ws + dir
formula[[4]] <- ~ (ws + I(ws^2))
formula[[5]] <- ~ ws + I(ws^2) + dir + ws:dir + I(ws^2):dir
formula[[6]] <- ~ sex	
formula[[7]] <- ~ lod	
formula[[8]] <- ~ ws + sex
formula[[9]] <- ~ ws + dir + sex
formula[[10]] <- ~ ws + I(ws^2) + sex	
formula[[11]] <- ~ ws + I(ws^2) + dir + sex	+ ws:dir + I(ws^2):dir
formula[[12]] <- ~ ws + lod	
formula[[13]] <- ~ ws+ ws:dir + dir + lod	
formula[[14]] <- ~ ws + I(ws^2) + lod
formula[[15]] <- ~ ws + I(ws^2) + dir + lod + ws:dir + I(ws^2):dir
formula[[16]] <- ~ lod + sex	
formula[[17]] <- ~ ws:sex + ws + sex	
formula[[18]] <- ~ ws:sex + ws:dir + ws+ dir + sex	+ ws:sex:dir
formula[[19]] <- ~ ws:sex + I(ws^2):sex + ws + I(ws^2) + sex	
formula[[20]] <- ~ ws:sex + I(ws^2):sex + ws + I(ws^2) + sex	+ dir + ws:dir + I(ws^2):dir + ws:sex:dir + I(ws^2):sex:dir
formula[[21]] <- ~ ws:lod + ws + lod	
formula[[22]] <- ~ ws:lod + ws + lod + ws:dir + dir + ws:dir:lod	
formula[[23]] <- ~ ws:lod + I(ws^2):lod + ws + I(ws^2) + lod	
formula[[24]] <- ~ ws:lod + I(ws^2):lod + ws + I(ws^2) + lod	+ dir + ws:lod:dir + I(ws^2):lod:dir + ws:dir + I(ws^2):dir
formula[[25]] <- ~ ws + lod + sex	
formula[[26]] <- ~ ws + lod + sex	+ dir + ws:dir
formula[[27]] <- ~ ws + I(ws^2) + lod + sex	
formula[[28]] <- ~ ws + I(ws^2) + lod + sex + dir + ws:dir + I(ws^2):dir
formula[[29]] <- ~ ws:sex + ws + sex + lod	
formula[[30]] <- ~ ws:sex + ws + sex + lod	+ dir + ws:sex:dir + ws:dir
formula[[31]] <- ~  ws:sex + I(ws^2):sex + ws + I(ws^2) + sex + lod	
formula[[32]] <- ~  ws:sex + I(ws^2):sex + ws + I(ws^2) + sex + lod	+ dir + ws:sex:dir + I(ws^2):sex:dir + ws:dir + I(ws^2):dir
formula[[33]] <- ~  ws:lod + ws + lod + sex	
formula[[34]] <- ~ ws:lod + ws + lod + sex	+ dir + ws:lod:dir + ws:dir
formula[[35]] <- ~ ws:lod + I(ws^2):lod + ws + I(ws^2) + lod + sex	
formula[[36]] <- ~  ws:lod + I(ws^2):lod + ws + I(ws^2) + lod + sex	+ dir + ws:lod:dir + I(ws^2):lod:dir + ws:dir + I(ws^2):dir
formula[[37]] <- ~  ws:lod + ws:sex + ws + sex + lod 	
formula[[38]] <- ~  ws:lod + ws:sex + ws + sex + lod + ws:lod:dir + ws:sex:dir + ws:dir
formula[[39]] <- ~  ws:lod + I(ws^2):lod + ws:sex + I(ws^2):sex + ws + I(ws^2) + lod + sex	
formula[[40]] <- ~  ws:lod + I(ws^2):lod + ws:sex + I(ws^2):sex + ws + I(ws^2) + lod + sex + dir + ws:dir + I(ws^2):dir + I(ws^2):sex:dir + ws:lod:dir + I(ws^2):lod:dir + ws:sex:dir 

# this function gets starting values for each model from existing null model fit for a specified covariate formula
Par <- list()
for (i in 2:length(formula)){
  Par[[i]] <- getPar0(model=m1, nbStates=3, formula = formula[[i]])
}




####### MODEL FORMULATION ########


stateNames<-c("travel","search", "rest")


##### RUNNING MODELS ###### 


# the following code iterates through and runs all 40 models, pasting each out in turn
m.list <- list()
for (i in 2:length(formula)) {
  print(i)
  m.list[[i]] <- fitHMM(data=dat, nbStates=3,
                        dist=list(step="gamma",angle="vm"),
                        Par0=list(step=Par[[i]]$Par$step, angle=Par[[i]]$Par$angle,delta0 = Par[[i]]$delta),
                        estAngleMean = list(angle=TRUE), beta0 = Par[[i]]$beta,
                        stateNames = stateNames, 
                        formula = formula[[i]])
  model <- m.list[[i]]
  file.out <- paste0("./Data_outputs/", paste0("Cro_mod_", i, ".RData"))
  save(model, file = file.out)
  
}







###########  2. OUTPUT MODEL COEFFICIENTS AND RESIDUAL AND AUTOCORRELATION PLOTS #######



############# 2A. SOUTH GEORGIA ##########


# iterate through each model set, load in models, output autocorrelation plots for step lenghts and turning angles, pseudo-residual and qq plots, extract model 
# coefficients for each transition and plot, and paste out AIC table

m.list <- list()
out.df <- list()
for (i in 1:length(formula)) {
  print(i)
  file.in <- paste0("./Data_outputs/", paste0("SG_mod_", i, ".RData"))
  load(file = file.in)
  if (i == 1) { m.list[[i]] <- m1} else { m.list[[i]] <- model}
  # outputting AICs into dataframe
  if (i == 1) { form_out <- 1} else { form_out <- as.character(formula[[i]])[2]}
  out.df[[i]] <- data.frame(Model = paste0("M", i),
                            Formula = form_out, AIC = AIC(m.list[[i]]))
  
  # outputting autocorrelation plots
  pr <- pseudoRes(m.list[[i]])
  # step length
  par(mfrow=c(1,1))
  acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)
  name.plot <- paste0("./Plots/", paste0("SG_mod_", i, "_acf_step.png"))
  dev.copy(png, name.plot,  width = 800, height = 500)
  dev.off()
  # turning angle
  par(mfrow=c(1,1))
  acf(pr$angleRes[!is.na(pr$angleRes)],lag.max = 300)
  name.plot <- paste0("./Plots/", paste0("SG_mod_", i, "_acf_angle.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
  
  # outputting pseudo residual plots
  res.df <- data.frame(stepres = pr$stepRes, angleres = pr$angleRes,
                       sex = m.list[[i]]$data$sex, ws = m.list[[i]]$data$ws, 
                       lod = m.list[[i]]$data$lod)
  par(mfrow=c(2,3))
  boxplot(stepres ~ sex, data = res.df)
  boxplot(stepres ~ lod, data = res.df)
  plot(stepres ~ ws, data = res.df)
  boxplot(angleres ~ sex, data = res.df)
  boxplot(angleres ~ lod, data = res.df)
  plot(angleres ~ ws, data = res.df)
  name.plot <- paste0("./Plots/", paste0("SG_mod_", i, "_pseudo-residuals.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
  
  # outputting qq plots
  plotPR(m.list[[i]])
  name.plot <- paste0("./Plots/", paste0("SG_mod_", i, "_acf_qq.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off() 
  
  # extracting and plotting model coefficients for each transition
  beta.full <- CIbeta(m.list[[i]])$beta
  beta.full.est <- as.data.frame(beta.full$est)
  beta.full.upr <- as.data.frame(beta.full$upper)
  beta.full.lwr <- as.data.frame(beta.full$lower)
  beta.df <- data.frame(Est = c(beta.full.est$`1 -> 2`, beta.full.est$`1 -> 3`,
                                beta.full.est$`2 -> 1`, beta.full.est$`2 -> 3`,
                                beta.full.est$`3 -> 1`, beta.full.est$`3 -> 2`), 
                        Upr = c(beta.full.upr$`1 -> 2`, beta.full.upr$`1 -> 3`,
                                beta.full.upr$`2 -> 1`, beta.full.upr$`2 -> 3`,
                                beta.full.upr$`3 -> 1`, beta.full.upr$`3 -> 2`), 
                        Lwr = c(beta.full.lwr$`1 -> 2`, beta.full.lwr$`1 -> 3`,
                                beta.full.lwr$`2 -> 1`, beta.full.lwr$`2 -> 3`,
                                beta.full.lwr$`3 -> 1`, beta.full.lwr$`3 -> 2`), 
                        Transitions = rep(colnames(beta.full.est), each = nrow(beta.full.est)),
                        Covariates = rep(rownames(beta.full.est), 3))
  beta.df$Covariates <- as.factor(as.character(beta.df$Covariates))
  beta.df$Transitions <- as.factor(as.character(beta.df$Transitions))
  pd <- position_dodge(width=0.7)
  # removing intercept to plot
  beta.df2 <- subset(beta.df, Covariates != "(Intercept)",)
  pl <- ggplot(beta.df2, aes(Covariates, Est)) + geom_hline(yintercept=0, linetype="dashed", size=1)+
    geom_point(aes(colour = Transitions),position=pd)+
    geom_errorbar(aes(ymin=Lwr, ymax=Upr, colour = Transitions), width=.8, position=pd) +theme_bw()
  print(pl)
  name.plot <- paste0("./Plots/", paste0("SG_mod_", i, "_coefficients.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
}
all.out <- do.call(rbind, out.df)
all.out <- all.out[order(all.out$AIC),]

# paste out AIC table
out.path <- "./Data_outputs/SG_AIC_table.csv"
write.csv(all.out, out.path, row.names=T)








############# 2B. CROZET ##########


# iterate through each model set, load in models, output autocorrelation plots for step lenghts and turning angles, pseudo-residual and qq plots, extract model 
# coefficients for each transition and plot, and paste out AIC table

m.list <- list()
out.df <- list()
for (i in 1:length(formula)) {
  print(i)
  file.in <- paste0("./Data_outputs/", paste0("Cro_mod_", i, ".RData"))
  load(file = file.in)
  if (i == 1) { m.list[[i]] <- m1} else { m.list[[i]] <- model}
  # outputting AICs into dataframe
  if (i == 1) { form_out <- 1} else { form_out <- as.character(formula[[i]])[2]}
  out.df[[i]] <- data.frame(Model = paste0("M", i),
                            Formula = form_out, AIC = AIC(m.list[[i]]))
  
  # outputting autocorrelation plots
  pr <- pseudoRes(m.list[[i]])
  # step length
  par(mfrow=c(1,1))
  acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)
  name.plot <- paste0("./Plots/", paste0("Cro_mod_", i, "_acf_step.png"))
  dev.copy(png, name.plot,  width = 800, height = 500)
  dev.off()
  # turning angle
  par(mfrow=c(1,1))
  acf(pr$angleRes[!is.na(pr$angleRes)],lag.max = 300)
  name.plot <- paste0("./Plots/", paste0("Cro_mod_", i, "_acf_angle.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
  
  # outputting pseudo residual plots
  res.df <- data.frame(stepres = pr$stepRes, angleres = pr$angleRes,
                       sex = m.list[[i]]$data$sex, ws = m.list[[i]]$data$ws, 
                       lod = m.list[[i]]$data$lod)
  par(mfrow=c(2,3))
  boxplot(stepres ~ sex, data = res.df)
  boxplot(stepres ~ lod, data = res.df)
  plot(stepres ~ ws, data = res.df)
  boxplot(angleres ~ sex, data = res.df)
  boxplot(angleres ~ lod, data = res.df)
  plot(angleres ~ ws, data = res.df)
  name.plot <- paste0("./Plots/", paste0("Cro_mod_", i, "_pseudo-residuals.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
  
  # outputting qq plots
  plotPR(m.list[[i]])
  name.plot <- paste0("./Plots/", paste0("Cro_mod_", i, "_acf_qq.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off() 
  
  # extracting and plotting model coefficients for each transition
  beta.full <- CIbeta(m.list[[i]])$beta
  beta.full.est <- as.data.frame(beta.full$est)
  beta.full.upr <- as.data.frame(beta.full$upper)
  beta.full.lwr <- as.data.frame(beta.full$lower)
  beta.df <- data.frame(Est = c(beta.full.est$`1 -> 2`, beta.full.est$`1 -> 3`,
                                beta.full.est$`2 -> 1`, beta.full.est$`2 -> 3`,
                                beta.full.est$`3 -> 1`, beta.full.est$`3 -> 2`), 
                        Upr = c(beta.full.upr$`1 -> 2`, beta.full.upr$`1 -> 3`,
                                beta.full.upr$`2 -> 1`, beta.full.upr$`2 -> 3`,
                                beta.full.upr$`3 -> 1`, beta.full.upr$`3 -> 2`), 
                        Lwr = c(beta.full.lwr$`1 -> 2`, beta.full.lwr$`1 -> 3`,
                                beta.full.lwr$`2 -> 1`, beta.full.lwr$`2 -> 3`,
                                beta.full.lwr$`3 -> 1`, beta.full.lwr$`3 -> 2`), 
                        Transitions = rep(colnames(beta.full.est), each = nrow(beta.full.est)),
                        Covariates = rep(rownames(beta.full.est), 3))
  beta.df$Covariates <- as.factor(as.character(beta.df$Covariates))
  beta.df$Transitions <- as.factor(as.character(beta.df$Transitions))
  pd <- position_dodge(width=0.7)
  # removing intercept to plot
  beta.df2 <- subset(beta.df, Covariates != "(Intercept)",)
  pl <- ggplot(beta.df2, aes(Covariates, Est)) + geom_hline(yintercept=0, linetype="dashed", size=1)+
    geom_point(aes(colour = Transitions),position=pd)+
    geom_errorbar(aes(ymin=Lwr, ymax=Upr, colour = Transitions), width=.8, position=pd) +theme_bw()
  print(pl)
  name.plot <- paste0("./Plots/", paste0("Cro_mod_", i, "_coefficients.png"))
  dev.copy(png, name.plot, width = 800, height = 500)
  dev.off()
}
all.out <- do.call(rbind, out.df)
all.out <- all.out[order(all.out$AIC),]

# paste out AIC table
out.path <- "./Data_outputs/Cro_AIC_table.csv"
write.csv(all.out, out.path, row.names=T)

