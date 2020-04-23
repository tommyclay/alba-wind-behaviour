library(momentuHMM)
library(plyr)


# AIM OF SCRIPT - to run the "best" model but with covariate values reshuffled to test for over-parameterization



############ RUNNING "BEST" MODELS WITH COVARIATE VALUES RESHUFFLED TO TEST FOR OVER-PARAMETERIZATION ########



####### 1. SOUTH GEORGIA ########


#####  1A. LOAD IN NULL AND "BEST" MODEL BASED ON MODEL SELECTION AND PREPARE DATA #####


# null model
file.in <- paste0("./Data_outputs/", paste0("SG_mod_", 1, ".RData"))
load(file = file.in) # null model is called m1
stateNames<-c("travel","search", "rest")

# best model
file.in <- paste0("./Data_outputs/", paste0("SG_mod_", 30, ".RData"))
load(file = file.in)
best.mod <- model
form <- ~ws:sex + ws + sex + lod	+ dir + ws:sex:dir + ws:dir


# resampling sex between (rather than within) individual trips #

# load in data from model
gps <- model$data
# subset male and female ids to later assign fake sexes
male.sub <- subset(gps, sex == "M",)
male.sub$ID <- as.factor(as.character(male.sub$ID))
female.sub <- subset(gps, sex == "F",)
female.sub$ID <- as.factor(as.character(female.sub$ID))
nlevels(male.sub$ID) # 2
nlevels(female.sub$ID) # 2

# for 50 iterations create randomized sex dataframe
n.iter <- 50
# summarize sex for each individual
sex.summ <- ddply(gps,.(ID), summarize, sex = sex[1])
# create fake data.frame with which to paste randomize sex values onto
full.sex.df <- data.frame(Id = rep(levels(gps$ID), n.iter), Sex_real = rep(sex.summ$sex, n.iter),
                          Sex = "NA",
                          Iter = rep(seq(1, 100, 1), each = nlevels(gps$ID)))
full.sex.df$Sex<- as.character(full.sex.df$Sex)

# run through and assign fake sexes to each iteration and paste out how many times fake sexes match real sexes
same_df <- data.frame(Niter = seq(1, n.iter, 1), Sexes_same = "FALSE")
same_df$Sexes_same <- as.character(same_df$Sexes_same)
for (i in 1:n.iter)  {
  print(i)
  # assign fake sex
  sexes <- sample(c(rep("M", nlevels(male.sub$ID)), rep("F", nlevels(female.sub$ID))), length(unique(gps$ID)), replace = FALSE)
  full.sex.df$Sex[full.sex.df$Iter == i] <- sexes
  # checking that every individual doesn't have to same sex as the observed dataset
  ToF <- full.sex.df$Sex[full.sex.df$Iter == i]==full.sex.df$Sex_real[full.sex.df$Iter == i]
  if (all(ToF == TRUE)) { same_df$Sexes_same[i] <- "TRUE"} else { }
}
# this shows how many cases the sexes were same as real dataset
table(same_df$Sexes_same)
# for reduced dataset (i.e. n = 4), this will be common, but not for full dataset 
# you want to iterate again so that all are different from real dataset - i.e. you want all same_df$Sexes_same to be "FALSE"




# RUN 50 MODELS TO RANDOMIZE EACH COVARIATE IN TURN - SO RUN SEPARATELY FOR THE 4 COVARIATES - WIND SPEED, DIRECTION, LOD AND SEX

# setting step lengths and turn angles outside of loop

# assigning step lengths
shape_0 <- c(18.73,5.59, 0.62) 
scale_0 <- c(5.59, 6.47, 0.32) 
stepPar0 <- c(shape_0,scale_0)

# assigning turning angles
location_0 <- c(0.0022,  -0.004, 0.035) 
concentration_0 <- c(52.97,  1.31,  51.27)
anglePar0 <- c(location_0,concentration_0)


######### 1B. RUNNING MODELS RANDOMIZING BY SEX #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # assign random sexes from fake sex dataframe above
  sex.df <- subset(full.sex.df, Iter == i,)
  random.f <- subset(sex.df, Sex == "F",)
  gps.l[[i]]$sex_r <- "M"
  gps.l[[i]]$sex_r[gps.l[[i]]$ID %in% random.f$Id] <- "F"
  gps.l[[i]]$sex_r <- as.factor(as.character(gps.l[[i]]$sex_r))
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex_r, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws, dir =  gps.l[[i]]$dir, lod = gps.l[[i]]$lod)
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
   # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form) # getting parameters from null model
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/SG_sex_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)



######### 1C. RUNNING MODELS RANDOMIZING BY WIND SPEED #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # randomize wind
  gps.l[[i]]$ws_r <- sample(gps.l[[i]]$ws)
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws_r, dir =  gps.l[[i]]$dir, lod = gps.l[[i]]$lod)
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
  # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form)
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/SG_ws_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)




######### 1D. RUNNING MODELS RANDOMIZING BY LOD #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # randomize lod
  gps.l[[i]]$lod_r <- sample(gps.l[[i]]$lod)
  gps.l[[i]]$lod_r <- as.factor(as.character(gps.l[[i]]$lod_r))
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws, dir =  gps.l[[i]]$dir, lod = gps.l[[i]]$lod_r)
  
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
    # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form)
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/SG_lod_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)


######### 1E. RUNNING MODELS RANDOMIZING BY WIND DIRECTION #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # randomize wind direction
  gps.l[[i]]$dir_r <- sample(gps.l[[i]]$dir)
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws, dir =  gps.l[[i]]$dir_r, lod = gps.l[[i]]$lod)
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
    # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form)
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/SG_dir_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)



########## 1F. PLOTTING IN RELATION TO REAL AIC VALUE #########


# AIC of "real" model
AIC.real <- AIC(best.mod)

# load in AIC tables for each covariate
path.sex <- "./Data_outputs/SG_sex_randomized_AIC_table.csv"
sex <- read.csv(file = path.sex, na.strings = "NA")
path.ws <- "./Data_outputs/SG_ws_randomized_AIC_table.csv"
sex <- read.csv(file = path.ws, na.strings = "NA")
path.dir <- "./Data_outputs/SG_dir_randomized_AIC_table.csv"
dir <- read.csv(file = path.dir, na.strings = "NA")
path.lod <- "./Data_outputs/SG_dir_randomized_AIC_table.csv"
lod <- read.csv(file = path.lod, na.strings = "NA")

sex$Covar <- "sex"
ws$Covar <- "ws"
lod$Covar <- "lod"
dir$Covar <- "dir"
all.merge <- rbind(sex, ws, lod,  dir)

# plotting randomized AICs in relation to AIC of real "best" model - this is the basis of Figure S4 in Appendix S3
real.df <- data.frame(AIC = rep(AIC.real), Covar <- c("sex", "ws", "lod", "dir"))
ggplot(all.merge, aes(Covar, AIC)) + geom_boxplot()+geom_hline(yintercept=AIC.real, linetype="dashed", 
                                                               color = "red", size = 1.1)+xlab("Dummy covariates")








####### 2. CROZET  ########


#####  2A. LOAD IN NULL AND "BEST" MODEL BASED ON MODEL SELECTION AND PREPARE DATA #####


# null model
file.in <- paste0("./Data_outputs/", paste0("Cro_mod_", 1, ".RData"))
load(file = file.in) # null model is called m1
stateNames<-c("travel","search", "rest")

# best model
file.in <- paste0("./Data_outputs/", paste0("Cro_mod_", 30, ".RData"))
load(file = file.in)
best.mod <- model
form <- ~ws:sex + ws + sex + lod	+ dir + ws:sex:dir + ws:dir


# resampling sex between (rather than within) individual trips #

# load in data from model
gps <- model$data
# subset male and female ids to later assign fake sexes
male.sub <- subset(gps, sex == "M",)
male.sub$ID <- as.factor(as.character(male.sub$ID))
female.sub <- subset(gps, sex == "F",)
female.sub$ID <- as.factor(as.character(female.sub$ID))
nlevels(male.sub$ID) # 2
nlevels(female.sub$ID) # 2

# for 50 iterations create randomized sex dataframe
n.iter <- 50
# summarize sex for each individual
sex.summ <- ddply(gps,.(ID), summarize, sex = sex[1])
# create fake data.frame with which to paste randomize sex values onto
full.sex.df <- data.frame(Id = rep(levels(gps$ID), n.iter), Sex_real = rep(sex.summ$sex, n.iter),
                          Sex = "NA",
                          Iter = rep(seq(1, 100, 1), each = nlevels(gps$ID)))
full.sex.df$Sex<- as.character(full.sex.df$Sex)

# run through and assign fake sexes to each iteration and paste out how many times fake sexes match real sexes
same_df <- data.frame(Niter = seq(1, n.iter, 1), Sexes_same = "FALSE")
same_df$Sexes_same <- as.character(same_df$Sexes_same)
for (i in 1:n.iter)  {
  print(i)
  # assign fake sex
  sexes <- sample(c(rep("M", nlevels(male.sub$ID)), rep("F", nlevels(female.sub$ID))), length(unique(gps$ID)), replace = FALSE)
  full.sex.df$Sex[full.sex.df$Iter == i] <- sexes
  # checking that every individual doesn't have to same sex as the observed dataset
  ToF <- full.sex.df$Sex[full.sex.df$Iter == i]==full.sex.df$Sex_real[full.sex.df$Iter == i]
  if (all(ToF == TRUE)) { same_df$Sexes_same[i] <- "TRUE"} else { }
}
# this shows how many cases the sexes were same as real dataset
table(same_df$Sexes_same)
# for reduced dataset (i.e. n = 4), this will be common, but not for full dataset 
# you want to iterate again so that all are different from real dataset - i.e. you want all same_df$Sexes_same to be "FALSE"




# RUN 50 MODELS TO RANDOMIZE EACH COVARIATE IN TURN - SO RUN SEPARATELY FOR THE 4 COVARIATES - WIND SPEED, DIRECTION, LOD AND SEX

# setting step lengths and turn angles outside of loop

# assigning step lengths
shape_0 <- c(12.46,3.95, 0.34) 
scale_0 <- c(3.734, 4.44, 	0.19) 
stepPar0 <- c(shape_0,scale_0)

# assigning turning angles
location_0 <- c(0.0033,  -0.016, 0.03) 
concentration_0 <- c(47.15,  1.16,  39.00)
anglePar0 <- c(location_0,concentration_0)



######### 2B. RUNNING MODELS RANDOMIZING BY SEX #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # assign random sexes from fake sex dataframe above
  sex.df <- subset(full.sex.df, Iter == i,)
  random.f <- subset(sex.df, Sex == "F",)
  gps.l[[i]]$sex_r <- "M"
  gps.l[[i]]$sex_r[gps.l[[i]]$ID %in% random.f$Id] <- "F"
  gps.l[[i]]$sex_r <- as.factor(as.character(gps.l[[i]]$sex_r))
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex_r, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws, dir =  gps.l[[i]]$dir, lod = gps.l[[i]]$lod)
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
  # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form) # getting parameters from null model
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/Cro_sex_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)



######### 2C. RUNNING MODELS RANDOMIZING BY WIND SPEED #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # randomize wind
  gps.l[[i]]$ws_r <- sample(gps.l[[i]]$ws)
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws_r, dir =  gps.l[[i]]$dir, lod = gps.l[[i]]$lod)
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
  # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form)
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/Cro_ws_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)




######### 2D. RUNNING MODELS RANDOMIZING BY LOD #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # randomize lod
  gps.l[[i]]$lod_r <- sample(gps.l[[i]]$lod)
  gps.l[[i]]$lod_r <- as.factor(as.character(gps.l[[i]]$lod_r))
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws, dir =  gps.l[[i]]$dir, lod = gps.l[[i]]$lod_r)
  
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
  # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form)
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/Cro_lod_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)


######### 2E. RUNNING MODELS RANDOMIZING BY WIND DIRECTION #######


mod <- list()
gps.l <- list()
AIC.df <- list()
for (i in 1:n.iter) {
  print(i)
  gps.l[[i]] <- gps
  # randomize wind direction
  gps.l[[i]]$dir_r <- sample(gps.l[[i]]$dir)
  
  # setting up data.frame
  dat =  data.frame(ID = gps.l[[i]]$ID, sex = gps.l[[i]]$sex, x = gps.l[[i]]$x, y = gps.l[[i]]$y, 
                    ws = gps.l[[i]]$ws, dir =  gps.l[[i]]$dir_r, lod = gps.l[[i]]$lod)
  
  # Use prepData to create input object
  dat <-prepData(dat, type= "LL", # longs and lats
                 coordNames = c("x", "y")) ## these are our names
  # artifact for zero step lengths
  ind_zero <- which(dat$step == 0)
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  ind_zero <- which(dat$step == "NA")
  if (length(ind_zero)>0){
    dat$step[ind_zero] <- runif(length(ind_zero))/10000
  }
  
  # running model
  par1 <- getPar0(model=m1, nbStates=3, formula = form)
  mod[[i]] <- fitHMM(data=dat, nbStates=3,
                     dist=list(step="gamma",angle="vm"),
                     Par0=list(step=par1$Par$step, angle=par1$Par$angle,delta0 = par1$delta),
                     estAngleMean = list(angle=TRUE), beta0 = par1$beta,
                     stateNames = stateNames, 
                     formula = form)
  model <- mod[[i]]
  AIC.df[[i]] <- data.frame(Iter = i, Formula = as.character(form)[2], AIC = AIC(model))
}
all.AIC <- do.call(rbind, AIC.df)
all.AIC <- all.AIC[order(all.AIC$AIC),]
all.AIC

# paste out
out.path <- "./Data_outputs/Cro_dir_randomized_AIC_table.csv"
write.csv(all.AIC, out.path, row.names=T)



########## 2F. PLOTTING IN RELATION TO REAL AIC VALUE #########


# AIC of "real" model
AIC.real <- AIC(best.mod)

# load in AIC tables for each covariate
path.sex <- "./Data_outputs/Cro_sex_randomized_AIC_table.csv"
sex <- read.csv(file = path.sex, na.strings = "NA")
path.ws <- "./Data_outputs/Cro_ws_randomized_AIC_table.csv"
sex <- read.csv(file = path.ws, na.strings = "NA")
path.dir <- "./Data_outputs/Cro_dir_randomized_AIC_table.csv"
dir <- read.csv(file = path.dir, na.strings = "NA")
path.lod <- "./Data_outputs/Cro_dir_randomized_AIC_table.csv"
lod <- read.csv(file = path.lod, na.strings = "NA")

sex$Covar <- "sex"
ws$Covar <- "ws"
lod$Covar <- "lod"
dir$Covar <- "dir"
all.merge <- rbind(sex, ws, lod,  dir)

# plotting randomized AICs in relation to AIC of real "best" model - this is the basis of Figure S4 in Appendix S3
real.df <- data.frame(AIC = rep(AIC.real), Covar <- c("sex", "ws", "lod", "dir"))
ggplot(all.merge, aes(Covar, AIC)) + geom_boxplot()+geom_hline(yintercept=AIC.real, linetype="dashed", 
                                                               color = "red", size = 1.1)+xlab("Dummy covariates")




