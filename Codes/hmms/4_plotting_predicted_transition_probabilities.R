library(ggplot2)
library(momentuHMM)
library(plyr)


# AIM OF SCRIPT - to plot predicted transition probabilites for a range of specified covariates and produce Figure 4 and outputs for Table 3


##### LOAD IN MODEL OUTPUTS #####


# load in best models for each site # remember these models have been trained on samples in the previous scripts
file.in <- paste0("./Data_outputs/", paste0("SG_mod_", 30, ".RData"))
load(file = file.in)
m.sg <- model

file.in <- paste0("./Data_outputs/", paste0("Cro_mod_", 30, ".RData"))
load(file = file.in)
m.cro <- model




##### 1. PLOTTING TRANSITION PROBABILITIES FOR CHANGES IN WIND SPEED - FOR FIGURE 4 ###############


# predict wind speed for average wind direction which is 75 degrees relative to bird direction, and for daylight hours when birds are most active. 


#### 1A. CREATE DATAFRAMES OF PREDICTED VALUES #####


# predict for a range of wind speeds separately for each population and sex 


# CROZET MALE #

cro.cov.m <- data.frame(ws=seq(from=min(m.cro$data$ws), max(m.cro$data$ws),by=0.2), lod = "L", 
                        sex = "M", dir = 75)
cro.cov.m$sex <- as.factor(as.character(cro.cov.m$sex))
cro.cov.m$lod <- as.factor(as.character(cro.cov.m$lod))
cro.cov.m$ws <- as.numeric(as.character(cro.cov.m$ws))
cro.cov.m$dir <- as.numeric(as.character(cro.cov.m$dir))
head(cro.cov.m)

# function CIreal only allows single row of covariates, so function iterates through range of covariate values and outputs list of predicted
# transition probabilities and upper and lower CIs 
ci.list <- lapply(1:nrow(cro.cov.m), function(x) {
  print(x)
  cov.sub.df <- cro.cov.m[x,]
  return(CIreal(m.cro,covs=cov.sub.df)$gamma)
})
# extract means, upper and lower bounds and put into separate lists
cro.means.m <- lapply(ci.list, '[[', 1) 
cro.lb.m <- lapply(ci.list,'[[',3)
cro.ub.m <- lapply(ci.list,'[[',4)


# CRO FEMALE #

cro.cov.f <- data.frame(ws=seq(from=min(m.cro$data$ws), max(m.cro$data$ws),by=0.2), lod = "L", 
                        sex = "F", dir = 75)
cro.cov.f$sex <- as.factor(as.character(cro.cov.f$sex))
cro.cov.f$lod <- as.factor(as.character(cro.cov.f$lod))
cro.cov.f$ws <- as.numeric(as.character(cro.cov.f$ws))
cro.cov.f$dir <- as.numeric(as.character(cro.cov.f$dir))

ci.list <- lapply(1:nrow(cro.cov.f), function(x) {
  print(x)
  cov.sub.df <- cro.cov.f[x,]
  return(CIreal(m.cro,covs=cov.sub.df)$gamma)
})

cro.means.f <- lapply(ci.list, '[[', 1) 
cro.lb.f <- lapply(ci.list,'[[',3)
cro.ub.f <- lapply(ci.list,'[[',4)


# SG MALE #

sg.cov.m <- data.frame(sex = "M", ws=seq(from=min(m.sg$data$ws), max(m.sg$data$ws),by=0.2), dir = 75,
                       lod = "L")
sg.cov.m$sex <- as.factor(as.character(sg.cov.m$sex))
sg.cov.m$lod <- as.factor(as.character(sg.cov.m$lod))
sg.cov.m$ws <- as.numeric(as.character(sg.cov.m$ws))
sg.cov.m$dir <- as.numeric(as.character(sg.cov.m$dir))

ci.list <- lapply(1:nrow(sg.cov.m), function(x) {
  print(x)
  cov.sub.df <- sg.cov.m[x,]
  return(CIreal(m.sg,covs=cov.sub.df)$gamma)
})

sg.means.m <- lapply(ci.list, '[[', 1) 
sg.lb.m <- lapply(ci.list,'[[',3)
sg.ub.m <- lapply(ci.list,'[[',4)


# SG FEMALE #

sg.cov.f <- data.frame(sex = "F", ws=seq(from=min(m.sg$data$ws), max(m.sg$data$ws),by=0.2), dir = 75,
                       lod = "L")
sg.cov.f$sex <- as.factor(as.character(sg.cov.f$sex))
sg.cov.f$lod <- as.factor(as.character(sg.cov.f$lod))
sg.cov.f$ws <- as.numeric(as.character(sg.cov.f$ws))
sg.cov.f$dir <- as.numeric(as.character(sg.cov.f$dir))

ci.list <- lapply(1:nrow(sg.cov.f), function(x) {
  print(x)
  cov.sub.df <- sg.cov.f[x,]
  return(CIreal(m.sg,covs=cov.sub.df)$gamma)
})

sg.means.f <- lapply(ci.list, '[[', 1) 
sg.lb.f <- lapply(ci.list,'[[',3)
sg.ub.f <- lapply(ci.list,'[[',4)




#### 1B. PLOTTING ALL TRANSITIONS INDIVIDUALLY #####


#### 1 -> 1, i.e. directed travel -> directed travel #####

state1 = 1
state2 = 1
nb_states = 3

cro.1_1_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.1_1_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.1_1_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.1_1_means.m,lower_bound=cro.1_1_lb.m,upper_bound=cro.1_1_ub.m)


cro.1_1_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.1_1_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.1_1_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.1_1_means.f,lower_bound=cro.1_1_lb.f,upper_bound=cro.1_1_ub.f)


sg.1_1_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.1_1_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.1_1_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.1_1_means.m,lower_bound=sg.1_1_lb.m,upper_bound=sg.1_1_ub.m)


sg.1_1_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.1_1_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.1_1_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.1_1_means.f,lower_bound=sg.1_1_lb.f,upper_bound=sg.1_1_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0.5, 1) +
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))



#### 1 -> 2, i.e. directed travel to search #####

state1 = 1
state2 = 2
nb_states = 3

cro.1_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.1_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.1_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.1_2_means.m,lower_bound=cro.1_2_lb.m,upper_bound=cro.1_2_ub.m)


cro.1_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.1_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.1_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.1_2_means.f,lower_bound=cro.1_2_lb.f,upper_bound=cro.1_2_ub.f)


sg.1_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.1_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.1_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.1_2_means.m,lower_bound=sg.1_2_lb.m,upper_bound=sg.1_2_ub.m)


sg.1_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.1_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.1_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.1_2_means.f,lower_bound=sg.1_2_lb.f,upper_bound=sg.1_2_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))



#### 1 -> 3, i.e. directed travel to rest #####

state1 = 1
state2 = 3
nb_states = 3

cro.1_3_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.1_3_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.1_3_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.1_3_means.m,lower_bound=cro.1_3_lb.m,upper_bound=cro.1_3_ub.m)


cro.1_3_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.1_3_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.1_3_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.1_3_means.f,lower_bound=cro.1_3_lb.f,upper_bound=cro.1_3_ub.f)


sg.1_3_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.1_3_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.1_3_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.1_3_means.m,lower_bound=sg.1_3_lb.m,upper_bound=sg.1_3_ub.m)


sg.1_3_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.1_3_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.1_3_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.1_3_means.f,lower_bound=sg.1_3_lb.f,upper_bound=sg.1_3_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))



#### 2 -> 1, i.e. search to directed travel #####

state1 = 2
state2 = 1
nb_states = 3

cro.2_1_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.2_1_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.2_1_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.2_1_means.m,lower_bound=cro.2_1_lb.m,upper_bound=cro.2_1_ub.m)


cro.2_1_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.2_1_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.2_1_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.2_1_means.f,lower_bound=cro.2_1_lb.f,upper_bound=cro.2_1_ub.f)


sg.2_1_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.2_1_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.2_1_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.2_1_means.m,lower_bound=sg.2_1_lb.m,upper_bound=sg.2_1_ub.m)


sg.2_1_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.2_1_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.2_1_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.2_1_means.f,lower_bound=sg.2_1_lb.f,upper_bound=sg.2_1_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 2 -> 2, i.e. search to search #####

state1 = 2
state2 = 2
nb_states = 3

cro.2_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.2_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.2_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.2_2_means.m,lower_bound=cro.2_2_lb.m,upper_bound=cro.2_2_ub.m)


cro.2_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.2_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.2_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.2_2_means.f,lower_bound=cro.2_2_lb.f,upper_bound=cro.2_2_ub.f)


sg.2_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.2_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.2_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.2_2_means.m,lower_bound=sg.2_2_lb.m,upper_bound=sg.2_2_ub.m)


sg.2_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.2_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.2_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.2_2_means.f,lower_bound=sg.2_2_lb.f,upper_bound=sg.2_2_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0.5, 1)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) +  #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 2 -> 3, i.e. search to rest #####

state1 = 2
state2 = 3
nb_states = 3

cro.2_3_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.2_3_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.2_3_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.2_3_means.m,lower_bound=cro.2_3_lb.m,upper_bound=cro.2_3_ub.m)


cro.2_3_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.2_3_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.2_3_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.2_3_means.f,lower_bound=cro.2_3_lb.f,upper_bound=cro.2_3_ub.f)


sg.2_3_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.2_3_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.2_3_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.2_3_means.m,lower_bound=sg.2_3_lb.m,upper_bound=sg.2_3_ub.m)


sg.2_3_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.2_3_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.2_3_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.2_3_means.f,lower_bound=sg.2_3_lb.f,upper_bound=sg.2_3_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))



#### 3 -> 1, i.e. search to directed travel #####

state1 = 3
state2 = 1
nb_states = 3

cro.3_1_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.3_1_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_1_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.3_1_means.m,lower_bound=cro.3_1_lb.m,upper_bound=cro.3_1_ub.m)


cro.3_1_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.3_1_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_1_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.3_1_means.f,lower_bound=cro.3_1_lb.f,upper_bound=cro.3_1_ub.f)


sg.3_1_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.3_1_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_1_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.3_1_means.m,lower_bound=sg.3_1_lb.m,upper_bound=sg.3_1_ub.m)


sg.3_1_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.3_1_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_1_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.3_1_means.f,lower_bound=sg.3_1_lb.f,upper_bound=sg.3_1_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  ylim(0, 0.5)+theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 23)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 3 -> 2, i.e. rest to search #####

state1 = 3
state2 = 2
nb_states = 3

cro.3_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.3_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.3_2_means.m,lower_bound=cro.3_2_lb.m,upper_bound=cro.3_2_ub.m)


cro.3_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.3_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.3_2_means.f,lower_bound=cro.3_2_lb.f,upper_bound=cro.3_2_ub.f)


sg.3_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.3_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.3_2_means.m,lower_bound=sg.3_2_lb.m,upper_bound=sg.3_2_ub.m)


sg.3_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.3_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.3_2_means.f,lower_bound=sg.3_2_lb.f,upper_bound=sg.3_2_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  theme_bw() + #ylab("Transition probability")+
  # scale_y_continuous(limits=c(0, 0.64), breaks = seq(0, 0.6, 0.1)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 3 -> 3, i.e. rest to rest #####

state1 = 3
state2 = 3
nb_states = 3

cro.3_3_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.3_3_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_3_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.3_3_means.m,lower_bound=cro.3_3_lb.m,upper_bound=cro.3_3_ub.m)


cro.3_3_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.3_3_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_3_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.3_3_means.f,lower_bound=cro.3_3_lb.f,upper_bound=cro.3_3_ub.f)


sg.3_3_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.3_3_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_3_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.3_3_means.m,lower_bound=sg.3_3_lb.m,upper_bound=sg.3_3_ub.m)


sg.3_3_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.3_3_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_3_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.3_3_means.f,lower_bound=sg.3_3_lb.f,upper_bound=sg.3_3_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = wind, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  theme_bw() + #ylab("Transition probability")+
  # scale_y_continuous(limits=c(0.36, 1), breaks = seq(0.4, 1, 0.1)) + #xlab("Wind speed (ms-1)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))





##### 2. PLOTTING TRANSITION PROBABILITIES FOR CHANGES IN WIND DIRECTION - FOR FIGURE 4 ###############


# predict wind direction for average wind speed which is 9 ms-1, and for daylight hours when birds are most active. 


#### 2A. CREATE DATAFRAMES OF PREDICTED VALUES #####


# CRO MALE #

cro.cov.m <- data.frame(dir=seq(from=min(m.cro$data$dir), max(m.cro$data$dir),by=1), lod = "L", 
                        sex = "M", ws = 9)
cro.cov.m$sex <- as.factor(as.character(cro.cov.m$sex))
cro.cov.m$lod <- as.factor(as.character(cro.cov.m$lod))
cro.cov.m$ws <- as.numeric(as.character(cro.cov.m$ws))
cro.cov.m$dir <- as.numeric(as.character(cro.cov.m$dir))

ci.list <- lapply(1:nrow(cro.cov.m), function(x) {
  print(x)
  cov.sub.df <- cro.cov.m[x,]
  return(CIreal(m.cro,covs=cov.sub.df)$gamma)
})

cro.means.m <- lapply(ci.list, '[[', 1) 
cro.lb.m <- lapply(ci.list,'[[',3)
cro.ub.m <- lapply(ci.list,'[[',4)


# CRO FEMALE #


cro.cov.f <- data.frame(dir=seq(from=min(m.cro$data$dir), max(m.cro$data$dir),by=1), lod = "L", 
                        sex = "F", ws = 9)
cro.cov.f$sex <- as.factor(as.character(cro.cov.f$sex))
cro.cov.f$lod <- as.factor(as.character(cro.cov.f$lod))
cro.cov.f$ws <- as.numeric(as.character(cro.cov.f$ws))
cro.cov.f$dir <- as.numeric(as.character(cro.cov.f$dir))

ci.list <- lapply(1:nrow(cro.cov.f), function(x) {
  print(x)
  cov.sub.df <- cro.cov.f[x,]
  return(CIreal(m.cro,covs=cov.sub.df)$gamma)
})

cro.means.f <- lapply(ci.list, '[[', 1) 
cro.lb.f <- lapply(ci.list,'[[',3)
cro.ub.f <- lapply(ci.list,'[[',4)



# SG MALE #

sg.cov.m <- data.frame(dir=seq(from=min(m.sg$data$dir), max(m.sg$data$dir),by=1), lod = "L", 
                       sex = "M", ws = 9)
sg.cov.m$sex <- as.factor(as.character(sg.cov.m$sex))
sg.cov.m$lod <- as.factor(as.character(sg.cov.m$lod))
sg.cov.m$ws <- as.numeric(as.character(sg.cov.m$ws))
sg.cov.m$dir <- as.numeric(as.character(sg.cov.m$dir))

ci.list <- lapply(1:nrow(sg.cov.m), function(x) {
  print(x)
  cov.sub.df <- sg.cov.m[x,]
  return(CIreal(m.sg,covs=cov.sub.df)$gamma)
})

sg.means.m <- lapply(ci.list, '[[', 1) 
sg.lb.m <- lapply(ci.list,'[[',3)
sg.ub.m <- lapply(ci.list,'[[',4)


# SG FEMALE #

sg.cov.f <- data.frame(dir=seq(from=min(m.sg$data$dir), max(m.sg$data$dir),by=1), lod = "L", 
                       sex = "F", ws = 9)
ci.list <- lapply(1:nrow(sg.cov.f), function(x) {
  print(x)
  cov.sub.df <- sg.cov.f[x,]
  return(CIreal(m.sg,covs=cov.sub.df)$gamma)
})

sg.means.f <- lapply(ci.list, '[[', 1) 
sg.lb.f <- lapply(ci.list,'[[',3)
sg.ub.f <- lapply(ci.list,'[[',4)





#### 2B. PLOTTING ALL TRANSITIONS INDIVIDUALLY #####



#### 1 -> 1 #####

state1 = 1
state2 = 1
nb_states = 3

cro.1_1_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.1_1_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.1_1_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.1_1_means.m,lower_bound=cro.1_1_lb.m,upper_bound=cro.1_1_ub.m)


cro.1_1_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.1_1_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.1_1_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.1_1_means.f,lower_bound=cro.1_1_lb.f,upper_bound=cro.1_1_ub.f)


sg.1_1_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.1_1_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.1_1_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.1_1_means.m,lower_bound=sg.1_1_lb.m,upper_bound=sg.1_1_ub.m)


sg.1_1_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.1_1_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.1_1_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.1_1_means.f,lower_bound=sg.1_1_lb.f,upper_bound=sg.1_1_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0.5, 1)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 1 -> 2 #####

state1 = 1
state2 = 2
nb_states = 3

cro.1_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.1_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.1_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.1_2_means.m,lower_bound=cro.1_2_lb.m,upper_bound=cro.1_2_ub.m)


cro.1_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.1_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.1_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.1_2_means.f,lower_bound=cro.1_2_lb.f,upper_bound=cro.1_2_ub.f)


sg.1_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.1_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.1_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.1_2_means.m,lower_bound=sg.1_2_lb.m,upper_bound=sg.1_2_ub.m)


sg.1_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.1_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.1_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.1_2_means.f,lower_bound=sg.1_2_lb.f,upper_bound=sg.1_2_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 1 -> 3 #####


state1 = 1
state2 = 3
nb_states = 3

cro.1_3_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.1_3_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.1_3_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.1_3_means.m,lower_bound=cro.1_3_lb.m,upper_bound=cro.1_3_ub.m)


cro.1_3_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.1_3_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.1_3_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.1_3_means.f,lower_bound=cro.1_3_lb.f,upper_bound=cro.1_3_ub.f)


sg.1_3_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.1_3_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.1_3_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.1_3_means.m,lower_bound=sg.1_3_lb.m,upper_bound=sg.1_3_ub.m)


sg.1_3_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.1_3_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.1_3_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.1_3_means.f,lower_bound=sg.1_3_lb.f,upper_bound=sg.1_3_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  ylim(0, 0.5)+theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 2 -> 1 #####


state1 = 2
state2 = 1
nb_states = 3

cro.2_1_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.2_1_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.2_1_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.2_1_means.m,lower_bound=cro.2_1_lb.m,upper_bound=cro.2_1_ub.m)


cro.2_1_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.2_1_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.2_1_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.2_1_means.f,lower_bound=cro.2_1_lb.f,upper_bound=cro.2_1_ub.f)


sg.2_1_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.2_1_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.2_1_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.2_1_means.m,lower_bound=sg.2_1_lb.m,upper_bound=sg.2_1_ub.m)


sg.2_1_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.2_1_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.2_1_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.2_1_means.f,lower_bound=sg.2_1_lb.f,upper_bound=sg.2_1_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 2 -> 2 #####


state1 = 2
state2 = 2
nb_states = 3

cro.2_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.2_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.2_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.2_2_means.m,lower_bound=cro.2_2_lb.m,upper_bound=cro.2_2_ub.m)


cro.2_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.2_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.2_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.2_2_means.f,lower_bound=cro.2_2_lb.f,upper_bound=cro.2_2_ub.f)


sg.2_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.2_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.2_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.2_2_means.m,lower_bound=sg.2_2_lb.m,upper_bound=sg.2_2_ub.m)


sg.2_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.2_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.2_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.2_2_means.f,lower_bound=sg.2_2_lb.f,upper_bound=sg.2_2_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0.5, 1)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 2 -> 3 #####


state1 = 2
state2 = 3
nb_states = 3

cro.2_3_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.2_3_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.2_3_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.2_3_means.m,lower_bound=cro.2_3_lb.m,upper_bound=cro.2_3_ub.m)


cro.2_3_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1)) 
cro.2_3_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.2_3_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.2_3_means.f,lower_bound=cro.2_3_lb.f,upper_bound=cro.2_3_ub.f)


sg.2_3_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.2_3_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.2_3_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.2_3_means.m,lower_bound=sg.2_3_lb.m,upper_bound=sg.2_3_ub.m)


sg.2_3_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.2_3_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.2_3_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.2_3_means.f,lower_bound=sg.2_3_lb.f,upper_bound=sg.2_3_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  # ylim(0, 0.5)+
  theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))


#### 3 -> 1 #####


state1 = 3
state2 = 1
nb_states = 3

cro.3_1_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.3_1_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_1_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.3_1_means.m,lower_bound=cro.3_1_lb.m,upper_bound=cro.3_1_ub.m)


cro.3_1_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.3_1_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_1_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.3_1_means.f,lower_bound=cro.3_1_lb.f,upper_bound=cro.3_1_ub.f)


sg.3_1_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.3_1_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_1_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.3_1_means.m,lower_bound=sg.3_1_lb.m,upper_bound=sg.3_1_ub.m)


sg.3_1_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.3_1_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_1_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.3_1_means.f,lower_bound=sg.3_1_lb.f,upper_bound=sg.3_1_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  ylim(0, 0.5)+theme_bw() + #ylab("Transition probability")+
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 3 -> 2 #####


state1 = 3
state2 = 2
nb_states = 3

cro.3_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.3_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.3_2_means.m,lower_bound=cro.3_2_lb.m,upper_bound=cro.3_2_ub.m)


cro.3_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.3_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.3_2_means.f,lower_bound=cro.3_2_lb.f,upper_bound=cro.3_2_ub.f)


sg.3_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.3_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.3_2_means.m,lower_bound=sg.3_2_lb.m,upper_bound=sg.3_2_ub.m)


sg.3_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1))
sg.3_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir,mean=sg.3_2_means.f,lower_bound=sg.3_2_lb.f,upper_bound=sg.3_2_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  theme_bw() + #ylab("Transition probability")+
  # scale_y_continuous(limits=c(0, 0.63), breaks = seq(0, 0.6, 0.1)) + 
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))




#### 3 -> 3 #####


state1 = 3
state2 = 3
nb_states = 3

cro.3_3_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1)) 
cro.3_3_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_3_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       dir=cro.cov.m$dir,mean=cro.3_3_means.m,lower_bound=cro.3_3_lb.m,upper_bound=cro.3_3_ub.m)


cro.3_3_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.3_3_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_3_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       dir=cro.cov.f$dir,mean=cro.3_3_means.f,lower_bound=cro.3_3_lb.f,upper_bound=cro.3_3_ub.f)


sg.3_3_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1))
sg.3_3_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_3_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      dir=sg.cov.m$dir,mean=sg.3_3_means.m,lower_bound=sg.3_3_lb.m,upper_bound=sg.3_3_ub.m)


sg.3_3_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.3_3_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_3_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      dir=sg.cov.f$dir, mean=sg.3_3_means.f,lower_bound=sg.3_3_lb.f,upper_bound=sg.3_3_ub.f)

all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))

# all together
ggplot(all.df, aes(x = dir, y=mean)) + facet_wrap(~site)+
  geom_ribbon(size = 1.3, aes(ymin=lower_bound, ymax=upper_bound, fill = sex), alpha=0.15) + 
  geom_line(size = 1.3, aes(linetype = sex)) + 
  theme_bw() + #ylab("Transition probability")+
  # scale_y_continuous(limits=c(0.35, 1), breaks = seq(0.4, 1, 0.1)) + 
  scale_x_continuous(limits=c(0, 180), breaks = seq(0, 180, 60)) + #xlab("Relative wind direction (degrees)")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.text.x = element_text(size = 18),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_fill_manual(values = c("grey30", "grey50", "grey50", "grey30", "grey30", "grey50"))





##########  3. PLOTTING PREDICTED TRANSITION PROBABILITIES FROM REST TO SEARCH (TAKE OFF) AT WIND SPEED INTERVALS - FOR TABLE 3 ######

########## RJ: PLOTTING? THERE'S NO PLOT

#### 3A. CREATE DATAFRAMES OF PREDICTED VALUES  - same as step 1A. #####


# CROZET MALE #

cro.cov.m <- data.frame(ws=seq(from=min(m.cro$data$ws), max(m.cro$data$ws),by=0.2), lod = "L", 
                        sex = "M", dir = 75)
cro.cov.m$sex <- as.factor(as.character(cro.cov.m$sex))
cro.cov.m$lod <- as.factor(as.character(cro.cov.m$lod))
cro.cov.m$ws <- as.numeric(as.character(cro.cov.m$ws))
cro.cov.m$dir <- as.numeric(as.character(cro.cov.m$dir))

ci.list <- lapply(1:nrow(cro.cov.m), function(x) {
  print(x)
  cov.sub.df <- cro.cov.m[x,]
  return(CIreal(m.cro,covs=cov.sub.df)$gamma)
})
# extract means, upper and lower bounds from predicted outputs
cro.means.m <- lapply(ci.list, '[[', 1) 
cro.lb.m <- lapply(ci.list,'[[',3)
cro.ub.m <- lapply(ci.list,'[[',4)


# CRO FEMALE #

cro.cov.f <- data.frame(ws=seq(from=min(m.cro$data$ws), max(m.cro$data$ws),by=0.2), lod = "L", 
                        sex = "F", dir = 75)
cro.cov.f$sex <- as.factor(as.character(cro.cov.f$sex))
cro.cov.f$lod <- as.factor(as.character(cro.cov.f$lod))
cro.cov.f$ws <- as.numeric(as.character(cro.cov.f$ws))
cro.cov.f$dir <- as.numeric(as.character(cro.cov.f$dir))

ci.list <- lapply(1:nrow(cro.cov.f), function(x) {
  print(x)
  cov.sub.df <- cro.cov.f[x,]
  return(CIreal(m.cro,covs=cov.sub.df)$gamma)
})

cro.means.f <- lapply(ci.list, '[[', 1) 
cro.lb.f <- lapply(ci.list,'[[',3)
cro.ub.f <- lapply(ci.list,'[[',4)


# SG MALE #

sg.cov.m <- data.frame(sex = "M", ws=seq(from=min(m.sg$data$ws), max(m.sg$data$ws),by=0.2), dir = 75,
                       lod = "L")
sg.cov.m$sex <- as.factor(as.character(sg.cov.m$sex))
sg.cov.m$lod <- as.factor(as.character(sg.cov.m$lod))
sg.cov.m$ws <- as.numeric(as.character(sg.cov.m$ws))
sg.cov.m$dir <- as.numeric(as.character(sg.cov.m$dir))

ci.list <- lapply(1:nrow(sg.cov.m), function(x) {
  print(x)
  cov.sub.df <- sg.cov.m[x,]
  return(CIreal(m.sg,covs=cov.sub.df)$gamma)
})

sg.means.m <- lapply(ci.list, '[[', 1) 
sg.lb.m <- lapply(ci.list,'[[',3)
sg.ub.m <- lapply(ci.list,'[[',4)


# SG FEMALE #

sg.cov.f <- data.frame(sex = "F", ws=seq(from=min(m.sg$data$ws), max(m.sg$data$ws),by=0.2), dir = 75,
                       lod = "L")
sg.cov.f$sex <- as.factor(as.character(sg.cov.f$sex))
sg.cov.f$lod <- as.factor(as.character(sg.cov.f$lod))
sg.cov.f$ws <- as.numeric(as.character(sg.cov.f$ws))
sg.cov.f$dir <- as.numeric(as.character(sg.cov.f$dir))

ci.list <- lapply(1:nrow(sg.cov.f), function(x) {
  print(x)
  cov.sub.df <- sg.cov.f[x,]
  return(CIreal(m.sg,covs=cov.sub.df)$gamma)
})

sg.means.f <- lapply(ci.list, '[[', 1) 
sg.lb.f <- lapply(ci.list,'[[',3)
sg.ub.f <- lapply(ci.list,'[[',4)



#### 3B. EXTRACTING PROBABILITIES FOR TRANSITION FROM REST TO SEARCH #####


state1 = 3
state2 = 2
nb_states = 3

cro.3_2_means.m <- unlist(lapply(cro.means.m,'[[',nb_states*(state2-1)+state1))
cro.3_2_lb.m <- unlist(lapply(cro.lb.m,'[[',nb_states*(state2-1)+state1))
cro.3_2_ub.m <- unlist(lapply(cro.ub.m,'[[',nb_states*(state2-1)+state1))
cro.df.m <- data.frame(sex = "M", site = "Crozet", 
                       wind=cro.cov.m$ws,mean=cro.3_2_means.m,lower_bound=cro.3_2_lb.m,upper_bound=cro.3_2_ub.m)
# summarizing into wind speed categories
cro.df.m$wind_cat <- 1
cro.df.m$wind_cat[cro.df.m$wind > 5 & cro.df.m$wind <= 15] <- 2
cro.df.m$wind_cat[cro.df.m$wind > 15] <- 3


cro.3_2_means.f <- unlist(lapply(cro.means.f,'[[',nb_states*(state2-1)+state1))
cro.3_2_lb.f <- unlist(lapply(cro.lb.f,'[[',nb_states*(state2-1)+state1))
cro.3_2_ub.f <- unlist(lapply(cro.ub.f,'[[',nb_states*(state2-1)+state1))
cro.df.f <- data.frame(sex = "F", site = "Crozet", 
                       wind=cro.cov.f$ws,mean=cro.3_2_means.f,lower_bound=cro.3_2_lb.f,upper_bound=cro.3_2_ub.f)
# summarizing into wind speed categories
cro.df.f$wind_cat <- 1
cro.df.f$wind_cat[cro.df.f$wind > 5 & cro.df.f$wind <= 15] <- 2
cro.df.f$wind_cat[cro.df.f$wind > 15] <- 3


sg.3_2_means.m <- unlist(lapply(sg.means.m,'[[',nb_states*(state2-1)+state1)) 
sg.3_2_lb.m <- unlist(lapply(sg.lb.m,'[[',nb_states*(state2-1)+state1))
sg.3_2_ub.m <- unlist(lapply(sg.ub.m,'[[',nb_states*(state2-1)+state1))
sg.df.m <- data.frame(sex = "M", site = "South Georgia", 
                      wind=sg.cov.m$ws,mean=sg.3_2_means.m,lower_bound=sg.3_2_lb.m,upper_bound=sg.3_2_ub.m)
# summarizing into wind speed categories
sg.df.m$wind_cat <- 1
sg.df.m$wind_cat[sg.df.m$wind > 5 & sg.df.m$wind <= 15] <- 2
sg.df.m$wind_cat[sg.df.m$wind > 15] <- 3


sg.3_2_means.f <- unlist(lapply(sg.means.f,'[[',nb_states*(state2-1)+state1)) 
sg.3_2_lb.f <- unlist(lapply(sg.lb.f,'[[',nb_states*(state2-1)+state1))
sg.3_2_ub.f <- unlist(lapply(sg.ub.f,'[[',nb_states*(state2-1)+state1))
sg.df.f <- data.frame(sex = "F", site = "South Georgia", 
                      wind=sg.cov.f$ws,mean=sg.3_2_means.f,lower_bound=sg.3_2_lb.f,upper_bound=sg.3_2_ub.f)
# summarizing into wind speed categories
sg.df.f$wind_cat <- 1
sg.df.f$wind_cat[sg.df.f$wind > 5 & sg.df.f$wind <= 15] <- 2
sg.df.f$wind_cat[sg.df.f$wind > 15] <- 3


all.df <- rbind(cro.df.m, cro.df.f, sg.df.m, sg.df.f)
all.df$sex <- as.factor(as.character(all.df$sex))
all.df$wind_cat <- as.factor(as.character(all.df$wind_cat))


#### 3C. SUMMARIZING PROBABILITIES FOR EACH WIND SPEED CATEGORY #####

wind_sum <- ddply(all.df,.(sex, site, wind_cat), summarize, w_mean = mean(mean), w_low = mean(lower_bound), w_upp = mean(upper_bound))
wind_sum
# We should have similar results
#   sex          site wind_cat     w_mean      w_low     w_upp
# 1    F        Crozet        1 0.09512124 0.05853763 0.1514885
# 2    F        Crozet        2 0.06858004 0.04511832 0.1042383
# 3    F        Crozet        3 0.05138183 0.02321813 0.1100571
# 4    F South Georgia        1 0.19577900 0.12125148 0.3012170
# 5    F South Georgia        2 0.20551201 0.12925069 0.3122516
# 6    F South Georgia        3 0.21364747 0.09878546 0.4026569
# 7    M        Crozet        1 0.14922044 0.08291460 0.2549658
# 8    M        Crozet        2 0.45953293 0.22189540 0.6809509
# 9    M        Crozet        3 0.75519568 0.32045571 0.9522864
# 10   M South Georgia        1 0.16367153 0.08747671 0.2880735
# 11   M South Georgia        2 0.34839603 0.21928099 0.4979386
# 12   M South Georgia        3 0.54017690 0.29296831 0.7682557


