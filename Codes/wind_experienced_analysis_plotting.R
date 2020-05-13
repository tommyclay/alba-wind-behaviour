
library(lme4)
library(MuMIn)
library(ggplot2)
library(visreg)
library(plyr)
library(circular)


# AIM OF SCRIPT - to load in wind values associated with locations for both populations and model how wind speeds and directions encountered 
# vary according to population and sex, and create Figure 2.



#### 1. LOAD IN DATA #####


path <- "./Data_inputs/BothSites_forwindexperienced.csv"
dat <- read.csv(file = path, na.strings = "NA")


#### 2A. MODELLING WIND SPEEDS ENCOUNTERED IN RELATION TO SITE AND SEX #####


hist(dat$WindSp)
# looks great, except for a bit of skewness on the right, Gaussian distribution could be reasonable

# transforming into factor variables
dat$Year <- as.factor(as.character(dat$Year))
table(dat$Year)
dat$ID <- as.factor(as.character(dat$ID))
dat$TripID <- as.factor(as.character(dat$TripID))
dat$Site <- as.factor(as.character(dat$Site))
dat$Sex <- as.factor(as.character(dat$Sex))

# check for na values
summary(dat)

wind.lm <- lmer(WindSp ~ Site*Sex + Year + (1|ID/TripID), data = dat)
wind.lm2 <- update(wind.lm, REML = FALSE)
options(na.action = "na.fail") # without this we get error "global model na action arguement not set" - this function returns the object only if it 
# contains no missing values. Setting na.action = "na.pass" doesn't change the results at all. 
m_set <- dredge(wind.lm2)
m_set
# Model selection table 
#     (Int) Sex Sit Yer Sex:Sit df   logLik    AICc delta weight
# 16 8.364   +   +   +       + 13 -1018031 2036088  0.00   0.91
# 8  8.413   +   +   +         12 -1018035 2036093  4.62   0.09
# 6  8.421   +       +         11 -1018044 2036109 20.74   0.00
# 7  8.702       +   +         11 -1018047 2036116 27.78   0.00
# 5  8.702           +         10 -1018055 2036130 41.48   0.00
# 12 8.474   +   +           +  7 -1018061 2036136 47.95   0.00
# all variables important

# best model is full model
visreg(wind.lm)
# from looking, values of wind speed for crozet higher than south georgia
# males experience faster wind speeds than females

# plotting residuals
plot(wind.lm)
plot(residuals(wind.lm, type = "pearson") ~ dat$Site)
plot(residuals(wind.lm, type = "pearson") ~ dat$Sex)
plot(residuals(wind.lm, type = "pearson") ~ dat$Year)


summary(wind.lm)
# Fixed effects:
#                         Estimate Std. Error t value
# (Intercept)              8.3655     0.4116  20.322
# SiteSouth Georgia       -0.8036     0.4259  -1.887
# SexM                     0.9017     0.1535   5.873
# Year2011                -1.0823     0.4551  -2.378
# Year2012                 1.3849     0.4834   2.865
# Year2013                 0.2590     0.4267   0.607
# Year2014                 0.2139     0.5103   0.419
# Year2015                -0.4569     0.5407  -0.845
# Year2016                 0.1356     0.4309   0.315
# SiteSouth Georgia:SexM  -1.1939     0.4684  -2.549




####### 2B. PLOTTING FOR FIGURE 2  ######



# plot modelled mean +- 95% CI for each category

pred.df <- expand.grid(Sex = levels(dat$Sex), Site = levels(dat$Site), ID = levels(dat$ID),
                       TripID = levels(dat$TripID), Year = levels(dat$Year)) 
pred.df2 <- predict(wind.lm, newdata = pred.df, re.form = NA, se.fit = TRUE)
pred.df$fit <- pred.df2$fit
pred.df$se.lwr <- pred.df$fit -pred.df2$se.fit
pred.df$se.upr <- pred.df$fit +pred.df2$se.fit
pred.df$ci.lwr <- pred.df$fit - pred.df2$se.fit*1.96
pred.df$ci.upr <- pred.df$fit + pred.df2$se.fit*1.96

# average for all years
pred.avg <- ddply(pred.df,.(Sex, Site), summarize, ci.upr = mean(ci.upr), ci.lwr = mean(ci.lwr), se.upr = mean(se.upr), se.lwr = mean(se.lwr),
                  fit = mean(fit))

pd <- position_dodge(0.9)
p <- ggplot(pred.avg, aes(x=Site, y=fit)) 
p + geom_errorbar(aes(ymin=ci.lwr, ymax=ci.upr, fill = Sex), width=.45, position=pd) +
  geom_point(size = 4,shape = 21, aes(fill = Sex), position=pd) + theme_bw() +
  ylab(expression(paste("Wind speed (", ms^-1, ")", sep="")))+xlab("Population")+
  scale_fill_manual(values=c("gray", "black")) +
  scale_colour_manual(values=c("gray50", "gray50")) +
  theme(axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18), 
        axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        strip.text.x = element_text(size=18),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        legend.title = element_text(size=18),
        legend.text  = element_text(size=16),
        legend.position="top", 
        legend.margin=margin(t = -0.2, unit='cm'), plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(colour=FALSE, fill=guide_legend(override.aes=list(size=4)))




####### 3. MODELLING WIND DIRECTIONS ENCOUNTERED IN RELATION TO SITE AND SEX   #########




#### 3A. SUMMARIZING MEDIAN WIND PER TRIP #####

dat$WindDir2 <- circular(dat$WindDir,
                         type = "angles",
                         units = "degrees", 
                         zero = 0, 
                         rotation = "counter")
trips.sum <- ddply(dat,.(Site, Sex, TripID), summarize, MedDir = median.circular(WindDir2))
hist(as.numeric(trips.sum$MedDir))
# checking that values are within -180 and 180
range(as.numeric(trips.sum$MedDir))


# wind direction for all individuals - mean and confidence intervals
mle.vonmises(trips.sum$MedDir)
# mu: 100.6  ( 2.197 )
# kappa: 2.359  ( 0.1459 )
mle.vonmises.bootstrap.ci(trips.sum$MedDir, alpha = 0.05)
# Mean Direction:            Low = 97.15   High = 103.58 
# Concentration Parameter:   Low = 2.01   High = 2.84 

# function to convert to where wind is coming FROM on 0-360 scale
dir_conv <- function(x) { if (x > 0) { dir <- 360-(180-x)} else { dir <- 180+x }
  return(dir)
}
dir_conv(100.6) # 280.6
dir_conv(103.58) # 283.58
dir_conv(97.15) # 277.15



##### 3B. COMPARING BY SITE ######


cro.df <- subset(trips.sum, Site == "Crozet",)
cro.dir <- cro.df$MedDir
sg.df <- subset(trips.sum, Site == "South Georgia",)
sg.dir <- sg.df$MedDir

# Watson test
watson.two.test(cro.dir,sg.dir, alpha=0)
# Test Statistic: 0.081 
# P-value > 0.10 
watson.two.test(cro.dir,sg.dir, alpha=0.1)
# Test Statistic: 0.081 
# Level 0.1 Critical Value: 0.152 
# Do Not Reject Null Hypothesis 


# getting mean confidence intervals around values for each site
mle.vonmises(cro.dir)
# mu: 100.4  ( 2.416 )
# kappa: 2.218  ( 0.145 )
mle.vonmises.bootstrap.ci(cro.dir, alpha = 0.05)
# Bootstrap Confidence Intervals for Mean Direction and Concentration 
# Confidence Level:   95 % 
# Mean Direction:            Low = 96.26   High = 104.08 
# Concentration Parameter:   Low = 1.87   High = 2.7 

# convert to where wind is coming FROM on 0-360 scale
dir_conv(100.4) # 280.4
dir_conv(104.08 ) # 284.08
dir_conv(96.26) # 276.26


mle.vonmises(sg.dir)
# mu: 102.4  ( 3.901 )
# kappa: 6.202  ( 1.349 )
mle.vonmises.bootstrap.ci(sg.dir, alpha = 0.05)
# Bootstrap Confidence Intervals for Mean Direction and Concentration 
# Confidence Level:   95 % 
# Mean Direction:            Low = 94.84   High = 110.48 
# Concentration Parameter:   Low = 4.04   High = 11.48

# to convert to where wind is coming FROM on 0-360 scale
dir_conv(102.4) # 282.4
dir_conv(110.48) # 290.48
dir_conv(94.84) # 274.84




##### 3C. CROZET BY SEX #####


cro.m.df <- subset(cro.df, Sex == "M",)
cro.m.dir <- cro.m.df$MedDir
cro.f.df <- subset(cro.df, Sex == "F",)
cro.f.dir <- cro.f.df$MedDir


watson.two.test(cro.m.dir,cro.f.dir, alpha=0)
# Test Statistic: 0.3466
# 0.001 < P-value < 0.01 
watson.two.test(cro.m.dir,cro.f.dir, alpha=0.01)
# Test Statistic: 0.3463 
# Level 0.01 Critical Value: 0.268 
# Reject Null Hypothesis 
watson.two.test(cro.m.dir,cro.f.dir, alpha=0.001)
# Test Statistic: 0.34636
# Level 0.001 Critical Value: 0.385 
# Do Not Reject Null Hypothesis 

# p value between 0.001 and 0.01


mle.vonmises(cro.m.dir)
# mu: 101  ( 2.63 )
# kappa: 3.258  ( 0.304 )
mle.vonmises.bootstrap.ci(cro.m.dir, alpha = 0.05)
# Bootstrap Confidence Intervals for Mean Direction and Concentration 
# Confidence Level:   95 % 
# Mean Direction:            Low = 96.57   High = 105.12 
# Concentration Parameter:   Low = 2.55   High = 4.42 

# to convert to where wind is coming FROM on 0-360 scale
dir_conv(101) # 281
dir_conv(105.12) # 285.12
dir_conv(96.57) # 276.57


mle.vonmises(cro.f.dir)
# mu: 99.52  ( 4.302 )
# kappa: 1.646  ( 0.1631 )
mle.vonmises.bootstrap.ci(cro.f.dir, alpha = 0.05)
# Bootstrap Confidence Intervals for Mean Direction and Concentration 
# Confidence Level:   95 % 
# Mean Direction:            Low = 93.13   High = 105.99 
# Concentration Parameter:   Low = 1.29   High = 2.1 

# to convert to where wind is coming FROM on 0-360 scale
dir_conv(99.52) # 279.52
dir_conv(105.99) # 285.99
dir_conv(93.13) # 273.13





##### 3D. SOUTH GEORGIA BY SEX #####



sg.m.df <- subset(sg.df, Sex == "M",)
sg.m.dir <- sg.m.df$MedDir
sg.f.df <- subset(sg.df, Sex == "F",)
sg.f.dir <- sg.f.df$MedDir

watson.two.test(sg.m.dir,sg.f.dir, alpha=0)
# Test Statistic: 0.1222 
# P-value > 0.10 
watson.two.test(sg.m.dir,sg.f.dir, alpha=0.1)
# Test Statistic: 0.1222 
# Level 0.1 Critical Value: 0.152 
# Do Not Reject Null Hypothesis 


mle.vonmises(sg.m.dir)
# mu: 99.94  ( 4.872 )
# kappa: 7.106  ( 2.098 )
mle.vonmises.bootstrap.ci(sg.m.dir, alpha = 0.05)
# Bootstrap Confidence Intervals for Mean Direction and Concentration 
# Confidence Level:   95 % 
# Mean Direction:            Low = 91.56   High = 109.47 
# Concentration Parameter:   Low = 4.04   High = 15 

# to convert to where wind is coming FROM on 0-360 scale
dir_conv(99.94) # 279.94
dir_conv(109.47) # 289.47
dir_conv(91.56) # 271.56


mle.vonmises(sg.f.dir)
#mu: 105.4  ( 6.237 )
#kappa: 5.493  ( 1.767 )
mle.vonmises.bootstrap.ci(sg.f.dir, alpha = 0.05)
# Bootstrap Confidence Intervals for Mean Direction and Concentration 
# Confidence Level:   95 % 
# Mean Direction:            Low = 91.11   High = 115.96 
# Concentration Parameter:   Low = 2.83   High = 15.54 

# to convert to where wind is coming FROM on 0-360 scale
dir_conv(105.4) # 285.4
dir_conv(115.96) # 295.96
dir_conv(91.11) # 271.11

