###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long September 2021

### Differences since June 2021 version:
##### Cubic fits removed
##### Total GM volume and ICV analyses in addition to region by region

### Differences since Sept 2021 version
#### Sex and MatEdu fits after cubic removed
#### LOESS fits for later added

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, nlme, tidyverse, methods, lmtest, gridExtra, grid, ggplot2, sjstats, patchwork)
install.packages("devtools")
devtools::install_github("cardiomoon/ggiraphExtra")
devtools::install_github("hoxo-m/magicfor")
library(magicfor)

install.packages('lmerTest')
library(lmerTest)
install.packages('pbkrtest')
library(pbkrtest)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Preschool_Volume/")

PSdata <- read.csv("Preschool_Volume_Longitudinal_July2022.csv")

mc_age <- scale(PSdata$Age, center=TRUE, scale=TRUE)
#recoding variables into differenty data types, for example string into numeric
PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$Age <- as.numeric(PSdata$Age)
PSdata$Subject <- as.factor(PSdata$Subject)
PSdata$SessionID_Split1 <- as.numeric(PSdata$SessionID_Split1)
PSdata$Sex.F1_M2. <- as.factor(PSdata$Sex.F1_M2.)
PSdata$OrderedByAgeAtSsn1 <- as.factor(PSdata$OrderedByAgeAtSsn1)
PSdata$income_3 <- as.factor(PSdata$income_3)
PSdata$schooling_ord <- as.factor(PSdata$schooling_ord)
PSdata$schooling_ord <- ordered(as.factor(PSdata$schooling_ord), levels=c("1","2","3","4"))
PSdata$schooling_ord <- as.numeric(PSdata$schooling_ord)

eta_sq(lmer_volume_quad)

##############################################################################################################
###################### trying various things related to getting peak age ############################################
##############################################################################################################

#predicted values versus raw volume values
#curve looks the same and peak is incorrect on both
lmer_volume_quad <- lmer(Total_Gray ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
eta_sq(lmer_volume_quad)

PSdata$lmerpred<-predict(lmer_volume_quad)

plotROI <- ggplot(PSdata, aes(x=Age, y=Left_TMP_temporal_pole))
plotROI3 <- plotROI + geom_point()+geom_smooth(aes(y=Left_TMP_temporal_pole), method='lm', se=FALSE, formula=y~x+I(x^2)) + scale_color_grey() +
  theme_classic() +
  theme (plot.title = element_text(size=40, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
  theme(legend.position = "none") +
  xlab("Age (Years)") + ylab(""~mm^3*"") +geom_vline(xintercept=7.23977) + scale_y_continuous(breaks=seq(2000, 15000, 2500)) + coord_cartesian(ylim=c(2000, 15000))+geom_vline(xintercept=7.2)
  

plotROI <- ggplot(PSdata, aes(x=Age, y=lmerpred))
plotROI4 <- plotROI + geom_point()+geom_smooth(aes(y=lmerpred), method='lm', se=FALSE, formula=y~x+I(x^2)) + scale_color_grey() +
  theme_classic() +
  theme (plot.title = element_text(size=40, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
  theme(legend.position = "none") +
  xlab("Age (Years)") + ylab(""~mm^3*"") +geom_vline(xintercept=7.23977) + scale_y_continuous(breaks=seq(2000, 15000, 2500)) + coord_cartesian(ylim=c(2000, 15000))+geom_vline(xintercept=7.2)


#plotROI4<- plotROI3 + geom_smooth(aes(x=Age, y=Left_TMP_temporal_pole), se=FALSE, method =lm, formula=y~x+I(x^2))

plotROI1 + 
plotROI + scale_y_continuous(breaks=seq(2000, 15000, 2500)) + coord_cartesian(ylim=c(2000, 15000))+geom_vline(xintercept=7.2)

#actual formula for the L_TMP trajectory, vertex is 7.2ish
curve(-77.53*x^2+1122.6*x+4623.86, from=2, to=9)
curve(-77.53*x^2+1500.6*x+4623.86, from=2, to=9)



curve(157723*log(x)-20585*x+581581, from=2, to=9)
curve(100000*log(x)-20585*x+581581, from=2, to=9)


curve(15.16*log(x)+742.41, from=2, to=9)

### July22 sex abs significant interactions
#Left MPoG male
curve(134.59+425.85*x-39.96*x^2, col=4, from=2, to=8)
#Left MPoG female
curve(1117.66+47.48*x-5.79*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=c(5.33, 4.09), lty=c(1,2))

#Left FO male
curve(2280.23+126.52*x, col=4, from=2, to=8)
#Left FO female
curve(2385.64+80.27*x, col=2, from=2, to=8, add=TRUE)

#Left basal forebrain male
curve(377.56+228.57*x-20.82*x^2, col=4, from=2, to=8)
#Left basal forebrain female
curve(537.12+123.4*x-11.27*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=c(5.49, 5.48), lty=c(1,2))

curve(1000*x+6000, from=2, to=8, col="white")
#Left sup temp gyrus male
curve(10979.68-131.27*x+19.59*x^2, col=4, from=2, to=8, add=TRUE)
#Left sup temporal gyrus female
curve(8023.81+613.29*x-44.67*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=c(3.35, 6.86), lty=c(1,2))

curve(.01*x+0.001*x^2, from=2, to=8, col="white")
#Left MPoG male
curve(0.0117+0.0293*x-0.0029*x^2, col=4, from=2, to=8, add=TRUE)
#Left MPoG female
curve(0.1056-0.0029*x-0.000017*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=5.06, lty=1)

#Right OFuG male
curve(0.4231+0.0054*x, col=4, from=2, to=8)
#Right OFuG female
curve(0.4676-0.0056*x, col=2, from=2, to=8, add=TRUE)


#Left STG male
curve(0.9597-0.0738*x+0.0059*x^2, col=4, from=2, to=8)
#Left STG female
curve(0.7504+0.0026*x-0.0007*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=6.26, lty=1)
abline(v=1.86, lty=2)

curve(.25+0.035*x, from=2, to=8, col="white")
#Left Putamen male
curve(0.3904-0.0051*x+0.0007*x^2, col=4, from=2, to=8, add=TRUE)
#Left Putamen female
curve(0.3662+0.0139*x-0.0013*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=3.47, lty=1)
abline(v=5.39, lty=2)

curve(.17+0.025*x, from=2, to=9, col="white")
#Left Caudate male
curve(0.3039-0.0124*x+0.0013*x^2, col=4, from=2, to=9, add=TRUE)
#Left Caudate female
curve(0.3109-0.0021*x+0.0001*x^2, col=2, from=2, to=9, add=TRUE)
curve(0.3-0*x, lty=2, from=2, to=9, add=TRUE)
abline(v=4.72, lty=1)
abline(v=7.52, lty=2)

curve(.52+0.035*x, from=2, to=9, col="white")
#Right FuG male
curve(0.7378-0.0353*x+0.0041*x^2, col=4, from=2, to=8, add=TRUE)
#Right FuG female
curve(0.6095+0.0139*x-0.0008*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=4.28, lty=1)
abline(v=9.01, lty=2)

curve(.3+0.015*x, from=2, to=8, col="white")
#Left AIns male
curve(0.4131-0.0087*x+0.0009*x^2, col=4, from=2, to=8, add=TRUE)
#Left AIns female
curve(0.3197+0.0186*x-0.0013*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=4.795, lty=1)
abline(v=7.05, lty=2)

curve(.3+0.015*x, from=2, to=9, col="white")
#Left CO male
curve(0.3668+0.0132*x-0.0017*x^2, col=4, from=2, to=8, add=TRUE)
#Left CO female
curve(0.4045-0.0064*x+0.0004*x^2, col=2, from=2, to=8, add=TRUE)
abline(v=3.98, lty=1)
abline(v=8.75, lty=2)


#exponential functions for fun
curve(2*exp(x/1), from=0, to=1*3)
curve(2*exp(x/3), from=0, to=3*3)
#log of those exponential functions
#y=m*log(x/k)
curve(1000*log(x/1), from=0, to=)
curve(2*log(x/1), from=0, to=31000)

#creating a function that we could then graph with geom_function()
LTMPfn <- function(x) {
  as.numeric(QuadCo*I(x^2) + LinCo*x + Interc) # outcome is y
}

lmer_volume_quad_MatEdu <- lmer(Right_Accumbens_Area ~ Age*schooling_ord + I(Age^2)*schooling_ord + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_quad_MatEdu)
emmeans(lmer_volume_quad_MatEdu, pairwise~schooling_ord)

lmer_volume_quad_Sex <- lmer(Right_Accumbens_Area ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_quad_Sex)
eta_sq(lmer_volume_quad_Sex)



#make a dataframe with the Intercept and coefficients from the model
CoefData <- as.data.frame(fixef(lmer_volume_quad))

#Access the 3 columns with the variables you want. These will get plugged into geom_function(fun=LTMPfn())
QuadCo <- as.numeric(CoefData[3,1])
LinCo <- as.numeric(CoefData[2,1])
Interc <- as.numeric(CoefData[1,1])

  
plotROI <- ggplot(PSdata, aes(x=Age, y=Right_Accumbens_Area))
plotROI3 <- plotROI + geom_point()+ geom_function(fun = LTMPfn, color="red") +
  theme(legend.position = "none") + geom_line()
  xlab("Age (Years)") + ylab(""~mm^3*"") +geom_vline(xintercept=7.23977)+
  theme_classic() +
  theme (plot.title = element_text(size=40, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
  theme(legend.position = "none") 

plotROI4<- plotROI3 + geom_smooth(aes(x=Age, y=Right_Accumbens_Area), se=FALSE, method =lm, formula=y~x+I(x^2))



lmer_volume_quad <- lmer(Right_Accumbens_Area ~ mc_age + I(mc_age^2) + (1|Subject), REML = FALSE, data = PSdata)
lmer_volume_quad2 <- lmer(Left_Accumbens_Area ~ mc_age + I(mc_age^2) + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_quad)
df.coef <- as.data.frame(coef(summary(lmer_volume_quad)))
df.coef2 <- as.data.frame(coef(summary(lmer_volume_quad2)))

df.coef['Region'] <- c(colnames(regions)[i])

if (as.logical(row.names(df.etasq))=TRUE) {
  print("It's true")
}else {
  print("...")
}

plotfilename <- paste0("NoCov_", colnames(regions)[i])
table_coef <- table(df.coef)


#centering age
mc_age = scale(PSdata$Age, center = TRUE, scale = TRUE)

######### Results Log File
sink("PSVolumeLog_Raw_NoCoV_May2022", append=FALSE, split=TRUE)

######### All MaCRUISE regions, not normalized 
regions <-as.data.frame(PSdata[,c(14:17)])
dim(regions)
magic_for(silent = TRUE)

lmer_tableinit.sex <- lmer(Right_Inf_Lat_Vent ~ Age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef.sex <- as.data.frame(coef(summary(lmer_tableinit.sex)))
df.master.coef.sex["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq.sex <- as.data.frame(eta_sq(lmer_tableinit.sex))
df.master.etasq.sex["Region"] <-c("Right_Inf_Lat_Vent")

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  #plotfilename <- paste0("NoCov_", colnames(regions)[i])
  finalmodel <- 0
  bestAIC <- 0 
  finalplot <-0
  #null model
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(ROI ~ Age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + Sex.F1_M2.+ I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  AICdif1 <- AIC_null - AIC_lin
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    
    #final_etasq <- eta_sq(finalmodel)
    
    lmer_volume_lin_Sex <- lmer(ROI ~ Age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
    summary(lmer_volume_lin_Sex)
    final_etasq_Sex <- eta_sq(lmer_volume_lin_Sex)
    
    df.etasq.sex <- as.data.frame(final_etasq_Sex)
    df.etasq.sex['Region'] <- c(colnames(regions)[i])
    df.master.etasq.sex <- rbind(df.master.etasq.sex, df.etasq.sex)
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI1 <- plotROI + geom_line(aes(y=predict(lmer_volume_lin), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      #geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      #theme_classic() +
      #theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      #theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      #theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      #theme(legend.position = "none") +
      #coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
   # finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
   # plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
   # plotROI2 <- plotROI + geom_line(aes(y=predict(lmer_volume_null), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    #  geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
     # theme_classic() +
    #  theme (plot.title = element_text(size=40, color="black",face="bold")) + 
     # theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
     # theme(legend.position = "none") +
    #  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
   # finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    #plotROI <- ggplot(PSdata, aes(x=Age, y=lmerpred))
    #plotROI3 <- plotROI + geom_point()+geom_smooth(aes(y=lmerpred), method='lm', se=FALSE, formula=y~x+I(x^2)) + scale_color_grey() +
      #geom_line(aes(y=lmerpred))+
      #theme_classic() +
      #theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      #theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      #theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      #theme(legend.position = "none") +
      #xlab("Age (Years)") + ylab(""~mm^3*"") +geom_vline(xintercept=7.23977)
    #plotROI4<- plotROI3 + geom_smooth(aes(x=Age, y=Left_TMP_temporal_pole), se=FALSE, method =lm, formula=y~x+I(x^2))
    
   # finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
 # print(finalmodel)
  write.csv(df.master.etasq.sex, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_NoCov_etasq_Sex.csv")



  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad)  
  #print(paste("AIC_null =", AIC_null))
  #print(paste("AIC_lin =", AIC_lin))
  #print(paste("AIC_quad =", AIC_quad))
  #final_etasq <- eta_sq(finalmodel)
  #print(finalmodel)
  #print(final_etasq)
  #print(coef(summary(finalmodel)))
  #print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  
}


Raw_NoCov <- magic_result_as_dataframe()
write.csv(Raw_NoCov, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Raw_NoCov_AIC.csv")
magic_free()
sink()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sink("PSVolumeLog_Norm_NoCoV_May2022", append=FALSE, split=TRUE)
#magic_for(silent = TRUE)
########## All MaCRUISE Regions, normalized by ICV
regions <-as.data.frame(PSdata[,c(149:281)])
dim(regions)

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  
  #plotfilename <- paste0("NoCov_", colnames(regions)[i])
  finalmodel <- 0
  bestAIC <- 0 
  finalplot <-0
  #null model
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  #linear model
  lmer_volume_lin <- lmer(ROI ~ mc_age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ mc_age + I(mc_age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  AICdif1 <- AIC_null - AIC_lin
  
  print(colnames(regions)[i]) 
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI1 <- plotROI + geom_line(aes(y=predict(lmer_volume_lin), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    # geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    # theme_classic() +
    # theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    # theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    # theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    # theme(legend.position = "none") +
    # coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    # finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    # plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    # plotROI2 <- plotROI + geom_line(aes(y=predict(lmer_volume_null), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    #  geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    # theme_classic() +
    #  theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    # theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    # theme(legend.position = "none") +
    #  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    # finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    # plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    # plotROI3 <- plotROI + geom_line(aes(y=predict(lmer_volume_quad), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    #  geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    # theme_classic() +
    #  theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    # theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    # theme(legend.position = "none") +
    #  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    # finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  print(finalmodel)
  #print(AIC_null)
  #print(AIC_lin)
  #print(AIC_quad)
  
  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad)
  print(paste("AIC_null =", AIC_null))
  print(paste("AIC_lin =", AIC_lin))
  print(paste("AIC_quad =", AIC_quad))

  
  final_etasq <- eta_sq(finalmodel)
  print(final_etasq)
  print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}
#Norm_NoCov <- magic_result_as_dataframe()
#write.csv(Norm_NoCov, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Norm_NoCov_AIC.csv")
#magic_free()
sink()



sink("PSVolumeLog_Raw_Sex_Oct2021", append=FALSE, split=TRUE)
  
############################# Sex Raw Lin Regions ################################
lmer_volume_lin <- lmer(Right_Pallidum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_Pallidum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_Ventral_DC ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_ACgG_anterior_cingulate_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_AOrG_anterior_orbital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_CO_central_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_CO_central_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_FO_frontal_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_FO_frontal_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_GRe_gyrus_rectus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_GRe_gyrus_rectus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_LiG_lingual_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_LOrG_lateral_orbital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MCgG_middle_cingulate_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MOrG_medial_orbital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_MOrG_medial_orbital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MPoG_postcentral_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_MPoG_postcentral_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MSFG_superior_frontal_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_OCP_occipital_pole ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_OpIFG_opercular_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_OrIFG_orbital_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_OrIFG_orbital_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_PHG_parahippocampal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_PO_parietal_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_PO_parietal_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_PrG_precentral_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_PT_planum_temporale ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_PT_planum_temporale ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_SCA_subcallosal_area ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SCA_subcallosal_area ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_SFG_superior_frontal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SFG_superior_frontal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_SMC_supplementary_motor_cortex ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SMC_supplementary_motor_cortex ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SMG_supramarginal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_TrIFG_triangular_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_TTG_transverse_temporal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

  
############################# Sex Raw Quad Regions ################################
lmer_volume_quad <- lmer(Total_Gray ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Accumbens_Area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Accumbens_Area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Amygdala ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Amygdala ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Caudate ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Caudate ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Hippocampus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Hippocampus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Putamen ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Putamen ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Thalamus_Proper ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Thalamus_Proper ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Ventral_DC ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Basal_Forebrain ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_ACgG_anterior_cingulate_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_AIns_anterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_AIns_anterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_AnG_angular_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_AnG_angular_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Calc_calcarine_cortex ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Calc_calcarine_cortex ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Ent_entorhinal_area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Ent_entorhinal_area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_FRP_frontal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_FRP_frontal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_FuG_fusiform_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_FuG_fusiform_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_IOG_inferior_occipital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_IOG_inferior_occipital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_ITG_inferior_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_ITG_inferior_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_LiG_lingual_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MFC_medial_frontal_cortex ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MFC_medial_frontal_cortex ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MFG_middle_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MFG_middle_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MOG_middle_occipital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MSFG_superior_frontal_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MTG_middle_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MTG_middle_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_OpIFG_opercular_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PCgG_posterior_cingulate_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PCgG_posterior_cingulate_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PCu_precuneus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PCu_precuneus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PHG_parahippocampal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PIns_posterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PIns_posterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_POrG_posterior_orbital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_POrG_posterior_orbital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PP_planum_polare ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PP_planum_polare ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_SMG_supramarginal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_SOG_superior_occipital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_STG_superior_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_STG_superior_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_TMP_temporal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_TMP_temporal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_TrIFG_triangular_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_agex^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)
  

sink()

sink("PSVolumeLog_Norm_Sex_Oct2021", append=FALSE, split=TRUE)

############################# Sex Norm Lin Regions ################################
lmer_volume_lin <- lmer(ICV_Norm_Right_Caudate ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Caudate ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_Ventral_DC ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Ventral_DC ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Basal_Forebrain ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_Basal_Forebrain ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_ACgG_anterior_cingulate_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_AOrG_anterior_orbital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_CO_central_operculum ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_Cun_cuneus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Cun_cuneus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_IOG_inferior_occipital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_IOG_inferior_occipital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_LOrG_lateral_orbital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_MCgG_middle_cingulate_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_MCgG_middle_cingulate_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_MOG_middle_occipital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_MOG_middle_occipital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_OCP_occipital_pole ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_OFuG_occipital_fusiform_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_OFuG_occipital_fusiform_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PCgG_posterior_cingulate_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_PCu_precuneus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PCu_precuneus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PoG_postcentral_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_PrG_precentral_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PrG_precentral_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_SMC_supplementary_motor_cortex ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_SOG_superior_occipital_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_SPL_superior_parietal_lobule ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_SPL_superior_parietal_lobule ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_STG_superior_temporal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_TTG_transverse_temporal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_TTG_transverse_temporal_gyrus ~ mc_age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

############################# Sex Norm Quad Regions ################################
lmer_volume_quad <- lmer(ICV_Norm_TotalGray ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Accumbens_Area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Amygdala ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Amygdala ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Hippocampus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Hippocampus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Putamen ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Putamen ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Thalamus_Proper ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Thalamus_Proper ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_AIns_anterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_AIns_anterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_AnG_angular_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Calc_calcarine_cortex ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_CO_central_operculum ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Ent_entorhinal_area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Ent_entorhinal_area ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_FRP_frontal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_FuG_fusiform_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_FuG_fusiform_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_ITG_inferior_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_LiG_lingual_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MFG_middle_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MPoG_postcentral_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_MPoG_postcentral_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MPrG_precentral_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_MPrG_precentral_gyrus_medial_segment ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MTG_middle_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_MTG_middle_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_PIns_posterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_PIns_posterior_insula ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_PoG_postcentral_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_SMG_supramarginal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_STG_superior_temporal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_TMP_temporal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_TMP_temporal_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_TrIFG_triangular_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)


sink()

PSdatacopy1 <- PSdata[sample(nrow(distinct(PSdata$SessionID_Split1)), 120),]
PSdatacopy1 <- PSdata %>% distinct(SessionID_Split1)
library(dplyr)

for (i in 1:length(PSdata)) {
  #print(i)
  #PSdata$SessionID_Split1[i]
  
  
}
PSdatacopy1 <- distinct(PSdata, SessionID_Split1, .keep_all = TRUE)



regions <-as.data.frame(PSdata[,c(14:130)])
dim(regions)

lmer_tableinit <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

for (i in 1:length(regions)){
  ROI <- regions[,i]
  lmer_volume_quad <- lmer(ROI ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  AIC_quad
  summary(lmer_volume_quad)
  final_etasq <- eta_sq(lmer_volume_quad)
  
  df.etasq <- as.data.frame(final_etasq)
  df.etasq['Region'] <- c(colnames(regions)[i])
  df.master.etasq <- rbind(df.master.etasq, df.etasq)
  
  #getting coefficients etc into a table
  df.coef <- as.data.frame(coef(summary(lmer_volume_quad)))
  df.coef['Region'] <- c(colnames(regions)[i])
  df.master.coef <- rbind(df.master.coef, df.coef)
  
  write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_Sex_allquad_Sumamry_trim.csv")
  write.csv(df.master.etasq, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_Sex_allquad_EtaSq_trim.csv")
}























plot_list <- list()
regions <-as.data.frame(PSdata[,c(15:130)])
finalplot <-0
for (i in 1:length(regions)) {
  #bigPlot <- ggplot() + ggplot() + ggplot() + ggplot() + plot_layout(ncol=4)
  i<-i
  ROI <- regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
  #finalplot <-0
  plotfilename <- paste0(colnames(regions)[i])
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  #null model
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  CoefDataNull <- as.data.frame(fixef(lmer_volume_null))
  IntercNull <- as.numeric(CoefDataNull[1,1])
  #linear model
  lmer_volume_lin <- lmer(ROI ~ Age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  CoefDataLin <- as.data.frame(fixef(lmer_volume_lin))
  LinCoLin <- as.numeric(CoefDataLin[2,1])
  IntercLin <- as.numeric(CoefDataLin[1,1])
  #quadratic model
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad))
  QuadCoQuad <- as.numeric(CoefDataQuad[3,1])
  LinCoQuad <- as.numeric(CoefDataQuad[2,1])
  IntercQuad <- as.numeric(CoefDataQuad[1,1])
  
  ROIfnlin <- function(x) {
    as.numeric(LinCoLin*x + IntercLin) # outcome is y
  }
  
  ROIfnquad <- function(x) {
    as.numeric(QuadCoQuad*I(x^2) + LinCoQuad*x + IntercQuad) # outcome is y
  }
  
  ROIfnnull <- function(x) {
    as.numeric(IntercNull) # outcome is y
  }
  
  
  print(colnames(regions)[i])
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    finalmodel<- lmer_volume_lin
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    finalplot <- plotROI + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
      geom_function(fun = ROIfnlin, color="red", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    finalplot <- plotROI + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
      geom_function(fun = ROIfnnull, color="black", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    finalplot <- plotROI + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
      geom_function(fun = ROIfnquad, color="orange", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
  }else {
    bestAIC <- bestAIC
    finalmodel<-finalmodel
    finalplot<-finalplot
  }
  #print(finalmodel)
  #print(finalplot)
  #bigPlot <- bigPlot + finalplot
  #print('----------------------------------------------------------------------------------------------------------------')
  #print('----------------------------------------------------------------------------------------------------------------')
  #bigPlot <- bigPlot + finalplot
  #ggsave(filename = plotfilename, device = "png", width = 15, height = 15)
  plot_list[[i]] <- ggplotGrob(finalplot)
}
bigPlot <- grid.arrange(grobs=plot_list, ncol=4)


ggsave("BigPlot.png", bigPlot, device="png", height=49, width=25, limitsize = FALSE)


###########
PSdata$Total_Gray <- PSdata$Total_Gray/10
lmer_volume_quad <- lmer(Total_Gray ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)

CoefData <- as.data.frame(fixef(lmer_volume_quad))

QuadCo <- as.numeric(CoefData[3,1])
LinCo <- as.numeric(CoefData[2,1])
Interc <- as.numeric(CoefData[1,1])


ROIfnquad <- function(x) {
  as.numeric(QuadCo*I(x^2) + LinCo*x + Interc) # outcome is y
}

plotTotGray <- ggplot(PSdata, aes(x=Age, y=Total_Gray, group=Subject)) 
finalplotcm3 <- plotTotGray + geom_point(data=PSdata, aes(y=predict(lmer_volume_quad), shape=Sex.F1_M2.),color='dark gray',size=1.7) + 
  geom_line(data=PSdata, aes(y=predict(lmer_volume_quad), linetype=Sex.F1_M2.), color="dark gray", size=.6)+
  geom_function(fun = ROIfnquad, color="orange", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Total Gray Matter ("~cm^3*")") +
  scale_linetype_discrete(name="Sex", labels=c("Female","Male"))+scale_shape_discrete(name="Sex", labels=c("Female","Male")) +
  theme(legend.text=element_text(size=20), legend.title=element_blank())

finalplotcm3

#######################

lmer_volume_quad1 <- lmer(Total_Gray_ICV_Norm ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)

CoefData1 <- as.data.frame(fixef(lmer_volume_quad1))

QuadCoa <- as.numeric(CoefData1[3,1])
LinCoa <- as.numeric(CoefData1[2,1])
Interca <- as.numeric(CoefData1[1,1])


ROIfnquad1 <- function(x) {
  as.numeric(QuadCoa*I(x^2) + LinCoa*x + Interca) # outcome is y
}

plotTotGray1 <- ggplot(PSdata, aes(x=Age, y=Total_Gray_ICV_Norm, group=Subject)) 
finalplotICV <- plotTotGray1 + geom_point(data=PSdata, aes(y=predict(lmer_volume_quad1), shape=Sex.F1_M2.),color='dark gray',size=1.7) + 
  geom_line(data=PSdata, aes(y=predict(lmer_volume_quad1), linetype=Sex.F1_M2.), color="dark gray", size=.6)+
  geom_function(fun = ROIfnquad1, color="orange", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Total Gray Matter (% ICV)") +
  scale_linetype_discrete(name="Sex", labels=c("Female","Male"))+scale_shape_discrete(name="Sex", labels=c("Female","Male")) +
  theme(legend.text=element_text(size=20), legend.title=element_blank())

finalplotICV

install.packages("ggpubr")
library(ggpubr)

ggarrange(finalplotcm3, finalplotICV, ncol=2, common.legend = TRUE, legend="right")


###########

lmer_volume_quad <- lmer(Right_AIns_anterior_insula ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)

CoefData <- as.data.frame(fixef(lmer_volume_quad))

QuadCo <- as.numeric(CoefData[3,1])
LinCo <- as.numeric(CoefData[2,1])
Interc <- as.numeric(CoefData[1,1])


ROIfnquad <- function(x) {
  as.numeric(QuadCo*I(x^2) + LinCo*x + Interc) # outcome is y
}

plotRAIns <- ggplot(PSdata, aes(x=Age, y=Right_AIns_anterior_insula, group=Subject)) 
plotRAIns1 <- plotRAIns + geom_point(data=PSdata, aes(y=predict(lmer_volume_quad), shape=Sex.F1_M2.),color='dark gray',size=1.7) + 
  geom_line(data=PSdata, aes(y=predict(lmer_volume_quad), linetype=Sex.F1_M2.), color="dark gray", size=.6)+
  geom_function(fun = ROIfnquad, color="plum 4", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Right Anterior Insula ("~mm^3*")") +
  scale_linetype_discrete(name="Sex", labels=c("Female","Male"))+scale_shape_discrete(name="Sex", labels=c("Female","Male")) +
  theme(legend.text=element_text(size=20), legend.title=element_blank())

plotRAIns1

#######################

lmer_volume_quad1 <- lmer(Left_SMC_supplementary_motor_cortex ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)

CoefData1 <- as.data.frame(fixef(lmer_volume_quad1))

QuadCoa <- as.numeric(CoefData1[3,1])
LinCoa <- as.numeric(CoefData1[2,1])
Interca <- as.numeric(CoefData1[1,1])


ROIfnquad1 <- function(x) {
  as.numeric(QuadCoa*I(x^2) + LinCoa*x + Interca) # outcome is y
}

plotsmc <- ggplot(PSdata, aes(x=Age, y=Left_SMC_supplementary_motor_cortex, group=Subject)) 
finalplotsmc <- plotsmc + geom_point(data=PSdata, aes(y=predict(lmer_volume_quad1), shape=Sex.F1_M2.),color='dark gray',size=1.7) + 
  geom_line(data=PSdata, aes(y=predict(lmer_volume_quad1), linetype=Sex.F1_M2.), color="dark gray", size=.6)+
  geom_function(fun = ROIfnquad1, color="plum 4", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Left Suppl. Motor Cortex ("~mm^3*")") +
  scale_linetype_discrete(name="Sex", labels=c("Female","Male"))+scale_shape_discrete(name="Sex", labels=c("Female","Male")) +
  theme(legend.text=element_text(size=20), legend.title=element_blank())

finalplotsmc

install.packages("ggpubr")
library(ggpubr)

ggarrange(plotRAIns1, finalplotsmc, ncol=2, common.legend = TRUE, legend="right")
