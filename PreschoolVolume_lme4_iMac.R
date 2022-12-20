###Analyses for Preschool Project, volume of all MaCRUISE output regions
###Madison Long September 2020
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, nlme, tidyverse, methods, lmtest, gridExtra, grid, ggplot2, sjstats)
install.packages("devtools")
devtools::install_github("cardiomoon/ggiraphExtra")
devtools::install_github("hoxo-m/magicfor")
library(magicfor)

install.packages('lmerTest')
library(lmerTest)
install.packages('pbkrtest')
library(pbkrtest)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Preschool_Volume")

PSdata <- read.csv("Preschool_Volume_Longitudinal_November2020.csv")
PSdataSex <- read.csv("Preschool_Volume_Longitudinal_November2020_Sex.csv")
PSdataSES <- read.csv("Preschool_Volume_Longitudinal_November2020_SES.csv")

PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$Subject <- as.factor(PSdata$Subject)
PSdata$Sex.F1_M2. <- as.factor(PSdata$Sex.F1_M2.)
PSdata$OrderedByAgeAtSsn1 <- as.factor(PSdata$OrderedByAgeAtSsn1)

PSdataSex$ICV <- as.numeric(PSdataSex$ICV)
PSdataSex$Subject <- as.factor(PSdataSex$Subject)
PSdataSex$Sex.F1_M2. <- as.factor(PSdataSex$Sex.F1_M2.)

PSdataSES$ICV <- as.numeric(PSdataSES$ICV)
PSdataSES$Subject <- as.factor(PSdataSES$Subject)
PSdataSES$Sex.F1_M2. <- as.factor(PSdataSES$Sex.F1_M2.)
PSdataSES$Income_3 <- as.factor(PSdataSES$Income_3)
PSdataSES$Race <- as.factor(PSdataSES$Race)


#centering age
mc_age = scale(PSdata$Age, center = TRUE, scale = TRUE)
mc_age_sex = scale(PSdataSex$Age, center = TRUE, scale = TRUE)
mc_age_ses = scale(PSdataSES$Age, center = TRUE, scale = TRUE)


###Mixed effects model example in lme4
lmerRCaudate <- lmer(Right_Caudate ~ mc_age + I(mc_age^2) + mc_ICV + (1|Subject), REML = FALSE, data = PSdata)
summary(lmerRCaudate)
parameters::p_value(lmerRCaudate)

######### Results Log File
sink("PSVolumeLog_Raw_NoCoV_check", append=FALSE, split=TRUE)

######### All MaCRUISE regions, not normalized 
regions <-as.data.frame(PSdata[,c(11:142)])
dim(regions)
magic_for(silent = TRUE)

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  plotfilename <- paste0("NoCov_", colnames(regions)[i])
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
  #cubic model
  lmer_volume_cub <- lmer(ROI ~ mc_age + I(mc_age^2) + I(mc_age^3) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_cub <- AIC(lmer_volume_cub)
  AICdif1 <- AIC_null - AIC_lin
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI1 <- plotROI + geom_line(aes(y=predict(lmer_volume_lin), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI2 <- plotROI + geom_line(aes(y=predict(lmer_volume_null), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI3 <- plotROI + geom_line(aes(y=predict(lmer_volume_quad), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  if ((bestAIC - AIC_cub) >= 5){
    bestAIC <- AIC_cub
    finalmodel <- lmer_volume_cub
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI4 <- plotROI + geom_line(aes(y=predict(lmer_volume_cub), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2) + I(x^3),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI4
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  

  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)  
  #print(paste("AIC_null =", AIC_null))
  #print(paste("AIC_lin =", AIC_lin))
  #print(paste("AIC_quad =", AIC_quad))
  #print(paste("AIC_cub =", AIC_cub))
  #final_etasq <- eta_sq(finalmodel)
  #print(finalmodel)
  #print(final_etasq)
  #print(coef(summary(finalmodel)))
  #print(summary(finalmodel))
  ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  
}
Raw_NoCov <- magic_result_as_dataframe()
write.csv(Raw_NoCov, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Raw_NoCov_AIC.csv")
magic_free()
sink()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sink("PSVolumeLog_Norm_NoCoV", append=FALSE, split=TRUE)
magic_for(silent = TRUE)
########## All MaCRUISE Regions, normalized by ICV
regions <-as.data.frame(PSdata[,c(143:274)])
dim(regions)
for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  plotfilename <- paste0("NoCov_", colnames(regions)[i])
  finalmodel <- 0
  bestAIC <- 0
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdata)
  AIC_null <- AIC(lmer_volume_null)
  lmer_volume_lin <- lmer(ROI ~ mc_age + (1|Subject), REML = FALSE, data = PSdata)
  AIC_lin <- AIC(lmer_volume_lin)
  lmer_volume_quad <- lmer(ROI ~ mc_age + I(mc_age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  lmer_volume_cub <- lmer(ROI ~ mc_age + I(mc_age^2) + I(mc_age^3) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_cub <- AIC(lmer_volume_cub)
  AICdif1 <- AIC_null - AIC_lin
  AICdif2 <- AIC_lin - AIC_quad
  AICdif3 <- AIC_quad - AIC_cub
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI1 <- plotROI + geom_line(aes(y=predict(lmer_volume_lin), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI2 <- plotROI + geom_line(aes(y=predict(lmer_volume_null), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI3 <- plotROI + geom_line(aes(y=predict(lmer_volume_quad), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  if ((bestAIC - AIC_cub) >= 5){
    bestAIC <- AIC_cub
    finalmodel <- lmer_volume_cub
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    plotROI4 <- plotROI + geom_line(aes(y=predict(lmer_volume_cub), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2) + I(x^3),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI4
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)
  #print(paste("AIC_null =", AIC_null))
  #print(paste("AIC_lin =", AIC_lin))
  #print(paste("AIC_quad =", AIC_quad))
  #print(paste("AIC_cub =", AIC_cub))
  #final_etasq <- eta_sq(finalmodel)
  #print(finalmodel)
  #print(final_etasq)
  #print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}
Norm_NoCov <- magic_result_as_dataframe()
write.csv(Norm_NoCov, "/Users/lebellab/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Norm_NoCov_AIC.csv")
magic_free()
sink()



sink("PSVolumeLog_Raw_Sex", append=FALSE, split=TRUE)
magic_for(silent = TRUE)
############ All MaCruise regions with Sex added, not normalized
regions <-as.data.frame(PSdataSex[,c(11:142)])
dim(regions)
for (i in 1:length(regions)) {
  ROI <- regions[,i]
  plotfilename <- paste0("Sex_", colnames(regions)[i])
  finalmodel <- 0
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_null <- AIC(lmer_volume_null)
  lmer_volume_lin <- lmer(ROI ~ mc_age_sex*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_lin <- AIC(lmer_volume_lin)
  lmer_volume_quad <- lmer(ROI ~ mc_age_sex*Sex.F1_M2. + I(mc_age_sex^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_quad <- AIC(lmer_volume_quad)
  lmer_volume_cub <- lmer(ROI ~ mc_age_sex*Sex.F1_M2. + I(mc_age_sex^2)*Sex.F1_M2. + I(mc_age_sex^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_cub <- AIC(lmer_volume_cub)
  AICdif1 <- AIC_null - AIC_lin
  AICdif2 <- AIC_lin - AIC_quad
  AICdif3 <- AIC_quad - AIC_cub
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point()
    plotROI1 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point(color="darkgrey")
    plotROI2 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point(color="darkgrey")
    plotROI3 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  if ((bestAIC - AIC_cub) >= 5){
    bestAIC <- AIC_cub
    finalmodel <- lmer_volume_cub
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point(color="darkgrey")
    plotROI4 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI4
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)
  print(paste("AIC_null =", AIC_null))
  print(paste("AIC_lin =", AIC_lin))
  print(paste("AIC_quad =", AIC_quad))
  print(paste("AIC_cub =", AIC_cub))
  #final_etasq <- eta_sq(finalmodel)
  print(finalmodel)
  #print(final_etasq)
  print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}
sink()
Raw_Sex <- magic_result_as_dataframe()
write.csv(Raw_Sex, "/Users/lebellab/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Raw_Sex_AIC.csv")
magic_free()

sink("PSVolumeLog_Norm_Sex", append=FALSE, split=TRUE)
magic_for(silent = TRUE)

################## All MaCRUISE regions normalized, sex interaction
regions <-as.data.frame(PSdataSex[,c(143:274)])
dim(regions)
for (i in 1:length(regions)) {
  ROI <- regions[,i]
  plotfilename <- paste0("Sex_", colnames(regions)[i])
  finalmodel <- 0
  bestAIC <- 0
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_null <- AIC(lmer_volume_null)
  lmer_volume_lin <- lmer(ROI ~ mc_age_sex*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_lin <- AIC(lmer_volume_lin)
  lmer_volume_quad <- lmer(ROI ~ mc_age_sex*Sex.F1_M2. + I(mc_age_sex^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_quad <- AIC(lmer_volume_quad)
  lmer_volume_cub <- lmer(ROI ~ mc_age_sex*Sex.F1_M2. + I(mc_age_sex^2)*Sex.F1_M2. + I(mc_age_sex^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdataSex)
  AIC_cub <- AIC(lmer_volume_cub)
  AICdif1 <- AIC_null - AIC_lin
  AICdif2 <- AIC_lin -AIC_quad
  AICdif3 <- AIC_quad - AIC_cub
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point()
    plotROI1 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point(color="darkgrey")
    plotROI2 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point(color="darkgrey")
    plotROI3 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  if ((bestAIC - AIC_cub) >= 5){
    bestAIC <- AIC_cub
    finalmodel <- lmer_volume_cub
    plotROI <- ggplot(PSdataSex, aes(x=Age, y=ROI, color=Sex.F1_M2.)) + geom_point(color="darkgrey")
    plotROI4 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI4
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)
  print(paste("AIC_null =", AIC_null))
  print(paste("AIC_lin =", AIC_lin))
  print(paste("AIC_quad =", AIC_quad))
  print(paste("AIC_cub =", AIC_cub))
  #final_etasq <- eta_sq(finalmodel)
  print(finalmodel)
  #print(final_etasq)
  print(summary(finalmodel))
  ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}
sink()
Norm_Sex <- magic_result_as_dataframe()
write.csv(Norm_Sex, "/Users/lebellab/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Norm_Sex_AIC.csv")
magic_free()

sink("PSVolumeLog_Raw_SES", append=FALSE, split=TRUE)
magic_for(silent = TRUE)
################## All MaCRUISE regions, not normalized, SES groups
regions <-as.data.frame(PSdataSES[,c(11:142)])
dim(regions)
for (i in 1:length(regions)) {
  ROI <- regions[,i]
  plotfilename <- paste0("SES_", colnames(regions)[i])
  finalmodel <- 0
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_null <- AIC(lmer_volume_null)
  lmer_volume_lin <- lmer(ROI ~ mc_age_ses*Income_3 + (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_lin <- AIC(lmer_volume_lin)
  lmer_volume_quad <- lmer(ROI ~ mc_age_ses*Income_3 + I(mc_age_ses^2)*Income_3 + (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_quad <- AIC(lmer_volume_quad)
  lmer_volume_cub <- lmer(ROI ~ mc_age_ses*Income_3 + I(mc_age_ses^2)*Income_3 + I(mc_age_ses^3)*Income_3 + (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_cub <- AIC(lmer_volume_cub)
  AICdif1 <- AIC_null - AIC_lin
  AICdif2 <- AIC_lin - AIC_quad
  AICdif3 <- AIC_quad - AIC_cub
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point()
    plotROI1 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point(color="darkgrey")
    plotROI2 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point(color="darkgrey")
    plotROI3 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  if ((bestAIC - AIC_cub) >= 5){
    bestAIC <- AIC_cub
    finalmodel <- lmer_volume_cub
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point(color="darkgrey")
    plotROI4 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    finalplot <- plotROI4
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)
  print(paste("AIC_null =", AIC_null))
  print(paste("AIC_lin =", AIC_lin))
  print(paste("AIC_quad =", AIC_quad))
  print(paste("AIC_cub =", AIC_cub))
  #final_etasq <- eta_sq(finalmodel)
  print(finalmodel)
  #print(final_etasq)
  print(summary(finalmodel))
  ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}

sink()
Raw_SES <- magic_result_as_dataframe()
write.csv(Raw_SES, "/Users/lebellab/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Raw_SES_AIC.csv")
magic_free()

sink("PSVolumeLog_Norm_SES", append=FALSE, split=TRUE)
magic_for(silent = TRUE)
################## All MaCRUISE regions normalized, SES groups
regions <-as.data.frame(PSdataSES[,c(143:274)])
dim(regions)
for (i in 1:length(regions)) {
  ROI <- regions[,i]
  plotfilename <- paste0("SES_", colnames(regions)[i])
  finalmodel <- 0
  lmer_volume_null <- lmer(ROI ~ (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_null <- AIC(lmer_volume_null)
  lmer_volume_lin <- lmer(ROI ~ mc_age_ses*Income_3 + (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_lin <- AIC(lmer_volume_lin)
  lmer_volume_quad <- lmer(ROI ~ mc_age_ses*Income_3 + I(mc_age_ses^2)*Income_3 + (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_quad <- AIC(lmer_volume_quad)
  lmer_volume_cub <- lmer(ROI ~ mc_age_ses*Income_3 + I(mc_age_ses^2)*Income_3 + I(mc_age_ses^3)*Income_3 + (1|Subject), REML = FALSE, data = PSdataSES)
  AIC_cub <- AIC(lmer_volume_cub)
  AICdif1 <- AIC_null - AIC_lin
  AICdif2 <- AIC_lin - AIC_quad
  AICdif3 <- AIC_quad - AIC_cub
  
  print(colnames(regions)[i])
  if (AICdif1 >= 5) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point()
    plotROI1 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point(color="darkgrey")
    plotROI2 <- plotROI + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI2
  }
  if ((bestAIC - AIC_quad) >= 5) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point(color="darkgrey")
    plotROI3 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  if ((bestAIC - AIC_cub) >= 5){
    bestAIC <- AIC_cub
    finalmodel <- lmer_volume_cub
    plotROI <- ggplot(PSdataSES, aes(x=Age, y=ROI, color=Income_3)) + geom_point(color="darkgrey")
    plotROI4 <- plotROI + geom_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), se = FALSE) +
      theme_classic() +
      theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
    finalplot <- plotROI4
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)
  print(paste("AIC_null =", AIC_null))
  print(paste("AIC_lin =", AIC_lin))
  print(paste("AIC_quad =", AIC_quad))
  print(paste("AIC_cub =", AIC_cub))
  #final_etasq <- eta_sq(finalmodel)
  print(finalmodel)
  #print(final_etasq)
  print(summary(finalmodel))
  ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}

###stop recording into log file
sink()
Norm_SES <- magic_result_as_dataframe()
write.csv(Norm_SES, "/Users/lebellab/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Norm_SES_AIC.csv")
magic_free()


#trying to loop file naming
regions <-as.data.frame(PSdata[,c(183:184)])
dim(regions)
for (i in 1:length(regions)) {
  ROI <- regions[,i]
  plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
  plotROI1 <- plotROI + geom_line(aes(y=predict(finalmodel), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    theme_classic() +
    theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    theme(legend.position = "none") +
    coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab("% ICV")
  plotfilename <- paste0(colnames(regions)[i], "_", i)
  ggsave(filename = plotfilename, plotROI1, device = "tiff", width = 10, height = 6)
}

for (i in 1:length(regions)) {
  ROI <- regions[,i]
  plotROI <- ggplot(data=PSdata, aes(x=Age,y=variablei, group=Subject))
  plotROI1 <- plotROI + geom_point() + geom_smooth(method="lm", se = FALSE)
  plotfilename <- paste0(colnames(regions)[i], "_", i)
  ggsave(filename = plotfilename, plotvar1, device = "tiff", width = 10, height = 6)
}





lmer_volume_quad <- lmer(ICV_Norm_Right_AnG_angular_gyrus ~ mc_age + I(mc_age^2) + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_quad)
print(lmer_volume_quad)
