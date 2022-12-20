###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long June 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, sjstats, pkbrtest)
devtools::install_github("cardiomoon/ggiraphExtra")

install.packages('lmerTest')
library(lmerTest)
install.packages('pbkrtest')
library(pbkrtest)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Preschool_Volume")

PSdata <- read.csv("Preschool_Volume_Longitudinal_July2022.csv")

#recoding variables into differenty data types, for example string into numeric
PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$Age <- as.numeric(PSdata$Age)
PSdata$Subject <- as.factor(PSdata$Subject)
PSdata$Sex.F1_M2. <- as.factor(PSdata$Sex.F1_M2.)
PSdata$OrderedByAgeAtSsn1 <- as.factor(PSdata$OrderedByAgeAtSsn1)
PSdata$income_3 <- as.factor(PSdata$income_3)
PSdata$schooling_5 <- as.factor(PSdata$schooling_5)
PSdata$schooling_3 <- as.factor(PSdata$schooling_3)


######### Results Log File
sink("PSVolumeLog_Raw_NoCoV_Null_Lin_Quad_NoAICcutoff_May2022", append=FALSE, split=TRUE)

######### All MaCRUISE regions, not normalized 
regions <-as.data.frame(PSdata[,c(14:130)])
dim(regions)
magic_for(silent = TRUE)

lmer_tableinit <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

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
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  
  print(colnames(regions)[i])
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    
    final_etasq <- eta_sq(finalmodel)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI1 <- plotROI + geom_line(aes(y=predict(lmer_volume_lin), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
     # geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
     # theme_classic() +
     # theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    #  theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
     # theme(legend.position = "none") +
      #coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    #finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI2 <- plotROI + geom_line(aes(y=predict(lmer_volume_null), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      #geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      #theme_classic() +
      #theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      #theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      #theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      #theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    #finalplot <- plotROI2
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    
    final_etasq <- eta_sq(finalmodel)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI3 <- plotROI + geom_line(aes(y=predict(lmer_volume_quad), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
      #geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
      #theme_classic() +
      #theme (plot.title = element_text(size=40, color="black",face="bold")) + 
      #theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
      #theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
      #theme(legend.position = "none") +
      #coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    #finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  #print(finalmodel)
  

  #magic getting AIC values into a table
  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad) 
  
  #getting coefficients etc into a table
  #df.coef <- as.data.frame(coef(summary(finalmodel)))
  #df.coef['Region'] <- c(colnames(regions)[i])
  #df.master.coef <- rbind(df.master.coef, df.coef)
  
  #print(paste("AIC_null =", AIC_null))
  #print(paste("AIC_lin =", AIC_lin))
  #print(paste("AIC_quad =", AIC_quad))
  #print(finalmodel)
  #print(final_etasq)
  #print(coef(summary(finalmodel)))
  #print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_NoCoV_finalmodelsummary_trim.csv")
  write.csv(df.master.etasq, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_NoCov_finalmodeletasq_trim.csv")
}

Raw_NoCov <- magic_result_as_dataframe()
write.csv(Raw_NoCov, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_NoCoV_AIC.csv")
magic_free()
sink()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sink("PSVolumeLog_Norm_NoCoV_", append=FALSE, split=TRUE)
#magic_for(silent = TRUE)
########## All MaCRUISE Regions, normalized by ICV
regions <-as.data.frame(PSdata[,c(131:247)])
dim(regions)
magic_for(silent = TRUE)

lmer_tableinit <- lmer(Right_TTG_transverse_temporal_gyrus ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  print(head(ROI))
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
  lmer_volume_quad <- lmer(ROI ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
  AIC_quad <- AIC(lmer_volume_quad)
  
  
  print(colnames(regions)[i])
  if (AIC_null >= AIC_lin) {
    bestAIC <- AIC_lin
    finalmodel <- lmer_volume_lin
    
    final_etasq <- eta_sq(finalmodel)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI1 <- plotROI + geom_line(aes(y=predict(lmer_volume_lin), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    # geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    # theme_classic() +
    # theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    #  theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #  theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    # theme(legend.position = "none") +
    #coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    #finalplot <- plotROI1
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI2 <- plotROI + geom_line(aes(y=predict(lmer_volume_null), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    #geom_smooth(data=PSdata, method=lm, formula = y ~ x,se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    #theme_classic() +
    #theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    #theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    #theme(legend.position = "none") +
    #coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    #finalplot <- plotROI2
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    finalmodel <- lmer_volume_quad
    
    final_etasq <- eta_sq(finalmodel)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    #plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    #plotROI3 <- plotROI + geom_line(aes(y=predict(lmer_volume_quad), group=Subject), linetype = "dashed", size=.9) + scale_color_grey() +
    #geom_smooth(data=PSdata, method=lm, formula = y ~ x + I(x^2),se=FALSE,color='plum4', size=3.0, fullrange = TRUE)+
    #theme_classic() +
    #theme (plot.title = element_text(size=40, color="black",face="bold")) + 
    #theme(axis.text.x=element_text(size=40, color="black")) + theme(axis.text.y=element_text(size=38, color="black")) + 
    #theme(axis.title.x=element_text(size=40, color="black")) + theme(axis.title.y=element_text(size=38, color="black")) + 
    #theme(legend.position = "none") +
    #coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + xlab("Age (Years)") + ylab(""~mm^3*"")
    #finalplot <- plotROI3
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  print(finalmodel)
  
  
  #magic getting AIC values into a table
  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad) 
  
  #getting coefficients etc into a table
  df.coef <- as.data.frame(coef(summary(finalmodel)))
  df.coef['Region'] <- c(colnames(regions)[i])
  df.master.coef <- rbind(df.master.coef, df.coef)
  
  #print(paste("AIC_null =", AIC_null))
  #print(paste("AIC_lin =", AIC_lin))
  #print(paste("AIC_quad =", AIC_quad))
  #print(finalmodel)
  #print(final_etasq)
  #print(coef(summary(finalmodel)))
  #print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Norm_NoCoV_finalmodelsummary_percentICV_trim.csv")
  write.csv(df.master.etasq, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Norm_NoCov_finalmodeletasq_percentICV_trim.csv")
}

Norm_NoCov <- magic_result_as_dataframe()
write.csv(Norm_NoCov, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Norm_NoCov_AIC.csv")
magic_free()
sink()



sink("PSVolumeLog_Raw_Sex_Oct2021", append=FALSE, split=TRUE)
  
############################# Sex Raw Lin Regions ################################
lmer_volume_lin <- lmer(Right_Pallidum ~ Age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)



  
############################# Sex Raw Quad Regions ################################
lmer_volume_quad <- lmer(Total_Gray ~ mc_Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)


sink()

sink("PSVolumeLog_Norm_Sex_Oct2021", append=FALSE, split=TRUE)

############################# Sex Norm Lin Regions ################################
lmer_volume_lin <- lmer(ICV_Norm_Right_Caudate ~ Age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)



############################# Sex Norm Quad Regions ################################
lmer_volume_quad <- lmer(ICV_Norm_TotalGray ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)


sink()

sink("PSVolumeLog_Raw_Sex_Oct2021", append=FALSE, split=TRUE)

############################# MatEdu Raw Lin Regions ################################
lmer_volume_lin <- lmer(Right_Pallidum ~ Age*schooling_3 + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)




############################# MatEdu Raw Quad Regions ################################
lmer_volume_quad <- lmer(Total_Gray ~ Age*schooling_3 + I(Age^2)*schooling_3 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)


sink()

############################# MatEdu Norm Lin Regions ################################
lmer_volume_lin <- lmer(Right_Pallidum ~ Age*schooling_3 + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)




############################# MatEdu Norm Quad Regions ################################
lmer_volume_quad <- lmer(Total_Gray ~ Age*schooling_3 + I(Age^2)*schooling_3 + (1|Subject), REML = FALSE, data = PSdata)
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)


