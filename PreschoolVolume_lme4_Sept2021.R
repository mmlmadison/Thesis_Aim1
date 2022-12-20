###Analyses for Preschool Project, volume of all gray matter MaCRUISE output regions
###Madison Long September 2021

### Differences since June 2021 version:
##### Cubic fits removed
##### Total GM volume and ICV analyses in addition to region by region

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

PSdata <- read.csv("Preschool_Volume_Longitudinal_Sept2021.csv")

PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$Subject <- as.factor(PSdata$Subject)
PSdata$Sex.F1_M2. <- as.factor(PSdata$Sex.F1_M2.)
PSdata$OrderedByAgeAtSsn1 <- as.factor(PSdata$OrderedByAgeAtSsn1)

#centering age
mc_age = scale(PSdata$Age, center = TRUE, scale = TRUE)

######### Results Log File
sink("PSVolumeLog_Raw_NoCoV_Sept2021", append=FALSE, split=TRUE)

######### All MaCRUISE regions, not normalized 
regions <-as.data.frame(PSdata[,c(14:146)])
dim(regions)
#magic_for(silent = TRUE)

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
  


  #put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)  
  print(paste("AIC_null =", AIC_null))
  print(paste("AIC_lin =", AIC_lin))
  print(paste("AIC_quad =", AIC_quad))
  final_etasq <- eta_sq(finalmodel)
  #print(finalmodel)
  print(final_etasq)
  #print(coef(summary(finalmodel)))
  #print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  
}


#Raw_NoCov <- magic_result_as_dataframe()
#write.csv(Raw_NoCov, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Raw_NoCov_AIC.csv")
#magic_free()
sink()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sink("PSVolumeLog_Norm_NoCoV_Sept2021", append=FALSE, split=TRUE)
#magic_for(silent = TRUE)
########## All MaCRUISE Regions, normalized by ICV
regions <-as.data.frame(PSdata[,c(147:279)])
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
  print(AIC_null)
  print(AIC_lin)
  print(AIC_quad)
  
 # put(colnames(regions)[i], AIC_null, AIC_lin, AIC_quad, AIC_cub)
  #print(paste("AIC_null =", AIC_null))
  #print(paste("AIC_lin =", AIC_lin))
  #print(paste("AIC_quad =", AIC_quad))
  #print(paste("AIC_cub =", AIC_cub))
  
  final_etasq <- eta_sq(finalmodel)
  print(final_etasq)
  #print(summary(finalmodel))
  #ggsave(filename = plotfilename, finalplot, device = "tiff", width = 10, height = 6)
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
}
#Norm_NoCov <- magic_result_as_dataframe()
#write.csv(Norm_NoCov, "/Users/lebellab/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Norm_NoCov_AIC.csv")
#magic_free()
sink()



sink("PSVolumeLog_Raw_Sex_Sept2021", append=FALSE, split=TRUE)

############################# Sex and Mat Edu Raw Null Regions ################################
lmer_volume_null <- lmer(Right_Basal_Forebrain ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)
AIC_null
lmer_volume_null
eta_sq(lmer_volume_null)

lmer_volume_null <- lmer(Right_Cun_cuneus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_Cun_cuneus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_LOrG_lateral_orbital_gyrusn ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_MCgG_middle_cingulate_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_MOG_middle_occipital_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_MPrG_precentral_gyrus_medial_segment ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_MPrG_precentral_gyrus_medial_segment ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_OCP_occipital_pole ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_OFuG_occipital_fusiform_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_OFuG_occipital_fusiform_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_PoG_postcentral_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_PoG_postcentral_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_PrG_precentral_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_SOG_superior_occipital_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_SPL_superior_parietal_lobule ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Left_SPL_superior_parietal_lobule ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)

lmer_volume_null <- lmer(Right_TTG_transverse_temporal_gyrus ~ (1|Subject), REML = FALSE, data = PSdata)
AIC_null <- AIC(lmer_volume_null)
  
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
  
############################# Sex and Raw Cub Regions ################################
lmer_volume_cub <- lmer(Left_Ventral_DC ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + I(mc_age^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)

lmer_volume_cub <- lmer(Right_AOrG_anterior_orbital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + I(mc_age^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)

lmer_volume_cub <- lmer(Right_OrIFG_orbital_part_of_the_inferior_frontal_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + I(mc_age^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)
  

sink()

sink("PSVolumeLog_Norm_Sex_QuadCub_Sept2021", append=FALSE, split=TRUE)

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


############################# Sex Norm Cub Regions ################################



lmer_volume_cub <- lmer(ICV_Norm_Left_OCP_occipital_pole ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + I(mc_age^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)

lmer_volume_cub <- lmer(ICV_Norm_Right_SOG_superior_occipital_gyrus ~ mc_age*Sex.F1_M2. + I(mc_age^2)*Sex.F1_M2. + I(mc_age^3)*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)


sink()


sink("PSVolumeLog_Raw_MatEdu_Sept2021", append=FALSE, split=TRUE)

############################# MatEdu Raw Lin Regions ################################

lmer_volume_lin <- lmer(Right_Pallidum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_Pallidum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_ACgG_anterior_cingulate_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_AOrG_anterior_orbital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_CO_central_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_CO_central_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_FO_frontal_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_FO_frontal_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_GRe_gyrus_rectus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_GRe_gyrus_rectus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_LiG_lingual_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_LOrG_lateral_orbital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MCgG_middle_cingulate_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MOrG_medial_orbital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_MOrG_medial_orbital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MPoG_postcentral_gyrus_medial_segment ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_MPoG_postcentral_gyrus_medial_segment ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_MSFG_superior_frontal_gyrus_medial_segment ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_OCP_occipital_pole ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_OpIFG_opercular_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_OrIFG_orbital_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_PHG_parahippocampal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_PO_parietal_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_PO_parietal_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_PrG_precentral_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_PT_planum_temporale ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_PT_planum_temporale ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_SCA_subcallosal_area ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SCA_subcallosal_area ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_SFG_superior_frontal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SFG_superior_frontal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Right_SMC_supplementary_motor_cortex ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SMC_supplementary_motor_cortex ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_SMG_supramarginal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_TrIFG_triangular_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(Left_TTG_transverse_temporal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)


############################# MatEdu Raw Quad Regions ################################

lmer_volume_quad <- lmer(Right_Accumbens_Area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Accumbens_Area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Amygdala ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Amygdala ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Caudate ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Caudate ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Hippocampus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Hippocampus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Putamen ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Putamen ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Thalamus_Proper ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Thalamus_Proper ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Ventral_DC ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Basal_Forebrain ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_ACgG_anterior_cingulate_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_AIns_anterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_AIns_anterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_AnG_angular_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_AnG_angular_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Calc_calcarine_cortex ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Calc_calcarine_cortex ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_Ent_entorhinal_area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_Ent_entorhinal_area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_FRP_frontal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_FRP_frontal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_FuG_fusiform_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_FuG_fusiform_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_IOG_inferior_occipital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_IOG_inferior_occipital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_ITG_inferior_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_ITG_inferior_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_LiG_lingual_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MFC_medial_frontal_cortex ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MFC_medial_frontal_cortex ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MFG_middle_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MFG_middle_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MOG_middle_occipital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MSFG_superior_frontal_gyrus_medial_segment ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_MTG_middle_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_MTG_middle_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_OpIFG_opercular_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PCgG_posterior_cingulate_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PCgG_posterior_cingulate_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PCu_precuneus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PCu_precuneus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PHG_parahippocampal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PIns_posterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PIns_posterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_POrG_posterior_orbital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_POrG_posterior_orbital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_PP_planum_polare ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_PP_planum_polare ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_SMG_supramarginal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_SOG_superior_occipital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_STG_superior_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_STG_superior_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_TMP_temporal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Left_TMP_temporal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(Right_TrIFG_triangular_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

############################# MatEdu Raw Cub Regions ################################




lmer_volume_cub <- lmer(Left_Ventral_DC ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + I(mc_age^3)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)

lmer_volume_cub <- lmer(Right_AOrG_anterior_orbital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + I(mc_age^3)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)

lmer_volume_cub <- lmer(Right_OrIFG_orbital_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + I(mc_age^3)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)
  
  
  

sink()


sink("PSVolumeLog_Norm_SES_Sept2021", append=FALSE, split=TRUE)
############################# MatEdu Norm Lin Regions ################################
lmer_volume_lin <- lmer(ICV_Norm_Right_Caudate ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Caudate ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_Ventral_DC ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Ventral_DC ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Basal_Forebrain ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_Basal_Forebrain ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_ACgG_anterior_cingulate_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_AOrG_anterior_orbital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_CO_central_operculum ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_Cun_cuneus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_Cun_cuneus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_IOG_inferior_occipital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_IOG_inferior_occipital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_LOrG_lateral_orbital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_MCgG_middle_cingulate_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_MCgG_middle_cingulate_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_MOG_middle_occipital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_MOG_middle_occipital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_OFuG_occipital_fusiform_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_OFuG_occipital_fusiform_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PCgG_posterior_cingulate_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_PCu_precuneus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PCu_precuneus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PoG_postcentral_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_PrG_precentral_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_PrG_precentral_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_SMC_supplementary_motor_cortex ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_SOG_superior_occipital_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_SPL_superior_parietal_lobule ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_SPL_superior_parietal_lobule ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_STG_superior_temporal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Right_TTG_transverse_temporal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)

lmer_volume_lin <- lmer(ICV_Norm_Left_TTG_transverse_temporal_gyrus ~ mc_age*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_lin <- AIC(lmer_volume_lin)
AIC_lin
summary(lmer_volume_lin)
eta_sq(lmer_volume_lin)


############################# MatEdu Norm Quad Regions ################################
lmer_volume_quad <- lmer(ICV_Norm_Left_Accumbens_Area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Accumbens_Area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Amygdala ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Amygdala ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Hippocampus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Hippocampus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Putamen ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Putamen ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Thalamus_Proper ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Thalamus_Proper ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_AIns_anterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_AIns_anterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_AnG_angular_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Calc_calcarine_cortex ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_CO_central_operculum ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_Ent_entorhinal_area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_Ent_entorhinal_area ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_FRP_frontal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_FuG_fusiform_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_FuG_fusiform_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_ITG_inferior_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_LiG_lingual_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MFG_middle_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MPoG_postcentral_gyrus_medial_segment ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_MPoG_postcentral_gyrus_medial_segment ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MPrG_precentral_gyrus_medial_segment ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_MPrG_precentral_gyrus_medial_segment ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_MTG_middle_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_MTG_middle_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_PIns_posterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_PIns_posterior_insula ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_PoG_postcentral_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_SMG_supramarginal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_STG_superior_temporal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_TMP_temporal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Left_TMP_temporal_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

lmer_volume_quad <- lmer(ICV_Norm_Right_TrIFG_triangular_part_of_the_inferior_frontal_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_quad <- AIC(lmer_volume_quad)
AIC_quad
summary(lmer_volume_quad)
eta_sq(lmer_volume_quad)

############################# MatEdu Norm Cub Regions ################################
lmer_volume_cub <- lmer(ICV_Norm_Left_OCP_occipital_pole ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + I(mc_age^3)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)

lmer_volume_cub <- lmer(ICV_Norm_Right_SOG_superior_occipital_gyrus ~ mc_age*schooling_5 + I(mc_age^2)*schooling_5 + I(mc_age^3)*schooling_5 + (1|Subject), REML = FALSE, data = PSdata)
AIC_cub <- AIC(lmer_volume_cub)
AIC_cub
summary(lmer_volume_cub)
eta_sq(lmer_volume_cub)


sink()

