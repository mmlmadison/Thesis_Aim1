if (!require("pacman")) install.packages("pacman")
pacman::p_load(mixedpower, foreign, plyr, lattice, lme4, tidyverse, methods, lmertest, gridExtra, grid, ggplot2, sjstats, pkbrtest, simr)

setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Preschool_Volume")

PSdata <- read.csv("Preschool_Volume_Longitudinal_July2022.csv")

#recoding variables into differenty data types, for example string into numeric
PSdata$ICV <- as.numeric(PSdata$ICV)
PSdata$Age <- as.numeric(PSdata$Age)
PSdata$Subject <- as.factor(PSdata$Subject)
PSdata$Sex.F1_M2. <- as.factor(PSdata$Sex.F1_M2.)

lmer_volume_quad <- lmer(Total_Gray ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)

powerSim(lmer_volume_quad, fixef=c("Age", "Sex.F1_M2."))
pc=powerCurve(lmer_volume_quad, breaks=c(100,200,300,393))
plot(pc)


