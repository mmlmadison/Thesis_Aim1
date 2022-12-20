setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Preschool_Volume")

lmer_tableinit <- lmer(Right_Accumbens_Area ~ Age + (1|Subject), REML = FALSE, data = PSdata)
df.master.coef <- as.data.frame(coef(summary(lmer_tableinit)))
df.master.coef["Region"] <-c("Right_Inf_Lat_Vent")
df.master.etasq <- as.data.frame(eta_sq(lmer_tableinit))
df.master.etasq["Region"] <-c("Right_Inf_Lat_Vent")

######################## Sex plots

regions <-as.data.frame(PSdata[,c(131:247)])
dim(regions)
plot_list <- list()

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
  #Fplot <-0
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  M_data <-subset(PSdata,Sex.F1_M2.=='2')
  F_data <-subset(PSdata,Sex.F1_M2.=='1')
  
  M_data_regions <- as.data.frame(M_data[,c(131:247)])
  ROIM <- M_data_regions[,i]
  
  F_data_regions <- as.data.frame(F_data[,c(131:247)])
  ROIF <- F_data_regions[,i]

  
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
    lmer_volume_lin_Sex <- lmer(ROI ~ Age + Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Sex))
    lin_Interc <- as.numeric(CoefDataLin[1,1])
    lin_beta_Age <- as.numeric(CoefDataLin[2,1])
    lin_beta_Sex <- as.numeric(CoefDataLin[3,1])
    
    final_etasq <- eta_sq(lmer_volume_lin_Sex)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(lmer_volume_lin_Sex)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    
    MaleIntLin <- lin_Interc+lin_beta_Sex*2
    
    FemIntLin <- lin_Interc+lin_beta_Sex*1
    
    ROIfnlin_Male <- function(x) {
      as.numeric(lin_beta_Age*x + MaleIntLin) # outcome is y
    }
    
    ROIfnlin_Fem <- function(x) {
      as.numeric(lin_beta_Age*x + FemIntLin) # outcome is y
    }
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, group=Subject)) 
    Mplot <- plotROI + geom_point(data=M_data, aes(x=Age,y=ROIM),color='purple',size=.6) + geom_line(data=M_data, aes(x=Age,y=ROIM), color="purple", size=.2)
    Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=ROIF),color='dark green',size=.6) + geom_line(data=F_data, aes(x=Age,y=ROIF), color="dark green", size=.2) +
      geom_function(fun = ROIfnlin_Male, color="purple", size=2) +
      geom_function(fun = ROIfnlin_Fem, color="dark green", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
      
    #print(paste(MaleLinLin, MaleIntLin))
    #print(paste(FemLinLin,FemIntLin))
    
  }else {
    bestAIC <- AIC_null
    
    Fplot <- ggplot()
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    lmer_volume_quad_Sex <- lmer(ROI ~ Age + Sex.F1_M2. + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad_Sex))
    quad_Interc <- as.numeric(CoefDataQuad[1,1])
    quad_beta_Age <- as.numeric(CoefDataQuad[2,1])
    quad_beta_Sex <- as.numeric(CoefDataQuad[3,1])
    quad_beta_Age2 <- as.numeric(CoefDataQuad[4,1])
    
    MaleIntQuad <- quad_Interc+quad_beta_Sex*2
    
    FemIntQuad <- quad_Interc+quad_beta_Sex*1
    
    ROIfnquad_Male <- function(x) {
      as.numeric(quad_beta_Age2*x^2 + quad_beta_Age*x + MaleIntQuad) # outcome is y
    }
    
    ROIfnquad_Fem <- function(x) {
      as.numeric(quad_beta_Age2*x^2 + quad_beta_Age*x + FemIntQuad) # outcome is y
    }

    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, group=Subject)) 
    Mplot <- plotROI + geom_point(data=M_data, aes(x=Age,y=ROIM),color='purple',size=.6) + geom_line(data=M_data, aes(x=Age,y=ROIM), color="purple", size=.2)
    Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=ROIF),color='dark green',size=.6) + geom_line(data=F_data, aes(x=Age,y=ROIF), color="dark green", size=.2) +
      geom_function(fun = ROIfnquad_Male, color="purple", size=2) +
      geom_function(fun = ROIfnquad_Fem, color="dark green", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
      
      final_etasq <- eta_sq(lmer_volume_quad_Sex)
    
    
      df.etasq <- as.data.frame(final_etasq)
      df.etasq['Region'] <- c(colnames(regions)[i])
      df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
      df.coef <- as.data.frame(coef(summary(lmer_volume_quad_Sex)))
      df.coef['Region'] <- c(colnames(regions)[i])
      df.master.coef <- rbind(df.master.coef, df.coef)
     # print(paste(MaleQuadQuad, MaleLinQuad, MaleIntQuad))
     # print(paste(FemQuadQuad, FemLinQuad, FemIntQuad))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  

  
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  #print(Fplot)
  plot_list[[i]] <- ggplotGrob(Fplot)
  
  write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_Sex_finalmodelsINTERACTIONSREMOVED.csv")
  write.csv(df.master.etasq, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Abs_Sex_etasqINTERACTIONSREMOVED.csv")
}

bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_Sex_Norm_INTERACTIONSREMOVED.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

#####Sex total Gray
PSdata$Total_Gray <- PSdata$Total_Gray/10
lmer_volume_quad_Sex <- lmer(Total_Gray ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2.+ (1|Subject), REML = FALSE, data = PSdata)
CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad_Sex))
IntercQuad <- as.numeric(CoefDataQuad[1,1])
LinCoQuad <- as.numeric(CoefDataQuad[2,1])
SexIntCoQuad <- as.numeric(CoefDataQuad[3,1])
QuadCoQuad <- as.numeric(CoefDataQuad[4,1])
SexLinCoQuad <- as.numeric(CoefDataQuad[5,1])
SexQuadCoQuad <- as.numeric(CoefDataQuad[6,1])

MaleQuadQuad <-QuadCoQuad+SexQuadCoQuad*2
MaleLinQuad <- LinCoQuad+SexLinCoQuad*2
MaleIntQuad <- IntercQuad+SexIntCoQuad*2

FemQuadQuad <-QuadCoQuad+SexQuadCoQuad*1
FemLinQuad <- LinCoQuad+SexLinCoQuad*1
FemIntQuad <- IntercQuad+SexIntCoQuad*1

ROIfnquad_Male <- function(x) {
  as.numeric(MaleQuadQuad*x^2 + MaleLinQuad*x + MaleIntQuad) # outcome is y
}

ROIfnquad_Fem <- function(x) {
  as.numeric(FemQuadQuad*x^2 + FemLinQuad*x + FemIntQuad) # outcome is y
}

M_data <-subset(PSdata,Sex.F1_M2.=='2')
F_data <-subset(PSdata,Sex.F1_M2.=='1')

plotTotGray <- ggplot(PSdata, aes(x=Age, y=Total_Gray, group=Subject)) 
Mplot <- plotTotGray + geom_point(data=M_data, aes(x=Age,y=Total_Gray),color='purple',size=.6) + geom_line(data=M_data, aes(x=Age,y=Total_Gray), color="purple", size=.2)
Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=Total_Gray),color='dark green',size=.6) + geom_line(data=F_data, aes(x=Age,y=Total_Gray), color="dark green", size=.2) +
  geom_function(fun = ROIfnquad_Male, color="purple", size=2) +
  geom_function(fun = ROIfnquad_Fem, color="dark green", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Total Gray ("~cm^3*")")
Fplot

######################## Sex Norm

regions <-as.data.frame(PSdata[,c(132:247)])
plot_list <- list()

for (i in 1:length(regions)) { 
  ROI <- regions[,i]
  finalmodel <- 0
  bestAIC <- 0 
  #Fplot <-0
  colNamesAll <- as.data.frame(colnames(regions))
  colName <- colNamesAll[i,1]
  
  
  
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
    lmer_volume_lin_Sex <- lmer(ROI ~ Age + Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Sex))
    lin_Interc <- as.numeric(CoefDataLin[1,1])
    lin_beta_Age <- as.numeric(CoefDataLin[2,1])
    lin_beta_Sex <- as.numeric(CoefDataLin[3,1])
    
    final_etasq <- eta_sq(lmer_volume_lin_Sex)
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(lmer_volume_lin_Sex)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    
    
  }else {
    bestAIC <- AIC_null
    
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    lmer_volume_quad_Sex <- lmer(ROI ~ Age + Sex.F1_M2. + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad_Sex))
    quad_Interc <- as.numeric(CoefDataQuad[1,1])
    quad_beta_Age <- as.numeric(CoefDataQuad[2,1])
    quad_beta_Sex <- as.numeric(CoefDataQuad[3,1])
    quad_beta_Age2 <- as.numeric(CoefDataQuad[4,1])
    
    
    final_etasq <- eta_sq(lmer_volume_quad_Sex)
    
    
    df.etasq <- as.data.frame(final_etasq)
    df.etasq['Region'] <- c(colnames(regions)[i])
    df.master.etasq <- rbind(df.master.etasq, df.etasq)
    
    df.coef <- as.data.frame(coef(summary(lmer_volume_quad_Sex)))
    df.coef['Region'] <- c(colnames(regions)[i])
    df.master.coef <- rbind(df.master.coef, df.coef)
    # print(paste(MaleQuadQuad, MaleLinQuad, MaleIntQuad))
    # print(paste(FemQuadQuad, FemLinQuad, FemIntQuad))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  
  
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  #print(Fplot)

  
  write.csv(df.master.coef, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Norm_Sex_finalmodelsINTERACTIONSREMOVED.csv")
  write.csv(df.master.etasq, "/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Table_Norm_Sex_etasqINTERACTIONSREMOVED.csv")
}

bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_Sex_Norm_INTERACTIONSREMOVED.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

