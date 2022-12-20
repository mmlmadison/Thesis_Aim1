setwd("/Users/madisonlong/Dropbox/UniversityofCalgary/LebelLab/LebelLab_Data/Preschool_Volume")

PSdata_toPivot <- read.csv("Preschool_Volume_Longitudinal_Sept2022_ToPivot.csv")
PSdata_wide <- pivot_wider(PSdata_toPivot, id_cols = SessionID_Split1, names_from = SessionID_Split2, values_from = c(Age, Sex.F1_M2., Race, marital.status., income, schooling))
PSdata_wide
write.csv(PSdata_wide, "Preschool_Volume_Longitudinal_Sept2022_wide.csv")

PSdata_wide_trim$schooling <- as.factor(PSdata_wide_trim$schooling)

PSdata_wide_trim <- read.csv("Preschool_Volume_Longitudinal_Sept2022_wide_trimmed.csv")

EducationBars <- ggplot(PSdata_wide_trim, aes(x=schooling)) + geom_bar(fill="#D2C60A")
EducationBars1 <- EducationBars + scale_x_discrete(limits=c("finished high school", "some post-secondary", "trade/technical diploma", "undergraduate degree", "some postgraduate", "postgraduate degree")) + 
  geom_text(stat="count", aes(label=..count..), vjust=-.7, size=10) + theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = -45, vjust = -.06, hjust = .08 , size = 25), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
EducationBars1

IncomeBars <- ggplot(PSdata_wide_trim, aes(x=income_Session01)) + geom_bar(fill="dark blue")
IncomeBars1 <- IncomeBars + scale_x_discrete(limits=c("$25 000 - $49 999", "$50 000 - $74 999", "$75 000 - $99 999", "$100 000 - $124 999", "$125 000 - $149 999", "$150 000 - $174 999", "over $175 000")) +
  geom_text(stat="count", aes(label=..count..), vjust=-.7, size=10) + theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = -45, vjust = -.06, hjust = .08 , size = 25), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
IncomeBars1

MaritalBars <- ggplot(PSdata_wide_trim, aes(x=marital.status.)) + geom_bar(fill="#D2C60A")
MaritalBars1 <- MaritalBars + scale_x_discrete(limits=c("single/never been married", "married/common law", "separated", "divorced")) +
  geom_text(stat="count", aes(label=..count..), vjust=-.7, size=10) + theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = -45, vjust = -.06, hjust = .08 , size = 25), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
MaritalBars1

RaceBars <- ggplot(PSdata_wide_trim, aes(x=Race)) + geom_bar(fill="dark green")
RaceBars1 <- RaceBars + scale_x_discrete(limits=c("White", "Asian/pacific islander", "$75 000 - $99 999", "$100 000 - $124 999", "$125 000 - $149 999", "$150 000 - $174 999", "over $175 000")) +
  geom_text(stat="count", aes(label=..count..), vjust=-.7, size=10) + theme_minimal() + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = -45, vjust = -.06, hjust = .08 , size = 25), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")
RaceBars1


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
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    finalplot <- plotROI + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
      geom_function(fun = ROIfnnull, color="black", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
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
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
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
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)


ggsave("BigPlot_Abs_All.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

#####Total Gray
PSdata$Total_Gray <- PSdata$Total_Gray/10
lmer_volume_quad <- lmer(Total_Gray ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad))
QuadCoQuad <- as.numeric(CoefDataQuad[3,1])
LinCoQuad <- as.numeric(CoefDataQuad[2,1])
IntercQuad <- as.numeric(CoefDataQuad[1,1])
plotTotGray <- ggplot(PSdata, aes(x=Age, y=Total_Gray, color=Subject)) + geom_point(color="darkgrey")
finalplot <- plotTotGray + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
  geom_function(fun = ROIfnquad, color="orange", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
  theme(legend.position = "none") +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Total Gray ("~cm^3*")")
finalplot




########ICV NoCov

plot_list <- list()
regions <-as.data.frame(PSdata[,c(132:247)])
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
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
    
  }else {
    bestAIC <- AIC_null
    finalmodel <- lmer_volume_null
    
    plotROI <- ggplot(PSdata, aes(x=Age, y=ROI, color=Subject)) + geom_point(color="darkgrey")
    finalplot <- plotROI + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
      geom_function(fun = ROIfnnull, color="black", size=2) +
      theme_classic() +
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
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
      theme (plot.title = element_text(size=20, color="black",face="bold")) + 
      theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
      theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
      theme(legend.position = "none") +
      coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab(paste(colName))
  }else {
    bestAIC <- bestAIC
    finalmodel<-finalmodel
    finalplot<-finalplot
  }
  print(colnames(regions)[i])
  print(colName)
  print(finalmodel)
  #print(finalplot)
  #bigPlot <- bigPlot + finalplot
  #print('----------------------------------------------------------------------------------------------------------------')
  #print('----------------------------------------------------------------------------------------------------------------')
  #bigPlot <- bigPlot + finalplot
  #ggsave(filename = plotfilename, device = "png", width = 15, height = 15)
  plot_list[[i]] <- ggplotGrob(finalplot)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_ICV_All.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

####ICV Total Gray
PSdata$Total_Gray_ICV_Norm <- PSdata$Total_Gray_ICV_Norm
lmer_volume_quad <- lmer(Total_Gray_ICV_Norm ~ Age + I(Age^2) + (1|Subject), REML = FALSE, data = PSdata)
CoefDataQuad <- as.data.frame(fixef(lmer_volume_quad))
QuadCoQuad <- as.numeric(CoefDataQuad[3,1])
LinCoQuad <- as.numeric(CoefDataQuad[2,1])
IntercQuad <- as.numeric(CoefDataQuad[1,1])
plotTotGray <- ggplot(PSdata, aes(x=Age, y=Total_Gray_ICV_Norm, color=Subject)) + geom_point(color="darkgrey")
finalplot <- plotTotGray + geom_line(aes(group=Subject), size=.4) + scale_color_grey(start = .6, end = .6) +
  geom_function(fun = ROIfnquad, color="orange", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) + 
  theme(legend.position = "none") +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Total Gray (%ICV)")
finalplot

######################## Sex ICV plots

regions <-as.data.frame(PSdata[,c(15:130)])
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
  
  M_data_regions <- as.data.frame(M_data[,c(15:130)])
  ROIM <- M_data_regions[,i]
  
  F_data_regions <- as.data.frame(F_data[,c(15:130)])
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
    lmer_volume_lin_Sex <- lmer(ROI ~ Age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Sex))
    IntercLin <- as.numeric(CoefDataLin[1,1])
    LinCoLin <- as.numeric(CoefDataLin[2,1])
    SexCoLin <- as.numeric(CoefDataLin[3,1])
    SexLinCoLin <- as.numeric(CoefDataLin[4,1])
    
    MaleLinLin <- LinCoLin+SexLinCoLin*2
    MaleIntLin <- IntercLin+SexCoLin*2
    
    FemLinLin <- LinCoLin+SexLinCoLin*1
    FemIntLin <- IntercLin+SexCoLin*1
    
    ROIfnlin_Male <- function(x) {
      as.numeric(MaleLinLin*x + MaleIntLin) # outcome is y
    }
    
    ROIfnlin_Fem <- function(x) {
      as.numeric(FemLinLin*x + FemIntLin) # outcome is y
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
      
    print(paste(MaleLinLin, MaleIntLin))
    print(paste(FemLinLin,FemIntLin))
    
  }else {
    bestAIC <- AIC_null
    
    Fplot <- ggplot()
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    lmer_volume_quad_Sex <- lmer(ROI ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2.+ (1|Subject), REML = FALSE, data = PSdata)
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
      
      print(paste(MaleQuadQuad, MaleLinQuad, MaleIntQuad))
      print(paste(FemQuadQuad, FemLinQuad, FemIntQuad))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  

  
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  #print(Fplot)
  plot_list[[i]] <- ggplotGrob(Fplot)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_Sex_Abs.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

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

######################## Sex Norm plots

regions <-as.data.frame(PSdata[,c(132:247)])
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
  
  M_data_regions <- as.data.frame(M_data[,c(132:247)])
  ROIM <- M_data_regions[,i]
  
  F_data_regions <- as.data.frame(F_data[,c(132:247)])
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
    lmer_volume_lin_Sex <- lmer(ROI ~ Age*Sex.F1_M2. + (1|Subject), REML = FALSE, data = PSdata)
    CoefDataLin <- as.data.frame(fixef(lmer_volume_lin_Sex))
    IntercLin <- as.numeric(CoefDataLin[1,1])
    LinCoLin <- as.numeric(CoefDataLin[2,1])
    SexCoLin <- as.numeric(CoefDataLin[3,1])
    SexLinCoLin <- as.numeric(CoefDataLin[4,1])
    
    MaleLinLin <- LinCoLin+SexLinCoLin*2
    MaleIntLin <- IntercLin+SexCoLin*2
    
    FemLinLin <- LinCoLin+SexLinCoLin*1
    FemIntLin <- IntercLin+SexCoLin*1
    
    ROIfnlin_Male <- function(x) {
      as.numeric(MaleLinLin*x + MaleIntLin) # outcome is y
    }
    
    ROIfnlin_Fem <- function(x) {
      as.numeric(FemLinLin*x + FemIntLin) # outcome is y
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
    
    print(paste(MaleLinLin, MaleIntLin))
    print(paste(FemLinLin,FemIntLin))
    
  }else {
    bestAIC <- AIC_null
    
    Fplot <- ggplot()
    
  }
  if (bestAIC >= AIC_quad) {
    bestAIC <- AIC_quad
    lmer_volume_quad_Sex <- lmer(ROI ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2.+ (1|Subject), REML = FALSE, data = PSdata)
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
    
    print(paste(MaleQuadQuad, MaleLinQuad, MaleIntQuad))
    print(paste(FemQuadQuad, FemLinQuad, FemIntQuad))
    
  }else {
    bestAIC <- bestAIC
    finalmodel <- finalmodel
  }
  
  
  
  
  print('----------------------------------------------------------------------------------------------------------------')
  print('----------------------------------------------------------------------------------------------------------------')
  #print(Fplot)
  plot_list[[i]] <- ggplotGrob(Fplot)
}
bigPlot <- marrangeGrob(grobs=plot_list, ncol=4, nrow=2)
ggsave("BigPlot_Sex_Norm.pdf", bigPlot, width=50, height=25, limitsize = FALSE)

#####Sex Norm total Gray

lmer_volume_quad_Sex <- lmer(Total_Gray_ICV_Norm ~ Age*Sex.F1_M2. + I(Age^2)*Sex.F1_M2.+ (1|Subject), REML = FALSE, data = PSdata)
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

plotTotGray <- ggplot(PSdata, aes(x=Age, y=Total_Gray_ICV_Norm, group=Subject)) 
Mplot <- plotTotGray + geom_point(data=M_data, aes(x=Age,y=Total_Gray_ICV_Norm),color='purple',size=.6) + geom_line(data=M_data, aes(x=Age,y=Total_Gray_ICV_Norm), color="purple", size=.2)
Fplot <- Mplot + geom_point(data=F_data, aes(x=Age,y=Total_Gray_ICV_Norm),color='dark green',size=.6) + geom_line(data=F_data, aes(x=Age,y=Total_Gray_ICV_Norm), color="dark green", size=.2) +
  geom_function(fun = ROIfnquad_Male, color="purple", size=2) +
  geom_function(fun = ROIfnquad_Fem, color="dark green", size=2) +
  theme_classic() +
  theme (plot.title = element_text(size=20, color="black",face="bold")) + 
  theme(axis.text.x=element_text(size=20, color="black")) + theme(axis.text.y=element_text(size=20, color="black")) + 
  theme(axis.title.x=element_blank()) + theme(axis.title.y=element_text(size=30, color="black")) +
  coord_cartesian(xlim=c(2, 8)) + scale_x_continuous(breaks=seq(2, 8, 1)) + ylab("Total Gray (%ICV)")
Fplot
