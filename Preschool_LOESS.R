### Code adapted for the Preschool data set by Madison Long from the Actuarial Data Scientist's GitHub
### Last updated November 2021


### Load Packages
require(dplyr)
require(datasets)
require(stats)
require(ggplot2)

set.seed(5)

### Find Best Span
# Metric of Interest
calc.RMSE <- function(pred,act){
  sqrt(mean((pred-act)^2)/length(pred))
}

### Split Data into Train/Valid
train.split.pct <- 0.7
train.row <- runif(n=nrow(PSdata))  <  train.split.pct
train <- PSdata[train.row, ]
valid <- PSdata[!train.row,]
nrow(train)
nrow(valid)

# Hyperparameter Search
results <- matrix(NA,0,4)
colnames(results) <- c("span","degree","train.rmse","valid.rmse")

for(degree in seq(0,2,1)){
  for(span in seq(0.15,1,0.01)){
    mod <- loess( formula=Total_Gray~Age,
                  data=train,
                  span = span,
                  degree = degree,
                  normalize = TRUE,
                  control = loess.control(surface = "direct"))
    train.rmse <- calc.RMSE(mod$fitted,train$Total_Gray)
    valid.rmse <- calc.RMSE(predict(mod,newdata=valid) %>% 
                              as.numeric(),valid$Total_Gray)
    
    results <- rbind(results,c(span,degree,train.rmse,valid.rmse))
  }
}

results <- as.data.frame(results)
best    <- results[which.min(results$valid.rmse),]

plot <- ggplot(results, aes(span, degree, fill= valid.rmse)) +
  geom_tile()+
  ggtitle(paste("Span:",round(best$span,2),"   Degree:",best$degree)) +
  scale_fill_gradient(low="red",high="white")

final.mod <- loess( formula=Total_Gray~
                      Age,
                    data=PSdata,
                    span = best$span,
                    degree = best$degree,
                    normalize = TRUE,
                    control = loess.control(surface = "direct"))

plot_Tot_Gray <- ggplot(PSdata, aes(x=Age, y=Total_Gray)) + geom_point(color="darkgrey")
plot_Tot_Gray + geom_line(aes(y=predict(final.mod)))
print(summary(final.mod))

#########################################################

results <- matrix(NA,0,4)
colnames(results) <- c("span","degree","train.rmse","valid.rmse")

for(degree in seq(0,2,1)){
  for(span in seq(0.15,1,0.01)){
    mod <- loess( formula=GM.ICV~Age,
                  data=train,
                  span = span,
                  degree = degree,
                  normalize = TRUE,
                  control = loess.control(surface = "direct"))
    train.rmse <- calc.RMSE(mod$fitted,train$GM.ICV)
    valid.rmse <- calc.RMSE(predict(mod,newdata=valid) %>% 
                              as.numeric(),valid$GM.ICV)
    
    results <- rbind(results,c(span,degree,train.rmse,valid.rmse))
  }
}

results <- as.data.frame(results)
best    <- results[which.min(results$valid.rmse),]

plot <- ggplot(results, aes(span, degree, fill= valid.rmse)) +
  geom_tile()+
  ggtitle(paste("Span:",round(best$span,2),"   Degree:",best$degree)) +
  scale_fill_gradient(low="red",high="white")

final.mod_GM.ICV <- loess( formula=GM.ICV~
                      Age,
                    data=PSdata,
                    span = best$span,
                    degree = best$degree,
                    normalize = TRUE,
                    control = loess.control(surface = "direct"))

plot_GM.ICV <- ggplot(PSdata, aes(x=Age, y=GM.ICV)) + geom_point(color="darkgrey")
plot_GM.ICV + geom_line(aes(y=predict(final.mod_GM.ICV)))
print(summary(final.mod_GM.ICV))


#########################################################
results <- matrix(NA,0,4)
colnames(results) <- c("span","degree","train.rmse","valid.rmse")

for(degree in seq(0,2,1)){
  for(span in seq(0.15,1,0.01)){
    mod <- loess( formula=Right_Pallidum~Age,
                  data=train,
                  span = span,
                  degree = degree,
                  normalize = TRUE,
                  control = loess.control(surface = "direct"))
    train.rmse <- calc.RMSE(mod$fitted,train$Right_Pallidum)
    valid.rmse <- calc.RMSE(predict(mod,newdata=valid) %>% 
                              as.numeric(),valid$Right_Pallidum)
    
    results <- rbind(results,c(span,degree,train.rmse,valid.rmse))
  }
}

results <- as.data.frame(results)
best    <- results[which.min(results$valid.rmse),]

plot <- ggplot(results, aes(span, degree, fill= valid.rmse)) +
  geom_tile()+
  ggtitle(paste("Span:",round(best$span,2),"   Degree:",best$degree)) +
  scale_fill_gradient(low="red",high="white")

final.mod_Right_Pallidum <- loess( formula=Right_Pallidum~
                             Age,
                           data=PSdata,
                           span = best$span,
                           degree = best$degree,
                           normalize = TRUE,
                           control = loess.control(surface = "direct"))

plot_Right_Pallidum <- ggplot(PSdata, aes(x=Age, y=Right_Pallidum)) + geom_point(color="darkgrey")
plot_Right_Pallidum + geom_line(aes(y=predict(final.mod_Right_Pallidum)))
print(summary(final.mod_Right_Pallidum))


#########################################################

results <- matrix(NA,0,4)
colnames(results) <- c("span","degree","train.rmse","valid.rmse")

for(degree in seq(0,2,1)){
  for(span in seq(0.15,1,0.01)){
    mod <- loess( formula=ICV_Norm_Right_Pallidum~Age,
                  data=train,
                  span = span,
                  degree = degree,
                  normalize = TRUE,
                  control = loess.control(surface = "direct"))
    train.rmse <- calc.RMSE(mod$fitted,train$ICV_Norm_Right_Pallidum)
    valid.rmse <- calc.RMSE(predict(mod,newdata=valid) %>% 
                              as.numeric(),valid$ICV_Norm_Right_Pallidum)
    
    results <- rbind(results,c(span,degree,train.rmse,valid.rmse))
  }
}

results <- as.data.frame(results)
best    <- results[which.min(results$valid.rmse),]

plot <- ggplot(results, aes(span, degree, fill= valid.rmse)) +
  geom_tile()+
  ggtitle(paste("Span:",round(best$span,2),"   Degree:",best$degree)) +
  scale_fill_gradient(low="red",high="white")

final.mod_ICV_Norm_Right_Pallidum <- loess( formula=ICV_Norm_Right_Pallidum~
                                     Age,
                                   data=PSdata,
                                   span = best$span,
                                   degree = best$degree,
                                   normalize = TRUE,
                                   control = loess.control(surface = "direct"))

plot_ICV_Norm_Right_Pallidum <- ggplot(PSdata, aes(x=Age, y=ICV_Norm_Right_Pallidum)) + geom_point(color="darkgrey")
plot_ICV_Norm_Right_Pallidum + geom_line(aes(y=predict(final.mod_ICV_Norm_Right_Pallidum)))
print(summary(final.mod_ICV_Norm_Right_Pallidum))

#########################################################

results <- matrix(NA,0,4)
colnames(results) <- c("span","degree","train.rmse","valid.rmse")

for(degree in seq(0,2,1)){
  for(span in seq(0.15,1,0.01)){
    mod <- loess( formula=Left_Accumbens_Area~Age,
                  data=train,
                  span = span,
                  degree = degree,
                  normalize = TRUE,
                  control = loess.control(surface = "direct"))
    train.rmse <- calc.RMSE(mod$fitted,train$Left_Accumbens_Area)
    valid.rmse <- calc.RMSE(predict(mod,newdata=valid) %>% 
                              as.numeric(),valid$Left_Accumbens_Area)
    
    results <- rbind(results,c(span,degree,train.rmse,valid.rmse))
  }
}

results <- as.data.frame(results)
best    <- results[which.min(results$valid.rmse),]

plot <- ggplot(results, aes(span, degree, fill= valid.rmse)) +
  geom_tile()+
  ggtitle(paste("Span:",round(best$span,2),"   Degree:",best$degree)) +
  scale_fill_gradient(low="red",high="white")

final.mod_Left_Accumbens_Area <- loess( formula=Left_Accumbens_Area~
                                              Age,
                                            data=PSdata,
                                            span = best$span,
                                            degree = best$degree,
                                            normalize = TRUE,
                                            control = loess.control(surface = "direct"))

plot_Left_Accumbens_Area <- ggplot(PSdata, aes(x=Age, y=Left_Accumbens_Area)) + geom_point(color="darkgrey")
plot_Left_Accumbens_Area + geom_line(aes(y=predict(final.mod_Left_Accumbens_Area)))
print(summary(final.mod_Left_Accumbens_Area))

#########################################################


#########################################################





                          