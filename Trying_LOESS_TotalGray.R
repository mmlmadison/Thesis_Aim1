### Code adapted for the Preschool data set by Madison Long from the Actuarial Data Scientist's GitHub
### Last updated November 2021


### Load Packages
require(dplyr)
require(datasets)
require(stats)
require(ggplot2)


### Split Data into Train/Valid
train.split.pct <- 0.7
train.row <- runif(n=nrow(PSdata))  <  train.split.pct
train <- PSdata[train.row, ]
valid <- PSdata[!train.row,]
nrow(train)
nrow(valid)



### LOESS function (testing it out)

mod <- loess( formula=Total_Gray ~ Age,
              data=train,
              span = 0.3,
              degree = 1,
              normalize = TRUE)


### Find Best Span
# Metric of Interest
calc.RMSE <- function(pred,act){
  sqrt(mean((pred-act)^2)/length(pred))
}
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

                          