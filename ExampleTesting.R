###########################################
# load the data set and functions
###########################################
load('hourly_data.RData')
attach(hourly_data)

source("methods.R")
source("Core_functions.R")

###########################################
# Produce plots of fitted values
###########################################
red <- rgb(170/255, 39/255, 47/255)
yellow <- rgb(253/255, 200/255, 47/255,1)
fush <- rgb(199/255, 67/255, 117/255)
#yellow <- rgb(180/255, 130/255, 40/255,1)

series <- GetSeries(24*365*1,24*7*1)

par(mfrow=c(1,1))
prediction1 <- method1(series,1)
score <- log(ErrorFunction(prediction$real, prediction$forecast))
#mean(abs(prediction$real-prediction$forecast))
plot(cr(prediction1$real), pch=20,col=yellow,ylab = "", xlab ="", main = "Traditional forecasting",type='l',lwd=3,cex.main=2)
points(cr(prediction1$forecast), pch=20, col=red, lwd =3)
lines(cr(series[[4]]), pch=20, lwd =3,ylab = "", xlab ="", main = "Temperature")
points(cr(prediction3$forecast), pch=20, col='blue', lwd =3)


prediction3 <- method3(series,1,1)
score <- log(ErrorFunction(prediction$real, prediction$forecast))
#mean(abs(prediction$real-prediction$forecast))
plot(cr(prediction$real), pch=20,col=yellow,ylab = "", xlab ="", main = "Conditional forecasting",type='l',lwd=3,cex.main=2)
points(cr(prediction$forecast), pch=20, col=red, lwd =3)
lines(cr(series[[4]]), pch=20, lwd =3,ylab = "", xlab ="", main = "Temperature")





###########################################
# Test Loop (box plots with errors)
###########################################
ntest <- 200
ntech <- 2
ErrorMatrix <- matrix(data=rep(0,(ntest*ntech)),nrow=ntest)

for (i in 1:ntest){
  
  series <- GetSeries(24*365*3,24*7*1)
  
  prediction <- method1(series,1)
  ErrorMatrix[i,1] <- mean(abs(prediction$real-prediction$forecast))
  
  
  #prediction <- method2(series,1,0.5)
  #ErrorMatrix[i,2] <-mean(abs(prediction$real-prediction$forecast))
  

  prediction <- method3(series,1,1)
  ErrorMatrix[i,2] <-mean(abs(prediction$real-prediction$forecast))
  

  #prediction <- method4(series,1,0.5)
  #ErrorMatrix[i,3] <-mean(abs(prediction$real-prediction$forecast))
  
}

Y <- log(as.vector(ErrorMatrix))
Group <- c(rep('Traditional forecasting',ntest),rep('Conditional forecasting',ntest))
#NAVV <- rep(NAV,Tech)
Group <- factor(Group,levels=unique(Group))

ToPlot <- data.frame(as.factor(Group),Y)

Summary <- ddply(.data=ToPlot, c("Group"), summarise,
                 Y.mean = mean(Y),
                 Y.sd   = sd(Y))

Plot <- ggplot(ToPlot, aes(Group, Y)) + geom_boxplot(color=red,lwd=1.5)  +  geom_sina(shape=21, fill=yellow, color='black', size=3) 

Plot <- Plot +scale_y_continuous(name="Mean absolute value error (logarithm)",limits=c(6,9))
Plot <- Plot + labs(x="Techniques")
Plot <- Plot + coord_flip()
Plot <- Plot  + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
#        panel.border = element_blank(),
 #       panel.background = element_blank(),
        text = element_text(size = 17)) 


Plot