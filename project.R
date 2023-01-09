library(data.table)
library(zoo)
library(xts)
library(reshape)
library(mgcv)
#data daily:
df_daily_GBPJPY<-fread(file = 'STA202/projet/data.csv')  #importation de la base de donnée.
daily_Date<- strptime(df_daily_GBPJPY$Date,"%Y-%m-%d")    #
n_daily<-length(daily_Date)
t_daily<-c(1:n_daily)
daily_GBPJPY <- xts(df_daily_GBPJPY$Close,order.by = daily_Date )   #création de la série temporelle.
plot(daily_GBPJPY)   # affichage de l'évolution temporelle de la serie.
mean(daily_GBPJPY) # calcule de la moyenne.
sd(daily_GBPJPY)   # calcule de l'écart-type.
boxplot(daily_GBPJPY, ylab="GBPJPY",main="GBPJPY boxplot.")  
hist(daily_GBPJPY)   # affichage de l'histogramme de la serie temporelle.
pacf(daily_GBPJPY, lag.max = 100, plot = T) # calcule de l'auto-correlation jusqu'à l'ordre 300.

mean.year<-tapply(daily_GBPJPY, as.factor(strftime(daily_Date,"%y-%m")),mean)
mean.year<-xts(mean.year,order.by=strptime(rownames(mean.year),"%Y-%m"))
plot(mean.year, type='l' )

#################regression

reg<-lm(daily_GBPJPY~t_daily)  #
summary(reg)
ychap.lm <- reg$fitted
plot(as.numeric(reg$residuals), type='l')
ychap.lm<-xts(as.numeric(reg$fitted),order.by=daily_Date)
# ychap.lm<-xts(reg$fitted,order.by=Date)
plot(daily_GBPJPY,type='l')
lines(ychap.lm,col='red')
summary(ychap.lm)
#Moving average:
l <- 19
daily_Moving_average<-filter(daily_GBPJPY, filter=array(1/l,dim=l),
                               method = c("convolution"),
                               sides = 2, circular = F)
mb<-xts(daily_Moving_average,order.by=daily_Date)
lines(mb,col='blue')

#################noyau Gaussien
h=1000
x<-seq(1,max(t_daily),length=n_daily)
# test<-lapply(x,function(x){x^2})
# test[[1]]
# test
noyau <- function(x){dnorm(x-t_daily,0,sd=sqrt(h/2))/sum(dnorm(x-t_daily,0,sd=sqrt(h/2)))}

W<-matrix(unlist(lapply(x,noyau)),ncol=n_daily,nrow=n_daily,byrow=F)

plot(W[,50], type='l')
ychap.kernel<-colSums(as.numeric(daily_GBPJPY)*W)
ychap.kernel<-xts(ychap.kernel,order.by=daily_Date)
plot(daily_GBPJPY,type='l')
lines(ychap.kernel,col='red')
#################polynomes locaux
daily_lo<-loess(daily_GBPJPY~t_daily, degree=1, span=0.9)
ychap.daily_lo<-xts(daily_lo$fitted,order.by=daily_Date)
plot(daily_GBPJPY,type='l')
lines(ychap.daily_lo,col='red')
#################regression sur bases de splines
g<-gam(daily_GBPJPY~s(t_daily, k=3))
summary(g)
ychap.gam<-xts(g$fitted,order.by=daily_Date)
plot(daily_GBPJPY,type='l')
lines(ychap.gam,col='red')

###################################################estimation de la partie saisonniere
#################regression
daily_GBPJPY.detrend<-daily_GBPJPY-ychap.lm
plot(daily_GBPJPY.detrend)

acf(as.numeric(daily_GBPJPY.detrend), lag.max=100)

w=2*pi/50
fourier<-cbind(cos(w*t_daily), sin(w*t_daily))
K<-20
for(i in c(2:K))
{
  fourier<-cbind(fourier,cos(i*w*t_daily), sin(i*w*t_daily))
}
matplot(fourier[,1:10],type='l')
dim(fourier)


reg<-lm(daily_GBPJPY.detrend~fourier[,1:2]-1)
summary(reg)
ychap.lm.season<-xts(as.numeric(reg$fitted),order.by=daily_Date)
plot(daily_GBPJPY.detrend,type='l')
lines(ychap.lm.season,col='red', lwd=2)

lines(xts(S,,order.by=daily_Date), col='blue', lwd=2)


ExpSmooth = function(x,alpha)
{
  n = length(x)
  xsmooth = array(NA, dim=n) 
  
  # Initialization
  # We have to start somewhere: this would be optimized in R, using the Holt-Winters function.
  # Here, we arbitrarily choose a starting point, that is not absurd (at the beginning we do not
  # smooth the time series since we have no past).
  xsmooth[1] = x[1]
  
  for(i in 2:n)
  {
    xsmooth[i] = alpha*x[i] + (1-alpha)*xsmooth[i-1] 
  }
  
  return(xsmooth)
}

DoubleExpSmooth = function(x,alpha)
{
  n = length(x)
  xsmooth = array(NA, dim=n) 
  
  # Initialization
  # Again, it is here arbitrary but not absurd. l_t represents the level, thus it is in the same
  # rough size. b_t, instead, represents the trend: it is the same rough size and concept than
  # x[2]-x[1]. For l_1, I put x[2] to be at the "same point" than b_1, but I could have chosen
  # x[1] for example.
  # Anyway, usually these parameters are optimized, and we hope here that their influence will 
  # quickly be insignificant, after a warm up period.
  xsmooth[1] = x[1]
  l = array(NA, dim=n) 
  l[1] = x[2]
  b = array(NA, dim=n)
  b[1] = x[2]-x[1]
  
  for(i in 2:n)
  {
    l[i] = xsmooth[i-1] + (1-(1-alpha)^2)*(x[i]-xsmooth[i-1]) 
    # when h=1, xsmooth[i-1] = l[i-1]+b[i-1]
    # xsmooth[i-1] represents xchap_i/(i-1)
    b[i] = b[i-1] + alpha^2*(x[i]-xsmooth[i-1])
    xsmooth[i] = l[i] + b[i] # same comment than 2 lines before
  }
  
  # the following lines of code allow to return multiple objects, stored in a list
  # return() does not allow to return multiple object (e.g. return(b,l) is not possible)
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  return(res)
}

DoubleExpSmoothHW=function(x,alpha,beta)
{
  n = length(x)
  xsmooth = array(NA, dim=n)
  
  xsmooth[1] = x[1]
  l = array(NA, dim=n)
  l[1] = x[2]
  b = array(NA, dim=n)
  b[1] = x[2]-x[1]
  
  for(i in 2:n)
  {
    l[i] = alpha*x[i] + (1-alpha)*xsmooth[i-1]
    b[i] = beta*(l[i]-l[i-1]) + (1-beta)*b[i-1]
    xsmooth[i] = l[i] + b[i]
  }
  
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  return(res)
}

SeasonalDoubleExpSmooth = function(x,alpha,beta,delta,period)
{
  n = length(x)
  xsmooth = array(NA, dim=n)
  
  xsmooth[1] = x[1]
  l = array(NA, dim=n)
  l[1] = x[2] # x[1] ?
  b = array(NA, dim=n)
  b[1] = x[2]-x[1]
  s = array(NA, dim=n)
  s[1] = x[1]
  
  for(i in 2:n)
  {
    l[i] = alpha*(x[i]-s[max(i-period,1)]) + (1-alpha)*(l[i-1]+b[i-1])
    b[i] = beta*(l[i]-l[i-1]) + (1-beta)*b[i-1]
    s[i] = delta*(x[i]-l[i]) + (1-delta)*s[max(i-period,1)]
    xsmooth[i] = l[i] + b[i] + s[i]
  }
  
  res = list()
  res$smooth = xsmooth
  res$l = l
  res$b = b
  res$s = s
  return(res)
}

#### Error functions

mse = function(x,y){
  return(mean((x-y)^2))
}

# Since we are going to forecast, we want to optimize our choice of alpha in this framework,
# and not in a smoothing framework. Thus, we build functions to evaluate the MSE of a 
# 1-day ahead forecast, given a smoothing => we have to shift the smoothing vector.

mse_smoothing_simple = function(true_values, smoothed_values){
  n = length(true_values)
  # we compare what we should have forecasted to what we smooth!
  return(mse(tail(true_values,n-1), head(smoothed_values,n-1))) 
}

mse_smoothing_methods = function(true_values, smoothed_object){
  n = length(true_values)
  return(mse(tail(true_values,n-1), head(smoothed_object$smooth,n-1)))
}

#### Application to our time series: exponential
alpha = 0.2
daily_GBPJPY.smooth = xts(ExpSmooth(daily_GBPJPY,alpha),order.by = daily_Date)
plot(daily_GBPJPY, type='l', xlab="Time", ylab=expression(X^1))
lines(daily_GBPJPY.smooth, col='red', lwd=2)
legend("bottomleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, lwd=c(1,2), col=c("black","red"), cex=0.6)

mse_smoothing_simple(daily_GBPJPY, daily_GBPJPY.smooth)

alpha = seq(0.05,0.95,length=100)

forecast = lapply(alpha, ExpSmooth, x=daily_GBPJPY)

erreur = sapply(forecast, mse_smoothing_simple, true_values=df_daily_GBPJPY$Close)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

daily_GBPJPY.smooth = xts(ExpSmooth(daily_GBPJPY, alpha[which.min(erreur)]),order.by = daily_Date)
plot(daily_GBPJPY, type='l', xlab="Time", ylab=expression(X^1))
lines(daily_GBPJPY.smooth, col='red')
legend("topleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, col=c("black","red"), cex=0.6)


#### Application to our detrended time series: exponential
alpha = 0.2
daily_GBPJPY.detrend.smooth = xts(ExpSmooth(daily_GBPJPY.detrend,alpha),order.by = daily_Date)
plot(daily_GBPJPY.detrend, type='l', xlab="Time", ylab=expression(X^1))
lines(daily_GBPJPY.detrend.smooth, col='red', lwd=2)
legend("bottomleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, lwd=c(1,2), col=c("black","red"), cex=0.6)

mse_smoothing_simple(daily_GBPJPY.detrend, daily_GBPJPY.detrend.smooth)

alpha = seq(0.05,0.95,length=100)

forecast = lapply(alpha, ExpSmooth, x=daily_GBPJPY.detrend)

erreur = sapply(forecast, mse_smoothing_simple, true_values=df_daily_GBPJPY$Close)
plot(alpha,erreur,type='l',xlab=expression(alpha),ylab="Mean squared error")

daily_GBPJPY.smooth = xts(ExpSmooth(daily_GBPJPY.detrend, alpha[which.min(erreur)]),order.by = daily_Date)
plot(daily_GBPJPY.detrend, type='l', xlab="Time", ylab=expression(X^1))
lines(daily_GBPJPY.detrend.smooth, col='red')
legend("topleft", legend=c("Observed data", "Smoothed time series"),
       lty=1, col=c("black","red"), cex=0.6)


########################################################
####Double
########################################################
DoubleExpSmooth=function(x,alpha)
{
  xsmooth=x
  l<-array(x[1],dim=length(x))
  b<-array(x[2]-x[1],dim=length(x))
  
  for(i in c(2:length(x)))
  {
    l[i]<-xsmooth[i-1]+(1-(1-alpha)^2)*(x[i]-xsmooth[i-1])
    b[i]<-b[i-1]+alpha^2*(x[i]-xsmooth[i-1])
    xsmooth[i]<-l[i]+b[i]
  }
  
  res<-list()
  res$smooth<-xsmooth
  res$l=l
  res$b<-b
  return(res)
}


alpha<-seq(0.05,0.95,length=100)
forecast<-lapply(alpha,DoubleExpSmooth,x=daily_GBPJPY.detrend)
erreur<-unlist(
  lapply(forecast,
         function(x){mean((tail(daily_GBPJPY.detrend,n_daily-1)-head(x$smooth,n-1))^2)}))
plot(alpha,erreur,type='l')

daily_GBPJPY.detrend.smooth<-DoubleExpSmooth(daily_GBPJPY.detrend,alpha[which.min(erreur)])
plot(daily_GBPJPY.detrend,type='l')
lines(daily_GBPJPY.detrend.smooth$smooth,col='red')

plot(daily_GBPJPY.detrend.smooth$l,type='l',ylim=range(daily_GBPJPY.detrend.smooth$l,daily_GBPJPY.detrend.smooth$b),col='blue')
lines(daily_GBPJPY.detrend.smooth$b,col='red')



alpha<-seq(0.05,0.95,length=100)
forecast<-lapply(alpha,DoubleExpSmooth,x=X2)
erreur<-unlist(
  lapply(forecast,
         function(x){mean((tail(X2,n-1)-head(x$smooth,n-1))^2)}))
plot(alpha,erreur,type='l')

X2.smooth<-DoubleExpSmooth(X2,alpha[which.min(erreur)])
plot(X2,type='l')
lines(X2.smooth$smooth,col='red')

plot(X2.smooth$l,type='l',ylim=range(X2.smooth$l,X2.smooth$b),col='blue')
plot(X2.smooth$b,col='red', type='l')





#data monthly:
df_monthly_GBPJPY<-fread(file = 'STA202/projet/data1.csv')
monthly_Date<- strptime(df_monthly_GBPJPY$Date,"%Y-%m-%d")
n_monthly<-length(monthly_Date)
t_monthly<-c(1:n_monthly)
monthly_GBPJPY <- xts(df_monthly_GBPJPY$Close,order.by = monthly_Date )
plot(monthly_GBPJPY)
mean(monthly_GBPJPY)
sd(monthly_GBPJPY)
boxplot(monthly_GBPJPY)
hist(monthly_GBPJPY)
acf(monthly_GBPJPY, lag.max = 300, plot = T)
#Moving average:
l <- 50
monthly_Moving_average<-filter(monthly_GBPJPY, filter=array(1/l,dim=l),
           method = c("convolution"),
           sides = 2, circular = F)
mb<-xts(monthly_Moving_average,order.by=monthly_Date)
lines(mb,col='blue')

#################noyau Gaussien
h=1000
x<-seq(1,max(t_monthly),length=n_monthly)
# test<-lapply(x,function(x){x^2})
# test[[1]]
# test
noyau <- function(x){dnorm(x-t_monthly,0,sd=sqrt(h/2))/sum(dnorm(x-t_monthly,0,sd=sqrt(h/2)))}

W<-matrix(unlist(lapply(x,noyau)),ncol=n_monthly,nrow=n_monthly,byrow=F)

plot(W[,50], type='l')
ychap.kernel<-colSums(as.numeric(monthly_GBPJPY)*W)
ychap.kernel<-xts(ychap.kernel,order.by=monthly_Date)
plot(monthly_GBPJPY,type='l')
lines(ychap.kernel,col='red')

#################polynomes locaux
monthly_lo<-loess(monthly_GBPJPY~t_monthly, degree=1, span=0.9)
ychap.monthly_lo<-xts(monthly_lo$fitted,order.by=monthly_Date)
plot(monthly_GBPJPY,type='l')
lines(ychap.monthly_lo,col='red')

#################regression sur bases de splines:
g<-gam(monthly_GBPJPY~s(t_monthly, k=3))
summary(g)
ychap.gam<-xts(g$fitted,order.by=monthly_Date)
plot(monthly_GBPJPY,type='l')
lines(ychap.gam,col='red')






