library(psych)
library(mvtnorm)
source("HW_3_functions.R")
#utilities <- read.csv(file="/home/jf37828/Documents/StatModeling2/utilities.csv", header=TRUE, sep=",")
utilities <- read.csv(file="c:/Users/15202/OneDrive/Documents/stat modeling 2/Ch 3/utilities.txt", header=TRUE, sep=",")


###############Kernel Smoothing
#test function creation------
x <- seq(-2*pi,2*pi, by=0.1)
y <- cos(x)
y <- AddUnifNoise(y,.3)
#----------------------------

x_new <- seq(-6.5,6.5,length.out = length(x)) #defining new x grid for fitting values
h <- .02 #the bandwidth

x <- x-mean(x) 
y<-y-mean(y)

y_new <- KernPredict(x,x_new,y,GaussKernel,h) #fitting the new values using kernel smoothing


plot(x,y, main = 'h=0.2') #plots
lines(x_new, y_new, lw=2, col = "red")

#######Cross Validation
#test function creation------------------------------  
x_obs <- matrix(seq(-6,2, length.out = 500))
x_new = matrix(x_obs)
y <-  cos(4*x_obs)# x_obs^2 #
y_obs <- matrix(AddUnifNoise(y,.2),nrow = length(y))
#----------------------------------------------------
h <- seq(0.01,2,by=0.01) #bandwidth candidate vector

h_opt = crossVal(x_obs, y_obs, h, GaussKernel) #returns optimal h value and plots the optimal reconstruction

hat <- makeHat(x_obs, x_obs, GaussKernel,h_opt) #creating the hat matrix
y_fit = hat %*% y_obs #calculating fitted values using the hat matrix

plot(x_obs,y_obs)
lines(x_obs, y_fit)

######Local linear regression with dataset

#creating data structure--------------------------
y_obs <- utilities$gasbill/utilities$billingdays
x_obs <- matrix(utilities$temp)
data <- data.frame(x_obs, y_obs)
data <- data[order(data$x), ]
#-------------------------------------------------

h <- seq(1,10,by=0.01) #bandwidth candidate vector
d = 1 #degree of polynomial regression

h_opt = LOOVC(matrix(data$x_obs), data$y_obs, d, GaussKernel, h) #finding optimal bandwidth using LOOCV

hat <- makeHat(data$x_obs, matrix(data$x_obs), GaussKernel, h_opt) #creating the hat matrix using the optimal bandwidth
y_fit = hat %*% data$y_obs #finding fitted values using the hat matrix

sig_estimate = norm(matrix(data$y_obs-y_fit), "2")/(length(x_obs) - 2*tr(hat) + tr(hat%*% hat))#estimating the error variance
est_variance = sig_estimate * apply(hat, 1, function(x) norm(matrix(x),"2")) #calculating the variance of each estimate
conf_int = matrix(NA, nrow =2, ncol = length(y_fit)) #creating the pointwise 95% confidence interval
conf_int[1,] = y_fit + 1.96*sqrt(est_variance)
conf_int[2,] = y_fit - 1.96*sqrt(est_variance)
#Plots--------------------------------------------
plot(data$x_obs,data$y_obs, xlab = "Temperature", ylab = "Average Daily Gas Bill", main = "Local Polynomial Regression Fit")
lines(data$x_obs, y_fit)

plot(data$x_obs, data$y_obs-y_fit, main = "Residuals", xlab = "temp", ylab = "residuals")
plot(data$x_obs, log((data$y_obs-y_fit)^2), main = "Log squared residuals", ylab = "log squared residual", xlab = "temp")


plot(data$x_obs,data$y_obs, xlab = "Temperature", ylab = "Average Daily Gas Bill", main = "Local Polynomial Regression Fit with 95% C.I.")
lines(data$x_obs, y_fit, col = "red", type = "o")
lines(data$x_obs,conf_int[1,],lty =3)
lines(data$x_obs, conf_int[2,], lty = 3)

###Fitting the model using log(y) to resolve heteroskedascity
y_obs <- utilities$gasbill/utilities$billingdays
y_obs = log(y_obs)
x_obs <- matrix(utilities$temp)
data <- data.frame(x_obs, y_obs)
data <- data[order(data$x), ]

h <- seq(1,10,by=0.01)
d = 1

h_opt = LOOVC(matrix(data$x_obs), data$y_obs, d, GaussKernel, h)

hat <- makeHat(data$x_obs, matrix(data$x_obs), GaussKernel, h_opt)
y_fit = hat %*% data$y_obs

plot(data$x_obs,data$y_obs, xlab = "Temperature", ylab = "Average Daily Gas Bill", main = "Local Polynomial Regression Log Fit")
lines(data$x_obs, y_fit)

plot(data$x_obs, data$y_obs-y_fit, main = "Residuals for Log Fit", xlab = "temp", ylab = "residuals")


####Fitting reweighted model using estimates of sigma via regression
# y_obs <- utilities$gasbill/utilities$billingdays
# x_obs <- matrix(utilities$temp)
# data <- data.frame(x_obs, y_obs)
# data <- data[order(data$x), ]
# 
# h <- seq(1,10,by=0.01)
# d = 1
# 
# h_opt = LOOVC(matrix(data$x_obs), data$y_obs, d, GaussKernel, h)
# 
# hat <- makeHat(data$x_obs, matrix(data$x_obs), GaussKernel, h_opt)
# y_fit = hat %*% data$y_obs
# 
# res = data$y_obs - y_fit
# logres = log(res^2)
# res_reg = lm(logres ~ log(data$x_obs))
# res_fit = exp(cbind(rep(1,length(data$x_obs)),log(data$x_obs)) %*% matrix(res_reg$coefficients, nrow = 2))
# 
# y_obs = y_obs/sqrt(res_fit)
# x_obs = x_obs/sqrt(res_fit)
# data <- data.frame(x_obs, y_obs)
# data <- data[order(data$x), ]
# 
# h <- seq(1,10,by=0.01)
# d = 1
# 
# h_opt = LOOVC(matrix(data$x_obs), data$y_obs, d, GaussKernel, h)
# 
# hat <- makeHat(data$x_obs, matrix(data$x_obs), GaussKernel, h_opt)
# y_fit = hat %*% data$y_obs
# 
# plot(data$x_obs,data$y_obs, xlab = "Temperature", ylab = "Average Daily Gas Bill", main = "Local Polynomial Regression Fit")
# lines(data$x_obs, y_fit)
# 
# plot(data$x_obs, data$y_obs-y_fit, main = "Residuals", xlab = "temp", ylab = "residuals")

##########Gaussian Process sampling ##################################################
par(mfrow=c(1,3))
#creating grid
x = seq(0,1,length.out = 100)

#initializing hyperparameters
mu = rep(0,length(x))
b=1
t1 = matrix(c(0.1, 1, 3), nrow = 1)
t2= 0.01

#calculating covariance matrices for SE and Matern 52 kernels
for(j in 1:3){
  covSE = getCov(x,maternSE,b,t1[j],t2)
  cov52 = getCov(x,matern52,b,t1[j],t2)
  colors= rainbow(3)
  samples = rmvnorm(3,mean = mu, cov52)
  plot(x,samples[1,], type = 'n', ylim = c(-5,5), main = paste('Matern52, b=',b, ', t1 = ', t1[j], ', t2 =', t2))
  
  for(i in 1:3){
    lines(x,samples[i,], col = colors[i])
  }
}


######## Nonparametric Regression and Spatial Smoothing###############################
#initializing data
y_obs <- utilities$gasbill/utilities$billingdays
x_obs <- matrix(utilities$temp)
data <- data.frame(x_obs, y_obs)
data <- data[order(data$x), ]
#initializing hyperparameters
sigma = 1
t1 = 10
t2 = 0.001
b = 20

#calculating posterior expected values at new points
postExp = apply(matrix(data$x_obs, nrow = 1), 2, function(x) postMean(x, data$x_obs, data$y_obs, sigma, maternSE, b, t1, t2))
#calculating posterior pointwise confidence interval
postVariance = apply(matrix(data$x_obs, nrow = 1), 2, function(x) postVar(x, data$x_obs, sigma, maternSE, b, t1, t2))

conf_int = matrix(NA, nrow =2, ncol = length(postExp)) #creating the pointwise 95% confidence interval
conf_int[1,] = postExp + 1.96*sqrt(postVariance)
conf_int[2,] = postExp - 1.96*sqrt(postVariance)

#plot of posterior expected value with 95% confidence interval
plot(data$x_obs, data$y_obs, main = "Posterior Expected Values with 95% Pointwise CI", ylab = "Average Gas Bill", xlab = 'Temperature')
lines(data$x_obs, postExp)
lines(data$x_obs,conf_int[1,],lty =3)
lines(data$x_obs, conf_int[2,], lty = 3)