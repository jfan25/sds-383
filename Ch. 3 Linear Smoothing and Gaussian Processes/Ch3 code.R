library(psych)
source("HW_3_functions.R")
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
