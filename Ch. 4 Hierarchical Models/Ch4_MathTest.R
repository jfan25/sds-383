library(ggplot2)
library(plyr)
library(mvtnorm)
library(LaplacesDemon)
library(purrr)
#import the data
mathscores <- read.csv(file="C:/Users/15202/OneDrive/Documents/stat modeling 2/Ch 4/mathtest.csv", header=TRUE, sep=",")
y = matrix(0,nrow = 100, ncol = max(count(mathscores$school,1)$freq)) 
count = count(mathscores$school, 1)$freq
track = 0
y = list()
for(i in 1:100)
{
  y[[i]] = c(mathscores$mathscore[(track+1):(track+count[i])])
  track = track+count[i]
}
y_avg = lapply(y, mean)

##########Math test hierachical gibbs model
#set the hyper parameters
a = 1
b = 1
c = 1
d = 1
gam = 0.01

#initialize the other parameters
N = 1000
mu = rep(0,N)
mu_o = mean(unlist(y_avg))
mu[1] = mean(unlist(y_avg))
lambda = rep(0,N)
lambda[1] = rgamma(1,1,1)
nu = rep(0,N)
nu[1] = rgamma(1,1,1)
theta = matrix(0,nrow = N, ncol = 100)
theta[1,] = rnorm(100, mean = mu[1], sd = sqrt(1/(nu[1]*lambda[1])))



#The gibbs sampling portion
for(i in 2:N){
  #update theta
  theta_sd = sqrt(1/(count*lambda[i-1]+nu[i-1]*lambda[i-1]))
  #theta_sd
  theta_mean = (lambda[i-1] * matrix(unlist(lapply(y,sum)), nrow = 100, ncol = 1) + nu[i-1]*lambda[i-1]*mu[i-1])*theta_sd^2
  #theta_mean
  theta[i,] = rnorm(100, mean = theta_mean, sd = theta_sd)
  # print(theta[i,])
  #update lambda
  lam_shape = (sum(count)+length(y))/2 + a
  total = 0
  for(j in 1:100){total = total + sum((y[[j]] - theta[i,j])^2)}
  lam_rate = 1/2*total + nu[i-1]/2*sum((theta[i,]-mu[i-1])^2) + b
  lambda[i] = rgamma(1, shape = lam_shape, rate = lam_rate)
  #update nu
  nu_shape = length(y)/2 + c
  nu_rate = lambda[i]/2 * sum((theta[i,]-mu[i-1])^2) + d
  nu[i] = rgamma(1, shape = nu_shape, rate = nu_rate)
  #nu[i]
  #update mu
  mu_sd = sqrt(1/(nu[i]*lambda[i]*100 + gam))
  mu_mean = (nu[i]*lambda[i]*sum(theta[i,])+gam*mu_o)*mu_sd^2
  mu[i]=rnorm(1, mu_mean, mu_sd)
}

######plot the shrinkage coefficient
post_theta = colMeans(theta[100:N,])#t(apply(theta[100:N, ],1, function(x) x + unlist(y_avg)))

kappa = abs((unlist(y_avg) - post_theta)/unlist(y_avg))

shrinkage = data.frame(kappa, y_avg, count)
plot(count,kappa, main = "Shrinkage as a function of sample size", xlab = "sample size",  ylab = "shrinkage")

###Trace Plots
par(mfrow = c(2,2))
plot(seq(1,N, by = 1),mu,type = 'n', main = "Trace plot of Mu", xlab = 'iterations')
lines(seq(1,N, by = 1), mu, col= 'sky blue')

plot(seq(500,N, by = 1),lambda[500:N],type = 'n', main = "Trace plot of Lambda", xlab = 'iterations')
lines(seq(500,N, by = 1), lambda[500:N], col= 'sky blue')

plot(seq(1,N, by = 1),nu,type = 'n', main = "Trace plot of Nu", xlab = 'iterations')
lines(seq(1,N, by = 1), nu, col= 'sky blue')


plot(seq(1,N, by = 1),theta[,1],type = 'n', main = "Trace plot of theta[1]", xlab = 'iterations')
lines(seq(1,N, by = 1), theta[,1], col= 'sky blue')

#compute the mean of the posterior draws
mean(mu[100:N])
mean(lambda[100:N])
mean(nu[100:N])