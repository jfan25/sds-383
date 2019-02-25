install.packages('mvtnorm')
library(mvtnorm)
data <- read.csv(file="/home/jf37828/Documents/StatModeling2/gdpgrowth.csv", header=TRUE, sep=",")
###Conjugate Gaussian Error model, basics, part d
#Constructing the data vectors
x = cbind(rep(1,length(data$DEF60)),data$DEF60)
y = data$GR6096

#Initializing the hyperparameters
lambda = diag(1,length(y))
k = diag(c(.01,.01))
m = c(0.02, 0.1)

#solving for the posterior parameters of p(beta|y)
A = t(x)%*% x + k
b = y %*% lambda %*% x + m%*%k
postBetaMean = solve(A)%*% t(b)

#fitting a line using posterior mean of Beta
y_hat = x%*%postBetaMean
  
plot(x[,2],y,main = "Bayesian Linear Regression plot")
lines(x[,2],y_hat)

############ Ch 2 Heavy Tail Error sampling, Gibbs Sampling, part c
#Setting the data
x = cbind(rep(1,length(data$DEF60)),data$DEF60)
y = data$GR6096

# Initializing hyperparameters
n = 1000 #number of samples to keep
h = 2
d = 2
eta = 2
m = c(0.02, 0.1)
k = diag(c(.01,.01))

#Creating data structures for each parameter
beta_samples = array(0,dim = c(2,1,n,1))
lambda_samples = array(0,dim = c(length(y), length(y),n,1))
omega_samples = rep(0,n)

#creating temporary variables
temp_beta = matrix(0, nrow = 2, ncol = 1)
temp_lambda = rep(0, length(y))
temp_omega = NULL

# Initializing values for each parameter

temp_omega = rgamma(1,d/2, eta/2)
temp_lambda = rgamma(length(y),h/2,h/2)
temp_beta = rmvnorm(1,mean = m, sigma = solve(temp_omega * k))
count = 0
#solving for the posterior parameters of p(beta|y)
for(i in 1:4000)
{
  #solving for the posterior parameters of p(beta|y, omega, lambda)  
    A = t(x)%*% x + k
    b = t(y %*% diag(temp_lambda) %*% x + m%*%k)
    postBetaMean = solve(A)%*% b
    postBetaVar = solve(temp_omega * A)
  #draw new Beta
    temp_beta = rmvnorm(1, mean = postBetaMean, sigma = postBetaVar)
    temp_beta = t(temp_beta)
  #calculate parameters for new omega
    postRate = 1/2*(t(temp_beta-solve(A)%*%b)%*% A %*% (temp_beta-solve(A)%*%b) + 
                      eta + t(b)%*% solve(A) %*% b + t(y)%*% diag(temp_lambda)%*%y + 
                      t(m)%*% k %*%m) 
  #draw new Omega
    temp_omega = rgamma(1,(d+length(y)+2)/2, postRate)
  #calculate parameters for new lambda
    lamRate = 1/2*(temp_omega*(y-x%*%temp_beta)^2 + h)
  #draw new lambdas
    temp_lambda = rgamma(lamRate, (h+1)/2, rate = lamRate)
  #saving the samples after burn in and with thinning
  if(i>1000 && i%%3==0){
    count = count + 1
    beta_samples[,,count,] = temp_beta
    lambda_samples[,,count,] = diag(temp_lambda)
    omega_samples[count] = temp_omega
  }
}

#Computing the mean of the beta samples
meanBeta0 = mean(beta_samples[1,,,])
meanBeta1 = mean(beta_samples[2,,,])

meanBeta = c(meanBeta0, meanBeta1)

y_hat_mean = x%*%meanBeta #calculting a prediction using the beta means

#Finding the mode of the beta samples
h <- hist(beta_samples[1,,,], breaks=200)
b_0mode <- max(h$mids[h$counts == max(h$counts)])

h <- hist(beta_samples[2,,,], breaks=200)
b_1mode <- h$mids[h$counts == max(h$counts)]

modeBeta = c(b_0mode,b_1mode)
y_hat_mode = x%*%modeBeta #calculating a prediction using the beta MAP

#Plot the original estimate against the posterior mean and the mode
#from the Gibbs sampling
plot(x[,2],y,main = "Bayesian Linear Regression plot")
lines(x[,2],y_hat)
lines(x[,2],y_hat_mean, col = 'red')
lines(x[,2], y_hat_mode, col = 'blue')
legend(0.075,0.005, legend = c("Posterior Mean", "Gibbs Post Mean", "Gibbs Post MAP"), 
       col = c("black", "red", "blue"), lty = 1)
