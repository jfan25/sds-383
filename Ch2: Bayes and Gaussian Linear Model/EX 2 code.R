install.packages('mvtnorm')
library(mvtnorm)

data <- read.csv(file="/home/jf37828/Documents/StatModeling2/gdpgrowth.csv", header=TRUE, sep=",")
###Conjugate Gaussian Error model, basics, part d
#Constructing the data vectors
x = cbind(rep(1,length(data$DEF60)),data$DEF60)
y = data$GR6096

#Initializing the hyperparameters
lambda = diag(1,length(y))
k_uninformative = diag(c(0.00001,0.00001)) #an uninformative prior precision
k_informative = diag(c(.1,.1))#diag(c(0.1, 0.1)) #an informative prior precision
m = c(.1,.1)


#OLS solution for beta
beta_hat = solve(t(x)%*% x)%*% t(x)%*%y
y_hat_ols = x%*% beta_hat

#solving for the posterior parameters of p(beta|y) using uninformative k
A = t(x)%*% x + k_uninformative
b = t(x) %*% lambda %*% y + k_uninformative%*%m
uninfBetaMean = solve(A)%*% b

#solving for the posterior parameters of p(beta|y) using informative k
A = t(x)%*% x + k_informative
b = t(x) %*% lambda %*% y + k_informative%*%m
infBetaMean = solve(A)%*% b

#fitting a line using posterior mean of Beta
y_hat_uninf = x%*%uninfBetaMean
y_hat_inf = x%*% infBetaMean

plot(x[,2],y,main = "Bayesian Linear Regression plot", ylab = "GDP", xlab = "Defense spending")
lines(x[,2],y_hat_ols, col = 'green', type = "o")
lines(x[,2], y_hat_uninf, col = 'red')
lines(x[,2], y_hat_inf, col = 'blue')

legend(0.075,0.005, legend = c("OLS", "Uninformative posterior", "Informative posterior"), 
       col = c("green", "red", "blue"), lty = 1)


############ Ch 2 Heavy Tail Error sampling, Gibbs Sampling, part c
#Setting the data
x = cbind(rep(1,length(data$DEF60)),data$DEF60)
y = data$GR6096

# Initializing hyperparameters
n = 2000 #number of samples to keep
N = 2000 + n*3
h = 10
d = 10
eta = 10
m = c(.4,.4)
k = diag(c(.01,.01))#diag(c(.2,.2))



#Creating data structures for each parameter
beta_samples = array(0,dim = c(2,1,n,1))
lambda_samples = array(0,dim = c(1, length(y),n,1))
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
set.seed(100)
for(i in 1:N)
{
  #solving for the posterior parameters of p(beta|y, omega, lambda)  
  A = t(x)%*% x + k
  #b = t(y %*% diag(temp_lambda) %*% x + m%*%k)
  b = t(x) %*% diag(temp_lambda) %*% y + k%*%m
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
  temp_omega = rgamma(1,(d+length(y))/2, rate = postRate)
  #calculate parameters for new lambda
  lamRate = 1/2*(temp_omega*(y-x%*%temp_beta)^2 + h)
  #draw new lambdas
  temp_lambda = rgamma(lamRate, (h+1)/2, rate = lamRate)
  #saving the samples after burn in and with thinning
  if(i>2000 && i%%3==0){
    count = count + 1
    beta_samples[,,count,] = temp_beta
    lambda_samples[,,count,] = temp_lambda
    omega_samples[count] = temp_omega
  }
}

#Computing the mean of the beta samples
meanBeta0 = mean(beta_samples[1,,,])
meanBeta1 = mean(beta_samples[2,,,])

meanBeta = c(meanBeta0, meanBeta1)

y_hat_mean = x%*%meanBeta #calculting a prediction using the beta means


#Plot the original estimate against the posterior mean 
#from the Gibbs sampling
plot(x[,2],y,main = "Bayesian Linear Regression plot", ylab = "GDP", xlab = "Defense spending")
lines(x[,2],y_hat_inf, col = "black") 
lines(x[,2],y_hat_mean, col = 'red')
#lines(x[,2], y_hat_mode, col = 'blue')
legend(0.075,0.005, legend = c("Posterior Mean", "Heavy Tail Error"), 
       col = c("black", "red"), lty = 1)
