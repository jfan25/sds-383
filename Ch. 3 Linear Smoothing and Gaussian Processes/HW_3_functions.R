AddUnifNoise <- function(y, nl){
  #-------------------------------------
  ##Usage: This function adds uniform noise to a given input vector
  #-------------------------------------
  ##Inputs: y - vector - the vector you'd like to add noise to 
  #         nl - float - noise level. The noise is generated from uniform(-nl,nl)
  #-------------------------------------
  ##Output: Noisy data vector
  #-------------------------------------
  
  out = y + runif(y,-nl,nl)
  return(out)
}
GaussKernel<- function(distance, h){
  #-------------------------------------
  ##Usage: a guassian kernel smoother
  #-------------------------------------
  ##Inputs: distance - vector - the distance between the observed x values and 
  #                             target x values
  #         h - float - the bandwidth for the kernel
  #-------------------------------------
  ##Output: kernel weights
  #-------------------------------------
  #output = 1/h * (1/sqrt(2*pi))*exp(-(distance/h)^2/2)
  output = 1/h * dnorm(distance/h,0,1)
  return(output)
}

KernWeights <- function(x_obs, x_target, kernel.fun, bandwidth){
  #-------------------------------------
  ##Usage: calculates the weights for kernel regression
  #-------------------------------------
  ##Inputs: x_obs - vector- x values for the observed data
  #         x_target - float- x values for which you would like to fit y_hat values
  #         kernel.fun - a kernel function. e.g. GaussKernel
  #         bandwidth - the bandwidth for the kernel function
  #-------------------------------------
  ##Output: weights for kernel regression
  #-------------------------------------
  diff = x_obs - x_target
  weight = kernel.fun(diff,bandwidth)
  weight = weight/sum(weight)
  return(weight)
}

KernPredict <- function(x_obs,x_target,y_obs,kernel.fun,bandwidth){
  #-------------------------------------
  ##Usage: produces fitted values for for a kernel regression
  #-------------------------------------
  ##Inputs: y - the vector you'd like to add noise to 
  #         nl - noise level. The noise is generated from uniform(-nl,nl)
  #-------------------------------------
  ##Output: Noisy data vector
  #-------------------------------------
  pred_y = rep(0,length(x_target))
  for(i in 1:length(x_target)){
    w = KernWeights(x_obs,x_target[i],kernel.fun,bandwidth)
   # w = w/sum(w)
    pred_y[i] = sum(w*y_obs)
  }
  return(pred_y)
}


crossVal<- function(x_obs, y_obs, h_test, kernel.fun){
  #-------------------------------------
  ##Usage: Finds the optimal bandwidth for a given kernal regressor using 2 group cross validation
  #-------------------------------------
  ##Inputs: x_obs - vector - observed x values
  #         y_obs - vector - observed y values
  #         h_test - vector- a vector of candidate bandwidth values
  #         kernel.fun - FUNCTION- a kernel smoother
  #-------------------------------------
  ##Output: h_opt, the optimal bandwidth
  #-------------------------------------
  ind <- seq(1,length(x_obs), by=1)
  train <- sample(ind,length(ind)/2,replace = FALSE)
  x_train <- x_obs[train]
  y_train <- y_obs[train]
  test <- setdiff(ind, train)
  x_test <- x_obs[test]
  y_test <- y_obs[test]
  
  
  y_fit <- matrix(NA,ncol = length(y_test), nrow = length(h))
  mse <- rep(NA,length(h))
  for (i in 1: length(h))
  {
    y_fit[i,] <- KernPredict(x_train,x_test,y_train,kernel.fun,h_test[i])
    mse[i] <- 1/length(h) * sum((y_test-y_fit[i,])^2)
  }
  
  h_opt <- h_test[which(mse == min(mse))]
  plot(x_obs,y_obs, main = paste('h_opt =', h_opt))
  lines(x_test,y_fit[which(mse==min(mse)),], col = 'red',lw = 2)
  return(h_opt)
}
# 
# MakeHat<- function(x_obs,x_new,d,kernel.fun, h){
#   #-------------------------------------
#   ##Usage: Creates the hat, or smoothing, matrix for local polynomial regression 
#   #-------------------------------------
#   ##Inputs: x_obs - vector - observed x values
#   #         x_new - float - the x value you wish to find the fitted y of
#   #         d - natural number - maximum polynomial degree
#   #         kernel.fun - FUNCTION- a kernel smoother  
#   #         h - float - bandwidth value for the kernel smoother
#   #-------------------------------------
#   ##Output: the hat matrix
#   #-------------------------------------
#   R = matrix(NA, ncol = d+1, nrow = length(x_obs))
#   R[,1] = 1
#   for(i in 1:d)
#   {
#     tempcol = (x_obs - x_new)^i
#     R[,i+1] = tempcol
#   }
#   
#   weights = KernWeights(x_obs,x_new,kernel.fun,h)
#   w = diag(weights)
#   
#   Hat = solve(t(R)%*%w%*%R)%*%t(R)%*%w
#   out = Hat[1,]
#   return(out)
# }
# 
# getHat <- function(x_obs, x_new, d, kernel.fun, h){
#  hat <-  apply(x_new, 1, function(x) MakeHat(x_obs = x_obs, x, d=d, kernel.fun= kernel.fun, h=h))
#  return(hat)
# }

makeHat <- function(x_obs, x_new, kernel.fun, h){ 
  #-------------------------------------
  #   ##Usage: Creates the hat, or smoothing, matrix for local polynomial regression with d=1
  #   #-------------------------------------
  #   ##Inputs: x_obs - vector - observed x values
  #   #         x_new - matrix vector - the vector of x values you wish to find the fitted y of
  #   #         kernel.fun - FUNCTION- a kernel smoother  
  #   #         h - float - bandwidth value for the kernel smoother
  #   #-------------------------------------
  #   ##Output: the hat matrix
  #   #-------------------------------------
  kern_weights <- apply(x_new, 1, function(x) KernWeights(x_obs = x_obs, x_target = x, kernel.fun, bandwidth = h))
  kern_weights <- t(kern_weights) * h
  diff <- apply(x_new, 1, function(x) x_obs - x)
  s_1 = diag(kern_weights %*% diff)
  s_2 = diag(kern_weights %*% diff^2)
  weights <- matrix(NA, nrow = length(x_new), ncol = length(x_obs))
  for(i in 1:length(x_new)){weights[i,] = kern_weights[i,]*(s_2[i]-diff[,i]*s_1[i])}
  weight_norm = apply(weights,1,sum)
  hat = weights/weight_norm
  return(hat)
}

LOOVC <- function(x_obs, y_obs, d, kernel.fun, h){
  #-------------------------------------
  ##Usage: Finds the optimal bandwidth for a given linear regressor using leave one out cross validation
  #-------------------------------------
  ##Inputs: x_obs - matrix vector - observed x values
  #         x_new - vector
  #         y_obs - vector - observed y values
  #         h_test - vector- a vector of candidate bandwidth values
  #         kernel.fun - FUNCTION- a kernel smoother
  #-------------------------------------
  ##Output: h_opt, the optimal bandwidth
  #-------------------------------------
  loovc = rep(NA, length(h))
  for(i in 1:length(h)){
  hat = makeHat(x_obs, x_obs, kernel.fun, h[i])
  y_hat = hat %*% y_obs
  loovc[i] = sum(((y_obs-y_hat)/(1-diag(hat)))^2)}
  h_opt = h[which(loovc == min(loovc))]
  return(h_opt)
}

expsqr<- function(x1,x2){
  exp(-1/2*abs(x1-x2)^2)
}
maternSE<- function(x1,x2,b,t1,t2){
  #-------------------------------------
  ##Usage: calculates the Matern squared exponential covariance for two inputs x1 and x2
  #-------------------------------------
  ##Inputs: x1 - scalar - input 1 for covariance kernel
  #         x2 - scalar - input 2 for covariance kernel
  #         b - scalar - bandwidth parameter. Controls smoothness/frequency of oscillations
  #         t1 - scalar - amplitude parameter, controls size of oscillations
  #         t2 - scalar - noise parameter
  #-------------------------------------
  ##Output: the covariance between x1 and x2
  #-------------------------------------
  out = t1^2*exp(-1/2*(abs(x1-x2)/b)^2) + t2^2*(x1==x2)
  return(out)
}

matern52<- function(x1,x2,b,t1,t2){
  #-------------------------------------
  ##Usage: calculates the Matern 52 covariance for two inputs x1 and x2
  #-------------------------------------
  ##Inputs: x1 - scalar - input 1 for covariance kernel
  #         x2 - scalar - input 2 for covariance kernel
  #         b - scalar - bandwidth parameter. Controls smoothness/frequency of oscillations
  #         t1 - scalar - amplitude parameter, controls size of oscillations
  #         t2 - scalar - noise parameter
  #-------------------------------------
  ##Output: the covariance between x1 and x2
  #-------------------------------------
  out = t1^2*(1+sqrt(5)*abs(x1-x2)/b + 5*(abs(x1-x2))^2/(3*b^2))*exp(-sqrt(5)*abs(x1-x2)/b) + t2*(x1==x2)
  return(out)
}

getCov <- function(x, cov.fun, b,t1,t2){
  #-------------------------------------
  ##Usage: calculates the covariance matrix using one of the Matern covariance kernels
  #-------------------------------------
  ##Inputs: x - vector - all x inputs
  #         cov.fun - function - which matern kernel to use
  #         b - scalar - bandwidth parameter. Controls smoothness/frequency of oscillations
  #         t1 - scalar - amplitude parameter, controls size of oscillations
  #         t2 - scalar - noise parameter
  #-------------------------------------
  ##Output: the covariance matrix of x
  #-------------------------------------
  cov = matrix(NA, length(x), length(x))
  for(i in 1:length(x)){
    for(j in i:length(x)){
      cov[i,j] = cov.fun(x[i],x[j],b,t1,t2)
      cov[j,i] = cov[i,j]
    }
  }
  return(cov)
}

postMean <- function(x_new, x, y, sigma, cov.fun, b, t1, t2){
  #-------------------------------------
  ##Usage: calculates the posterior mean of a gaussian process for nonlinear smoothing
  #-------------------------------------
  ##Inputs: xnew - scalar - new x value to compute posterior expected values at
  #         x - vector- x grid for observed data
  #         y - vector - observed values
  #         cov.fun - function - which matern kernel to use
  #         b - scalar - bandwidth parameter. Controls smoothness/frequency of oscillations
  #         t1 - scalar - amplitude parameter, controls size of oscillations
  #         t2 - scalar - noise parameter
  #-------------------------------------
  ##Output: posterior mean at new predicted point
  #-------------------------------------
  covxnew = matrix(NA, 1, length(x))
  for (i in 1:length(x))
    {
    covxnew[i] = cov.fun(x_new, x[i], b, t1, t2)
  }
  covxx = getCov(x, cov.fun, b, t1, t2)
  expectation = covxnew %*% solve(covxx + sigma*diag(length(covxnew))) %*% y
  return(expectation)
}

postVar <- function(x_new, x, sigma, cov.fun, b, t1, t2){
  #-------------------------------------
  ##Usage: calculates the posterior covariance of a gaussian process for nonlinear smoothing
  #-------------------------------------
  ##Inputs: xnew - scalar - new x value to compute posterior expected values at
  #         x - vector- x grid for observed data
  #         cov.fun - function - which matern kernel to use
  #         b - scalar - bandwidth parameter. Controls smoothness/frequency of oscillations
  #         t1 - scalar - amplitude parameter, controls size of oscillations
  #         t2 - scalar - noise parameter
  #-------------------------------------
  ##Output: posterior variance at new predicted point
  #-------------------------------------
  covx_new = getCov(x_new, cov.fun, b, t1, t2)
  covxnew = matrix(NA, 1, length(x))
  for (i in 1:length(x))
  {
    covxnew[i] = cov.fun(x_new, x[i], b, t1, t2)
  }
  covxx = getCov(x, cov.fun, b, t1, t2)
  
  variance = covx_new - covxnew %*% solve(covxx + sigma*diag(length(covxnew))) %*% t(covxnew)
  return(variance)
}
