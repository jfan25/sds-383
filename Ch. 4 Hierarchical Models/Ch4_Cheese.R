#-----------------Cheese Problem------------------------#
library(mvtnorm)
library(LaplacesDemon)
library(purrr)

#####Functions for this problem
make_sig_S<- function(beta, mu,psi_0){
  #-------------------------------------
  ##Usage: Creates the conditional posterior matrix for inverse wishart draw of Sigma
  #-------------------------------------
  ##Inputs: beta - matrix - Current beta values for all stores.
  #         mu - vector - the current mean for the betas. For the top level of the hierachy
  #         psi_0 - matrix - hyper parameter matrix of inverse wishart
  #-------------------------------------
  ##Output: updated parameter matrix for inverse wishart
  #-------------------------------------
  stuff = apply(beta, 1, function(x) x-mu)
  morestuff = array(0,dim = c(4,4,88,1))
  for(i in 1:88){
    morestuff[,,i,] = stuff[,i]%*%t(stuff[,i])
  }
  total = apply(morestuff, c(1,2), sum)
  sig_s = psi_0 + total
  return(sig_s)
}

make_lam_rate <- function(temp_beta, x,y,b){
  #-------------------------------------
  ##Usage: Creates rate parameter for updating the between group variance from a gamma
  #-------------------------------------
  ##Inputs: temp_beta - matrix - Current beta values for all stores.
  #         x - list of matrices - x values for each store: log price, display, interaction
  #         y - list of vectors - log volume for each store
  #         b - scalar - hyperparameter rate from prior
  #-------------------------------------
  ##Output: updated rate parameter
  #-------------------------------------
  beta_list = apply(temp_beta, 1, function(t) as.list(t))
  xb = map2(x,beta_list, function(x,y) unlist(x)%*% unlist(y))
  yxb_sqr = map2(y,xb,function(y,x) (unlist(y)-unlist(x))^2)
  lamp_rate = b + 1/2*sum(sapply(yxb_sqr,sum))
  return(lamp_rate)
}
############ Read in the data
#cheese <- read.csv('/home/jf37828/Documents/StatModeling2/cheese.csv',header = TRUE, sep = ',')
cheese <- read.csv(file="C:/Users/15202/OneDrive/Documents/stat modeling 2/Ch 4/cheese.csv", header=TRUE, sep=",")

###Creating data structures
y <- list()
x <- list()
for(i in 1:length(unique(cheese$store))){ #creating list of matrices for X and Y
  y[i] = list(log(c(cheese[cheese$store == unique(cheese$store)[i],]$vol)))
  price = log(c(cheese[cheese$store == unique(cheese$store)[i],]$price))
  ind = c(cheese[cheese$store == unique(cheese$store)[i],]$disp)
  interact = price * ind
  temp_matrix = cbind(rep(1,length(ind)), price, ind, interact)
  x[i] = list(temp_matrix)
}

#hyperparameters
a = 1/2 #scale for gamma prior for sigma^2, the error variance
b = 1/2 #rate for gamma prior for sigma^2, the error variance
v_0 = 2 # dgf for inverse wishart prior
psi_0 = diag(rep(1,4)) #S matrix for inverse wishart prior
#Setting up storage data containers

##########chain with inverse wishart
N=2000 #iteration number for Gibbs


mu = matrix(rep(NA,4), ncol = 4, nrow= N)
sigma = array(0, dim = c(4,4,N,1))
beta = array(0, dim = c(88,4,N,1))
lambda = c(rep(0,N))

temp_mu = c(10,-2,0,0)
temp_sig = diag(rep(1,4))
temp_beta = rmvn(n = 88, mu = temp_mu, Sigma = temp_sig)
temp_lam = 1

mn = nrow(cheese)
n=88
xtx = lapply(x,function(x) t(unlist(x))%*% unlist(x))
xty = map2(x,y, function(x,y) t(unlist(x))%*% unlist(y))

for(i in 1:N){
  #update mu
  mu_mean = colMeans(temp_beta)
  mu_sig = solve(n*solve(temp_sig))
  temp_mu = rmvnorm(1,mean = mu_mean, sigma = mu_sig)
  #update Sigma
  sig_nu = n + v_0
  sig_S = make_sig_S(temp_beta, temp_mu, psi_0)
  temp_sig = rinvwishart(nu = sig_nu, S = sig_S)
  #update lambda
  lam_shape = mn/2+a
  lam_rate = make_lam_rate(temp_beta,x,y,b)
  temp_lam = rgamma(1, lam_shape, rate = lam_rate)
  #update betas
  beta_var = lapply(xtx,function(x) solve(temp_lam*x + temp_sig))
  beta_mean = map2(beta_var,xty, function(x,y) x%*%(temp_lam* y + temp_sig%*%t(temp_mu)))
  temp_beta = matrix(unlist(map2(beta_mean,beta_var, function(m,v) rmvnorm(1,mean = m, sigma = v))), ncol = 4, byrow = TRUE)
  #store values
  mu[i,] = temp_mu
  sigma[,,i,] = temp_sig
  beta[,,i,] = temp_beta
  lambda[i] = temp_lam
  print(i)
}

par(mfrow = c(4,1))
plot(mu[500:N,1], type='l', col = 'sky blue', main = "Trace Plot of Mu_1")
plot(mu[500:N,2], type='l', col = 'sky blue', main = "Trace Plot of Mu_2")
plot(mu[500:N,3], type='l', col = 'sky blue', main = "Trace Plot of Mu_2")
plot(mu[500:N,4], type='l', col = 'sky blue', main = "Trace Plot of Mu_2")


plot(1/lambda[500:N], type='l', col = 'sky blue', main = "Trace Plot of Sigma^2",ylab = "sigma^2")
###Calculating the MSE
beta = beta[,,500:N,] #burning the first 500 iterations 
avgbeta = apply(beta,c(1,2), mean) #averaging the betas over all iterations
avgbetalist = apply(avgbeta, 1, as.list) #converting to a list to pass into map2 function
y_hat = map2(x, avgbetalist, function(x,y) unlist(x) %*% unlist(y)) #calculating the predicted values
tempmse = map2(y, y_hat, function(x,y) sum((x-y)^2)) #calculating sum((y-y_hat)^2) for each store
mse = lapply(tempmse, mean) #calculating the MSE for each store

####Calculating the posterior fitted lines
x_grid = seq(0, 1.5, length.out = 100)
x_grid = cbind(rep(1,100), x_grid)

beta_no_disp = cbind(avgbeta[,1], avgbeta[,2])
beta_disp = cbind(avgbeta[,1] + avgbeta[,3], avgbeta[,2] + avgbeta[,4])
y_lines_nd = apply(beta_no_disp, 1, function(x) x_grid %*% x) #calculating the predicted values of no display
y_lines_d = apply(beta_disp, 1, function(x) x_grid %*% x) #calculating the predicted values of display

###Plotting the posterior fitted lines with the log data
par(mar=c(2,2,1,1), mfrow = c(5,5), family = 'Palatino')

for(i in 1:25){
  ind_np <- which(x[[i]][,3]==0)
  ind_disp <- which(x[[i]][,3]==1)
  if(length(ind_np) < length(ind_disp)){
    plot(x[[i]][ind_disp, 2], y[[i]][ind_disp], col = 'blue',  main = i, xlab = "log price", ylab = 'log volume', ylim = c(6,11))
    points(x[[i]][ind_np,2],y[[i]][ind_np], col = 'red')
  }else{
    plot(x[[i]][ind_np,2],y[[i]][ind_np], col = 'red', main = i, xlab = "log price", ylab = 'log volume', ylim = c(6,11))
    points(x[[i]][ind_disp, 2], y[[i]][ind_disp], col = 'blue')
  }
  lines(x_grid[,2], y_lines_nd[,i], col = 'red')
  lines(x_grid[,2], y_lines_d[,i], col = 'blue')
}
#####price recommendation
opt_price = matrix(NA, ncol = 1501, nrow = 88)
for(i in 1:88)
{
  opt_price[i,] = beta[i,2,]/(beta[i,2,]+1)
}

