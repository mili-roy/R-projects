#The idea behind this monte carlo simulation is to perform the MLE for simple linear regression and 
# and to calculate the bias for intercept and slope from their true values. 
rm(list = ls())
set.seed(2024) # for reproducibility

library(corpcor)
library("numDeriv")

true.alpha= 1; true.beta = 2; true.sigma = 1
n<- 1000 # sample size
nsim<- 300 # number of repeats.

######## to store the results #########
alpha <- numeric(nsim)
beta <- numeric(nsim)
sigma <- numeric(nsim)
estimate<- matrix(NA, nsim, 3)
se<-  matrix(NA, nsim, 3)


#### Data generation ########
  
for (i in 1:nsim){
  
  x<- runif(n, -1, 1)
  # print(paste0("x: ", x))
  ep<- rnorm(n, 0, true.sigma)
  # print(paste0("ep: ", ep))
  y<- true.alpha + true.beta*x + ep
  

###################### MLE function ########

myfun<- function(pars, y, x){
  
  alpha<- pars[1]
  beta<- pars[2]
  sigma<- pars[3]
  n<- length(y)
  
    ep<- y - (alpha + beta*x)
    logl<- -0.5*n*log(2*pi)-0.5*n*log(sigma**2)-sum(ep^2)/(2*sigma**2)
    
    return(-logl)

}


fit<-optim(c(0.90,1.95,0.92), myfun, method= "BFGS", y=y, x=x, hessian = TRUE)

             
estimate[i,]<- fit$par
se[i,]<- sqrt(diag(solve(fit$hessian)))
alpha[i]<- fit$par[1]
beta[i]<- fit$par[2]
sigma[i]<- fit$par[3]

}




alpha_bias<- true.alpha - mean(alpha)
beta_bias<- true.beta - mean(beta)
s_bias<- true.sigma - mean(sigma)


estimate<- data.frame(alpha_bias, beta_bias, s_bias)
print(estimate)


