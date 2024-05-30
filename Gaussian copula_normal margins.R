rm(list = ls())
set.seed(2078)
libname <- c("MASS","mnormt","psych", "foreign","mvtnorm","stats4","miscTools",
             "matrixcalc","corpcor","Matrix")
lapply(libname, require, character.only=T)
library(openxlsx)

ptm <- proc.time()
cat("The program started:", date(), "\n")
#library(pracma)


alpha0_true<- 1
alpha1_true<- 1.5
beta0_true<- 2
beta1_true<- 2.5; sd_true1<- 1
sd_true2<-1; rho_true<- 0.5
sigma<- matrix(c(1,rho_true,rho_true,1),2,2)
n<- 500 #sample size
r<- 100 # no of repeats

x<- matrix(NA,n,1)
z1<- matrix(NA,n,1)
z2<- matrix(NA,n,1)
w1<- matrix(NA,n,1)
w2<- matrix(NA,n,1)
mu1<- matrix(NA,n,1)
mu2<- matrix(NA,n,1)
y1<- matrix(NA,n,1)
y2<- matrix(NA,n,1)
CDF_y1<- matrix(NA,n,1)

estimate<- matrix(NA,r,7)
SE<- matrix(NA,r,7)
alpha0<- matrix(NA, r, 1)
alpha1<- matrix(NA, r, 1)
beta0<- matrix(NA, r, 1)
beta1<- matrix(NA, r, 1)
sd1<- matrix(NA,r,1)
sd2<- matrix(NA, r, 1)
rho<- matrix(NA, r, 1)

for (j in 1:r){
  print(paste0(" j: ", j))
  z<- rmvnorm(n,mean=c(0,0),sigma) #To generate data from bivariate normal distribution
  u<- pnorm(z)
  #z1<- matrix(z[,1], n, 1) # to consider 1st column from z
  #z2<- matrix(z[,2], n, 1) # to consider 2nd column from z
  w1<- sort(runif(n,1,2)) # to generate covariate for 1st outcome
  w2<- sort(runif(n,1,2)) # to generate covariate for 2nd outcome
  mu1<- alpha0_true+(alpha1_true*w1)
  mu2<- beta0_true+(beta1_true*w2) # mean for second outcome
  y1<- qnorm(u[,1],mean= mu1, sd_true1,lower.tail=TRUE, log.p=FALSE) # to generate second outcome from logistic density
  y2<- qnorm(u[,2],mean= mu2, sd_true2,lower.tail=TRUE, log.p=FALSE) # to generate second outcome from logistic density  
  
  
  data<-data.frame(y1=y1, y2=y2, mu1=mu1, mu2=mu2, w1=w1, w2=w2)
  #head(data)
  
  ################# Analyze the data, it analyzes only first data set, need to include repeated data in the analysis#########
  myloglik <- function(pars,data){
    data<- data
    alpha0<- pars[1]
    alpha1<- pars[2]
    beta0<- pars[3]
    beta1<- pars[4]
    sd1<- pars[5]
    sd2<- pars[6]
    eta<- pars[7]
    rho<- atanh(eta)# Fisher z transformation
    #print(paste0("rho: ", rho))
    w1<- w1
    w2<- w2
    mu1<- alpha0+(alpha1*w1)
    mu2<- beta0+(beta1*w2) # mean for second outcome
    sd1[sd1<=0]<- 0.001
    sd2[sd2<=0]<- 0.001
    y1<- y1
    y2<- y2
    CDF_y1<- pnorm(y1, mean= mu1, sd1,lower.tail=TRUE, log.p=FALSE)# cdf of second outcome
    v1<- qnorm(CDF_y1)
    CDF_y2<- pnorm(y2, mean= mu2, sd2,lower.tail=TRUE, log.p=FALSE)# cdf of second outcome
    v2<- qnorm(CDF_y2)
    log_phi_z1<-  (-0.5)*log(2*pi)-0.5*v1^2
    log_phi_z2<-  (-0.5)*log(2*pi)-0.5*v2^2
    comp1<-  (-1)*log(2*pi)-(0.5*log((1-rho^2)))-((v1^2+v2^2-2*rho*v1*v2)/(2*(1-rho^2))) #bi variate standard normal for copula density
    comp2<-  (-0.5)*log(2*pi)-log(sd1)-0.5*((y1-mu1)/sd1)**2
    comp3<-  (-0.5)*log(2*pi)-log(sd2)-0.5*((y2-mu2)/sd2)**2
    comp4<-  log_phi_z1
    comp5<-  log_phi_z2
    ll<-  comp1+comp2+comp3-comp4-comp5
    loglik<- sum(ll)
    return(-loglik)
  }
  
  #true_value= c(1, 1.5, 2, 2.5, 1,1, 0.5)  
  
  
#fit <- optim(par= as.numeric(c(0.93, 1.45, 1.97, 2.48, 0.94, 0.97, 0.47)), myloglik, method = "BFGS", data=data,
 # hessian=TRUE,control = list(trace = 3, maxit = 500, ndeps = rep(1e-5, 7)))


fit <- try(optim(par = as.numeric(c(0.97, 1.48, 1.95, 2.45, 0.95, 0.98, 0.47)), myloglik, method="L-BFGS-B", data=data,
              lower= c(-Inf,-Inf,-Inf,-Inf, 0.001, 0.001, -1), upper= c(Inf, Inf, Inf, Inf, Inf, Inf, 1), 
              hessian = TRUE, control=list(trace=1,  fnscale=1e4)), silent=TRUE)                                                                                                               
  
  mat<- fit$hessian
  #is.positive.definite(mat)
  mat1<- make.positive.definite(mat)
  OI<- solve(mat1,tol = 1e-27)
  estimate[j,]<- fit$par
  SE[j,]<- sqrt(diag(OI))
  alpha0[j]<- fit$par[1]
  alpha1[j]<- fit$par[2]
  beta0[j]<- fit$par[3]
  beta1[j]<- fit$par[4]
  sd1[j]<- fit$par[5]
  sd2[j]<- fit$par[6]
  rho[j]<- tanh(fit$par[7])
  
}

#estimate
#SE

alpha0 <- as.numeric(alpha0)
alpha1 <- as.numeric(alpha1)
beta0 <- as.numeric(beta0)
beta1 <- as.numeric(beta1)
sd1<- as.numeric(sd1)
sd2 <- as.numeric(sd2)
rho <- as.numeric(rho)

# ###### mean and bias #########

alpha0_bias<- (mean(alpha0)-alpha0_true)
alpha1_bias<-(mean(alpha1)-alpha1_true)
beta0_bias<- (mean(beta0)-beta0_true)
beta1_bias<-(mean(beta1)-beta1_true)
sd1_bias<-(mean(sd1)-sd_true1)
sd2_bias<-(mean(sd2)-sd_true2)
rho_bias<-(mean(rho)-rho_true)
# # ###MSE
alpha0_mse<-(sum(alpha0-alpha0_true)^2)/r
alpha1_mse<-(sum(alpha1-alpha1_true)^2)/r
beta0_mse<-(sum(beta0-beta0_true)^2)/r
beta1_mse<-(sum(beta1-beta1_true)^2)/r
sd1_mse<-(sum(sd1-sd_true1)^2)/r
sd2_mse<-(sum(sd2-sd_true2)^2)/r
rho_mse<-(sum(rho-rho_true)^2)/r
# 
########### output ##############
bias<- data.frame(alpha0_bias,alpha1_bias,beta0_bias,beta1_bias,sd1_bias,sd2_bias,rho_bias)
bias<- t(bias)

MSE<- data.frame(alpha0_mse,alpha1_mse,beta0_mse,beta1_mse,sd1_mse,sd2_mse,rho_mse)
MSE<- t(MSE)

output<- round(cbind(bias, MSE),5)
output

library(openxlsx)
Estimate<- data.frame(alpha0= estimate[,1],alpha1= estimate[,2],beta0= estimate[,3],beta1= estimate[,4],sd1= estimate[,5],
                sd2= estimate[,6], rho=estimate[,7])
SE<- data.frame(alpha0= SE[,1],alpha1= SE[,2],beta0= SE[,3],beta1=SE[,4],sd1= SE[,5],sd2= SE[,6],rho= SE[,7])

dataset_names <- list('Sheet1' = Estimate, 'Sheet2' = SE)
openxlsx::write.xlsx(dataset_names, file = 'optim.xlsx')

