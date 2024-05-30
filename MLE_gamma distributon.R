library(insuranceData)
data(dataCar)
head(dataCar)

# take only non-zeros responses

claims<-subset(dataCar[,5],(dataCar[,5]!=0))
head(claims)
plot(density(claims),main="Density estimate of data")

# log likelihood equations
gammaeq1 = function(u,v,y) {
  n = length(y)
  -n*v/u+sum(y)*v/u^2
}

gammaeq2 = function(u,v,y) {
  n = length(y)
  -n*digamma(v)+n*(log(v)+1)-n*log(u)+sum(log(y))-sum(y)/u
}

# construct objective function
gammamle = function(theta,y) {
  gammaeq1(theta[1],theta[2],y)^2+gammaeq2(theta[1],theta[2],y)^2 
}

claimmod = nlm(f= function(theta) {
  gammamle(theta,claims)
},p=c(mean(claims),1))
# find the mle estimation

claimmod$estimate
# make histogram
gammahist = hist(claims,col="red",xlim=c(0,30e3),breaks=20,freq=FALSE)
curve(dgamma(x,scale=claimmod$estimate[1],
             shape=claimmod$estimate[2]),add=TRUE,col="green",lwd=2,lty=2)