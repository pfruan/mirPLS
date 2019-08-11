# derivative of SCAD penalty function
dp<-function(theta,lambda,a=3.7){
  p<-length(theta)
  b1<-rep(0,p)
  b1[theta>lambda]<-1
  b2<-rep(0,p)
  b2[theta<(lambda*a)]<-1
  lambda*(1-b1)+((lambda*a)-theta)*b2/(a-1)*b1
}
# deviance residuals of a fit
devresd<-function(y,fit,family=c("binomial","poisson","gaussian")){
  family <- match.arg(family)
  z<-exp(fit)
  r <- switch(family, poisson = {
    u<-z
    r<-y
    id1<-(1:length(y))[y>0]
    id0<-(1:length(y))[y==0]
    if (length(id1)==length(y)){
      r<-sqrt(2*(y*log(y)-y*fit-y+u))*sign(y-u)
    }
    else{
      if(length(id0)==length(y)){ r<-sqrt(2*(u))*sign(-u) }
      else{
        r[id1]<-sqrt(2*(y[id1]*log(y[id1])-y[id1]*fit[id1]-y[id1]+u[id1]))*sign(y[id1]-u[id1])
        r[id0]<-sqrt(2*(u[id0]))*sign(-u[id0])
      }
    }
  }, binomial = {
    u<-z/(1+z)
    r<-sqrt(-y*fit+log(1+z))*sign(y-u)
  }, gaussian = {
    u<-fit
    r<-y-u})
  r
}
# negative log likelihood
loglik<-function(eta,y,family=c("binomial","poisson","gaussian")){
  family <- match.arg(family)
  n<-length(y)
  r <- switch(family,binomial = {
    sum(-y*eta+log(1+exp(eta)))
  }, poisson = {
    sum(-y*eta+exp(eta))
  }, gaussian = {
    sum((-y+eta)^2/2)
  })
  r
}
