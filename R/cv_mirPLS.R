cv_mirPLS<-function(x,y,family=c("binomial","poisson","gaussian"),s1,s2,nfolds=5){
  family <- match.arg(family)
  all.folds <- cv.folds(length(y), nfolds)
  residmat <- matrix(Inf, length(s1)*length(s2), nfolds)
  for (i in seq(nfolds)) {
    omit <- all.folds[[i]]
    for (j in 1:(length(s1)*length(s2))){
      j1 <- (j-1)%/%length(s1)+1
      j2 <- j-(j1-1)*length(s1)
      try({
        cat("lambda1=",s1[j1],"lambda2=",s2[j2],"\n")
        fwk <- mirPLS(x=x[-omit, ],y=y[-omit],family=family,lambda1=s1[j1],lambda2=s2[j2], selection=T)
        fit <- predict_mirPLS(fwk,x[omit, , drop = FALSE])
        ##sum of squares of deviance residuals
        residmat[j, i] <- mean((devresd(y[omit],fit,family=family))^2)
      })
    }
  }
  cv <- apply(residmat, 1, mean)
  jmin <- which.min(cv)
  j1min <- (jmin-1)%/%length(s1)+1
  j2min <- jmin-(j1min-1)*length(s1)
  s1min=s1[j1min]
  s2min=s2[j2min]
  list(s1min = s1min, s2min=s2min)
}
