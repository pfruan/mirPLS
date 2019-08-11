mirPLS<-function(x,y,family=c("binomial","poisson","gaussian"),lambda1,lambda2,selection=T,nknots=1,beta0=NULL,eps=1e-6,btol=1e-3){
  call <- match.call()
  family <- match.arg(family)
  nobs <- length(y)
  nx <- dim(x)[2]
  ## constants for cubic B-splines: change here for B-splines of other orders
  ## assemble Z matrix with dim=nobs*(nbasis*nx)
  knots <- c(rep(0,3),seq(0,1,length=nknots+2),rep(1,3))
  ## Note: need to get rid of one basis for model identifiability
  nbasis <- nknots+3
  zmat <- NULL
  for (j in 1:nx) {
    wk <- splineDesign(knots,x[,j])
    zmat <- cbind(zmat,wk[,1:nbasis])
  }
  #browser()
  ## assemble integral matrices A and D
  nmesh <- 100
  quad <- gauss.quad(nmesh,c(0,1))
  wk0 <- splineDesign(knots,quad$pt)
  wk2 <- splineDesign(knots,quad$pt,derivs=rep(2,nmesh))
  amat <- dmat <- NULL
  for (i in 1:nbasis) {
    amat <- cbind(amat,apply(quad$wt*wk0[,i]*wk0[,1:nbasis],2,sum))
    dmat <- cbind(dmat,apply(quad$wt*wk2[,i]*wk2[,1:nbasis],2,sum))
  }

  ## use the minimum lambda in the lambda sequence
  if (is.null(beta0)) {
    fwk <- glmnet(zmat,y,family=family)
    lamlasso <- sort(cv.glmnet(zmat,y,family=family)$lambda)
    lambda <- lamlasso[1]
    beta0 <- drop(Matrix(coef(fwk,s=lambda),sparse=FALSE))
  }

  betatnew <- beta0; betac <- beta0[-1]
  beta <- beta0[-1]
  ## xpow: vector of x-powers, taking values 0,1,2.
  ## xnzind: vector of indices for nonzero x_j's.
  ## indmat: matrix of starting and ending indices in zmat and betac for
  ##         predictors, c(indmat[1,j],indmat[2,j]) are indices for x_j.
  xpow <- rep(2,nx)
  xnzind <- 1:nx
  indmat <- rbind((0:(nx-1))*nbasis+1,(1:nx)*nbasis)
  ## start RWLS loop
  bool.rwls <- nit.rwls <- 1
  npred <- nx
  while(bool.rwls) {
    cat("nit.rwls=",nit.rwls,"\n")
    ## compute w_i, u_i, tildeY_i, tildeY_i*, muy
    betat <- betatnew
    etatx <- betat[1] + drop(zmat%*%betat[-1])

    switch(family, binomial = {
      uwk <- -y + 1/(1+exp(-etatx))+.Machine$double.eps
      wt <- 1/(1+exp(-etatx))/(1+exp(etatx))+.Machine$double.eps
    },
    poisson = {
      uwk <- -y + exp(etatx)
      wt <- exp(etatx)
    },
    gaussian = {
      uwk <- -y + etatx
      wt <- 1
    } )
    ywk <- etatx - uwk/wt
    muywk <- sum(ywk*wt)/sum(wt)
    ystar <- sqrt(wt)*(ywk - muywk)
    ## start loop for the LQA
    bool.pen <- nit.pen <- 1
    while (bool.pen) {
      cat("nit.pen=",nit.pen,"\n")
      norma <- normd <- rep(0,npred)
      for (j in 1:npred) {
        jwk <- xnzind[j]
        ## compute norma and normd differently according to linear or nonlinear
        ## for linear term, norma=abs(b)*sqrt(int(x*x)dx) = abs(b)/sqrt(3)
        if(xpow[jwk]==1) norma[j] <- abs(betac[indmat[1,jwk]])/sqrt(3)
        else {
          wk <- betac[indmat[1,jwk]:indmat[2,jwk]]
          norma[j] <- sqrt(max(t(wk)%*%amat%*%wk,0))
          normd[j] <- sqrt(max(t(wk)%*%dmat%*%wk,0))
        }

        ## further actions for new zero and linear terms
        if (normd[j] <=eps) {
          if (norma[j]<=eps) { # zero component
            norma[j] <- 0
            normd[j] <- 0
            ## if all remaining terms are zeros, exit the for loop of 1:npred
            if (length(indmat[1,jwk]:indmat[2,jwk])==dim(zmat)[2]) {
              zmat <- betac <- NULL
              indmat <- matrix(0,2,nx)
              break  # out of the for(j in 1:npred) loop
            }
            zmat <- zmat[,-(indmat[1,jwk]:indmat[2,jwk]),drop=FALSE]
            betac <- betac[-(indmat[1,jwk]:indmat[2,jwk])]
            #cat("jwk=",jwk,"indmat[,jwk]=",c(indmat[1,jwk],indmat[2,jwk]),"\n")
            if(jwk<nx) {
              wk <- ((jwk+1):nx)[indmat[1,(jwk+1):nx]>0]
              indmat[,wk] <- indmat[,wk]-(indmat[2,jwk]-indmat[1,jwk]+1)
            }
            indmat[,jwk] <- rep(0,2)
            xpow[jwk] <- 0
          }
          ## if x_jwk is linear and was nonlinear,
          ## then truncate zmat and update indmat, xpow.
          else if (xpow[jwk]>1) {
            zmat[,indmat[1,jwk]] <- x[,jwk]
            betac <- betac[-((indmat[1,jwk]+1):indmat[2,jwk])]
            zmat <- zmat[,-((indmat[1,jwk]+1):indmat[2,jwk]),drop=FALSE]
            if(jwk<nx) {
              wk <- ((jwk+1):nx)[indmat[1,(jwk+1):nx]>0]
              indmat[,wk] <- indmat[,wk]-(indmat[2,jwk]-indmat[1,jwk])
            }
            indmat[,jwk] <- rep(indmat[1,jwk],2)
            xpow[jwk] <- 1
          }
        }

      }

      wk <- (1:npred)[norma==0 & normd==0]
      if (length(wk)>0) {
        if (length(wk)==npred) {
          norma <- normd <- xnzind <- NULL
          npred <- 0
        }
        else {
          norma <- norma[-wk]
          normd <- normd[-wk]
          xnzind <- xnzind[-wk]
          npred <- npred-length(wk)
        }
      }

      if (npred==0) {
        muzwk <- 0
        beta <- NULL
        bool.pen <- 0
        break
      }
      ## compute Z*, muz
      muzwk <- apply(wt*zmat,2,sum)/sum(wt)
      zmatstar <- sqrt(wt)*t(t(zmat) - muzwk)
      ## compute penalty matrix
      dwk <- sum(xpow==1)+nbasis*sum(xpow==2)
      pmat <- matrix(0,dwk,dwk)
      for (j in 1:npred) {
        jwk <- xnzind[j]
        iwk <- indmat[1,jwk]
        #cat("npred=",npred,"j=",j,"jwk=",jwk,"\n")
        if (xpow[jwk]==1) pmat[iwk,iwk] <- dp(norma[j],lambda1)/norma[j]/3
        else if (xpow[jwk]==2) {
          c1 <- dp(norma[j],lambda1)/norma[j]
          c2 <- dp(normd[j],lambda2)/normd[j]
          pmat[(iwk:(iwk+nbasis-1)),(iwk:(iwk+nbasis-1))] <- c1*amat+c2*dmat
        }
      }

      hmat <- t(zmatstar)%*%zmatstar+pmat
      hdim <- dim(hmat)[1]
      try({
        zz <- chol(hmat,pivot=TRUE)
        hcho <- zz
        rkv <- hdim
        while(hcho[rkv,rkv]<hcho[1,1]*eps) rkv <- rkv-1
        if (rkv<hdim) hcho[(rkv+1):hdim,(rkv+1):hdim] <- hcho[1,1]*diag(hdim-rkv)
        betanew <- t(zmatstar)%*%ystar
        betanew <- betanew[attr(zz,"pivot")]
        betanew <- backsolve(hcho,betanew,tran=TRUE)
        if (rkv<hdim) betanew[(rkv+1):hdim] <- 0
        betanew <- backsolve(hcho,betanew)
        betanew[attr(zz,"pivot")] <- betanew
      })
      ## check convergence of LQA

      if (length(beta)==length(betanew))
        if (max(abs(beta-betanew))<btol) bool.pen <- 0
      betac <- beta <- betanew
      nit.pen <- nit.pen + 1
      if (nit.pen > 20) bool.pen <- 0
    } # end of the LQA loop
    ## new betat for approximating lkhd
    if (is.null(beta)) {
      betanew <- muywk
      break
    }
    betatnew <- c(muywk-sum(muzwk*beta),beta)
    ## check convergence of IRLS
    if (length(betat)==length(betatnew))
      if (max(abs(betat-betatnew))<btol) bool.rwls <- 0
    nit.rwls <- nit.rwls+1
    if (nit.rwls>20) bool.rwls <- 0
  } # end of RWLS loop
  #cat("xpow=",xpow,"\n")

  if (selection== T){
    idx=seq(dim(x)[2])[xpow==2]
    res_gam=rep(999,length(idx))
    for (i in 1:length(idx)){
      fit3=gam(y~s(x[,idx[i]]),family = "binomial")
      res_gam[i]=summary(fit3)$anova$`P(Chi)`[2]
    }
    res_gam=p.adjust(res_gam,"fdr")
    type=xpow
    type[idx[res_gam>0.05]]=1
  }
  list(type=type,beta=beta,tbeta=betat,zmat=zmat,xpow=xpow,indmat=indmat,xnzind=xnzind,knots=knots)
}
