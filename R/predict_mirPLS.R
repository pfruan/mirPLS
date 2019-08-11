predict_mirPLS <- function(object,newx,index=NULL)
{
  beta <- object$tbeta
  xpow <- object$xpow
  knots <- object$knots
  ## assuming uniform cubic B-splines are used
  nbasis <- length(knots)-5
  nx <- length(xpow)
  if (is.null(index)) {
    if (nx!=dim(newx)[2]) stop("Wrong dimension of new x!")
    zmat <- NULL
    for (j in 1:nx) {
      if (xpow[j]==0) next
      if (xpow[j]==1) zmat <- cbind(zmat,newx[,j])
      else {
        wk <- splineDesign(knots,newx[,j])
        zmat <- cbind(zmat,wk[,1:nbasis,drop=FALSE])
        ##zmat <- cbind(zmat,t(t(wk)-apply(wk,2,mean))[,1:nbasis])
      }
    }
    val <- drop(zmat%*%beta[-1])+beta[1]
  }
  else {
    if (length(index)>1 | index<1 | index>nx) stop("Wrong index of covariate!")
    if(xpow[index]==0) val <- rep(0,length(newx))
    else {
      if (xpow[index]==1) val <- beta[object$indmat[1,index]]*newx
      else val <- drop(splineDesign(knots,newx)[,1:nbasis]%*%beta[object$indmat[1,index]:object$indmat[2,index]])
    }
  }
  val
}
