\name{cv_mirPLS}
\alias{cv_mirPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validation in mirPLS 
}
\description{
Determine regularization parameters in mirPLS
}
\usage{
cv_mirPLS(x,y,family=c("binomial","poisson","gaussian"),s1,s2,nfolds=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of miRNA data, should be normalized into 0-1 range. Columns of x represent patients and rows of x represent miRNAs.
}
  \item{y}{
A vector of responses.
}
  \item{family}{
One of "binomial","poisson","gaussian".
}
  \item{s1}{
The range of the regularization parameter lambda1 which determines whether miRNAs have zero associations or not.
}
  \item{s2}{
The range of the regularization parameter lambda2 which determines whether miRNAs have non-linear or linear associations.
}
  \item{nfolds}{
Number of folders.
}
}
\value{
Regularization parameters suggested by cross validation.
}
\references{
mirPLS: a partial linear structure identifier algorithm for miRNA data. Peifeng Ruan, Shuang Wang, Hua Liang.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

# load data
data(x)
data(y)


# determine regularization parameters using cross validation
lamseq1 <- seq(10)/10
lamseq2 <- seq(10)
cvfit <- cv_mirPLS(x,y,"binomial",lamseq1,lamseq2, nfolds=5)

# fit model 
fit <- mirPLS(x,y,"binomial",cvfit$s1min,cvfit$s2min)
fit$type


}
