\name{predict_mirPLS}
\alias{predict_mirPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions of new data for mirPLS.
}
\description{
predicting the regression mean of new data for mirPLS.
}
\usage{
predict_mirPLS(object,newx,index=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A matrix of miRNA data, should be normalized into 0-1 range. Columns of x represent patients and rows of x represent miRNAs.
}
  \item{newx}{
a data frame of new x if index is NULL, or a vector of new x if index is an integer.
}
  \item{index}{
index - index of the covariate whose component function estimate is to be returned; if null, this gives predicted value for eta.
}
}
\value{
Regression means of new data.
}
\references{
mirPLS: a partial linear structure identifier method for cancer subtyping using MicroRNAs. Peifeng Ruan, Shuang Wang, Hua Liang.
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
