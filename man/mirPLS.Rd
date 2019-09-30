\name{mirPLS}
\alias{mirPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial linear structure identifier for miRNA data
}
\description{
Select the miRNAs linearly or non-linearly associated to outcomes 
}
\usage{
mirPLS(x,y,family=c("binomial","poisson","gaussian"),lambda1,lambda2,selection=T,nknots=1,beta0=NULL,eps=1e-6,btol=1e-3)
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
  \item{lambda1}{
A regularization parameter determines whether miRNAs have zero associations or not.
}
  \item{lambda2}{
A regularization parameter determines whether miRNAs have non-linear or linear associations.
}
  \item{selection}{
Whether perform the proposed further selection step.
}
}
\value{
Returns a mirPLS object. In the returned object, type represents a vector of association types of miRNAs. 2 represents a non-linear miRNA, 1 represents a linear association and 0 represents no association.
}
\references{
Nonlinear molecular omics features for cancers subtyping. Peifeng Ruan, Shuang Wang, Hua Liang.
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
