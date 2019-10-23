\name{sp_clustering}
\alias{sp_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectral clustering
}
\description{
Perform the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
}
\usage{
sp_clustering(data, k = 20, l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Input data, a n*p matrix, n is the number of patients and p is the number of miRNAs.
}
  \item{k}{
A parameter in constructing similarity matrix. Default is 20. It can also be set as [n/10], n is the number of patients.
}
  \item{l}{
The number of clusters. It can be a signle number or a list which contains the possible choices of number of clusters. When l is a list, the number of clusters will be determined by eigengap
}
}
\value{
A vector consisting of cluster labels of each sample.
}
\references{
mirPLS: a partial linear structure identifier method for cancer subtyping using MicroRNAs. Peifeng Ruan, Shuang Wang, Hua Liang.
Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Bo Wang, et al.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

data(x)
data=x[1:200,1:6]
label=sp_clustering(data,l=4)
data(true_lable)
table(label,true_lable)

}
