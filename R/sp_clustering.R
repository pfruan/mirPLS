sp_clustering=function(data,k=20,alpha=0.5,l,Normalization=T){
  if (Normalization==T){
    Data1 = standardNormalization(data)
  }else{
    Data1 = data}
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
  W1 = affinityMatrix(Dist1, k, alpha)
  c=estimateNumberOfClustersGivenGraph(W1, NUMC=l)[[1]]
  label = spectralClustering(W1, c)
  return(label)
}
