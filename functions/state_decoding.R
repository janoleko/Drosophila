cond_stateProbs = function(Gamma, delta, eta, size, X, N=2, L=48, K=3){
  IDs = unique(X$ID)
  nbAnimals = length(IDs)
  PHIs = matrix(nrow = 0, ncol = N)
  for(i in 1:length(IDs)){
    X_id = X[which(X$ID == IDs[i]),]
    DD = which(X_id$condition == "DD")[1]
    
    phi = matrix(nrow = nrow(X_id), ncol = N)
    
    allprobs = matrix(1, nrow(X_id), N)
    ind = which(!is.na(X_id$activity))
    for (j in 1:N){ allprobs[ind,j] = dnbinom(X_id$activity[ind], size = size[j], mu = eta[j]) }
    
    foo = delta%*%diag(allprobs[1,])
    phi[1,] = foo/sum(foo)
    for(t in 2:(DD-1)){
      foo = phi[t-1,]%*%Gamma[,,X_id$tod[t],1]%*%diag(allprobs[t,])
      phi[t,] = foo/sum(foo) }
    for(t in DD:nrow(X_id)){
      foo = phi[t-1,]%*%Gamma[,,X_id$tod[t],2]%*%diag(allprobs[t,])
      phi[t,] = foo/sum(foo) }
    PHIs = rbind(PHIs, phi)
  }
  return(PHIs)
}

cond_stateProbs_hom = function(Gamma, delta, eta, size, X, N=2){
  IDs = unique(X$ID)
  nbAnimals = length(IDs)
  PHIs = matrix(nrow = 0, ncol = N)
  for(i in 1:length(IDs)){
    X_id = X[which(X$ID == IDs[i]),]
    DD = which(X_id$condition == "DD")[1]
    phi = matrix(nrow = nrow(X_id), ncol = N)
    allprobs = matrix(1, nrow(X_id), N)
    ind = which(!is.na(X_id$activity))
    for (j in 1:N){ allprobs[ind,j] = dnbinom(X_id$activity[ind], size = size[j], mu = eta[j]) }
    foo = delta%*%diag(allprobs[1,])
    phi[1,] = foo/sum(foo)
    for(t in 2:(DD-1)){
      foo = phi[t-1,]%*%Gamma[,,1]%*%diag(allprobs[t,])
      phi[t,] = foo/sum(foo) }
    for(t in DD:nrow(X_id)){
      foo = phi[t-1,]%*%Gamma[,,2]%*%diag(allprobs[t,])
      phi[t,] = foo/sum(foo) }
    PHIs = rbind(PHIs, phi)
  }
  return(PHIs)
}