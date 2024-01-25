ddwell_t = function(
    r, # vector of dwell times to compute
    t, # time point at which the state is entered
    state, # which state to compute
    Gamma # array of dim c(N,N,L)
){
  L = dim(Gamma)[3] # cycle length
  ind = (t+(1:max(r))-1)%%L
  ind[which(ind==0)] = L
  gamma_ii = Gamma[state,state,ind]
  pmf = c(1, cumprod(gamma_ii)[-length(ind)])*(1-gamma_ii)
  return(pmf[r])
}

ddwell = function(
    r, # vector of dwell times to compute
    state, # which state to compute
    Gamma # array of dim c(N,N,L)
){
  L = dim(Gamma)[3]
  N = dim(Gamma)[1]
  delta = matrix(NA, L, N) # calculate p-stationary
  GammaT=Gamma[,,1]
  for (k in 2:L){ GammaT=GammaT%*%Gamma[,,k] }
  delta[1,] = solve(t(diag(N)-GammaT+1), rep(1,N))
  for(k in 2:L){ delta[k,] = delta[k-1,]%*%Gamma[,,k-1] }
  weights=numeric(L) # calculate weights
  weights[1] = sum(delta[L,-state] * Gamma[-state,state, L])
  for (k in 2:L){ weights[k] = sum(delta[k-1,-state] * Gamma[-state,state, k-1]) }
  weights = weights / sum(weights)
  pmfs_weighted = matrix(NA, L, max(r))# calculate all weighted d_i^t's
  for(k in 1:L){ pmfs_weighted[k,] = weights[k] * ddwell_t(1:max(r), k, state, Gamma) }
  pmf = apply(pmfs_weighted, 2, sum)
  return(pmf[r])
}