# negative log-likelihood function with random effects in state-dependent process

mllk = function(theta.star, X, L = 48, K = 3, M = 25, low, up){
  N = 2
  source("./functions/periodic_variation.R")
  # state process
  coef = array(NA, dim = c(N*(N-1), 1+2*K, 2))
  coef[,,1] = matrix(theta.star[1:((1+2*K)*N*(N-1))], 
                     nrow = N*(N-1), ncol = (1+2*K)) # LD
  coef[,,2] = matrix(theta.star[((1+2*K)*(N-1)*N) + 1:((1+2*K)*N*(N-1))], 
                     nrow = N*(N-1), ncol = (1+2*K)) # DD
  # state-dependent process
  # random effects: params for gamma distribution
  mu = exp(theta.star[(2*(1+2*K)*(N-1)*N) + 1:N]) 
  sigma = exp(theta.star[(2*(1+2*K)*(N-1)*N) + N + 1:N])
  # fixed dispersion parameter for negative binomial distribution
  size = exp(theta.star[(2*(1+2*K)*(N-1)*N) + 2*N + 1:N])
  # building 48 tpms
  Gamma = array(NA, dim = c(N,N,L,2))
  for(cond in 1:2){
    for(k in 1:L){
      G = diag(N)
      G[!G] = exp(pv(coef[,,cond], time = (k-1), degree = K, L = L))
      G = G/rowSums(G)
      Gamma[,,k,cond] = G }}
  # periodic stationary start
  GammaT=Gamma[,,42,1]
  for(k in 42+1:(L-1)){
    k=ifelse(k%%L==0,L,k%%L)
    GammaT = GammaT%*%Gamma[,,k,1] }
  delta = solve(t(diag(nrow(GammaT))-GammaT+1), rep(1,nrow(GammaT)))
  # grid points for numerical intergration
  a1 = seq(from=low[1], to=up[1], length.out = M) # Intervals for first RE
  a2 = seq(from=low[2], to=up[2], length.out = M) # Intervals for second RE
  a.m1 = (a1[-1] + a1[-M])*0.5 # midpoints interval first RE
  a.m2 = (a2[-1] + a2[-M])*0.5 # midpoints interval second RE 
  am = cbind(a.m1, a.m2)
  # correcting standardization factor
  stan.fac = (stats::pgamma(max(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(min(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
    (stats::pgamma(max(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(min(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
  probs = numeric((M-1)^2)
  for(j1 in 1:(M-1)){
    for(j2 in 1:(M-1)){
      r = (j1-1)*(M-1)+1+(j2-1)
      probs[r] = (stats::pgamma(a1[j1+1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(a1[j1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
        (stats::pgamma(a2[j2+1], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(a2[j2], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
    }}
  probs = probs/stan.fac
  IDs = unique(X$ID)
  nbAnimals = length(IDs)
  l = 0
  for(k in 1:nbAnimals){
    X_k = X[which(X$ID == IDs[k]),]
    startDD = which(X_k$condition == "DD")[1]
    nObs = nrow(X_k)
    ind = which(!is.na(X_k$activity))
    one.animal = numeric((M-1)^2)
    for(j1 in 1:(M-1)){
      for(j2 in 1:(M-1)){
        allprobs = matrix(1, nObs, N)
        for (j in 1:N){ 
          eta = c(am[j1,1], am[j2,2])
          allprobs[ind,j] = stats::dnbinom(X_k$activity[ind], size = size[j], mu = eta[j]) }
        # forward algorithm in CPP
        llk = Lcpp:::forward_cpp_flies(allprobs, delta, Gamma[,,,1], Gamma[,,,2], startDD, X_k$tod-1)
        r = (j1-1)*(M-1)+1+(j2-1)
        one.animal[r] = llk
      }}   
    ma = max(one.animal)
    logL = numeric(nObs)
    notNull = which(ma - one.animal < 700)
    logL[notNull] = exp(one.animal[notNull] - ma)
    l = l + (log(sum(logL*probs)) + ma) 
  }
  return(-l)
}

mllk_hom = function(theta.star, X, L = 48, M, low, up){
  source("./functions/periodic_variation.R")
  N = 2
  # state parameters
  Gamma_hom = array(dim = c(N, N, 2))
  for(j in 1:2){
    G = diag(N)
    G[!G] = exp(theta.star[(j-1)*N*(N-1) + 1:N*(N-1)])
    G = G / rowSums(G)
    Gamma_hom[,,j] = G }
  Gamma = array(dim = c(N, N, L, 2))
  for(j in 1:2){ # I know this is stupid, but I don't want to rewrite the CPP code
    for(t in 1:L){
      Gamma[,,t,j] = Gamma_hom[,,j] } }
  # state-dependent parameters
  ## random effects
  mu = exp(theta.star[2*(N-1)*N + 1:N])
  sigma = exp(theta.star[2*(N-1)*N + N + 1:N])
  # non-random dispersion
  size = exp(theta.star[2*(N-1)*N + 2*N + 1:N])
  # initial distribution
  delta = solve(t(diag(N)-Gamma_hom[,,1]+1), rep(1,N))
  # stuff for numerical intergration
  a1 = seq(from=low[1], to=up[1], length.out = M) # Intervals for first RE
  a2 = seq(from=low[2], to=up[2], length.out = M) # Intervals for second RE
  a.m1 = (a1[-1] + a1[-M])*0.5 # midpoints interval first RE
  a.m2 = (a2[-1] + a2[-M])*0.5 # midpoints interval second RE 
  am = cbind(a.m1, a.m2)
  stan.fac = (stats::pgamma(max(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(min(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
    (stats::pgamma(max(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(min(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
  probs = numeric((M-1)^2)
  for(j1 in 1:(M-1)){
    for(j2 in 1:(M-1)){
      r = (j1-1)*(M-1)+1+(j2-1)
      probs[r] = (stats::pgamma(a1[j1+1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(a1[j1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
        (stats::pgamma(a2[j2+1], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(a2[j2], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
    }}
  probs = probs/stan.fac
  IDs = unique(X$ID)
  nbAnimals = length(IDs)
  l = 0
  for(k in 1:nbAnimals){
    X_k = X[which(X$ID == IDs[k]),]
    startDD = which(X_k$condition == "DD")[1]
    nObs = nrow(X_k)
    ind = which(!is.na(X_k$activity))
    one.animal = numeric((M-1)^2)
    for(j1 in 1:(M-1)){
      for(j2 in 1:(M-1)){
        allprobs = matrix(1, nObs, N)
        for (j in 1:N){ 
          eta = c(am[j1,1], am[j2,2])
          allprobs[ind,j] = stats::dnbinom(X_k$activity[ind], size = size[j], mu = eta[j]) }
        llk = Lcpp:::forward_cpp_flies(allprobs, delta, Gamma[,,,1], Gamma[,,,2], startDD, X_k$tod-1)
        r = (j1-1)*(M-1)+1+(j2-1)
        one.animal[r] = llk
      }}   
    ma = max(one.animal)
    logL = numeric(nObs)
    notNull = which(ma - one.animal < 700)
    logL[notNull] = exp(one.animal[notNull] - ma)
    l = l + (log(sum(logL*probs)) + ma) 
  }
  return(-l)
}

