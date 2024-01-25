get_delta = function(theta.star, t, state, N=2, K=3, cond){
  coef = array(theta.star[1:(2*(1+2*K)*N*(N-1))], dim = c(N*(N-1), 1+2*K, 2))
  coef_star = coef[,c(1,2,5,3,6,4,7),] # reorder columns
  Gamma = Lcpp::tpm_p(0:47, L = 48, beta = coef_star[,,cond], degree = K)
  delta = Lcpp::stationary_p(Gamma, t)
  return(delta[state])
}
delta_method = function(theta.star, t, state, N=2, I, K=3, cond){
  gradient = numDeriv::grad(get_delta, theta.star, t=t, state=state, N=N, K=K, cond=cond)
  # method = "Richardson", method.args=list(eps=1e-6, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=12, v=2, show.details=FALSE))
  Var = t(gradient)%*%I%*%gradient
  return(Var)
}
get_delta_cont = function(theta.star, t, state, N=2, L=48, K=3, cond){
  coef = array(NA, dim = c(N*(N-1), 1+2*K, 2))
  coef[,,1] = matrix(theta.star[1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  coef[,,2] = matrix(theta.star[((1+2*K)*(N-1)*N) + 1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  timeseq = (t + (0:47)/2) %% 24
  Gamma = array(NA, dim = c(N,N,L))
  for (k in 1:L){
    eta = pv(coef[,,cond], time = timeseq[k], degree = K, L = 24)
    G = diag(N)
    G[!G] = exp(eta)
    G = G/rowSums(G)
    Gamma[,,k] = G }
  GammaT=Gamma[,,1]
  for(k in 2:(L-1)){
    GammaT = GammaT%*%Gamma[,,k] }
  delta = solve(t(diag(N)-GammaT+1), rep(1,N))
  return(delta[state])
}
delta_method_cont = function(theta.star, t, state, N=2, L=48, I, K=3, cond){
  gradient = numDeriv::grad(get_delta_cont, theta.star, t=t, state=state, N=N, L=L, K=K, cond=cond,
                            method = "Richardson", method.args=list(eps=1e-10, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=12, v=2, show.details=FALSE))
  Var = t(gradient)%*%I%*%gradient
  return(Var)
}

get_rho = function(theta.star, t, state, N=2, K=3, cond){
  coef = array(theta.star[1:(2*(1+2*K)*N*(N-1))], dim = c(N*(N-1), 1+2*K, 2))
  coef_star = coef[,c(1,2,5,3,6,4,7),] # reorder columns
  Gamma = Lcpp::tpm_p(0:47, L = 48, beta = coef_star[,,cond], degree = K)
  rho = solve(t(diag(N)-Gamma[,,t]+1), rep(1,N))
  return(rho[state])
}
delta_method_rho = function(theta.star, t, state, N=2, I, K=3, cond){
  gradient = numDeriv::grad(get_rho, theta.star, t=t, state=state, N=N, K=K, cond=cond)
  # method = "Richardson", method.args=list(eps=1e-6, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=12, v=2, show.details=FALSE))
  Var = t(gradient)%*%I%*%gradient
  return(Var)
}
get_rho_cont = function(theta.star, t, state, N=2, K=3, cond){
  coef = array(NA, dim = c(N*(N-1), 1+2*K, 2))
  coef[,,1] = matrix(theta.star[1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  coef[,,2] = matrix(theta.star[((1+2*K)*(N-1)*N) + 1:((1+2*K)*N*(N-1))], nrow = N*(N-1), ncol = (1+2*K))
  eta = pv(coef[,,cond], time = t, degree = K, L = 24)
  Gamma = diag(N)
  Gamma[!Gamma] = exp(eta)
  Gamma = Gamma / rowSums(Gamma)
  rho = solve(t(diag(N)-Gamma+1), rep(1,N))
  return(rho[state])
}
delta_method_rho_cont = function(theta.star, t, state, N=2, I, K=3, cond){
  gradient = numDeriv::grad(get_rho_cont, theta.star, t=t, state=state, N=N, K=K, cond=cond)
  # method = "Richardson", method.args=list(eps=1e-6, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=12, v=2, show.details=FALSE))
  Var = t(gradient)%*%I%*%gradient
  return(Var)
}

get_mean = function(theta.star, N=2, K=3, t, state, cond){
  coef = array(theta.star[1:(2*(1+2*K)*N*(N-1))], dim = c(N*(N-1), 1+2*K, 2))
  coef_star = coef[,c(1,2,5,3,6,4,7),] # reorder columns
  Gamma = Lcpp::tpm_p(0:47, L = 48, beta = coef_star[,,cond], degree = K)
  # calculating the mean as explained in paper
  pmf = ddwell_t(1:L, t = t, state = state, Gamma = Gamma)
  E = (L + sum(1:L * pmf)) / sum(pmf) - L 
  return(E)
}
delta_method_mean = function(theta.star, t, state, N=2, L=48, I, K=3, cond){
  gradient = numDeriv::grad(get_mean, theta.star, t=t, state=state, N=N, K=K, cond=cond)
  Var = t(gradient)%*%I%*%gradient
  return(Var)
}

get_dwell_times = function(states, k){
  ind = which(states == k)
  dwell = 0
  counter = 1
  for (i in 1:(length(ind)-1)){
    if(ind[i+1] == ind[i]+1){
      counter = counter+1
    } else{
      dwell = c(dwell, counter)
      counter = 1
    }
  }
  return(dwell[2:length(dwell)])
}