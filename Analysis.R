
# Packages ----------------------------------------------------------------

library(optimParallel)
library(scales)
library(dplyr)


# install.packages("devools")
# to install Lcpp, you need to have a functioning C++ compiler. See
# https://teuder.github.io/rcpp4everyone_en/020_install.html for details
# devtools::install_github("janoleko/Lcpp")
library(Lcpp)


# Data loading and colors -------------------------------------------------

load("fruitflies.RData")

## defining colors
color = c("orange", "deepskyblue")


# Exploratory data analysis -----------------------------------------------

dataLD = data[which(data$condition=="LD"),]
dataDD = data[which(data$condition=="DD"),]

medianLD = dataLD %>% group_by(tod) %>% summarise(median = median(activity, na.rm=T))
medianDD = dataDD %>% group_by(tod) %>% summarise(median = median(activity, na.rm=T))

colbp = "deepskyblue"

# pdf("./figures/boxplot_tod.pdf", width = 8, height = 4)
par(mfrow = c(1,2), mar = c(4.5, 4, 2.5, 0.5) + 0.1)
boxplot(dataLD$activity ~ dataLD$tod, lwd = .7, ylab = "activity count", 
        xlab = "time of day", xaxt = "none", frame = F, col = "gray97", outline = F, 
        main = "LD", ylim = c(0,200))
axis(side = 1, at = c(0,4,8,12,16,20,24)*2, labels = c(0,4,8,12,16,20,24))
lines(medianLD, lwd = 3, col = alpha(colbp, 0.7))
par(xpd = T)
polygon(c(0:48, 48:0), c(rep(-8, 49), rep(-4, 49)),col = "black", border = "black")
polygon(c(16:40, 40:16), c(rep(-8, 25), rep(-4, 25)),col = "white", border = "black")
par(xpd = F)

boxplot(dataDD$activity ~ dataDD$tod, lwd = .7, ylab = "activity count", 
        xlab = "time of day", xaxt = "none", frame = F, col = "gray97", outline = F, 
        main = "DD", ylim = c(0,200))
axis(side = 1, at = c(0,4,8,12,16,20,24)*2, labels = c(0,4,8,12,16,20,24))
lines(medianDD, lwd = 3, col = alpha(colbp, 0.7))
par(xpd = T)
polygon(c(0:48, 48:0), c(rep(-8, 49), rep(-4, 49)),col = "black", border = "black")
polygon(c(16:40, 40:16), c(rep(-8, 25), rep(-4, 25)),col = "gray45", border = "black")
par(xpd = F)
# dev.off()

# summary statistics

round(median(data$activity))
round(mean(data$activity))
max(data$activity)


# Fitting homogeneous random effects model --------------------------------

# likelihood functions and helper function for trigonometric link
source("./functions/periodic_variation.R")
source("./functions/likelihoods_RE.R")

# initial parameter value
theta.star = c(rep(-1,4), # Gamma
               log(c(4, 55)), log(c(2,15)), # mu and sigma for random effects
               log(c(0.1, 1.7))) # dispersion parameter for negbin

# practical range for numerical integration
low = c(0.001, 15)
up = c(20, 130)

# homogeneous model fit
# cl = makeCluster(detectCores()-1); setDefaultCluster(cl=cl); t1 = Sys.time()
# mod_hom = optimParallel(theta.star, mllk_hom, X = data, M = 25, low = low, up = up,
#                             parallel = list(forward = T),
#                             control = list(trace = 4, maxit = 1e4, ndeps=rep(1e-5, length(theta.star))))
# Sys.time()-t1; stopCluster(cl)
# saveRDS(mod_hom, "./models/mod_hom.rds")

mod_hom = readRDS("./models/mod_hom.rds")

# AIC and BIC
(AIC_hom = 2*mod_hom$value + 2*length(mod_hom$par))
(BIC_hom = 2*mod_hom$value + log(nrow(data))*length(mod_hom$par))

theta.hat.hom = mod_hom$par

# Fitting inhomogeneous random effects model ------------------------------

# initial parameter value
theta.star = c(5, -6, 3, -5, 8, -6, 5, -2, 12, -7, 8, -2, 1, 0, -2, -5, 0, -9, 
               0, -8, -1, -2, -1, -7, 2, -4, 1, -1, 2, 4, 1, 3, -3, 1)

# model fit
# cl = makeCluster(detectCores()-1); setDefaultCluster(cl=cl); t1 = Sys.time()
# mod = optimParallel(theta.star, mllk, X = data, K = 3, M = 27, low = low, up = up,
#                          parallel = list(forward=T), 
#                     control = list(trace=4, maxit=1e4, ndeps=rep(1e-7, length(theta.star))))
# Sys.time()-t1; stopCluster(cl)

# nlm() is more accurate with respect to the stopping criterium
# therefore we finish estimation by nlm()
# mod_final = nlm(mllk, mod$par, X = data, K = 3, M = 27, low = low, up = up, 
#                 print.level = 2, iterlim = 1000)
# saveRDS(mod_final, "./models/mod_RE3.rds")

mod_final = readRDS("./models/mod_RE3.rds")

# calculating the approximate Hessian at parameter estimate
# H3 = numDeriv::hessian(mllk, mod_final$estimate, X = data, M = 27, low = low, up = up,
#                              method = "Richardson", method.args = list(
#                                eps=1e-7, d=0.0001, zero.tol = sqrt(.Machine$double.eps/7e-7), 
#                                r=12, v=2, show.details=TRUE))
# saveRDS(H3, "./models/H_RE3.rds") 

H3 = readRDS("./models/H_RE3.rds")

# inverting the Hessian to get Fisher information
I3 = solve(H3, tol = 1e-30)

# AIC and BIC
(AIC = 2*mod_final$minimum + 2*length(mod_final$estimate))
(BIC = 2*mod_final$minimum + log(nrow(data))*length(mod_final$estimate))

theta.hat = mod_final$estimate

# backtransform estimated parameters
N = 2; K = 3; L = 48
coef = array(theta.hat[1:(2*(1+2*K)*N*(N-1))], dim = c(N*(N-1), 1+2*K, 2))
coef_hat = coef[,c(1,2,5,3,6,4,7),] # reorder columns: sin1, cos1, sin2, cos2, ...

## state process
Gamma_hat = array(dim = c(N,N,L,2))
Gamma_hat[,,,1] = tpm_p(0:47, L = 48, beta = coef_hat[,,1], degree = 3)
Gamma_hat[,,,2] = tpm_p(0:47, L = 48, beta = coef_hat[,,2], degree = 3)
delta_hat = stationary_p(Gamma_hat[,,,1], t = data$tod[1])

## state-dependent process
### random effects
mu_hat = exp(theta.hat[(2*(1+2*K)*(N-1)*N) + 1:N])
sigma_hat = exp(theta.hat[(2*(1+2*K)*(N-1)*N) + N + 1:N])
### fixed effects
size_hat = exp(theta.hat[(2*(1+2*K)*(N-1)*N) + 2*N + 1:N])


# Results -----------------------------------------------------------------

# State-dependent parameter distributions

par(mfrow = c(1,1))
curve(dgamma(x, shape = mu_hat[1]^2/sigma_hat[1]^2, scale = sigma_hat[1]^2/mu_hat[1]), 
      xlim = c(0, 140), col = color[1], bty = "n", n = 500, lwd = 2, xlab = "eta", ylab = "density")
curve(dgamma(x, shape = mu_hat[2]^2/sigma_hat[2]^2, scale = sigma_hat[2]^2/mu_hat[2]), 
      add = T, col = color[2], n = 500, lwd = 2)


# Periodic stationary distribution (delta) --------------------------------

source("./functions/auxiliary_functions.R")

# # calculating deltas
Delta = Rho = array(NA, dim = c(L,N,2))
for(t in 1:L){
  for(cond in 1:2){
    Delta[t,,cond] = get_delta(theta.hat, t=t, state=1:2, cond=cond)
    Rho[t,,cond] = get_rho(theta.hat, t=t, state=1:2, cond=cond)
  }
}
n = 200 # smoothness of line in plot
timeseq = seq(0,24, length=n)
Delta_cont = Rho_cont = array(NA, dim = c(n,N,2))
for(k in 1:n){
  for(cond in 1:2){
    Delta_cont[k,,cond] = get_delta_cont(theta.hat, t=timeseq[k], state=1:2, cond=cond)
    Rho_cont[k,,cond] = get_rho_cont(theta.hat, t=timeseq[k], state=1:2, cond=cond)
  }
}

# calculate confidence intervals via Monte Carlo, based on approximate normal distribution of MLE
# you can also just load the confidence intervals:

CIs = readRDS("./models/confidence_intervals.rds")
DeltaCI = CIs$DeltaCI
RhoCI = CIs$RhoCI

n = 2000 # this takes some time
timeseq = seq(0,24, length=200)
set.seed(123)
thetas = mvtnorm::rmvnorm(n, mean = theta.hat, sigma = I3)
Deltas = Rhos = array(dim = c(200,N,2,n))
for(t in 1:200){
  print(t)
  for(cond in 1:2){
    for(i in 1:n){
      Deltas[t,,cond,i] = get_delta_cont(thetas[i,], t=timeseq[t], state = 1:2, cond = cond)
      Rhos[t,,cond,i] = get_rho_cont(thetas[i,], t=timeseq[t], state = 1:2, cond = cond)
    }
  }
}
DeltaCI = RhoCI = array(dim = c(200,N,2,2))
for(t in 1:200){
  for(cond in 1:2){
    for(state in 1:2){
      DeltaCI[t,state,cond,] = quantile(Deltas[t,state,cond,], probs = c(0.025, 0.975))
      RhoCI[t,state,cond,] = quantile(Rhos[t,state,cond,], probs = c(0.025, 0.975))
    }
  }
}
# CIs = list(DeltaCI = DeltaCI, RhoCI = RhoCI)
# saveRDS(CIs, "./models/confidence_intervals.rds")


# pdf("./figures/p_stationary.pdf", width=7, height=3.8)

colrho = "gold2"
colbox = c("white", "gray45")
conditions = c("LD", "DD")
a = 0.3 # alpha value

m = matrix(c(1,1,2,3),nrow = 2, ncol = 2,byrow = TRUE)
layout(mat = m, heights = c(0.13, 1, 1))
par(mar = c(0,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
state = 2; colrho = "gold2"
legend("top", inset=0, 
       legend = c(expression(delta), expression(rho)), pch = 16, col = c(color[state], colrho), lwd = 2,
       horiz=T, bty = "n", text.width = c(0, 0.03))

par(mar = c(5.5,4,2,1), xpd = F)

state = 2 # only for high-activity state
cond = 1 # LD condition

for(cond in 1:2){
  ## rho
  plot(timeseq, Rho_cont[,state,cond], type = "l", lwd = 1, col = colrho, bty = "n", ylim = c(0,1), 
       ylab = "Pr(high-activity state)", xlab = "time of day", xaxt = "n", main = conditions[cond])
  points(0:47/2, Rho[,state,cond], pch = 16, col = colrho)
  polygon(c(timeseq, rev(timeseq)), c(RhoCI[,state,cond,1], rev(RhoCI[,state,cond,2])), col = alpha(colrho, a), border = F)
  
  ## delta
  lines(timeseq, Delta_cont[,state,cond], col = color[state], lwd = 1)
  points(1:48/2, Delta[,state,cond], pch = 16, col = color[state])
  polygon(c(timeseq, rev(timeseq)), c(DeltaCI[,state,cond,1], rev(DeltaCI[,state,cond,2])), col = alpha(color[2], a), border = F)
  
  ## axes
  axis(1, at=seq(0, 24, by = 4), label = seq(0,24,by = 4))
  par(xpd = T)
  polygon(c(0:24, 24:0), c(rep(-0.025, 25), rep(-0.045, 25)), col = "black", border = "black")
  polygon(c(8:20, 20:8), c(rep(-0.025, 13), rep(-0.045, 13)), col = colbox[cond], border = "black")
  par(xpd = F)
}

# dev.off()



# Time-varying dwell-time distribution ------------------------------------

source("./functions/dwell_distr.R")

timepoints = (1:8)*3 + 20 # select only a few time points for time-varying plot

### calculating all time-varying dwell-time distributions
PMF = array(dim = c(48, 25, 2, 2))
for(j in 1:2){
  for(cond in 1:2){
    for(t in 1:48){ PMF[t,,j,cond] = ddwell_t(1:25, t = t, state = j, Gamma = Gamma_hat[,,,cond]) }
  }
}


# pdf("./figures/time-varying_distr.pdf", width=8.5, height=4)

## plotting 8 time varying distributions
par(mfrow = c(2,4), mar = c(4.5, 4, 1, 0.5) + 0.1)
state = 2 # high-activity state
cond = 1 # LD condition
for(i in ((1:8)*2 + 14)){
  if(i == 16 | i == 24){
    yLab = "probabilities"
  } else{ yLab = "" }
  if(i %in% c(16,18,20,22)){
    xLab = ""
  } else{ xLab = "dwell time (in hours)" }
  plot(PMF[i,,state,cond][1:24], type = "h", bty = "n", col = color[state], lwd = 2, xaxt = "n", xlim = c(0,24), ylim = c(0, 0.35),
       ylab = yLab, xlab = xLab)
  axis(1, at = seq(0, 24, by = 2), label = seq(0,12, by = 1)) # dwell time in hours
  legend("top", bty = "n", legend = paste0(i/2, ":00"))
}

# dev.off()

## plot dynamically
timelab = paste0(rep(0:23,each=2), ":", c("00", "30"))
par(mfrow = c(1,1))
state = 2 # high-activity state
cond = 2 # DD condition
for(t in 1:48){
  plot(1:24, PMF[t,1:24,state,cond], type = "h", lwd = 2, col = color[state], bty = "n", 
       ylim = c(0,0.7), xlim = c(0,24), xaxt = "n", 
       ylab = "time-varying distribution", xlab = "dwell time", main = timelab[t])
  axis(1, at = seq(0, 24, by = 2), label = seq(0,12, by = 1)) # dwell time in hours
  Sys.sleep(0.15)
}

# time-varying mean dwell times

# calculate the time-varying mean dwell times (and standard errors)
ET = array(dim = c(L, 2, 2))
for(j in 1:2){
  for(cond in 1:2){
    for(t in 1:L){ 
      ET[t,j,cond] = get_mean(theta.hat, t=t, state = j, cond = cond) 
    }
  }
}

# Monte carlo based confidence intervals
n = 2000 # this takes some time
ETs = array(dim = c(L,N,2,n))
for(t in 1:L){
  print(t)
  for(cond in 1:2){
    for(state in 1:2){
      for(i in 1:n){
        ETs[t,state,cond,i] = get_mean(thetas[i,], t=t, state = state, cond = cond) 
      }
    }
  }
}
ETCI = array(dim = c(L,N,2,2))
for(t in 1:L){
  for(cond in 1:2){
    for(state in 1:2){
      ETCI[t,state,cond,] = quantile(ETs[t,state,cond,], probs = c(0.025, 0.975))
    }
  }
}

# pdf("./figures/mean_dwell_times.pdf", width=8.5, height=4)

m = matrix(c(1,1,2,3),nrow = 2, ncol = 2,byrow = TRUE)
layout(mat = m, heights = c(0.12, 1, 1))
par(mar = c(0,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top", inset=0, 
       legend = c("low-activity state", "high-activity state"), pch = 16, lwd = 2, col = color,
       horiz=T, bty = "n", text.width = c(0, 0.15))

par(mar = c(4.5,4,2.5,1))

for (cond in 1:2){
  plot(1:L, ET[,1,cond], type = "l", bty = "n", xaxt = "n", ylim = c(0,22), lwd = 2, col = color[1],
       ylab = "mean dwell time", xlab = "time of day")
  points(1:L, ET[,1,cond], pch = 16, col = color[1])
  polygon(c(1:L,rev(1:L)),c(ETCI[,1,cond,1], rev(ETCI[,1,cond,2])), col = alpha(color[1], 0.3), border = FALSE)
  
  lines(1:L, ET[,2,cond], lwd = 2, col = color[2])
  points(1:L, ET[,2,cond], pch = 20, col = color[2])
  polygon(c(1:L,rev(1:L)),c(ETCI[,2,cond,1], rev(ETCI[,2,cond,2])), col = alpha(color[2], 0.3), border = FALSE)
  
  legend("top", bty = "n", legend = conditions[cond])
  axis(1, at=seq(0, 48, by = 8), label = seq(0,24,by = 4))
}

# dev.off()



# State decoding ----------------------------------------------------------

# viterbi decoding is not possible in random effects model.
# here we use local (forward) decoding

source("./functions/state_decoding.R")

M = 40
probs = array(dim = c(nrow(data), 2, M, M))
low = c(0.001, 15)
up = c(20, 130)
a1 = seq(from=low[1], to=up[1], length.out = M) # Intervals for first RE
a2 = seq(from=low[2], to=up[2], length.out = M) # Intervals for second RE
a.m1 = (a1[-1] + a1[-M])*0.5 # midpoints interval first RE
a.m2 = (a2[-1] + a2[-M])*0.5 # midpoints interval second RE 
mu = mu_hat; sigma = sigma_hat
am = cbind(a.m1, a.m2)
stan.fac = (stats::pgamma(max(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(min(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
  (stats::pgamma(max(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(min(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))

for(j1 in 1:M){
  print(j1)
  for(j2 in 1:M){
    prob = (stats::pgamma(a1[j1+1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(a1[j1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
      (stats::pgamma(a2[j2+1], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(a2[j2], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
    eta = c(a.m1[j1], a.m2[j2])
    phis = cond_stateProbs(Gamma_hat, delta_hat, eta, size_hat, X = data)
    probs[,,j1,j2] = phis
  }
}
overallprobs = apply(probs, MARGIN = c(1,2), FUN = mean, na.rm = T)
states = apply(overallprobs, MARGIN = 1, FUN = which.max)


# plotting the time series colored by decoded states

# pdf("./figures/timeseries_decoded.pdf", width=9, height=4)

par(mfrow = c(1,1))
ind = c(1:400)
plot(data$activity[ind], bty = "n", type = "h", col = alpha(color[states[ind]],0.3), lwd = 1, xlab = "time", ylab = "activtiy count")
points(data$activity[ind], pch = 20, col = alpha(color[states[ind]],0.7))

# dev.off()


# decoding states for homogeneous model

theta.hat.hom = mod_hom$par
N = 2
Gamma_hom = array(dim = c(N, N, 2))
for(j in 1:2){
  G = diag(N)
  G[!G] = exp(theta.hat.hom[(j-1)*N*(N-1) + 1:N*(N-1)])
  G = G / rowSums(G)
  Gamma_hom[,,j] = G }
mu_hom = exp(theta.hat.hom[2*(N-1)*N + 1:N])
sigma_hom = exp(theta.hat.hom[2*(N-1)*N + N + 1:N])
size_hom = exp(theta.hat.hom[2*(N-1)*N + 2*N + 1:N])
delta_hom = solve(t(diag(N)-Gamma_hom[,,1]+1), rep(1,N))

M = 40
probs_hom = array(dim = c(nrow(data), 2, M, M))
low = c(0.001, 15)
up = c(20, 130)
a1 = seq(from=low[1], to=up[1], length.out = M) # Intervals for first RE
a2 = seq(from=low[2], to=up[2], length.out = M) # Intervals for second RE
a.m1 = (a1[-1] + a1[-M])*0.5 # midpoints interval first RE
a.m2 = (a2[-1] + a2[-M])*0.5 # midpoints interval second RE 
mu = mu_hom; sigma = sigma_hom
am = cbind(a.m1, a.m2)
stan.fac = (stats::pgamma(max(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(min(a1), shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
  (stats::pgamma(max(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(min(a2), shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))

for(j1 in 1:M){
  print(j1)
  for(j2 in 1:M){
    prob = (stats::pgamma(a1[j1+1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1])-stats::pgamma(a1[j1], shape=mu[1]^2/sigma[1]^2,scale=sigma[1]^2/mu[1]))*
      (stats::pgamma(a2[j2+1], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2])-stats::pgamma(a2[j2], shape=mu[2]^2/sigma[2]^2,scale=sigma[2]^2/mu[2]))
    eta = c(a.m1[j1], a.m2[j2])
    phis = cond_stateProbs_hom(Gamma_hom, delta_hom, eta, size_hom, X = data)
    probs_hom[,,j1,j2] = phis
  }
}
overallprobs_hom = apply(probs_hom, MARGIN = c(1,2), FUN = mean, na.rm = T)
states_hom = apply(overallprobs_hom, MARGIN = 1, FUN = which.max)

# decoded = list(state_probs = overallprobs, dec_states = states, 
#                state_probs_hom = overallprobs_hom, dec_states_hom = states_hom)
# saveRDS(decoded, "./models/decoded.rds")

# Overall dwell-time distribution -----------------------------------------

PMFO = array(dim = c(50, 2, 2)) # support up to 50 is enough
for(state in 1:2){
  for(cond in 1:2){ PMFO[,state,cond] = ddwell(1:50, state = state, Gamma = Gamma_hat[,,,cond]) }
}

# overall mean dwell times
m1 = sum((1:50) * PMFO[,2,1]) / 2 # in hours
m2 = sum((1:50) * PMFO[,2,2]) / 2 # in hours
round(m1,3)
round(m2,3)

# pdf("./figures/overall_distr.pdf", width=8.5, height=4)

m = matrix(c(1,1,2,3),nrow = 2, ncol = 2,byrow = TRUE)
layout(mat = m, heights = c(0.15, 1, 1))
par(mar = c(0,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

state = 2
legend("top", inset=0, 
       legend = c("model-implied", "empirical"), pch = c(NA, 19), lwd = c(4, NA), col = c(color[state], alpha("black",0.3)),
       horiz=T, bty = "n", text.width = c(0, 0.2))

par(mar = c(4.5,4,2,1))

cond = 1
lim = 24
plot(1:lim, PMFO[,state,cond][1:lim], type = "h", lwd = 4, bty = "n", col = color[state], xaxt = "n",
     xlab = "dwell time (in hours)", ylab = "probabilities", ylim = c(0,0.4), xlim = c(0,lim))
axis(1, at = seq(0, lim, by = 2), label = seq(0, lim/2, by = 1)) # dwell time in hours

# empirical
states_sub = states[which(data$condition == conditions[cond])]
emp_pmf = prop.table(table(get_dwell_times(states_sub, state)))
points(as.numeric(names(emp_pmf)), as.numeric(emp_pmf), pch = 19, col = alpha(1, 0.3))
# problem: Don't know how to even decode the states in RE model
legend("top", legend = c("LD"), bty = "n")

cond = 2
lim = 24
plot(1:lim, PMFO[,state,cond][1:lim], type = "h", lwd = 4, bty = "n", col = color[state], xaxt = "n",
     xlab = "dwell time (in hours)", ylab = "probabilities", ylim = c(0,0.4), xlim = c(0,lim))
axis(1, at = seq(0, lim, by = 2), label = seq(0, lim/2, by = 1)) # dwell time in hours

# empirical
states_sub = states[which(data$condition == conditions[cond])]
emp_pmf = prop.table(table(get_dwell_times(states_sub, state)))
points(as.numeric(names(emp_pmf)), as.numeric(emp_pmf), pch = 19, col = alpha(1, 0.3))
# problem: Don't know how to even decode the states in RE model
legend("top", legend = c("DD"), bty = "n")
# legend("topright", bty = "n", legend = c("model-implied", "empirical"), pch = c(NA, 19), lwd = c(4, NA), col = c(color[state], alpha("black",0.3)))

# dev.off()


# Overall dwell-time distribution homogeneous model -----------------------

# pdf("./figures/overall_distr_hom.pdf", width=8.5, height=4)

m = matrix(c(1,1,2,3),nrow = 2, ncol = 2,byrow = TRUE)
layout(mat = m, heights = c(0.12, 1, 1))
par(mar = c(0,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

state = 2
legend("top", inset=0, 
       legend = c("model-implied distribution", "empirical distribution"), pch = c(NA, 19), lwd = c(4, NA), col = c(color[state], alpha("black",0.3)),
       horiz=T, bty = "n", text.width = c(0, 0.2))

par(mar = c(4.5,4,2,1))

cond = 1
lim = 24
plot(1:lim, dgeom((1:lim)-1, prob = 1-diag(Gamma_hom[,,cond])[state]), type = "h", lwd = 4, bty = "n", col = color[state], xaxt = "n",
     xlab = "dwell time (in hours)", ylab = "probabilities", ylim = c(0,0.4), xlim = c(0,lim))
axis(1, at = seq(0, lim, by = 2), label = seq(0, lim/2, by = 1)) # dwell time in hours

# empirical
states_sub_hom = states_hom[which(data$condition == conditions[cond])]
emp_pmf_hom = prop.table(table(get_dwell_times(states_sub_hom, state)))
points(as.numeric(names(emp_pmf_hom)), as.numeric(emp_pmf_hom), pch = 19, col = alpha(1, 0.3))
# problem: Don't know how to even decode the states in RE model
legend("top", legend = c("LD"), bty = "n")

cond = 2
lim = 24
plot(1:lim, dgeom((1:lim)-1, prob = 1-diag(Gamma_hom[,,cond])[state]), type = "h", lwd = 4, bty = "n", col = color[state], xaxt = "n",
     xlab = "dwell time (in hours)", ylab = "probabilities", ylim = c(0,0.4), xlim = c(0,lim))
axis(1, at = seq(0, lim, by = 2), label = seq(0, lim/2, by = 1)) # dwell time in hours

# empirical
states_sub_hom = states_hom[which(data$condition == conditions[cond])]
emp_pmf_hom = prop.table(table(get_dwell_times(states_sub_hom, state)))
points(as.numeric(names(emp_pmf_hom)), as.numeric(emp_pmf_hom), pch = 19, col = alpha(1, 0.3))
# problem: Don't know how to even decode the states in RE model
legend("top", legend = c("DD"), bty = "n")
# legend("topright", bty = "n", legend = c("model-implied", "empirical"), pch = c(NA, 19), lwd = c(4, NA), col = c(color[state], alpha("black",0.3)))

# dev.off()



# State-dependent distribution --------------------------------------------

# pdf("./figures/state_dependent.pdf", width=8, height=4)

delta_overall = apply(Delta, 2, mean)
nSample = 200 # this can be increased
eta = matrix(NA, nrow = nSample, ncol = 2)
eta[,1] = rgamma(nSample, shape=mu_hat[1]^2/sigma_hat[1]^2,scale=sigma_hat[1]^2/mu_hat[1])
eta[,2] = rgamma(nSample, shape=mu_hat[2]^2/sigma_hat[2]^2,scale=sigma_hat[2]^2/mu_hat[2])
lw = 1

par(mar = c(4.5,4,2,1), mfrow = c(1,1))
plot(delta_overall[1]*dnbinom(0:200, size = size_hat[1], mu = eta[1,1]), type = "l", 
     lwd = lw, col = alpha(color[1], 0.01), bty = "n", ylab = "probabilities", xlab = "activity count", ylim = c(0,0.02))
lines(delta_overall[2]*dnbinom(0:200, size = size_hat[2], mu = eta[1,2]), lwd = lw, col = alpha(color[2], 0.01))

for(k in 2:nSample){
  points(delta_overall[1]*dnbinom(0:200, size = size_hat[1], mu = eta[k,1]), pch = 20, col = alpha(color[1], 0.01))
  points(delta_overall[2]*dnbinom(0:200, size = size_hat[2], mu = eta[k,2]), pch = 20, col = alpha(color[2], 0.01))
}
legend("topright", bty = "n", lwd = 1, col = color, legend = c("low-activity state", "high activity state"))

# dev.off()

