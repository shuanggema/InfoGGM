#rm(list=ls())
#dev.off()

# Step 1: load functions
source("info_incorporatedGGM.R")
source("networks.R")
library(MASS)
library(mvtnorm)
library(glasso)
library(ggcorrplot)

# Step 2: simulation
# ---------- step 2a: set parameters ---------- #
lambda1 <- seq(0.01, 1, length.out = 10)
lambda2 <- seq(0.01,0.5, length.out = 10)
tau1 <- seq(0,1,length.out = 5)

p.1 <- 50; p <- 100; n <- 300 # these are p's and n for 100*100 matrices
tol1 <- 0.005; tol2 <- 0.005; maxIter <- 500

# ---------- step 2b: generate cov matrix ---------- #
Truth   <- networkBD(p=p) # code for banded network - positive structure, add u=0.5

ggcorrplot(Truth)
weak_rate <- 0.2 # weaken signal
Truth_weak <- Truth - Truth*weak_rate*(1-diag(1,p))
Network <- abs(sign(Truth_weak))
sum(Network)-p
sigma <- chol2inv(chol(Truth_weak))

# ---------- step 2c: generate prior info ---------- #
P1 <- 1-abs(Network)  ### 100% T prior, number of zeros
P2 <- prior.generator2(p=p, n_total=700, T_prob = 0.7, network = Network)
P3 <- prior.generator2(p=p, n_total=700, T_prob = 0.5, network = Network)
P4 <- prior.generator2(p=p, n_total=700, T_prob = 0.2, network = Network)

# replicates begin:
# set.seed(4321)

replicate <- 100
rep <- 1
Est.1 = Est.2 = Est.3 = Est.4 <- list() 
prior1 = prior2 = prior3 = prior4 <- list() 

while(rep <= replicate){
  skip_to_next <- FALSE
  print(paste("Replicate =", rep))
  X <- mvrnorm(n,rep(0,p),sigma)
  S <- var(X)
  
  print("Prior #1")
  result <- list()
  prior <- list()
  
  tryCatch(
    expr = {
      for (i in 1:length(lambda1)){
        prior[[i]] <- est.prior(S=S, P=P1, lambda=lambda1[i], tol=tol1, maxIter=maxIter)
        
        S.prior <- solve(prior[[i]])
        #S.prior <- round(S.prior, 4) #DY: 2020.02.07
        for (k in 1:length(tau1)){
          for (j in 1:length(lambda2)){
            result[[length(result)+1]] <- est.composite(S=S, S.prior=S.prior, u=0.5, tau=tau1[k], lambda=lambda2[j], tol=tol2, maxIter=maxIter)
          }
        }
      }
    },
    error = function(e){
      message('Do not converge!')
      print(e)
      skip_to_next <<- TRUE
    }
  )
  if(skip_to_next) {
    next
  }
  prior1[[rep]] <- prior
  Est.1[[rep]] <- result
  
  print("Prior #2")
  result <- list()
  prior <- list()
  
  tryCatch(
    expr = {
      for (i in 1:length(lambda1)){
        prior[[i]] <- est.prior(S=S, P=P2, lambda=lambda1[i], tol=tol1, maxIter=maxIter)   
        S.prior <- solve(prior[[i]])
        
        for (k in 1:length(tau1)){
          for (j in 1:length(lambda2)){
            result[[length(result)+1]] <- est.composite(S=S, S.prior=S.prior, u=0.5, tau=tau1[k], lambda=lambda2[j], tol=tol2, maxIter=maxIter)
          }
        }
      }
    },
    error = function(e){
      message('Do not converge!')
      print(e)
      skip_to_next <<- TRUE
    }
  )
  if(skip_to_next) {
    prior1[[rep]] <- NULL
    Est.1[[rep]] <- NULL
    next
  }
  prior2[[rep]] <- prior
  Est.2[[rep]] <- result
  
  print("Prior #3")
  result <- list()
  prior <- list()
  
  tryCatch(
    expr = {
      for (i in 1:length(lambda1)){
        prior[[i]] <- est.prior(S=S, P=P3, lambda=lambda1[i], tol=tol1, maxIter=maxIter)   
        S.prior <- solve(prior[[i]])
        
        for (k in 1:length(tau1)){
          for (j in 1:length(lambda2)){
            result[[length(result)+1]] <- est.composite(S=S, S.prior=S.prior, u=0.5, tau=tau1[k], lambda=lambda2[j], tol=tol2, maxIter=maxIter)
          }
        }
      }
    },
    error = function(e){
      message('Do not converge!')
      print(e)
      skip_to_next <<- TRUE
    }
  )
  if(skip_to_next) {
    prior1[[rep]] = prior2[[rep]] <- NULL
    Est.1[[rep]] = Est.2[[rep]] <- NULL
    next
  }
  prior3[[rep]] <- prior
  Est.3[[rep]] <- result
  
  print("Prior #4")
  result <- list()
  prior <- list()
  
  tryCatch(
    expr = {
      for (i in 1:length(lambda1)){
        prior[[i]] <- est.prior(S=S, P=P4, lambda=lambda1[i], tol=tol1, maxIter=maxIter)   
        S.prior <- solve(prior[[i]])
        
        for (k in 1:length(tau1)){
          for (j in 1:length(lambda2)){
            result[[length(result)+1]] <- est.composite(S=S, S.prior=S.prior, u=0.5, tau=tau1[k], lambda=lambda2[j], tol=tol2, maxIter=maxIter)
          }
        }
      }
    },
    error = function(e){
      message('Do not converge!')
      print(e)
      skip_to_next <<- TRUE
    }
  )
  if(skip_to_next) {
    prior1[[rep]] = prior2[[rep]] = prior3[[rep]] <- NULL
    Est.1[[rep]] = Est.2[[rep]] = Est.3[[rep]] <- NULL
    next
  }
  prior4[[rep]] <- prior
  Est.4[[rep]] <- result
  
  rep <- rep + 1
}
