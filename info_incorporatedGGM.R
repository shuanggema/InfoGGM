penSpMcp <- function(x,lambda,gamma){
  x        <- abs(x)
  sign.1   <- 1*(x<(lambda*gamma))
  val.1    <- lambda*x -x^2/(2*gamma)
  sign.2   <- 1-sign.1
  val.2    <- (lambda*gamma)^2 /(2*gamma)
  rlt      <- val.1 * sign.1 + val.2* sign.2
  return(rlt)
}

penSpMcp.Deriv <- function(x,lambda,gamma){
  x         <- abs(x)
  sign.1    <- 1*(x<(lambda*gamma))
  val.1     <- lambda -x/gamma
  rlt       <- val.1 * sign.1 
  return(rlt)
}


NormDiff <- function(xnew,xold,measure="Fnorm2",type="R"){
  if (measure=="Fnorm2"){
    num1 <- sum((xnew-xold)^2)
    num2 <- max(c(1, sum(xnew^2), sum(xold^2) ) )
  }
  if (measure=="Fnorm"){
    num1 <-  sqrt(sum((xnew-xold)^2))
    num2 <- max(c(1, sqrt(sum(xnew^2)), sqrt(sum(xold^2))))
  }
  if (measure=="abs"){
    num1 <-  sum(abs(xnew-xold))
    num2 <- sum(abs(xold))
  }
  val <- ifelse(type=="R",num1/num2,num1)
  return(val)
}

est.prior <- function(S, P, lambda, a=3, u=1, tol=1E-9, maxIter=500){
  # est.prior() is a function to compute the precision matrix given 
  # a empirical sigma matrix and prior info.
  # S      : the empirical sigma matrix
  # P      : Prior information
  # lambda : tunining parameter lambda_1
  # a, u   : parameters
  # tol    : convergent criterion
  # maxIter: the max number of iterations
  p = dim(S)[2]

  #----------- * initialize * -------------------------#
  Olist = lapply(1:2, function(m) matrix(0, p, p))
  V     = matrix(0, p, p)
  iIter = 0 
  diffVal.O = list()
  diffVal = 1000
  
  while((iIter <= (maxIter-1)) & (diffVal > tol)){
    iIter = iIter + 1
    Olist.pre = Olist
    V.pre = V
    
    ### Update O_1
    eigenComp  = eigen(u*(Olist[[2]]-V)-S)
    Di         = eigenComp$values
    vec        = eigenComp$vectors
    Di.til     = (Di + sqrt(Di^2 + 4*u))/(2*u)
    Olist[[1]] = round(vec %*% diag(Di.til) %*% t(vec), digits = 10)
    
    ### Update O_2
    temp = u*abs(Olist[[1]]+V) - penSpMcp.Deriv(abs(Olist[[2]]),lambda, a)*P
    Olist[[2]] = round(sign(Olist[[1]]+V) *  temp * (temp>0), digits = 10)
    
    ### Update V
    V = round(V.pre + Olist[[1]]-Olist[[2]], digits = 10)
    
    ## record Wlist and calculate the difference btw iterations
    diffVal.O = c(sapply(1:2, function(m) NormDiff(Olist.pre[[m]], Olist[[m]])), 
                  NormDiff(Olist[[1]], Olist[[2]]), max(Olist[[1]]-Olist[[2]]))
    diffVal           = max(diffVal.O)
  }
  if (iIter == maxIter){
    stop(paste("Lambda1 =", lambda))
  }
  else return(Olist[[2]])
}

est.composite <- function(S, S.prior, tau, lambda, a=3, u=1, tol=1E-3, maxIter=500){
  # est.composite() is a function to compute the precision matrix given 
  # a empirical sigma matrix and sigma matrix computed using prior info.
  # S      : the empirical sigma matrix
  # S.prior: the sigma matrix computed using prior info
  # tau    : tuning parameter controling the proportion of the effect of prior info
  #          when tau=0, it will not consider the prior info and equivalent to using
  #          observed data only.
  # lambda : tuning parameter lambda_2
  # a, u   : parameters
  # tol    : convergent criterion
  # maxIter: the max number of iterations
  p = dim(S)[2]
  
  #----------- * initialize * -------------------------#
  Olist = lapply(1:2, function(m) matrix(0, p, p))
  V     = matrix(0, p, p)
  iIter = 0 
  diffVal.O = list()
  diffVal = 1000
  
  while((iIter <= (maxIter-1)) & (diffVal > tol)){
    iIter = iIter + 1
    Olist.pre = Olist
    V.pre = V
    
    ### Update O_1
    eigenComp  = eigen(u*(Olist[[2]]-V)- (1-tau)*S - tau*S.prior)
    Di         = eigenComp$values
    vec        = eigenComp$vectors
    Di.til     = (Di + sqrt(Di^2 + 4*u))/(2*u)
    Olist[[1]] = round(vec %*% diag(Di.til) %*% t(vec), digits = 10)
    
    ### Update O_2
    temp = u*abs(Olist[[1]]+V) - penSpMcp.Deriv(abs(Olist[[2]]),lambda, a)*(1-diag(1,p))
    Olist[[2]] = round(sign(Olist[[1]]+V) *  temp * (temp>0), digits = 10)
    
    ### Update V
    V = round(V.pre + Olist[[1]]-Olist[[2]], digits = 10)
    
    ## record Wlist and calculate the difference btw iterations
    diffVal.O = c(sapply(1:2, function(m) NormDiff(Olist.pre[[m]], Olist[[m]])), 
                  NormDiff(Olist[[1]], Olist[[2]]), max(Olist[[1]]-Olist[[2]]))
    diffVal   = max(diffVal.O)
  }
  if (iIter == maxIter){
    stop(paste("Lambda2 =", lambda))
  }
  else return(Olist[[2]])
}

library(MASS)
library(glasso)
library(ggcorrplot)