# Blocked Network with two Erdos-Renyi subnetwork
networkER <- function(p, p1){
  # input: p - total dim; p1 - dim of block 1
  A <- matrix(0,p,p)
  for (i in 1:p1){
    for (j in 1:p1) A[i,j] <- ifelse(rbinom(1, size=1, 0.05), ifelse(rbinom(1, size=1, 0.5), 
                                                                     runif(1, min = -0.4, max = -0.1), 
                                                                     runif(1, min = 0.1, max = 0.4)),0)
  }
  
  for (i in (1+p1):p){
    for (j in (1+p1):p) A[i,j] <-  ifelse(rbinom(1, size=1, 0.07), ifelse(rbinom(1, size=1, 0.5), 
                                                                          runif(1, min = -0.4, max = -0.1), 
                                                                          runif(1, min = 0.1, max = 0.4)),0)
  }
  
  A <- (A+t(A))
  for (i in 1:p) A[i,i] <- sum(abs(A[i,]))+0.1
  
  return(A)
}

# Scale-free network
networkSF <- function(p, m0=5, m=2){
  # 1) Add m0 nodes to G.
  # 2) Connect every node in G to every other node in G, i.e. create a complete graph.
  G <- diag(1, p, p)
  for (i in 1:m0){
    for (j in 1:m0){
      if (i !=j) G[i,j] = G[j,i] <- 1
    }
  }
  
  # 3) Create a new node i.
  d <- m0
  for (i in (m0+1):p){
    k <- 0
    edge <- NULL
    while(k<m){
      # 4) Pick a node j uniformly at random from the graph G. Set P = (k(j)/k_tot)^a.
      j <- sample(1:d, 1)
      while (j %in% edge) j <- sample(1:d, 1)
      pj <- sum(G[j,])/(sum(G)-p)
      if (pj > runif(1)) {
        G[i,j] = G[j,i] <- 1
        edge <- c(edge, j)
        k <- k+1
      }
    }
    d <- d+1
  }
  
  # input values
  for (i in 2:p){
    for (j in 1:(i-1)) G[i,j] = G[j,i] <- ifelse(G[i,j] == 1, ifelse(rbinom(1, size=1, 0.5), 
                                                                     runif(1, min = -0.4, max = -0.1), 
                                                                     runif(1, min = 0.1, max = 0.4)),0)
  }
  
  for (i in 1:p) G[i,i] <- sum(abs(G[i,]))-1+0.1
  
  return(G)
}

# Nearest-neighbor network
networkNN <- function(p, k=4){
  x <- runif(p); y <- runif(p)
  coord <- cbind(x,y)
  dist <- matrix(0, p, p)
  for (i in 2:p){
    for (j in 1:(i-1)){
      dist[i,j] <- sqrt((coord[i,1]-coord[j,1])^2 + (coord[i,2]-coord[j,2])^2)
    }
  }  
  dist <- dist + t(dist)
  B <- matrix(0, p, p)
  for (i in 1:p){
    ind <- order(dist[i,])
    B[i, ind[1:k]] <- 1
  }  
  A <- matrix(0, p, p)
  for (i in 2:p){
    for (j in 1:(i-1)) A[i,j] = A[j,i] <- ifelse(B[i,j] == 1 | B[j,i] == 1, ifelse(rbinom(1, size=1, 0.5), 
                                                                                   runif(1, min = -0.4, max = -0.1), 
                                                                                   runif(1, min = 0.1, max = 0.4)),0)
  }  
  for (i in 1:p) A[i,i] <- sum(abs(A[i,]))+0.1
  
  return(A)
}

# Banded network
networkBD <- function(p=50){ #this function needs to be updated according to the dim change
  A1 <- diag(0.5, 5)
  for (i in 2:5){
    for (j in 1:(i-1)){
      if (i-j==1) A1[i,j] <- 0.66
      if (i-j==2) A1[i,j] <- 0.3
      if (i-j==3) A1[i,j] <- 0.1
    }
  }
  A1 <- A1 + t(A1)
  
  A2 <- diag(0.5, 10)
  for (i in 2:10){
    for (j in 1:(i-1)){
      if (i-j==1) A2[i,j] <- 0.8
      if (i-j==2) A2[i,j] <- 0.5
      if (i-j==3) A2[i,j] <- 0.3
      if (i-j==4) A2[i,j] <- 0.1
    }
  }
  A2 <- A2 + t(A2)
  A3 <- diag(0.5, 5)
  for (i in 2:5){
    for (j in 1:(i-1)){
      if (i-j==1) A3[i,j] <- 0.5
      if (i-j==2) A3[i,j] <- 0.2
    }
  }
  A3 <- A3 + t(A3)
  
  B <- matrix(0, 50, 50)
  B[1:5, 1:5] = B[6:10, 6:10] <- A1
  B[11:20, 11:20] = B[21:30, 21:30] = B[31:40, 31:40] <- A2
  B[41:45, 41:45] = B[46:50, 46:50] <- A3
  
  if (p==100){
    C <- matrix(0, 100, 100)
    C[1:50, 1:50] = C[51:100, 51:100] = B
    C[46:55, 46:55] = A2
    return(C)
  } else{
    return(B)
  }
}

networkBD2 <- function(p=50){ #this function needs to be updated according to the dim change
  A1 <- diag(0.5, 5)
  for (i in 2:5){
    for (j in 1:(i-1)){
      if (i-j==1) A1[i,j] <- -0.66
      if (i-j==2) A1[i,j] <- 0.3
      if (i-j==3) A1[i,j] <- 0.1
    }
  }
  A1 <- A1 + t(A1)
  
  A2 <- diag(0.5, 10)
  for (i in 2:10){
    for (j in 1:(i-1)){
      if (i-j==1) A2[i,j] <- -0.8
      if (i-j==2) A2[i,j] <- 0.5
      if (i-j==3) A2[i,j] <- -0.3
      if (i-j==4) A2[i,j] <- 0.1
    }
  }
  A2 <- A2 + t(A2)
  
  A3 <- diag(0.5, 5)
  for (i in 2:5){
    for (j in 1:(i-1)){
      if (i-j==1) A3[i,j] <- -0.5
      if (i-j==2) A3[i,j] <- 0.2
    }
  }
  A3 <- A3 + t(A3)
  
  B <- matrix(0, 50, 50)
  B[1:5, 1:5] = B[6:10, 6:10] <- A1
  B[11:20, 11:20] = B[21:30, 21:30] = B[31:40, 31:40] <- A2
  B[41:45, 41:45] = B[46:50, 46:50] <- A3
  
  if (p==100){
    C <- matrix(0, 100, 100)
    C[1:50, 1:50] = C[51:100, 51:100] = B
    C[46:55, 46:55] = A2
    return(C)
  } else{
    return(B)
  }
}

prior.generator <- function(p, T_prob, F_prob, Network){
  # p      : is the number of nodes, i.e., the sigma matrix dim = p*p
  # T_prob : the prop of correct signals
  # F_prob : the prop of false non-signals
  # Network: the true network structure
  
  P = temp <- matrix(0, p, p)
  P <- 1 - (Network*rbinom(p*p, 1, T_prob) + (1-Network)*rbinom(p*p,1,F_prob))
  temp[upper.tri(temp)] <- P[upper.tri(P)]
  P <- temp + t(temp)
  print(isSymmetric(P))
  return(P)
}

prior.generator2 <- function(p, n_total, T_prob, network){
  P <- 1-abs(network)
  temp1 = temp2 <- matrix(0, p, p)
  n_signal = sum(abs(network)) - p
  fn <- (n_signal - n_total*T_prob)/2
  fp <- n_total*T_prob/2
  
  temp1[upper.tri(temp1)] <- P[upper.tri(P)]
  temp1[sample(which(temp1==1), fn)] <- 0
  
  temp2[upper.tri(temp2, diag = F)] <- network[upper.tri(network, diag = F)]
  temp2[sample(which(temp2==1), fp)] <- 0
  
  P <- temp1+temp2 + t(temp1+temp2)
  print(isSymmetric(P))
  return(P)
}
