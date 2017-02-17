Forward <- function(Obs, Ini, A, E){
  # Takes a sequence of observations, calculate forward probability matrix
  #
  # Args: 
  #   Obs: a numeric vector
  #   Ini: q_start to q_1 prob distribution
  #   A: Trans matrix
  #   E: Emission matrix, e_q\sigma
  #
  # Output: Forward matrix: F_qi, row is state, column is time index
  nQ <- dim(A)[1]
  n <- length(Obs)
  a <- matrix(Ini * E[, Obs[1]], nQ, 1)
  for (i in 2:n){
    temp <- a[, i-1] * A * (rep(1, nQ) %*% t(E[, Obs[i]]))
    a <- cbind(a, colSums(temp))
  }
  return(a)
}

Backward <- function(Obs, A, E){
  # Takes a sequence of observations, calculate backward probability matrix
  #
  # Args: 
  #   Obs: a numeric vector
  #   A: Trans matrix
  #   E: Emission matrix, e_q\sigma
  #
  # Output: Backward matrix: F_qi, row is state, column is time index
  nQ <- dim(A)[1]
  n <- length(Obs)
  b <- matrix(rep(1, nQ), nQ, 1)
  for (i in (n - 1):1){
    temp <- A * (rep(1, nQ) %*% t(E[, Obs[i + 1]])) * (rep(1, nQ) %*% t(b[, 1]))
    b <- cbind(rowSums(temp), b)
  }
  return(b)
}

SoftCoding <- function(Obs, Ini, A, E){
  # Takes a sequence of observations, calculate posterior distribution
  #
  # Args: 
  #   Obs: a numeric vector
  #   Ini: q_start to q_1 prob distribution
  #   A: Trans matrix
  #   E: Emission matrix, e_q\sigma
  #
  # Output:
  #   Probability matrix of posterior distribution at each index(time)
  #
  b <- Backward(Obs, A, E)
  f <- Forward(Obs, Ini, A, E)
  fb <- f * b
  output <- apply(fb, 2, function(x) x / sum(x))
  return(t(output))
}

BW <- function(nQ, nS, Obs, Ini = NULL, epsilon = 1e-09, maxIter = 100, IniA = NULL, IniE = NULL){
  # Take a list of observations, return HMM parameters
  #
  # Args:
  #   nQ: number of states
  #   nS: number of alphabet
  #
  # Output:
  #   Ini: q_start to q_1 distribution   
  #   A: transition matrix
  #   E: emission matrix
  #
  # Initialization
  n <- length(Obs)
  if(is.null(IniA)){
    A <- matrix(abs(runif(nQ * nQ)), nQ, nQ)
    A <- t(apply(A, 1, function(x) x / sum(x)))
  } else {
    A <- IniA
  }
  if(is.null(IniE)){
    E <- matrix(abs(runif(nQ * nS)), nQ, nS)
    E <- t(apply(E, 1, function(x) x / sum(x)))
  } else {
    E <- IniE
  }
  
  if(is.null(Ini)){
    Ini <- abs(runif(nQ)) 
    Ini <- Ini / sum(Ini)
  }
  iter <- 1
  diff <- 1
  # EM algorithm
  while (diff > epsilon && iter <= maxIter){
    alpha <- Forward(Obs, Ini, A, E)
    beta <- Backward(Obs, A, E)
    alphabeta <- alpha * beta
    gamma <- apply(alphabeta, 2, function(x) x / sum(x))
    A.nu <- 0
    for (i in 1:(n - 1)){
      A.nu <- A.nu + alpha[, i] * A * (rep(1, nQ) %*% t(E[, Obs[i+1]])) * (rep(1, nQ) %*% t(beta[, i+1]))
    }
    newA <- t(apply(A.nu, 1, function(x) x / sum(x))) 
    E.nu <- NULL
    for (i in 1:nS){
      E.nu <- cbind(E.nu, rowSums(alphabeta[, which(Obs == i)]))
    }
    newE <- t(apply(E.nu, 1, function(x) x / sum(x)))
    diff <- sqrt(sum((newA - A) ^ 2)) + sqrt(sum((newE - E) ^ 2))
    A <- newA
    E <- newE
    Ini <- gamma[, 1]
    iter <- iter + 1
  }
  if(iter == maxIter) print("Maximum iteration reached, algorithm didn't converge.")
  output <- list(Ini = Ini, A = A, E = E, Difference = diff)
  return(output)
}

library(HMM)

Obstest <- c("1", "2", "1", "1", "2", "3", "2", "3", "3", "3")
IniA <- matrix(c(.75, .25, .25, .75), 2, 2)
IniE <- matrix(rep(1/3, 6), 2, 3)
Ini <- rep(1/2, 2)
Obs <- as.numeric(Obstest)
hmm <- initHMM(c("1", "2"), c("1", "2", "3"))

library(microbenchmark)
microbenchmark(BW(2, 3, Obs, Ini, IniA = IniA, IniE = IniE), baumWelch(hmm, Obstest))
