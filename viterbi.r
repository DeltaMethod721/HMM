Viterbi <- function(nQ, nS, A, E, Ini, Obs){
  # Viterbi's Algorithm: given HMM model & observed sequence, infer state path
  #
  # Args: 
  #   mQ: number of states
  #   nS: number of different observed values
  #   A: transition matrix
  #   E: Emission matrix
  #   Ini: initial state prob
  #   Obs: observed sequence, an integer vector
  #
  # Output: A list of path(integer vector) and maximum likelihood(numerical)
  
  # Initialize vq
  vq <- Ini * E[, Obs[1]]
  nl <- length(Obs)
  # Update vq and save path
  i <- 2
  pathlist <- matrix(0, nQ, nl)
  pathlist[, 1] <- 1:nQ
  while (i <= nl){
    mat <- vq * A * (rep(1, nQ) %*% t(E[, Obs[i]]))
    vq <- apply(mat, 2, max)
    nextState <- apply(mat, 2, which.max)
    pathlist[, i] <- nextState
    i <- i + 1
  }
  ml <- max(vq)
  endState <- which.max(vq)
  path <- rep(0, nl)
  path[nl] <- endState
  for (i in (nl - 1):1){
    path[i] <- pathlist[path[i + 1], i + 1]
  }
  output <- list(Path = path, Likelihood = ml)
  return(output)
}



E <- matrix(c(.7, .9, .6, .2, .05, .3, .1, .05, .1), 3, 3)
A <- matrix(c(.1, .9, .7, .5, .05, .2, .4, .05, .1), 3, 3)
Obs <- c(3, 1, 2, 3, 1)
Ini <- c(.2, .4, .4)
Viterbi(3, 3, A, E, Ini, Obs)
#"3" "1" "3" "1" "2"

#Comparison between "HMM" and this function. Viterbi is 30% faster.
library(HMM)
hmm <- initHMM(c("1", "2", "3"), c("a", "b", "c"), Ini, A, E)
Obsletter <- c("c", "a", "b", "c", "a")

library(microbenchmark)
microbenchmark(viterbi(hmm, Obsletter), Viterbi(3, 3, A, E, Ini, Obs))
