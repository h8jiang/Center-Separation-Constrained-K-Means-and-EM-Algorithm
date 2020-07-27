####################################################################
#### Center Separation Constrained EM Algorithm in Dimension 1  ####
####################################################################

#' cons.em.sep
#' Constrained EM Algorithm with Arbitrary Center Separation in Gaussian Mixture Models
#' The function fits a Gaussian Mixture Model (GMM) on the input data with the desired number of components, given certain separation between adjacent component mean values

#' Inputs:
#' x: the vector containing the data used to fit the GMM, vector of length N
#' K: the desired number of components of the GMM, number
#' init_mu, init_pi, init_sigsq: if provided, the initialization of values of the mean, weight, variance parameters used in the EM algorithm, each vector of length $K$; notice given_init needs to be set to TRUE when these values are provided
#' delta: if given as a number, the lower bound of separation between adjacent center means; if given as a vector of length $K-1$, entry $k$ corresponds to the lower bound of separation of mean values between component $k$ and $k+1$, $k=1,...,K-1$
#' delta_2 if given as a number, the upper bound of separation between adjacent center means; if given as a vector of length $K-1$, entry $k$ corresponds to the upper bound of separation of mean values between component $k$ and $k+1$, $k=1,...,K-1$; note that if we are only interested in lower bound separation, this value could be set to a really large number, for example $max(x) - min(x) + 100$ in the case of same lower bound separation, or to a vector of length $K-1$ of large numbers, for example $rep(max(x) - min(x) + 100, K-1)$
#' tol: tolerance value, where if the change of each of mu, pi, sigsq values are all less than this value, we stop the EM iterations, number 
#' given_init: whether the initial parameter values of the EM iteration are given, boolean; if TRUE init_mu, init_pi, init_sigsq has to be provided; if FALSE the initialization would be carried out by Ckmeans.1d.dp
#' prog, whether to print the result of each step, boolean
#' Returns:
#' A list of containing pi, mu, sig_sq, s, l_l


# required packages
require(Ckmeans.1d.dp) 
require(quadprog) 

cons.em.sep = function(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, given_init, prog){
  
  # 1 initializations
  if (given_init == TRUE){
  }
  else{
    init_obj = Ckmeans.1d.dp(x, K)
    init_mu_uo = as.numeric(init_obj$centers)
    init_pi_uo = as.numeric(init_obj$size / N)
    init_sigsq_uo = init_obj$withinss / init_obj$size
    ord = order(init_mu_uo)
    r_m = cbind(ord, init_pi_uo, init_mu_uo, init_sigsq_uo)
    sort_m =r_m[order(r_m[,1]), ]
    init_pi = as.numeric(sort_m[,2])
    init_mu = as.numeric(sort_m[,3])
    init_sigsq = as.numeric(sort_m[,4])
  }
  mu = init_mu
  sig_sq = init_sigsq
  pi = init_pi
  N = length(x)
  s = 1 
  sat = FALSE 
  w_mat = matrix(NA, nrow = N, ncol = K) 
  dens_mat = matrix(NA, nrow = N, ncol = K) 
  
  # 2 EM
  while (sat == FALSE){
    l_l = 0 
    
    ## 2.1 E Step: 
    for (k in 1:K){
      dens_mat[ , k] = dnorm(x , mu[k] , sqrt(sig_sq[k]) , log=FALSE) 
    }  
    for (i in 1:N){
      for (k in 1:K){
        w_mat[i,k] = pi[k] * dens_mat[i,k]
      }
      denom = sum(w_mat[i,])        
      w_mat[i,] = w_mat[i,] / denom
    } 
    
    ## 2.2 M Step: 
    N_k = numeric(K) 
    mu_new = numeric(K)
    sig_sq_new = numeric(K)
    for (k in 1:K){
      N_k[k] = sum(w_mat[,k])
    }
    pi_new = N_k / N 
    
    ### update mu: constraint optimization
    GIK = matrix(0, nrow = N, ncol = K)
    G = matrix(0, nrow = K, ncol = K)
    d_trans = numeric(K) 
    vec_1 = rep(1, K)
    for (i in 1:N){
      for (k in 1:K){
        GIK[i, k] = w_mat[i, k] / (2 * sig_sq[k])
      }
      G_I = diag(GIK[i, ])
      d_trans = d_trans + 2 * (x[i] * t(vec_1) %*% G_I)
      G = G + G_I
    }
    A = matrix(0, nrow=2*K-2, ncol=K)
    for (k in 1:(K-1)){
      A[k,k] =  - 1
      A[k, k+1] = 1
      A[k+K-1,k] =  1
      A[k+K-1, k+1] = - 1
    }
    if (length(delta) == 1){ # delta and delta2 are numbers in this case
      bvec = c(rep(delta, K-1), rep(-delta2, K-1))
    }
    else{ # delta and delta2 are vectors in this case
      bvec = length(2*K - 2)
      for (k in 1:(K-1)){
        bvec[k] = delta[k]
      }
      for (k in K:(2*K-2)){
        bvec[k] = -delta2[k - K + 1]
      }
    }
    Amat = t(A)
    Dmat = 2 * G 
    dvec = t(d_trans)
    cons_op = solve.QP(Dmat, dvec, Amat, bvec, meq=0)
    mu_new = cons_op$solution
    
    for (k in 1:K){
      diff_sq = (x - mu_new[k])^2
      sig_sq_new[k] = as.numeric(diff_sq %*% w_mat[,k])
    }
    sig_sq_new = sig_sq_new / N_k 
    
    # 3 check for convergence and update parameters
    sat = all(abs(pi_new - pi) <= tol) && all(abs(mu_new - mu) <= tol) && all(abs(sig_sq_new - sig_sq) <= tol)
    if ( any(pi_new < 1e-16) ){
      sat = TRUE
    }
    pi = pi_new
    mu = mu_new
    sig_sq = sig_sq_new
    for (id in 1:N){
      dens_id = 0
      for (k in 1:K){
        dens_id = dens_id + dens_mat[id, k] * pi[k]
      }
      l_l = l_l + log(dens_id)
    }
    if (prog == TRUE){
      print(paste0("iteration: ", s))
      print(paste0("weights: ", pi, " means: ", mu, " sigsqs: ", sig_sq))
      print(paste0("log likelihood: ", l_l))
    }
    s = s + 1
    if (sat == TRUE){
      return(list(pi, mu, sig_sq, s, l_l))
    }
  }
}

