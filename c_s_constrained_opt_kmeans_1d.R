#######################################################################
#### Optimal Center Separation Constrained K Means in Dimension 1  ####
#######################################################################


#' # cons.kmeans.1d.main
#' Main function, with relevant quantity computations
#' The function computes relevant matrices and quantities of the optimal constrained K Means algorithm in 1D
#' Inputs:
#' x: a length N sorted vector
#' N: length of x
#' K: desired number of clusters
#' delta: minimum separation between cluster centers
#' Returns:
#' A list of containing D, I, B, C, W, beta, U

cons.kmeans.1d.main = function(x, N, K, delta){
  
  ## 1 Preliminary and Create Necassary Matrices 
  W = matrix(0, nrow = N, ncol = N)
  C = matrix(0, nrow = N, ncol = N)
  for (j in 1:(N-1)){
    x_comp = x[j:N]
    d = numeric(N-j+1)
    mu = numeric(N-j+1)
    d[1] = 0
    mu[1] = x_comp[1] 
    for (i in 2:(N-j+1)){
      d[i] = d[i-1] + ((i-1)/i) * (x_comp[i] - mu[i-1])^2 
      mu[i] = (x_comp[i] + (i-1)*mu[i-1]) / i
    }
    W[j, j:N] = d
    C[j, j:N] = mu
  }
  C[N, N] = x[N]
  beta = W[1, N] + 100 
  D = matrix(beta, nrow = N, ncol = K)
  I = matrix(0, nrow = N, ncol = K) 
  B = matrix(0, nrow = N, ncol = K)
  U = array(beta, c(N,N,K))
  D[1,1] = 0
  for (i in 2:N){
    D[i,1] = (i-1)*var(x[1:i]) # within cluster some of squres
  }
  for (i in 1:N){
    I[i,1] = 1
    B[i,1] = 1
  }
  U[,,1] = matrix(0, N, N)
  
  ## 2 Update D, I, B, U
  
  # 2.1 K=1 case
  if (K==1){
    # end the D,I,B function
    return(list(D, I, B, C, W, beta, U))
  }
  
  # 2.2 K=2 case
  for (i in 2:N){ 
    p = rep(beta, i - 1)    
    for (j in 2:i){ 
      if ( (C[j,i] - C[1,(j-1)] ) >= delta ){ 
        p[j - 1] = W[j,i] + W[1,(j-1)]
        U[j,i,2] = W[1,(j-1)]
      }
      else{next}
    } 
    D[i, 2] = min(p) 
    if (D[i,2] == beta){
      I[i, 2] = 0
      B[i, 2] = 0
    }
    else{
      I[i, 2] = min( which(p == min(p)) ) + 1 
      B[i, 2] = min( which(p < beta) ) + 1 
    }
  } 
  if (K==2){
    return(list(D, I, B, C, W, beta, U))
  }
  
  ## 2.3 K>=3 general case
  for (m in 3:K){ 
    for (i in m:N){ 
      if (D[i-1, m-1] < beta){ 
        p = rep(beta, i - m + 1)    
        for (j in m:i){ 
          if( B[j-1, m-1] != 0 ){ 
            ### 2.3.1 when there is only one situation starting point of t to consider
            if ( C[j,i] - C[ I[j-1,m-1] , j-1 ] >= delta){ 
              p[j - m + 1] = D[j-1, m-1] + W[j,i]
              U[j,i,m] = D[j-1, m-1] 
            }
            ### 2.3.2 when there are multiple situations for starting point of t to consider 
            else{ 
              for (t in I[j-1,m-1]:B[j-1,m-1]){ 
                if ( (C[j,i] - C[t,j-1] >= delta)  ){ 
                  p[j - m + 1] = min(p[j - m + 1], U[t,j-1,m-1] + W[t,j-1] + W[j,i])
                  U[j,i,m] = U[t,j-1,m-1] + W[t, j-1]
                  break
                }
              }
            }
          } 
        } 
        D[i, m] = min(p) 
        if (D[i,m] == beta){
          I[i, m] = 0
          B[i, m] = 0
        }
        else{
          I[i, m] = min( which(p == min(p)) ) + m - 1 
          B[i, m] = min( which(p < beta) ) + m - 1 
        }
      } 
    } 
  } 
  return(list(D, I, B, C, W, beta, U))
} 




#' # cons.kmeans.1d.si
#' Cluster Starting Index Computation
#' The function computes the optimal starting index of each cluster in the optimal constrained K Means algorithm in 1D
#' Inputs:
#' D: D from cons.kmeans.1d.main
#' I: I from cons.kmeans.1d.main
#' B: B from cons.kmeans.1d.main
#' C: C from cons.kmeans.1d.main
#' W: W from cons.kmeans.1d.main
#' beta: beta from cons.kmeans.1d.main
#' delta: minimum separation between cluster centers
#' N: length of x
#' K: desired number of clusters
#' Returns:
#' A vector b of the starting indices of each cluster 


cons.kmeans.1d.si = function(D, I, B, C, W, beta, delta, N, K){
  b = numeric(K) 
  b[K] = I[N,K]
  l = b[K] - 1 
  if (K==1){
    return(b)
  }
  for (j in I[l,K-1]:B[l,K-1]){ 
    if ( C[ l+1, N ] - C[j,l] >= delta){
      b[K-1] = j
      l = j-1
      break
    }
  }
  if (K==2){
    return(b)
  }
  for (k in (K-2):1){
    for (j in I[l,k]:B[l,k]){ 
      if ( (C[ l+1, b[k+2]-1 ] - C[j,l]) >= delta ){
        b[k] = j
        l = j-1
        break
      }
    }
  }
  return(b)
}



#' # Result reporting function
#' The function computes results of the optimal constrained K Means algorithm in 1D
#' Inputs:
#' x_input: the data vector used for clustering
#' K: desired number of clusters
#' delta: minimum separation between cluster centers
#' Returns: 
#' A list of containing cluster.labels, s.d.start.ind, cluster.start.value, centers, withinss, size, tot.withinss, betweenss, totss 


cons.kmeans.1d.res = function(x_input, K, delta){
  
  x_o = order(x_input)
  x = sort(x_input)
  N = length(x)
  cent = numeric(K) 
  ws = numeric(K) 
  size = numeric(K) 
  bs = 0 
  g_cent = mean(x) 
  l_clus = numeric(N) 
  
  d_comp = cons.kmeans.1d.main(x, N, K, delta)
  D = d_comp[[1]]
  I = d_comp[[2]]
  B = d_comp[[3]]
  C = d_comp[[4]]
  W = d_comp[[5]]
  beta = d_comp[[6]]
  if (D[N,K] == beta){
    print("Putting x_input into K clusters with center separation delta is not possible")
    return()
  }
  
  b = cons.kmeans.1d.si(D, I, B, C, W, beta, delta, N, K)
  for (k in 1:(K-1)){
    cent[k] = C[b[k], b[k+1]-1]
    ws[k] = W[b[k], b[k+1]-1]
    size[k] = b[k+1] - b[k]
  }
  cent[K] = C[b[K], N]
  ws[K] = W[b[K], N]
  size[K] = N - b[K] + 1
  for (k in 1:K){
    bs = bs + size[k] * (cent[k] - g_cent)^2
  }
  b_e = c(b,N+1)
  
  for (i in 1:N){
    l_o_i = which(x_o == i) 
    for (k in 1:K){
      if (l_o_i < b_e[k+1] && l_o_i >= b_e[k]){
        l_clus[i] = k
      }
    }
  }
  
  return( list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) )
}

