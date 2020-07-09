#### K Means and GMM with a Separation Constraint Experiments ####
#### He Jiang, Ery Arias-Castro ####

## Note: please load the attached datasets used by the authors for exact results

#### Require packages
require(MASS)
require(Ckmeans.1d.dp) # optimal k means Wang 2011
require(quadprog) # primal dual interior point method
require(haven) # load sps files
require(fossil) # rand.index function
require(ggplot2) # for producing plots with confidence bands


#### 4.2 NYPANS High School Height Dataset 

# Data Processing
fpath = paste0(getwd(), "/nypans2010.sav")
dataset = read_sav(fpath)


# Finding most appropriate mean age of each group
dat_inco_age = cbind(dataset$e_q1, dataset$e_htwt_height, dataset$e_q3) # age, height, grade
dat = dat_inco_age[complete.cases(dat_inco_age),]
any(is.na(dat)) 
grade = as.numeric(dat[,3])
na_grade_row = which(is.na(grade))
wrong_row = which( (as.numeric(grade)) == 5)
dat = dat[-c(na_grade_row,wrong_row),]
dim(dat) # 9980*3 
g9 = dat[dat[,3]==1,]
g9_age = mean(as.numeric(g9[,1])) + 11
g10 = dat[dat[,3]==2,]
g10_age = mean(as.numeric(g10[,1])) + 11
g11 = dat[dat[,3]==3,]
g11_age = mean(as.numeric(g11[,1])) + 11
g12 = dat[dat[,3]==4,]
g12_age = mean(as.numeric(g12[,1])) + 11
g9_age # [1] 14.70245 
g10_age # [1] 15.66468
g11_age # [1] 16.64928
g12_age # [1] 17.54117
delta=c(0.019023, 0.00895, 0.00154) # read from 2001 tables
delta2=c(0.029023, 0.01895, 0.01154)


# Getting dataframe with only grade and height
dat_inco = cbind(dataset$e_htwt_height, dataset$e_q3) # height and grade(9-12)
dat = dat_inco[complete.cases(dat_inco),]
grade = as.numeric(dat[,2]) 
na_grade_row = which(is.na(grade))
wrong_row = which( (as.numeric(grade)) == 5)
dat = dat[-c(na_grade_row, wrong_row), ]
any(is.na(dat)) # no NA terms
dim(dat) # [1] 9980    2
K = 4 
tol = 10^(-3)
given_init = TRUE
prog = FALSE
grade = as.numeric(dat[,2])
x = as.numeric(dat[,1])
dat = cbind(x, grade) 
lab_h = as.numeric(dat[,2]) 


# Actual Values
sum(grade==1)/length(grade) # [1] 0.2613226
sum(grade==2)/length(grade) # [1] 0.252505
sum(grade==3)/length(grade) # [1] 0.243988
sum(grade==4)/length(grade) # [1] 0.2421844
mean(dat[grade==1,1]) # 1.656931
mean(dat[grade==2,1]) # 1.674162
mean(dat[grade==3,1]) # 1.686936
mean(dat[grade==4,1]) # 1.694234
var(dat[grade==1,1]) # 0.008134112
var(dat[grade==2,1]) # 0.009035427
var(dat[grade==3,1]) # 0.01003764
var(dat[grade==4,1]) # 0.01020333


# Regular EM Function
em.reg = function(x, K, init_mu, init_sigsq, init_pi, tol){
  
  # initialization
  N = length(x)
  mu = init_mu
  sig_sq = init_sigsq
  pi = init_pi
  s = 1 
  sat = FALSE 
  w_mat = matrix(NA, nrow = N, ncol = K) 
  dens_mat = matrix(NA, nrow = N, ncol = K) 
  
  # begin EM
  while (sat == FALSE){
    
    # E Step: update density matrix in order to update weight matrix
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
    
    # M Step: 
    N_k = numeric(K) 
    mu_new = numeric(K)
    sig_sq_new = numeric(K)
    for (k in 1:K){
      N_k[k] = sum(w_mat[,k])
    }
    pi_new = N_k / N # update pi
    for (k in 1:K){
      mu_new[k] = as.numeric(x %*% w_mat[,k])
    }
    mu_new = mu_new / N_k # update mu
    for (k in 1:K){
      diff_sq = (x - mu_new[k])^2
      sig_sq_new[k] = as.numeric(diff_sq %*% w_mat[,k])
    }
    sig_sq_new = sig_sq_new / N_k # update sig_sq
    
    # check whether the tolerence is satisfied or not; or if one weight is super small, we stop
    sat = all(abs(pi_new - pi) <= tol) && all(abs(mu_new - mu) <= tol) && all(abs(sig_sq_new - sig_sq) <= tol)
    if ( any(pi_new < 1e-16) ){
      sat = TRUE
    }
    
    # update the parameters
    pi = pi_new
    mu = mu_new
    sig_sq = sig_sq_new
    l_l = 0
    for (id in 1:N){
      dens_id = 0
      for (k in 1:K){
        dens_id = dens_id + dens_mat[id, k] * pi[k]
      }
      l_l = l_l + log(dens_id)
    }
    s = s + 1
    if (sat == TRUE){
      return(list(pi, mu, sig_sq, s, l_l))
    }
  }
}


# EM Initialization
N = length(x)
init_obj = Ckmeans.1d.dp(x, K)
init_mu_uo = as.numeric(init_obj$centers)
init_pi_uo = as.numeric(init_obj$size / N)
init_sigsq_uo = init_obj$withinss / init_obj$size
ord = order(init_mu_uo)
r_m = cbind(ord, init_pi_uo, init_mu_uo, init_sigsq_uo)
sort_m =r_m[order(r_m[,1]), ]
init_mu = as.numeric(sort_m[,3])
init_pi = as.numeric(sort_m[,2])
init_sigsq = as.numeric(sort_m[,4])


# Run the EM Algorithms with labels
com_samp = cbind(lab_h, x)
lab_mat = com_samp[order(com_samp[,2]),] 
sort_comp = lab_mat[,1]
x = lab_mat[,2]
res_reg_em = em.reg(x, K, init_mu, init_sigsq, init_pi, tol) # Regular EM
res_cons_em = cons.em.sep(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, given_init, prog) # Constrained EM

pi_reg_em = res_reg_em[[1]] # regular em weight values
mu_reg_em = res_reg_em[[2]] # regular em mean values
sigsq_reg_em = res_reg_em[[3]] # regular em variance values
pi_cons_em = res_cons_em[[1]] # constrained em weight values
mu_cons_em = res_cons_em[[2]] # constrained em mean values
sigsq_cons_em = res_cons_em[[3]] # constrained em variance values

lab_reg_em = numeric(N)
lab_cons_em = numeric(N)
for (i in 1:N){
  l_cd_reg = numeric(K) 
  l_cd_cons = numeric(K)
  for (k in 1:K){
    l_cd_reg[k] = pi_reg_em[k] * dnorm(x[i], mu_reg_em[k], sqrt(sigsq_reg_em[k]) )
    l_cd_cons[k] = pi_cons_em[k] * dnorm(x[i], mu_cons_em[k], sqrt(sigsq_cons_em[k]) )
  }
  lab_reg_em[i] = which.max(l_cd_reg)
  lab_cons_em[i] = which.max(l_cd_cons)
}

# rand index comparisons
rand.index(sort_comp, lab_reg_em)
rand.index(sort_comp, lab_cons_em)

# mean absolute value comparisons
mus = c(1.656931, 1.674162, 1.686936, 1.694234)
sigsqs = c(0.008134112, 0.009035427, 0.01003764, 0.01020333)
pis = c(0.2613226,0.252505,0.243988,0.2421844)
mean( abs(mu_reg_em - mus) )
mean( abs(mu_cons_em - mus) )

# all parameters absolute value comparisons
mean( abs(mu_reg_em - mus) ) + mean( abs(pi_reg_em - pis) ) + mean( abs(sigsq_reg_em - (sigsqs) ) ) 
mean( abs(mu_cons_em - mus) ) + mean( abs(pi_cons_em - pis) ) + mean( abs(sigsq_cons_em - (sigsqs) ) ) 

######################## ends 4.2 experiment ##########################









#### 4.1.12 K Means Experiments

## Experiment 1 
N = 500
K = 5
mus = c(0, 2, 4, 6, 8)
sds = c(0.25, 0.75, 1.25, 0.75, 0.25)
delta = 1.95
num_exp = 1000

exp1_sam_mat = matrix(0, nrow = num_exp, ncol = N) 
exp1_lab_mat = matrix(0, nrow = num_exp, ncol = N) 
center_d_mat = matrix(0, nrow = num_exp, ncol = 2) 
weight_d_mat = matrix(0, nrow = num_exp, ncol = 2) 
rand_ind_mat = matrix(0, nrow = num_exp, ncol = 2) 

for (i_exp in 1:num_exp){
  # generate the data
  components = sample(1:K,prob=c(0.1,0.2,0.4,0.2,0.1),size=N,replace=TRUE)
  samples = rnorm(n=N,mean=mus[components],sd=sds[components])
  
  # record the data and labels
  exp1_sam_mat[i_exp,] = samples
  exp1_lab_mat[i_exp,] = components
  
  # computations
  cor_size = numeric(K)
  cor_center = numeric(K)
  for (k in 1:K){ 
    cor_size[k] = sum(components == k) 
    cor_center[k] = mean(samples[components == k])
  }
  
  com_samp = cbind(components, samples)
  lab_mat = com_samp[order(com_samp[,2]),]  
  sort_comp = lab_mat[,1]
  x = lab_mat[,2]  # x is sorted here
  
  res = cons.kmeans.1d.res(x, K, delta)
  # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1])  

  #res[[1]] # labels
  #res[[4]] # centers
  #res[[6]] # sizes
  
  # Wang 2011 method
  reg_obj = Ckmeans.1d.dp(x, K)
  
  # center absolute differences
  center_d_mat[i_exp, 1] = sum(abs(res[[4]] - cor_center))
  center_d_mat[i_exp, 2] = sum(abs(reg_obj$centers - cor_center))
  # size absolute differences
  weight_d_mat[i_exp, 1] = sum(abs(res[[6]] - cor_size))
  weight_d_mat[i_exp, 2] = sum(abs(reg_obj$size - cor_size))
  # rand index
  rand_ind_mat[i_exp,1] = rand.index(sort_comp, res[[1]])
  rand_ind_mat[i_exp,2] = rand.index(sort_comp, reg_obj$cluster)
  
  # see iteration
  print(paste0('iteration ', i_exp, ' is complete'))
}

mean(center_d_mat[,1])
mean(center_d_mat[,2])
sd(center_d_mat[,1])
sd(center_d_mat[,2])

mean(weight_d_mat[,1])
mean(weight_d_mat[,2])
sd(weight_d_mat[,1])
sd(weight_d_mat[,2])

mean(rand_ind_mat[,1])
mean(rand_ind_mat[,2])
sd(rand_ind_mat[,1])
sd(rand_ind_mat[,2])



## Experiment 2 
N = 500
K = 3
mus = c(0, 2, 4)
sds = c(0.75, 1.5, 0.75)
delta = 1.95
num_exp = 1000

exp2_sam_mat = matrix(0, nrow = num_exp, ncol = N) 
exp2_lab_mat = matrix(0, nrow = num_exp, ncol = N) 
center_d_mat_exp2 = matrix(0, nrow = num_exp, ncol = 2) 
weight_d_mat_exp2 = matrix(0, nrow = num_exp, ncol = 2) 
rand_ind_mat_exp2 = matrix(0, nrow = num_exp, ncol = 2) 

for (i_exp in 1:num_exp){
  components = sample(1:K,prob=c(0.45,0.1,0.45),size=N,replace=TRUE)
  samples = rnorm(n=N,mean=mus[components],sd=sds[components])
  
  # record the data and labels
  exp2_sam_mat[i_exp,] = samples
  exp2_lab_mat[i_exp,] = components
  
  # computations
  cor_size = numeric(K)
  cor_center = numeric(K)
  for (k in 1:K){ 
    cor_size[k] = sum(components == k) 
    cor_center[k] = mean(samples[components == k])
  }
  
  com_samp = cbind(components, samples)
  lab_mat = com_samp[order(com_samp[,2]),]  
  sort_comp = lab_mat[,1]
  x = lab_mat[,2]  # x is sorted here
  
  res = cons.kmeans.1d.res(x, K, delta)
  # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1])  
  #res[[1]] # labels
  #res[[4]] # centers
  #res[[6]] # sizes
  
  # Wang 2011 method
  reg_obj = Ckmeans.1d.dp(x, K)

  # center absolute differences
  center_d_mat_exp2[i_exp, 1] = sum(abs(res[[4]] - cor_center))
  center_d_mat_exp2[i_exp, 2] = sum(abs(reg_obj$centers - cor_center))
  
  # size absolute differences
  weight_d_mat_exp2[i_exp, 1] = sum(abs(res[[6]] - cor_size))
  weight_d_mat_exp2[i_exp, 2] = sum(abs(reg_obj$size - cor_size))
  
  # rand index
  rand_ind_mat_exp2[i_exp,1] = rand.index(sort_comp, res[[1]])
  rand_ind_mat_exp2[i_exp,2] = rand.index(sort_comp, reg_obj$cluster)
  
  # see iteration
  print(paste0('iteration ', i_exp, ' is complete'))
}


mean(center_d_mat_exp2[, 1])
sd(center_d_mat_exp2[, 1])
mean(center_d_mat_exp2[, 2])
sd(center_d_mat_exp2[, 2])

mean(weight_d_mat_exp2[, 1])
sd(weight_d_mat_exp2[, 1])
mean(weight_d_mat_exp2[, 2])
sd(weight_d_mat_exp2[, 2])

mean(rand_ind_mat_exp2[, 1])
sd(rand_ind_mat_exp2[, 1])
mean(rand_ind_mat_exp2[, 2])
sd(rand_ind_mat_exp2[, 2])





## 4.1.345 Constrained EM Experiments

## Experiment 3 
N = 500
K = 3
pis = c(0.45,0.1,0.45)
mus = c(0, 2, 4)
sds = c(0.75, 1.5, 0.75)
delta = 1.9
delta2 = 2.1
tol = 10^(-3)
num_exp = 1000


exp3_sam_mat = matrix(0, nrow = num_exp, ncol = N) 
exp3_lab_mat = matrix(0, nrow = num_exp, ncol = N) 
center_d_mat_exp3 = matrix(0, nrow = num_exp, ncol = 2) 
all_d_mat_exp3 = matrix(0, nrow = num_exp, ncol = 2) 
rand_ind_mat_exp3 = matrix(0, nrow = num_exp, ncol = 2) 

for (i_exp in 1:num_exp){
  # generate the data in exp3
  components = sample(1:K,prob=pis,size=N,replace=TRUE)
  samples = rnorm(n=N,mean=mus[components],sd=sds[components])
  exp3_sam_mat[i_exp,] = samples
  exp3_lab_mat[i_exp,] = components
  
  com_samp = cbind(components, samples)
  lab_mat = com_samp[order(com_samp[,2]),]  
  sort_comp = lab_mat[,1]
  x = lab_mat[,2]  # x is sorted here
  
  # constrained k means initialization
  res_cons_km = cons.kmeans.1d.res(x, K, delta)
  # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) 
  #      [[4]] # centers
  #      [[5]] # within sum of squares vector
  #      [[6]] # size
  
  init_mu = res_cons_km[[4]]
  init_sigsq = res_cons_km[[5]]/res_cons_km[[6]]
  init_pi = res_cons_km[[6]]/N
  
  # EM Algorithms
  res_reg_em = em.reg(x, K, init_mu, init_sigsq, init_pi, tol) 
  res_cons_em = cons.em.sep(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, TRUE, FALSE) # use same initialization
  
  pi_reg_em = res_reg_em[[1]] # regular em weight values
  mu_reg_em = res_reg_em[[2]] # regular em mean values
  sigsq_reg_em = res_reg_em[[3]] # regular em variance values
  pi_cons_em = res_cons_em[[1]] # constrained em weight values
  mu_cons_em = res_cons_em[[2]] # constrained em mean values
  sigsq_cons_em = res_cons_em[[3]] # constrained em variance values
  
  lab_reg_em = numeric(N)
  lab_cons_em = numeric(N)
  for (i in 1:N){
    l_cd_reg = numeric(K) 
    l_cd_cons = numeric(K)
    for (k in 1:K){
      l_cd_reg[k] = pi_reg_em[k] * dnorm(x[i], mu_reg_em[k], sqrt(sigsq_reg_em[k]) )
      l_cd_cons[k] = pi_cons_em[k] * dnorm(x[i], mu_cons_em[k], sqrt(sigsq_cons_em[k]) )
    }
    lab_reg_em[i] = which.max(l_cd_reg)
    lab_cons_em[i] = which.max(l_cd_cons)
  }

  # rand index comparisons
  rand_ind_mat_exp3[i_exp, 1] = rand.index(sort_comp, lab_reg_em)
  rand_ind_mat_exp3[i_exp, 2] = rand.index(sort_comp, lab_cons_em)
  # mean absolute value comparisons
  center_d_mat_exp3[i_exp, 1] = mean( abs(mu_reg_em - mus) )
  center_d_mat_exp3[i_exp, 2] = mean( abs(mu_cons_em - mus) )
  # all parameters absolute value comparisons
  all_d_mat_exp3[i_exp, 1] = mean( abs(mu_reg_em - mus) ) + mean( abs(pi_reg_em - pis) ) + mean( abs(sigsq_reg_em - (sds^2) ) ) 
  all_d_mat_exp3[i_exp, 2] = mean( abs(mu_cons_em - mus) ) + mean( abs(pi_cons_em - pis) ) + mean( abs(sigsq_cons_em - (sds^2) ) ) 
  
  # see iteration
  print(paste0('iteration ', i_exp, ' is complete'))
}


mean(center_d_mat_exp3[, 1])
sd(center_d_mat_exp3[, 1])
mean(center_d_mat_exp3[, 2])
sd(center_d_mat_exp3[, 2])

mean(all_d_mat_exp3[, 1])
sd(all_d_mat_exp3[, 1])
mean(all_d_mat_exp3[, 2])
sd(all_d_mat_exp3[, 2])

mean(rand_ind_mat_exp3[, 1])
sd(rand_ind_mat_exp3[, 1])
mean(rand_ind_mat_exp3[, 2])
sd(rand_ind_mat_exp3[, 2])




## Experiment 4
N = 500
K = 2
pis = c(1/3, 2/3)
mus = c(0, 2)
sds = c(1, 1)
delta = 1.9
delta2 = 2.1
tol = 10^(-3)

num_exp = 1000

exp4_sam_mat = matrix(0, nrow = num_exp, ncol = N) 
exp4_lab_mat = matrix(0, nrow = num_exp, ncol = N) 
center_d_mat_exp4 = matrix(0, nrow = num_exp, ncol = 2) 
all_d_mat_exp4 = matrix(0, nrow = num_exp, ncol = 2) 
rand_ind_mat_exp4 = matrix(0, nrow = num_exp, ncol = 2) 

for (i_exp in 1:num_exp){
  # generate the data in exp3
  components = sample(1:K,prob=pis,size=N,replace=TRUE)
  samples = rnorm(n=N,mean=mus[components],sd=sds[components])
  exp3_sam_mat[i_exp,] = samples
  exp3_lab_mat[i_exp,] = components
  
  com_samp = cbind(components, samples)
  lab_mat = com_samp[order(com_samp[,2]),]  
  sort_comp = lab_mat[,1]
  x = lab_mat[,2]  # x is sorted here
  
  # constrained k means initialization
  res_cons_km = cons.kmeans.1d.res(x, K, delta)
  # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) 
  #      [[4]] # centers
  #      [[5]] # within sum of squares vector
  #      [[6]] # size
  
  init_mu = res_cons_km[[4]]
  init_sigsq = res_cons_km[[5]]/res_cons_km[[6]]
  init_pi = res_cons_km[[6]]/N
  
  # EM Algorithms
  res_reg_em = em.reg(x, K, init_mu, init_sigsq, init_pi, tol) 
  res_cons_em = cons.em.sep(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, TRUE, FALSE) # use same initialization
  
  pi_reg_em = res_reg_em[[1]] # regular em weight values
  mu_reg_em = res_reg_em[[2]] # regular em mean values
  sigsq_reg_em = res_reg_em[[3]] # regular em variance values
  pi_cons_em = res_cons_em[[1]] # constrained em weight values
  mu_cons_em = res_cons_em[[2]] # constrained em mean values
  sigsq_cons_em = res_cons_em[[3]] # constrained em variance values
  
  lab_reg_em = numeric(N)
  lab_cons_em = numeric(N)
  for (i in 1:N){
    l_cd_reg = numeric(K) 
    l_cd_cons = numeric(K)
    for (k in 1:K){
      l_cd_reg[k] = pi_reg_em[k] * dnorm(x[i], mu_reg_em[k], sqrt(sigsq_reg_em[k]) )
      l_cd_cons[k] = pi_cons_em[k] * dnorm(x[i], mu_cons_em[k], sqrt(sigsq_cons_em[k]) )
    }
    lab_reg_em[i] = which.max(l_cd_reg)
    lab_cons_em[i] = which.max(l_cd_cons)
  }
  
  # rand index comparisons
  rand_ind_mat_exp4[i_exp, 1] = rand.index(sort_comp, lab_reg_em)
  rand_ind_mat_exp4[i_exp, 2] = rand.index(sort_comp, lab_cons_em)
  # mean absolute value comparisons
  center_d_mat_exp4[i_exp, 1] = mean( abs(mu_reg_em - mus) )
  center_d_mat_exp4[i_exp, 2] = mean( abs(mu_cons_em - mus) )
  # all parameters absolute value comparisons
  all_d_mat_exp4[i_exp, 1] = mean( abs(mu_reg_em - mus) ) + mean( abs(pi_reg_em - pis) ) + mean( abs(sigsq_reg_em - (sds^2) ) ) 
  all_d_mat_exp4[i_exp, 2] = mean( abs(mu_cons_em - mus) ) + mean( abs(pi_cons_em - pis) ) + mean( abs(sigsq_cons_em - (sds^2) ) ) 
  
  # see iteration
  print(paste0('iteration ', i_exp, ' is complete'))
}


mean(center_d_mat_exp4[, 1])
sd(center_d_mat_exp4[, 1])
mean(center_d_mat_exp4[, 2])
sd(center_d_mat_exp4[, 2])

mean(all_d_mat_exp4[, 1])
sd(all_d_mat_exp4[, 1])
mean(all_d_mat_exp4[, 2])
sd(all_d_mat_exp4[, 2])

mean(rand_ind_mat_exp4[, 1])
sd(rand_ind_mat_exp4[, 1])
mean(rand_ind_mat_exp4[, 2])
sd(rand_ind_mat_exp4[, 2])



## Experiment 5 
N = 500
K = 5
pis = c(1/5, 1/5, 1/5, 1/5, 1/5)
mus = c(0, 2, 4, 6, 8)
sds = c(1, 1, 1, 1, 1)
delta = 1.9
delta2 = 2.1
tol = 10^(-3)

num_exp = 1000

exp5_sam_mat = matrix(0, nrow = num_exp, ncol = N) 
exp5_lab_mat = matrix(0, nrow = num_exp, ncol = N) 
center_d_mat_exp5 = matrix(0, nrow = num_exp, ncol = 2) 
all_d_mat_exp5 = matrix(0, nrow = num_exp, ncol = 2) 
rand_ind_mat_exp5 = matrix(0, nrow = num_exp, ncol = 2) 

for (i_exp in 1:num_exp){
  # generate the data in exp3
  components = sample(1:K,prob=pis,size=N,replace=TRUE)
  samples = rnorm(n=N,mean=mus[components],sd=sds[components])
  exp3_sam_mat[i_exp,] = samples
  exp3_lab_mat[i_exp,] = components
  
  com_samp = cbind(components, samples)
  lab_mat = com_samp[order(com_samp[,2]),]  
  sort_comp = lab_mat[,1]
  x = lab_mat[,2]  # x is sorted here
  
  # constrained k means initialization
  res_cons_km = cons.kmeans.1d.res(x, K, delta)
  # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) 
  #      [[4]] # centers
  #      [[5]] # within sum of squares vector
  #      [[6]] # size
  
  init_mu = res_cons_km[[4]]
  init_sigsq = res_cons_km[[5]]/res_cons_km[[6]]
  init_pi = res_cons_km[[6]]/N
  
  # EM Algorithms
  res_reg_em = em.reg(x, K, init_mu, init_sigsq, init_pi, tol) 
  res_cons_em = cons.em.sep(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, TRUE, FALSE) # use same initialization
  
  pi_reg_em = res_reg_em[[1]] # regular em weight values
  mu_reg_em = res_reg_em[[2]] # regular em mean values
  sigsq_reg_em = res_reg_em[[3]] # regular em variance values
  pi_cons_em = res_cons_em[[1]] # constrained em weight values
  mu_cons_em = res_cons_em[[2]] # constrained em mean values
  sigsq_cons_em = res_cons_em[[3]] # constrained em variance values
  
  lab_reg_em = numeric(N)
  lab_cons_em = numeric(N)
  for (i in 1:N){
    l_cd_reg = numeric(K) 
    l_cd_cons = numeric(K)
    for (k in 1:K){
      l_cd_reg[k] = pi_reg_em[k] * dnorm(x[i], mu_reg_em[k], sqrt(sigsq_reg_em[k]) )
      l_cd_cons[k] = pi_cons_em[k] * dnorm(x[i], mu_cons_em[k], sqrt(sigsq_cons_em[k]) )
    }
    lab_reg_em[i] = which.max(l_cd_reg)
    lab_cons_em[i] = which.max(l_cd_cons)
  }
  
  # rand index comparisons
  rand_ind_mat_exp5[i_exp, 1] = rand.index(sort_comp, lab_reg_em)
  rand_ind_mat_exp5[i_exp, 2] = rand.index(sort_comp, lab_cons_em)
  # mean absolute value comparisons
  center_d_mat_exp5[i_exp, 1] = mean( abs(mu_reg_em - mus) )
  center_d_mat_exp5[i_exp, 2] = mean( abs(mu_cons_em - mus) )
  # all parameters absolute value comparisons
  all_d_mat_exp5[i_exp, 1] = mean( abs(mu_reg_em - mus) ) + mean( abs(pi_reg_em - pis) ) + mean( abs(sigsq_reg_em - (sds^2) ) ) 
  all_d_mat_exp5[i_exp, 2] = mean( abs(mu_cons_em - mus) ) + mean( abs(pi_cons_em - pis) ) + mean( abs(sigsq_cons_em - (sds^2) ) ) 
  
  # see iteration
  print(paste0('iteration ', i_exp, ' is complete'))
}


mean(center_d_mat_exp5[, 1])
sd(center_d_mat_exp5[, 1])
mean(center_d_mat_exp5[, 2])
sd(center_d_mat_exp5[, 2])

mean(all_d_mat_exp5[, 1])
sd(all_d_mat_exp5[, 1])
mean(all_d_mat_exp5[, 2])
sd(all_d_mat_exp5[, 2])

mean(rand_ind_mat_exp5[, 1])
sd(rand_ind_mat_exp5[, 1])
mean(rand_ind_mat_exp5[, 2])
sd(rand_ind_mat_exp5[, 2])




 
#### 4.1.67 EM Bound Experiments

## Experiment 6 gradual increase of lower bound on ONE-SIDED bound

N = 500
K = 2
pis = c(1/3, 2/3)
mus = c(0, 2)
sds = c(1, 1)
tol = 10^(-3)
pos_delta = c(0.85, 0.9, 1.1, 1.3, 1.5, seq(1.6,2.3,0.05), 2.4)
am_delta = length(pos_delta) 
# center separation, all parameters separation, rand index(with standard deviations and 90% confidence bands)
res_mat_delta = matrix(0, nrow = 16, ncol = am_delta) 


num_exp = 100


# once for regular version
reg_mean_d = numeric(num_exp)
reg_all_d = numeric(num_exp)
reg_rand_d = numeric(num_exp)

for (i_exp in 1:num_exp){
  samples = exp4_sam_mat[i_exp+100,] 
  components = exp4_lab_mat[i_exp+100,] 
  com_samp = cbind(components, samples)
  lab_mat = com_samp[order(com_samp[,2]),]  
  sort_comp = lab_mat[,1]
  x = lab_mat[,2]  # x is sorted here
  
  res_cons_km = cons.kmeans.1d.res(x, K, delta)
  # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) 
  #      [[4]] # centers
  #      [[5]] # within sum of squares vector
  #      [[6]] # size
  
  init_mu = res_cons_km[[4]]
  init_sigsq = res_cons_km[[5]]/res_cons_km[[6]]
  init_pi = res_cons_km[[6]]/N
  
  # Regular EM Algorithms
  res_reg_em = em.reg(x, K, init_mu, init_sigsq, init_pi, tol) # Regular EM
  pi_reg_em = res_reg_em[[1]] # regular em weight values
  mu_reg_em = res_reg_em[[2]] # regular em mean values
  sigsq_reg_em = res_reg_em[[3]] # regular em variance values
  
  lab_reg_em = numeric(N)
  for (i in 1:N){
    l_cd_reg = numeric(K) 
    for (k in 1:K){
      l_cd_reg[k] = pi_reg_em[k] * dnorm(x[i], mu_reg_em[k], sqrt(sigsq_reg_em[k]) )
    }
    lab_reg_em[i] = which.max(l_cd_reg)
  }
  
  # Column 1 is original
  # rand index comparisons
  reg_rand_d[i_exp] = rand.index(sort_comp, lab_reg_em)
  # mean absolute value comparisons
  reg_mean_d[i_exp] = mean( abs(mu_reg_em - mus) )
  # all parameters absolute value comparisons
  reg_all_d[i_exp] = mean( abs(mu_reg_em - mus) ) + mean( abs(pi_reg_em - pis) ) + mean( abs(sigsq_reg_em - (sds^2) ) ) 
  # see iteration
  print(paste0('iteration ', i_exp, ' is complete'))
}



## constrained versions
for (i_del in 1:am_delta){
  delta = pos_delta[i_del]
  # delta2 = delta + 0.2
  delta2 = 10000 # one-sided constraint
  
  # constrained versions
  cons_mean_d = numeric(num_exp)
  cons_all_d = numeric(num_exp)
  cons_mis_d = numeric(num_exp)
  cons_rand_d = numeric(num_exp)

  for (i_exp in 1:num_exp){

   # data has been selected as the 101th to 200th data of the randomly generated samples in Exp 4
   samples = exp4_sam_mat[i_exp+100,] 
   components = exp4_lab_mat[i_exp+100,] 
   
   # computations
   com_samp = cbind(components, samples)
   lab_mat = com_samp[order(com_samp[,2]),]  
   sort_comp = lab_mat[,1]
   x = lab_mat[,2]  # x is sorted here
   
   # constrained k means initialization
   res_cons_km = cons.kmeans.1d.res(x, K, delta)
   # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) 
   #      [[4]] # centers
   #      [[5]] # within sum of squares vector
   #      [[6]] # size
   init_mu = res_cons_km[[4]]
   init_sigsq = res_cons_km[[5]]/res_cons_km[[6]]
   init_pi = res_cons_km[[6]]/N
   
   # EM Algorithms
   res_cons_em = cons.em.sep(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, TRUE, FALSE) 
   pi_cons_em = res_cons_em[[1]] # constrained em weight values
   mu_cons_em = res_cons_em[[2]] # constrained em mean values
   sigsq_cons_em = res_cons_em[[3]] # constrained em variance values
   
   lab_cons_em = numeric(N)
   
   for (i in 1:N){
     l_cd_cons = numeric(K)
     for (k in 1:K){
       l_cd_cons[k] = pi_cons_em[k] * dnorm(x[i], mu_cons_em[k], sqrt(sigsq_cons_em[k]) )
     }
     lab_cons_em[i] = which.max(l_cd_cons)
   }
   
   # Column 1 is original, Column 2 is ours
   # rand index comparisons
   cons_rand_d[i_exp] = rand.index(sort_comp, lab_cons_em)
   # mean absolute value comparisons
   cons_mean_d[i_exp] = mean( abs(mu_cons_em - mus) )
   # all parameters absolute value comparisons
   cons_all_d[i_exp] = mean( abs(mu_cons_em - mus) ) + mean( abs(pi_cons_em - pis) ) + mean( abs(sigsq_cons_em - (sds^2) ) ) 
   
   res_mat_delta[,i_del] = c( mean(cons_mean_d),sd(cons_mean_d),  mean(cons_all_d),sd(cons_all_d),  mean(cons_rand_d),sd(cons_rand_d),  0,0, as.numeric(quantile(cons_mean_d, c(0.05, 0.95))),   as.numeric(quantile(cons_all_d, c(0.05, 0.95))),   as.numeric(quantile(cons_rand_d, c(0.05, 0.95))),0,0  )
   # see iteration
   print(paste0(i_del, '-th delta and iteration ', i_exp, ' is complete'))
  }
}
gradual_lower = res_mat_delta





## Experiment 7 gradual increase of lower bound on TWO-SIDED bound
# The data here, as well as regular versions, carry over from Exp 6
# center separation, all parameters separation, rand index(with standard deviations and 90% confidence bands)
res_mat_delta = matrix(0, nrow = 16, ncol = am_delta) 

## constrained versions
for (i_del in 1:am_delta){
  delta = pos_delta[i_del]
  delta2 = delta + 0.2 # two sided version
  
  # constrained versions
  cons_mean_d = numeric(num_exp)
  cons_all_d = numeric(num_exp)
  cons_mis_d = numeric(num_exp)
  cons_rand_d = numeric(num_exp)
  
  for (i_exp in 1:num_exp){
    # data has been selected as the 101th to 200th data of the randomly generated samples in Exp 4
    samples = exp4_sam_mat[i_exp+100,] 
    components = exp4_lab_mat[i_exp+100,] 
    
    # computations
    com_samp = cbind(components, samples)
    lab_mat = com_samp[order(com_samp[,2]),]  
    sort_comp = lab_mat[,1]
    x = lab_mat[,2]  # x is sorted here
    
    # constrained k means initialization
    res_cons_km = cons.kmeans.1d.res(x, K, delta)
    # list(cluster.labels = l_clus, s.d.start.ind = b, cluster.start.value = x[b], centers = cent, withinss = ws, size = size, tot.withinss = D[N,K], betweenss = bs, totss = D[N,1]) 
    #      [[4]] # centers
    #      [[5]] # within sum of squares vector
    #      [[6]] # size
    init_mu = res_cons_km[[4]]
    init_sigsq = res_cons_km[[5]]/res_cons_km[[6]]
    init_pi = res_cons_km[[6]]/N
    
    # Constrained EM Algorithms
    res_cons_em = cons.em.sep(x, K, init_mu, init_pi, init_sigsq, delta, delta2, tol, TRUE, FALSE) 
    pi_cons_em = res_cons_em[[1]] # constrained em weight values
    mu_cons_em = res_cons_em[[2]] # constrained em mean values
    sigsq_cons_em = res_cons_em[[3]] # constrained em variance values
    
    lab_cons_em = numeric(N)
    
    for (i in 1:N){
      l_cd_cons = numeric(K)
      for (k in 1:K){
        l_cd_cons[k] = pi_cons_em[k] * dnorm(x[i], mu_cons_em[k], sqrt(sigsq_cons_em[k]) )
      }
      lab_cons_em[i] = which.max(l_cd_cons)
    }
    
    # Column 1 is original, Column 2 is ours
    # rand index comparisons
    cons_rand_d[i_exp] = rand.index(sort_comp, lab_cons_em)
    # mean absolute value comparisons
    cons_mean_d[i_exp] = mean( abs(mu_cons_em - mus) )
    # all parameters absolute value comparisons
    cons_all_d[i_exp] = mean( abs(mu_cons_em - mus) ) + mean( abs(pi_cons_em - pis) ) + mean( abs(sigsq_cons_em - (sds^2) ) ) 
    
    res_mat_delta[,i_del] = c( mean(cons_mean_d),sd(cons_mean_d),  mean(cons_all_d),sd(cons_all_d),  mean(cons_rand_d),sd(cons_rand_d),  0,0, as.numeric(quantile(cons_mean_d, c(0.05, 0.95))),   as.numeric(quantile(cons_all_d, c(0.05, 0.95))),   as.numeric(quantile(cons_rand_d, c(0.05, 0.95))), 0,0 )
    # see iteration
    print(paste0(i_del, '-th delta and iteration ', i_exp, ' is complete'))
  }
}
gradual_both = res_mat_delta


## Exp 6,7 Plots
options(warn = -1)

## mean error
center_lower = gradual_lower[1,]
lower_m_err_lower = gradual_lower[9,]
upper_m_err_lower = gradual_lower[10,]

center_both = gradual_both[1,]
lower_m_err_both = gradual_both[9,]
upper_m_err_both = gradual_both[10,]

center_reg = mean(reg_mean_d) # 0.2442039
lower_m_err_reg = quantile(reg_mean_d,0.05)
upper_m_err_reg = quantile(reg_mean_d,0.95)

df = data.frame(delta = pos_delta, mean_error = center_reg, center_lower = center_lower, center_both = center_both, lb_r = lower_m_err_reg, ub_r = upper_m_err_reg,  lb_l = lower_m_err_lower, ub_l = upper_m_err_lower, lb_b = lower_m_err_both, ub_b = upper_m_err_both)
ggplot(df, aes(delta)) + 
  geom_vline(xintercept=2, size=0.75, linetype="dotted") + 
  geom_line(aes(y=mean_error), colour="black", size=1.5,alpha=0.85) + 
  geom_line(aes(y=center_lower), colour="red", size=1.5,alpha=0.85) + 
  geom_line(aes(y=center_both), colour="blue", size=1.5,alpha=0.85) + 
  geom_ribbon(aes(ymin=lb_r, ymax=ub_r), colour ='black',alpha=0.15) + 
  geom_ribbon(aes(ymin=lb_l, ymax=ub_l), colour ='red',alpha=0.15) +
  geom_ribbon(aes(ymin=lb_b, ymax=ub_b), colour ='blue',alpha=0.15) +
  scale_y_continuous(limits = c(0, 1)) + 
  theme(axis.text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()
        )



## all error

all_middle_os = gradual_lower[3,]
all_lower_os = gradual_lower[11,]
all_upper_os = gradual_lower[12,]

all_middle_ts = gradual_both[3,]
all_lower_ts = gradual_both[11,]
all_upper_ts = gradual_both[12,]

all_middle_reg = mean(reg_all_d) 
all_lower_reg = quantile(reg_all_d,0.05)
all_upper_reg = quantile(reg_all_d,0.95)


df = data.frame(delta = pos_delta, all_error = all_middle_reg, 
                all_middle_os = all_middle_os, all_middle_ts = all_middle_ts, 
                all_lower_os = all_lower_os, all_upper_os = all_upper_os,  
                all_lower_ts = all_lower_ts, all_upper_ts = all_upper_ts, 
                all_lower_reg = all_lower_reg, all_upper_reg = all_upper_reg, row.names = NULL)
df = data.frame(delta = pos_delta, all_error = all_middle_reg)
df$all_middle_os = all_middle_os
df$all_middle_ts = all_middle_ts
df$all_lower_os = all_lower_os
df$all_upper_os = all_upper_os
df$all_lower_ts = all_lower_ts
df$all_upper_ts = all_upper_ts
df$all_lower_reg = all_lower_reg
df$all_upper_reg = all_upper_reg

ggplot(df, aes(delta)) + 
  geom_vline(xintercept=2, size=0.75, linetype="dotted") + 
  geom_line(aes(y=all_error), colour="black", size=1.5,alpha=0.85) + 
  geom_line(aes(y=all_middle_os), colour="red", size=1.5,alpha=0.85) + 
  geom_line(aes(y=all_middle_ts), colour="blue", size=1.5,alpha=0.85) + 
  geom_ribbon(aes(ymin=all_lower_reg, ymax=all_upper_reg), colour ='black',alpha=0.15) + 
  geom_ribbon(aes(ymin=all_lower_os, ymax=all_upper_os), colour ='red',alpha=0.15) +
  geom_ribbon(aes(ymin=all_lower_ts, ymax=all_upper_ts), colour ='blue',alpha=0.15) +
  scale_y_continuous(limits = c(0, 2.5)) + 
  theme(axis.text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()
  )




## rand index

all_middle_os = gradual_lower[5,]
all_lower_os = gradual_lower[13,]
all_upper_os = gradual_lower[14,]

all_middle_ts = gradual_both[5,]
all_lower_ts = gradual_both[13,]
all_upper_ts = gradual_both[14,]

all_middle_reg = mean(reg_rand_d) 
all_lower_reg = quantile(reg_rand_d,0.05)
all_upper_reg = quantile(reg_rand_d,0.95)

df = data.frame(delta = pos_delta, rand_index = all_middle_reg, 
                all_middle_os = all_middle_os, all_middle_ts = all_middle_ts, 
                all_lower_os = all_lower_os, all_upper_os = all_upper_os,  
                all_lower_ts = all_lower_ts, all_upper_ts = all_upper_ts, 
                all_lower_reg = all_lower_reg, all_upper_reg = all_upper_reg, row.names = NULL)


ggplot(df, aes(delta)) + 
  geom_vline(xintercept=2, size=0.75, linetype="dotted") + 
  geom_line(aes(y=rand_index), colour="black", size=1.5,alpha=0.85) + 
  geom_line(aes(y=all_middle_os), colour="red", size=1.5,alpha=0.85) + 
  geom_line(aes(y=all_middle_ts), colour="blue", size=1.5,alpha=0.85) + 
  geom_ribbon(aes(ymin=all_lower_reg, ymax=all_upper_reg), colour ='black',alpha=0.15) + 
  geom_ribbon(aes(ymin=all_lower_os, ymax=all_upper_os), colour ='red',alpha=0.15) +
  geom_ribbon(aes(ymin=all_lower_ts, ymax=all_upper_ts), colour ='blue',alpha=0.15) +
  scale_y_continuous(limits = c(0.4, 0.9)) + 
  theme(axis.text=element_text(size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()
  )








