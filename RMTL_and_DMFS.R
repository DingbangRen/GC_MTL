
## RMTL methods
# create simulated data for regression and classification problem
library(RMTL)
{
nMSE <- function(X,Y,beta){
  estimate <- X %*% c(beta)
  nmse <- (Y-estimate)^2 %>% sum() / 
    (nrow(X) * var(Y) )
  return(nmse)
}
rMSE <- function(X,Y,beta){
  estimate <- X %*% c(beta)
  ori_mse <- (Y-estimate)^2 %>% sum() / 
    (nrow(X))
  
  return((ori_mse)^(1/2))
}
Rsquare <-function(X,Y,beta){
  estimate <- X %*% c(beta)
  denomi <- (Y - estimate)^2 %>% sum()
  nomina <- ((Y - mean(Y))^2 %>% sum())
  return(1 - denomi/nomina)
}
}
#####
## RMTL
### 此处数据输入的是 all tasks data
rmtl_type = "Regression"
#rmtl_type = "Classification"
{
  Xlist_new <- X_train
  Ylist_new <- Y_train
  ## L21
  cvfitr_L21_simulation<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="L21")
  
  # train
  modelRMTL_L21_simulation<-MTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="L21",Lam2 = 0,
                                Lam1=cvfitr_L21_simulation$Lam1.min, Lam1_seq=cvfitr_L21_simulation$Lam1_seq)
  
  ## Lasso
  cvfitr_Lasso_simulation<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Lasso")
  
  # train
  modelRMTL_Lasso_simulation<-MTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Lasso",
                                  Lam1=cvfitr_Lasso_simulation$Lam1.min, Lam1_seq=cvfitr_Lasso_simulation$Lam1_seq)
  
  ## Trace
  cvfitr_Trace_simulation<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Trace")
  
  # train
  modelRMTL_Trace_simulation<-MTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Trace",
                                  Lam1=cvfitr_Trace_simulation$Lam1.min, Lam1_seq=cvfitr_Trace_simulation$Lam1_seq)
}

###############
{
  G_network_mean_reg <- matrix(0,nrow = K, ncol = K)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      G_network_mean_reg[i,j] <- -1/K
    }
  }
  
  G_network_mean_reg <- t(G_network_mean_reg)+G_network_mean_reg+diag((K-1)/K,K)
  
  cvfitr_network_mean_reg<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Graph", G=G_network_mean_reg)
  #Train
  modelRMTL_network_mean_reg_simulation <- MTL(Xlist_new, Ylist_new,type=rmtl_type, Regularization="Graph", 
                                               Lam1=cvfitr_network_mean_reg$Lam1.min, Lam1_seq=cvfitr_network_mean_reg$Lam1_seq, G=G_network_mean_reg) 
}
### only use when simulation
{
  ## network structure
  
  G_network_mean_reg <- matrix(0,nrow = K, ncol = K)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      G_network_mean_reg[i,j] <- -1/K
    }
  }
  
  G_network_mean_reg <- t(G_network_mean_reg)+G_network_mean_reg+diag((K-1)/K,K)
  
  ##
  ## C_4^2 + 1 = 7
  G_network_priornet <- matrix(0,nrow = 7, ncol = K)
  G_network_priornet[1,] <- c(1,-1,0,0,0,0,0,0)
  G_network_priornet[2,] <- c(1,0,-1,0,0,0,0,0)
  G_network_priornet[3,] <- c(1,0,0,-1,0,0,0,0)
  G_network_priornet[4,] <- c(0,1,-1,0,0,0,0,0)
  G_network_priornet[5,] <- c(0,1,0,-1,0,0,0,0)
  G_network_priornet[6,] <- c(0,0,1,-1,0,0,0,0)
  G_network_priornet[7,] <- c(0,0,0,0,1,-1,0,0)
  
  G_network_priornet <- t(G_network_priornet)
  
  
  cvfitr_network_mean_reg<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Graph", G=G_network_mean_reg)
  #Train
  modelRMTL_network_mean_reg_simulation <- MTL(Xlist_new, Ylist_new,type=rmtl_type, Regularization="Graph", 
                                               Lam1=cvfitr_network_mean_reg$Lam1.min, Lam1_seq=cvfitr_network_mean_reg$Lam1_seq, G=G_network_mean_reg)
  
  ## 
  cvfitr_network_priornet<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="Graph", G=G_network_priornet)
  #Train
  modelRMTL_network_priornet_simulation <- MTL(Xlist_new, Ylist_new,type=rmtl_type, Regularization="Graph", 
                                               Lam1=cvfitr_network_priornet$Lam1.min, Lam1_seq=cvfitr_network_priornet$Lam1_seq, G=G_network_priornet)
  
  ###K-means clustering
  #Create datasets
  K_simulation = 4
  
  cvfit_clustering<-cvMTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="CMTL", k=K_simulation)
  #> Loading required namespace: MASS
  #> Loading required namespace: psych
  modelRMTL_clustering = MTL(Xlist_new, Ylist_new, type=rmtl_type, Regularization="CMTL", 
                             Lam1=cvfit_clustering$Lam1.min, Lam1_seq=cvfit_clustering$Lam1_seq, k=K_simulation)
}


#####
## Dirty model for feature selection:
#install.packages('doMC')
#library(NORMT3) not available
{
  library(doMC)
  source('C:\\Users\\scott\\Downloads\\R-DMFS\\beta.R')
  source('C:\\Users\\scott\\Downloads\\R-DMFS\\ep_robust-multi-task-natural_noise.R')
  source('C:\\Users\\scott\\Downloads\\R-DMFS\\probit.R')
  source('C:\\Users\\scott\\Downloads\\R-DMFS\\SB_prior.R')
  source('C:\\Users\\scott\\Downloads\\R-DMFS\\SB_prior_natural.R')
  source('C:\\Users\\scott\\Downloads\\R-DMFS\\ep_robust-multi-task-natural_noise.R')
  
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\myInitMessages.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\dnormt3.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\dst.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\erf.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\erfc.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\ic1.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\is1.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\normt3ip.R')
  source('C:\\Users\\scott\\Downloads\\NORMT3\\R\\wofz.R')
  
  ### 源代码里面的这个原函数用不了 需要简化成以下形式
  erfc <- function(z){
    ansx=as.double(Re(z))
    ansy=as.double(Im(z))
    return(complex(real=ansx, imaginary=ansy))
    
  }
  ### 增加以下ep的niter_dmfs, 原来的niter_dmfs = 200, 可以适当增一些
  ep_RM <- function(X, Y, class = FALSE, a0 = rep(1, 5), b0 = rep(1, 5), verbose = FALSE, 
                    alpha0 = 5, beta0 = 5, initial_damping = NULL, n_cores = 1, shared_noise = TRUE,
                    niter_dmfs = 2e2) {
    
    # We initialize the posterior approximation
    
    a <- initialize_approximation(X, Y, class, a0, b0, alpha0, beta0, n_cores, shared_noise)
    
    iters <- 0
    convergence <- FALSE
    
    if (is.null(initial_damping)){
      damping <- .5
      
    }else{
      damping <- initial_damping
      
    }
    
    q <- compute_titled_distribution(a)
    
    evidence_old <- +Inf
    
    while (iters < niter_dmfs && ! convergence) {
      
      aOld <- a
      
      valid_parameters <- FALSE
      factor <- 1
      
      while (! valid_parameters) {
        
        a <- aOld
        
        # We process the factors 
        
        if (class == TRUE){
          a <- process_likelihood_factors(q, a, factor * damping)
          
        }else{
          a <- process_likelihood_factors_regression(q, a, factor * damping)
          
        }
        
        a <- process_prior_model_coefficients(q, a, factor * damping)
        
        a <- process_prior_binary_latent_variables(q, a, factor * damping)
        
        qNew <- compute_titled_distribution(a)
        
        valid_parameters <- all(qNew$v_w > 0) & all(sapply(qNew$v_targets, function(x) all(x > 0))) & 
          qNew$a_z > -1 & qNew$b_z > -1 & qNew$a_omega > -1 & qNew$b_omega > -1 & qNew$a_gamma > -1 & qNew$b_gamma > -1 &
          qNew$a_tau > -1 & qNew$b_tau > -1 & qNew$a_eta > -1 & qNew$b_eta > -1
        
        if (valid_parameters == FALSE) {
          factor <- factor * 0.5
          if (verbose)
            cat("\tReducing damping size:", damping * factor, "\n")
        }
        
      }
      
      # We look for convergence in the EP algorithm
      
      cat("iter:", iters, " ")
      
      #		convergence <- check_convergence(a, aOld, verbose, damping)
      convergence <- check_convergence_q(qNew, q, TRUE, damping)
      
      iters <- iters + 1
      
      damping <- damping * 0.99
      
      q <- qNew
      
      #		evidence_new <- evaluate_evidence(a, class)
      #		cat("Itr:", iters, "logZ:", evaluate_evidence(a, class), "Change:", evidence_new - evidence_old, 
      #			"Percnt:", abs(evidence_new - evidence_old) / abs(evidence_new), "Damping:", damping, "\n")
      #		if (abs(evidence_new - evidence_old) / abs(evidence_new) < 1e-3) convergence <- TRUE
      #		evidence_old <- evidence_new
    }
    
    # We evaluate the models evidence
    
    logZ <- evaluate_evidence(a, class)
    
    q <- compute_titled_distribution(a)
    
    # We compute the probs that the coefficients are different from zero
    
    P1 <- sigmoid(matrix(q$logit_p_z, a$nTasks, a$d, byrow = TRUE)) * sigmoid(q$logit_p_eta)
    P2 <- sigmoid(matrix(-q$logit_p_z, a$nTasks, a$d, byrow = TRUE)) * sigmoid(matrix(q$logit_p_omega, a$nTasks, a$d)) * sigmoid(q$logit_p_tau)
    P3 <- sigmoid(matrix(-q$logit_p_z, a$nTasks, a$d, byrow = TRUE)) * sigmoid(matrix(-q$logit_p_omega, a$nTasks, a$d)) * 
      sigmoid(matrix(q$logit_p_gamma, a$nTasks, a$d, byrow = TRUE))
    
    P <- P1 + P2 + P3
    
    # We return the model evidence, the posterior approximation and the feature ranking
    # W and P have in each column a task and in each row a model coefficient
    
    list(logZ = logZ, a = a, q = q, P = t(P), W = t(q$m_w))
  }
}


compute_titled_distribution <- function(a) {
  
  # We set some useful constants. 
  
  d <- a$d
  
  v_w <- m_w <- matrix(0, a$nTasks, d)
  m_targets <- list()
  v_targets <- list()
  
  # These are the posterior of the targets one multiplied by the delta function
  
  m_targets_delta <- list()
  v_targets_delta <- list()
  
  # First we reconstruct the Gaussian approximation over the Targets
  
  for (k in 1 : a$nTasks) {
    v_targets[[ k ]] <- a$fTildes$vTilde[[ k ]]^-1
    m_targets[[ k ]] <- v_targets[[ k ]] * a$fTildes$mTilde[[ k ]]
  }
  
  # Second we reconstruct the Bernoulli approximation over the latent variables
  
  logit_p_z <- colSums(a$gTildes$zTilde) + a$hTildes_z$logit_zTilde
  logit_p_omega <- rowSums(a$gTildes$omegaTilde) + a$hTildes_omega$logit_omegaTilde
  logit_p_gamma <- colSums(a$gTildes$gammaTilde) + a$hTildes_gamma$logit_gammaTilde
  logit_p_tau <- a$gTildes$tauTilde + a$hTildes_tau$logit_tauTilde
  logit_p_eta <- a$gTildes$etaTilde + a$hTildes_eta$logit_etaTilde
  
  # If it is a classification problem we assume the first feature is an oulier feature
  # that is always activated
  
  if (a$class == TRUE) {
    logit_p_z[ 1 ] <- 1e4
    logit_p_eta[ ,1 ] <- 1e4
    logit_p_gamma[ 1 ] <- 0
    logit_p_tau[ ,1 ] <- 0
  }
  
  # Fourth we reconstruct the Beta approximation over the prior probabilities of the latent variables
  
  a_z <- sum(a$hTildes_z$a_zTilde) + a$lTildes_z$a_Tilde 
  b_z <- sum(a$hTildes_z$b_zTilde) + a$lTildes_z$b_Tilde 
  
  a_omega <- sum(a$hTildes_omega$a_omegaTilde) + a$lTildes_omega$a_Tilde 
  b_omega <- sum(a$hTildes_omega$b_omegaTilde) + a$lTildes_omega$b_Tilde 
  
  a_gamma <- sum(a$hTildes_gamma$a_gammaTilde) + a$lTildes_gamma$a_Tilde 
  b_gamma <- sum(a$hTildes_gamma$b_gammaTilde) + a$lTildes_gamma$b_Tilde 
  
  a_tau <- sum(a$hTildes_tau$a_tauTilde) + a$lTildes_tau$a_Tilde 
  b_tau <- sum(a$hTildes_tau$b_tauTilde) + a$lTildes_tau$b_Tilde 
  
  a_eta <- sum(a$hTildes_eta$a_etaTilde) + a$lTildes_eta$a_Tilde 
  b_eta <- sum(a$hTildes_eta$b_etaTilde) + a$lTildes_eta$b_Tilde 
  
  # Fourth we reconstruct the Gaussian approximation over the model coefficients 
  
  if (a$n_cores != 1)  {
    
    result <- foreach(k = 1 : a$nTasks) %dopar% {
      
      n <- nrow(a$X[[ k ]])
      
      # We extract the mean and the variance of the marginals
      # We compute first the lambda matrix and the upsilon vector
      
      tX <- t.default(a$X[[ k ]])
      
      dX <- matrix(a$gTildes$nuTilde[ k, ]^-1, n, d, byrow = TRUE) * a$X[[ k ]]
      #		U <- solve(diag(v_targets[[ k ]]) + a$X[[ k ]] %*% t.default(dX))
      U <- matrix(v_targets[[ k ]]^-1, n, n, byrow = TRUE) * solve(diag(length(v_targets[[ k ]])) + 
                                                                     (matrix(v_targets[[ k ]]^-1, n, d) * a$X[[ k ]]) %*% t.default(dX))
      UdX <- U %*% dX
      
      v_w <- a$gTildes$nuTilde[ k, ]^-1 - colSums(UdX * dX)
      
      v_aux <- c((m_targets[[ k ]] / v_targets[[ k ]]) %*% a$X[[ k ]] + a$gTildes$muTilde[ k, ])
      
      m_w <- c(a$gTildes$nuTilde[ k, ]^-1 * v_aux) - c(t.default(UdX) %*% (dX %*% v_aux))
      
      # Now we compute the posterior of the targets times the delta function 
      
      v_targets_delta <- rowSums(dX * a$X[[ k ]]) -  rowSums((a$X[[ k ]] %*% t(dX) %*% UdX) * a$X[[ k ]])
      m_targets_delta <- c(a$X[[ k ]] %*% m_w)
      
      list(m_w = m_w, v_w = v_w, v_targets_delta = v_targets_delta, m_targets_delta = m_targets_delta)
    }
    
    for (k in 1 : a$nTasks) {
      m_w[ k, ] <- result[[ k ]]$m_w
      v_w[ k, ] <- result[[ k ]]$v_w
      v_targets_delta[[ k ]] <- result[[ k ]]$v_targets_delta
      m_targets_delta[[ k ]] <- result[[ k ]]$m_targets_delta
    }
    
  } else {
    for (k in 1 : a$nTasks) {
      
      n <- nrow(a$X[[ k ]])
      
      # We extract the mean and the variance of the marginals
      # We compute first the lambda matrix and the upsilon vector
      
      tX <- t.default(a$X[[ k ]])
      
      dX <- matrix(a$gTildes$nuTilde[ k, ]^-1, n, d, byrow = TRUE) * a$X[[ k ]]
      #		U <- solve(diag(v_targets[[ k ]]) + a$X[[ k ]] %*% t.default(dX))
      U <- matrix(v_targets[[ k ]]^-1, n, n, byrow = TRUE) * solve(diag(length(v_targets[[ k ]])) + 
                                                                     (matrix(v_targets[[ k ]]^-1, n, d) * a$X[[ k ]]) %*% t.default(dX))
      UdX <- U %*% dX
      
      v_w[ k, ] <- a$gTildes$nuTilde[ k, ]^-1 - colSums(UdX * dX)
      
      v_aux <- c((m_targets[[ k ]] / v_targets[[ k ]]) %*% a$X[[ k ]] + a$gTildes$muTilde[ k, ])
      
      m_w[ k, ] <- c(a$gTildes$nuTilde[ k, ]^-1 * v_aux) - c(t.default(UdX) %*% (dX %*% v_aux))
      
      # Now we compute the posterior of the targets times the delta function 
      
      v_targets_delta[[ k ]] <- rowSums(dX * a$X[[ k ]]) -  rowSums((a$X[[ k ]] %*% t(dX) %*% UdX) * a$X[[ k ]])
      m_targets_delta[[ k ]] <- c(a$X[[ k ]] %*% m_w[ k, ])
    }
  }
  
  # Now the noise. We only work with natural parameters
  
  if (a$class == FALSE) {
    if (a$shared_noise == TRUE) {
      alpha <- rep(sum(sapply(a$fTildes$aTilde, function(x) sum(x))) + a$kTildes$aTilde, a$nTasks)
      beta <- rep(sum(sapply(a$fTildes$bTilde, function(x) sum(x))) + a$kTildes$bTilde, a$nTasks)
    } else {
      alpha <- sapply(a$fTildes$aTilde, function(x) sum(x)) + a$kTildes$aTilde
      beta <- sapply(a$fTildes$bTilde, function(x) sum(x)) + a$kTildes$bTilde
    }
    
    alphaNoise <- alpha 
    betaNoise <- beta
  } else {
    alphaNoise <- NULL
    betaNoise <- NULL
  }
  
  # We return a list with all the quantities computed
  
  list(m_w = m_w, v_w = v_w, m_targets = m_targets, v_targets = v_targets, logit_p_z = logit_p_z, 
       logit_p_omega = logit_p_omega, logit_p_gamma = logit_p_gamma, logit_p_tau = logit_p_tau, 
       logit_p_eta = logit_p_eta, a_z = a_z, b_z = b_z, a_omega = a_omega, b_omega = b_omega,
       a_gamma = a_gamma, b_gamma = b_gamma, a_tau = a_tau, b_tau = b_tau, a_eta = a_eta, b_eta = b_eta,
       m_targets_delta = m_targets_delta, v_targets_delta = v_targets_delta, alphaNoise = alphaNoise, betaNoise = betaNoise)
}


{
  Y_train_DMFS <- list()
  ## 需要将y response转为vector, 原数据中是matrix
  for(k in 1:K){
    Y_train_DMFS[[k]] <- Y_train[[k]]%>%c()
  }
  
  
  X = X_train
  Y = Y_train_DMFS
  ## alpha_0 beta_0是给sigma_k^2的hyperparameter 见原文P3
  ## a0 = b0 = 1 recommended in P3 but perform badly
  
  beta_const = 1
  #beta_const = 0.01
  ## for sacros
  #beta_const = 0.001
  #beta_const = 0.01
  ### for school beta should be
  #beta_const = 0.005
  ###
  ing_gamma_hyper = 10 ## 这是task response variation, 可以调整他来使得regression适应更high variation
  dmfs_results_simulation = ep_RM(X=X_train, Y=Y_train_DMFS, class = FALSE,
                                  a0 = rep(beta_const , 5), b0 = rep(beta_const, 5), 
                                  alpha0 = ing_gamma_hyper, beta0 = ing_gamma_hyper, shared_noise = FALSE,niter_dmfs = 300)
  ##
}
dmfs_results_simulation$W %>% round(2)

###########
{## bayesian kernelized MTL
  source('C:\\Users\\scott\\Downloads\\Bayesian Kernelized MTL R\\bayesian_multitask_multiple_kernel_learning_train.R')
  source('C:\\Users\\scott\\Downloads\\Bayesian Kernelized MTL R\\bayesian_multitask_multiple_kernel_learning_test.R')
  
  #initalize the parameters of the algorithm
  parameters_bmtmkl <- list()
  
  #set the hyperparameters of gamma prior used for sample weights
  parameters_bmtmkl$alpha_lambda <- 1
  parameters_bmtmkl$beta_lambda <- 1
  
  #set the hyperparameters of gamma prior used for intermediate noise
  parameters_bmtmkl$alpha_upsilon <- 1
  parameters_bmtmkl$beta_upsilon <- 1
  
  #set the hyperparameters of gamma prior used for bias
  parameters_bmtmkl$alpha_gamma <- 1
  parameters_bmtmkl$beta_gamma <- 1
  
  #set the hyperparameters of gamma prior used for kernel weights
  parameters_bmtmkl$alpha_omega <- 1
  parameters_bmtmkl$beta_omega <- 1
  
  #set the hyperparameters of gamma prior used for output noise
  parameters_bmtmkl$alpha_epsilon <- 1
  parameters_bmtmkl$beta_epsilon <- 1
  
  ### IMPORTANT ###
  #For gamma priors, you can experiment with three different (alpha, beta) values
  #(1, 1) => default priors
  #(1e-10, 1e+10) => good for obtaining sparsity
  #(1e-10, 1e-10) => good for small sample size problems (like in Nature Biotechnology paper)
  
  #set the number of iterations
  parameters_bmtmkl$iteration <- 400
  
  #determine whether you want to calculate and store the lower bound values
  parameters_bmtmkl$progress <- 0
  
  #set the seed for random number generator used to initalize random variables
  parameters_bmtmkl$seed <- 1606
  
  #set the number of tasks (e.g., the number of compounds in Nature Biotechnology paper)
  ## K
  #set the number of kernels (e.g., the number of views in Nature Biotechnology paper)
  ## 此处是one view
  P <- 1
  
  #initialize the kernels and outputs of each task for training
  Ktrain <- vector("list", K)
  ytrain <- vector("list", K)
  
  ## 关于kernel的设置 可以参考https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.2877/MediaObjects/41587_2014_BFnbt2877_MOESM1_ESM.pdf
  ## 该文章的P6
  ## 简单来说就是 将每个task的training sample进行Z normalization之后
  ## 计算基于sample的gaussian kernel as a similarity metric, 进而得到task-wise kernel\in R ^{N_k \times N_k}
  
  gaussian_kernel <- function(Xtrain_k){
    var_kernel = Xtrain_k %>% ncol() ## variance估计为 num of features
    N_k = nrow(Xtrain_k)
    Kernel_k = matrix(0,nrow = N_k,ncol = N_k)
    
    for(i in 1:(N_k-1)){
      for(j in (i+1):N_k)
        Kernel_k[i,j] = exp( - 0.5 * sum((Xtrain_k[i,] - Xtrain_k[j,])^2) / var_kernel)
    }
    final_kernel_K = Kernel_k + t(Kernel_k) + diag(rep(1,N_k), ncol = N_k)
    return(final_kernel_K)
  }
  
  ### training num \times testing num similarity
  gaussian_kernel_test <- function(Xtrain_k,Xtest_k){
    var_kernel = Xtrain_k %>% ncol() ## variance估计为 num of features
    nrow_k = nrow(Xtrain_k)
    ncol_k = nrow(Xtest_k)
    Kernel_k = matrix(0,nrow = nrow_k,ncol = ncol_k)
    
    for(i in 1:nrow_k){
      for(j in 1:ncol_k){
        Kernel_k[i,j] = exp( - 0.5 * sum((Xtrain_k[i,] - Xtest_k[j,])^2) / var_kernel)
      }
    }
    
    return(Kernel_k)
  }
  
  
  for (k in 1:K) {
    Nk = nrow(X_train[[k]])
    Ktrain[[k]] <- gaussian_kernel(Xtrain_k = X_train[[k]]) %>% array(dim = c(Nk,Nk,1))
    #should be an Ntra x Ntra x P matrix containing similarity values between training samples of task t
    ## denote the similarity metric as gaussian kernel with variance == J
    ytrain[[k]] <- Y_train[[k]] 
    
    #should be an Ntra x 1 matrix containing target outputs of task t
  }
  
  #perform training
  bmtmkl_train_result <- bayesian_multitask_multiple_kernel_learning_train(Ktrain, ytrain, parameters_bmtmkl)
  
  #display the kernel weights
  ## 由于P = 1 所以此项在此无效
  #print(bmtmkl_train_result$be$mu[(K+1):(K+P)])
  
  #initialize the kernels of each task for testing
  Ktest <- vector("list", K)
  for (k in 1:K) {
    Nk_test = nrow(X_test[[k]])
    Nk_train = nrow(X_train[[k]])
    Ktest[[k]] <- gaussian_kernel_test(Xtrain_k=X_train[[k]],
                                       Xtest_k = X_test[[k]])%>% array(dim = c(Nk_train,Nk_test,1))
    #should be an Ntra x Ntest x P matrix containing similarity values between training and test samples of task t
  }
  
  #perform prediction
  prediction <- bayesian_multitask_multiple_kernel_learning_test(Km = Ktest, state = bmtmkl_train_result)
  
  #display the predictions for each task
  for (k in 1:K) {
    print(prediction$f[[k]]$mean)
  }
}
