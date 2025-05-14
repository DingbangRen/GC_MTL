#############

{
  ini_value<-priorspecification_allnew(J=J,K=K,a_psi = 1,b_psi = 1)
  Z <- ini_value$Z
  Psi <- ini_value$Psi
  Beta <- ini_value$Beta
  Tau <- ini_value$Tau
  Lambda <- ini_value$Lambda
  ## fix the Sigma
  
  Sigma_needed = T
  #K=8
  Sigma <- diag(1,K)
  #Sigma <- Sigma_true
  #Sigma <- Sigma_block3
  
  ## for p > n, much larger alpha needs to be chosen
  alpha_const = 3
  alpha = matrix(alpha_const,nrow = J, ncol = K)
  eta <- 1
  #eta <- (1+alpha_const)^(1/2)
  
  mean_eps <- 0.1
  
  t_record <- c()
  resultlist <- list()
  errbound <- 1.5
  interval <- 1
}


{niter <- 300
  burnin <- 50
  iniburn <- 50
  Acc_prob <- rep(1,K)
  
  #regress_type = 'binary'
  regress_type = 'continuous'
}


{
########## 固定true Sigma 看看是否能帮助beta的estimate

first_Beta = Beta_alltasks_J_40
first_Beta[,4] <- abs(Beta_alltasks_J_40[,4])
preSigma = cor(matrix(rnorm(J*K, mean = 0, sd = 0.01),nrow = J) + first_Beta) %>% 
  as.matrix() %>% round(5)
colnames(preSigma) = rownames(preSigma) = NULL


first_Beta = W
preSigma = cor(matrix(rnorm(J*K, mean = 0, sd = 0.01),nrow = J) + first_Beta) %>% 
  as.matrix() %>% round(5)
colnames(preSigma) = rownames(preSigma) = NULL

}

### 可尝试fix Sigma = I 观察Copula结构的有效与否
i
for(i in 1:niter){
  if(i == 1){t_start <- Sys.time()}
  
  if(regress_type == 'binary'){
    PG <- PG_posterior_logistic(X_train,Beta)
    Beta <- Beta_posterior_binary_logistic(X_train,PG,Y_train,Psi,Tau)
    #Z_probit <- Z_posterior_binary(X_train, Y_train, Beta)
    #Beta <- Beta_posterior_binary_probit(X_train,Z_probit,Psi,Tau)
    
    Y_binary_pred <- Y_to_binary_compare0_list(X_train,Beta)
    Loss_binary <- Loss_binary_func(trueY_list = Y_train,predY_list = Y_binary_pred)
    Log_loss <- Logistic_loss(Y_train,X_train,Beta)
    #Loss_binary <- ...
  }else{
    Errvar <- Errvar_posterior(X_train, Y_train, Beta,a_err = 0.001,b_err = 0.001,K=K)
    Beta <- Beta_posterior(X_train,Y_train,Errvar,Psi,Tau) # how to improve instead of using forceSymmetric??
  }
  ### fix Beta and see if Lambda and Sigma will be properly estimated
  #Beta = Beta_alltasks_J_40 %>% as.matrix()
  #colnames(Beta) = NULL
  #Beta[which(Beta==0)] = 1e-03
  ## a_psi = b_psi = 1 尝试一下large psi samplings
  ###### for continuous
  #Psi <- Psi_posterior(Tau,Beta,a_psi=2,b_psi=2)
  #Tau <- Tau_posteriorGDP(Lambda,Beta,Psi)
  
  ###### for binary
  
  Psi <- Psi_posterior(Tau,Beta,a_psi=2,b_psi=2)
  Tau <- Tau_posteriorGDP(Lambda,Beta,Psi)
  
  
  ####### 固定Sigma 看是否能帮助Lambda的estimation 进而帮助Beta
  #Sigma = preSigma
  
  if(Sigma_needed == T){
    if(sum(Acc_prob<0.2)!=0){
      optrate = 0.85
      L = 10
      ini_eps = 0.005
      const_zero = max(floor(L/100),1)
      eps_const = max(floor(L/100),1)
    }else{
      
      #if(i %% 20 == 0){
      #  L = 600
      #  ini_eps = 0.005
      #  optrate = 0.7
      #  const_zero = 5
      #  eps_const = 5
      #}else{
      #L = max(floor(3 / mean(mean_eps)),80)
      L = max(floor(10 / mean(mean_eps)),100 )
      #optrate = runif(1,0.7,0.7)
      optrate = 0.7
      #optrate = 0.85
      #   if(i %% 15 == 0){
      #    L = 600
      #    optrate = 0.9
      #  }
      #}
      #   }
      ini_eps = 0.005
      const_zero = 10
      eps_const = 1
      #eps_const = max(floor(L/5),1)
    }
    ### adaptive simplified M-MALA
    #Lambda_and_Accprob <- adaptive_simplified_M_MALA_Lambda_GDP(Beta,Sigma,Psi,shape=alpha,rate = eta,#grad_U = -Gradient_lambda, U=-Lam_posterior,
    #                                                            epsilon=ini_eps, L=L,
    #                                                            loweps = 0.0001, uppereps = 5,
    #                                                            opt_rate = optrate, 
    #                                                            const_0 = const_zero, epsconst = eps_const,
    #                                                            inisd=0.2)
   # Lambda <- Lambda_and_Accprob[[1]]
  #  Acc_prob <- Lambda_and_Accprob[[2]]
  #  mean_eps<- Lambda_and_Accprob[[3]]

    #adaptive_full_M_MALA_Lambda_GDP
    ## full MMALA
    Lambda_and_Accprob <- adaptive_full_M_MALA_Lambda_GDP(Beta,Sigma,Psi,shape=alpha,rate = eta,
                                                            epsilon = ini_eps, L=L,
                                                            ## loweps不能太小，即使被rej也有更新的余地
                                                            ## 起码能在min step下 纠正bad initialization
                                                            loweps = 0.0001, uppereps = 5,
                                                            opt_rate = optrate, const_0 = const_zero,
                                                            epsconst = eps_const,
                                                            inisd = 0.01)
    ## 但其实acc rate 并不都是接近opt rate的，甚至还有压根没update的变量，可能与Initialization有关
    ## 但是第一种 有update但是acc rate没有接近opt rate是何原因
    Lambda <- Lambda_and_Accprob[[1]]
    Acc_prob <- Lambda_and_Accprob[[2]]
    mean_eps <- Lambda_and_Accprob[[3]]
    
    ## the convergence test??
    Sigmapre <- diag(0, K)
    #P <- floor(runif(1,2,10))
    P<-10
    for(p in 1:P){
      Sigmapre <- Sigmapre + Sigma_PX_posterior_GDP(Lambda = Lambda, Sigma=Sigma, shape = alpha,rate = eta)
    }
    Sigma <- (Sigmapre / P) %>% round(digits = 8)
    
  }else{
    Lambda <- Lambda_posteriornew(Beta = Beta, Psi = Psi,
                                  shape = alpha,rate = eta)
  }
  
  
  ### Assume fix Sigma == Identity 观察Copula结构是否有利于performance of estimation
  #Sigma = diag(1,K)
    #Sigma <- Sigma_posterior_GDP(Lambda,shape=alpha)## 为什么Sigma每次生成的值都能有很多变动
  
  if(i %% 100 == 0){
    t_record[i/100] <- Sys.time()
    print(paste(i/100,'00 iterations finished',sep=''))
    
  }
  #  if(i %% 50 == 0 & i >= 1500){
  #    save(resultlist_per50,
  #         file = paste('C:\\Users\\scott\\Resultdata_iter=',i,'.RData',sep=''))
  #  }
  if(i >= burnin & i %% interval == 0){
    resultlist_per10 <- list(
      'Psi'=Psi,
      'Beta'=Beta,
      'Tau'=Tau,
      'Lambda'=Lambda, 
      'Error Variance'=Errvar
    )
    if(Sigma_needed == T){
      resultlist_per10$`Sigma` <- Sigma
      resultlist_per10$`Acceptance prob` = Acc_prob
      resultlist_per10$`Mean step size` = mean_eps
    }
    resultlist[[i/interval-burnin/interval+1]] = resultlist_per10
  }
  if(i==niter){t_end <- Sys.time()}
  
}

#save(resultlist,file = 'C:\\Users\\scott\\Sim_data_result_DMFS.RData')
#save(resultlist,file = 'C:\\Users\\scott\\Limited_data_J=40_K=80.RData')
#load(file='C:\\Users\\scott\\Limited_data_J=40_K=80.RData')

### 
#save(resultlist,file = 'C:\\Users\\scott\\Limited_data_J=40_(100,30).RData')
#load(file = 'C:\\Users\\scott\\Limited_data_J=40_(100,30).RData')

#save(resultlist,file = 'C:\\Users\\scott\\Limited_data_J=80_(80,20).RData')
#load(file = 'C:\\Users\\scott\\Limited_data_J=80_(80,20).RData')

#save(resultlist,file = 'C:\\Users\\scott\\schooldata_result.RData')
#load(file = 'C:\\Users\\scott\\schooldata_result.RData')

#save(resultlist,file = 'C:\\Users\\scott\\isoletdata_result.RData')
#load(file = 'C:\\Users\\scott\\isoletdata_result.RData')


#save(resultlist,file = 'C:\\Users\\scott\\Limited_binary_data_J=80_(80,20).RData')

#save(resultlist,file = 'C:\\Users\\scott\\sacrosdata_result.RData')
#load(file = 'C:\\Users\\scott\\sacrosdata_result.RData')

Errvar%>%round(3) # 当errvar 太小 而covbeta 非正定问题如何解决

i#Loss_binary %>% round(3)
# need the diagostics and think about the inference for sigma
Beta%>%round(digits = 2)
Sigma
i
nsamples
########################
  
mcmcSigma <- function(n, b, task,iniburn = burnin){
  sigma_mcmc <-  c()
  
  for(j in ((b-iniburn)/interval + 1):((n-iniburn)/interval) ){
    sigma_mcmc<- rbind(sigma_mcmc,resultlist[[j]]$Sigma[,task]%>%t())
  } 
  return(mcmc(sigma_mcmc))
}

tasknames<- c()
for(k in 1:K){
  tasknames[k] <- paste('task',as.character(k),sep = ' ')
}

library(coda)
d<-1
d<-2
d<-3
d<-4
d<-5
d<-6
d<-7
d<-8

niter = i-1
burnin  = burnin + 50
mcmcSigma(n=niter,b=burnin,iniburn=burnin,task = d)%>%plot()


{
library(modeest)
## credible I
lowquan <- 0.025
upperquan <- 0.975
apply(mcmcSigma(n=niter,b=burnin,task = d)[,-d], MARGIN = 2, function(x){return(quantile(x,probs= c(lowquan, upperquan)))})%>%round(4)%>%
  as.data.frame()%>%setNames(tasknames[-d]) -> corr_CredI
rownames(corr_CredI) <-  c(paste(tasknames[d],as.character(scales::percent(lowquan,0.01)),sep=' '),
                           paste(tasknames[d],as.character(scales::percent(upperquan,0.01)),sep=' '))
corr_CredI

## produce median and posterior mean
## including 0 is not a sign of rejecting in the bayesian credible interval 

apply(mcmcSigma(n=niter,b=burnin,task = d)[,-d], MARGIN = 2, function(x){return(mean(x))})%>%round(4)%>%
  setNames(tasknames[-d]) %>% as.matrix(nrow = 1) -> corr_post_mean
colnames(corr_post_mean) = tasknames[d]
t(corr_post_mean) -> corr_post_mean
corr_post_mean


apply(mcmcSigma(n=niter,b=burnin,task = d)[,-d], MARGIN = 2, function(x){return(median(x))})%>%round(4)%>%
  setNames(tasknames[-d]) %>% as.matrix(nrow = 1) -> corr_post_median
colnames(corr_post_median) = tasknames[d]
t(corr_post_median) -> corr_post_median
corr_post_median


apply(mcmcSigma(n=niter,b=burnin,task = d)[,-d], MARGIN = 2, function(x){
  dens <- density(x)  # Estimate the density
  mode_value <- dens$x[which.max(dens$y)]  
  return(mode_value)})%>%round(4)%>%
  setNames(tasknames[-d])%>% as.matrix(nrow = 1)-> corr_post_mode

colnames(corr_post_mode) = tasknames[d]
t(corr_post_mode) -> corr_post_mode
corr_post_mode

dens <- density(x)  # Estimate the density
mode_value <- dens$x[which.max(dens$y)]  

}

########
## generating the posterior median of estimated Sigma
library(modeest)

tasknames<- c()
for(k in 1:K){
  tasknames[k] <- paste('task',as.character(k),sep = ' ')
}
{
  niter = 300
  burnin= 150
  covnmat_post_median = matrix(0,nrow = K, ncol = K)
  for(d in 1:K){
    apply(mcmcSigma(n=niter,b=burnin,task = d)[,-d], MARGIN = 2, function(x){return(median(x))})%>%round(4)%>%
     as.matrix(nrow = 1) -> corr_post_median
    t(corr_post_median) -> corr_post_median
    corr_post_median_full = append(corr_post_median, 1, after = d-1)
    covnmat_post_median[d,] = corr_post_median_full
  }
  rownames(covnmat_post_median) = tasknames
  colnames(covnmat_post_median) = tasknames
  covnmat_post_median
}
## the posterior mean
{
  niter = 300
  burnin= 150
  covnmat_post_mean = matrix(0,nrow = K, ncol = K)
  for(d in 1:K){
    apply(mcmcSigma(n=niter,b=burnin,task = d)[,-d], MARGIN = 2, function(x){return(mean(x))})%>%round(4)%>%
      as.matrix(nrow = 1) -> corr_post_mean
    t(corr_post_mean) -> corr_post_mean
    corr_post_mean_full = append(corr_post_mean, 1, after = d-1)
    covnmat_post_mean[d,] = corr_post_mean_full
  }
  rownames(covnmat_post_mean) = tasknames
  colnames(covnmat_post_mean) = tasknames
  covnmat_post_mean
}

######


#save(resultlist, file = 'C:\\Users\\scott\\GCVS_isolet_result.RData')

