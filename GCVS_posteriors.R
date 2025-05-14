library(dplyr)
library(LaplacesDemon)
require(psych)
library('truncnorm')
library('BayesLogit')
################ 
## prior specification
{
priorspecification_all <- function(J=J,K=K,
                                   alpha = 1,eta = 1,
                                   a_psi = 1,b_psi = 1){
  
  v_0 = K+1
  V_0 = diag(rep(1, K))
  
  Lambda_ini <- matrix(0, nrow = J, ncol = K)
  Tau_ini <- matrix(0, nrow = J, ncol = K)
  
  V <- rinvwishart(nu = v_0, S=v_0 * V_0 )
  Sigma_ini <- psych::cor.smooth(cov2cor(V)) %>% round(digits = 8)
  
  Z_ini <- matrix(0, J,K)
  Psi_ini <- c()
  Beta_ini <- matrix(0,J,K)
  
  # initialization 
  for(j in 1:J){
    Z_ini[j,] <- rmvn(n=1, mu = rep(0,K),Sigma=Sigma_ini)[1,] ## drop the rownames
    Psi_ini[j] <- rinvgamma(1,a_psi/2,b_psi/2)
    for(k in 1:K){
      Lambda_ini[j,k] <- qgamma(pnorm(Z_ini[j,k],0,1),shape = alpha,rate = eta) ## F^{-1}(\Phi())
      Tau_ini[j,k] <- rexp(1,rate = Lambda_ini[j,k])
      Beta_ini[j,k] <- rnorm(1,mean = 0, sd = sqrt(Psi_ini[j]*Tau_ini[j,k])) # note that the sd requires the standard deviation
    }
  }
  resultlist <- list('Z' = Z_ini,
                     'Psi'=Psi_ini,
                     'Beta'=Beta_ini,
                     'Tau'=Tau_ini,
                     'Lambda'=Lambda_ini,
                     'Sigma'=Sigma_ini)
  return(resultlist)
}

priorspecification_allnew <- function(J=J, K=K,
                                      alpha = 1,eta = 1,
                                      a_psi = 1,b_psi = 1){
  
  v_0 = K+1
  V_0 = diag(rep(1, K))
  Lambda_ini <- matrix(0, nrow = J, ncol = K)
  Tau_ini <- matrix(0, nrow = J, ncol = K)
  
  V <- rinvwishart(nu = v_0, S=v_0 * V_0 )
  Sigma_ini <- psych::cor.smooth(cov2cor(V)) %>% round(digits = 8)
  
  Z_ini <- matrix(0, J,K)
  Psi_ini <- c()
  Beta_ini <- matrix(0,J,K)
  
  # initialization 
  for(j in 1:J){
    Z_ini[j,] <- rmvn(n=1, mu = rep(0,K),Sigma=Sigma_ini)[1,] ## drop the rownames
    Psi_ini[j] <- rinvgamma(1,a_psi/2,b_psi/2)
    for(k in 1:K){
      Lambda_ini[j,k] <- rgamma(1,alpha,eta)
      Tau_ini[j,k] <- qexp(pnorm(Z_ini[j,k],0,1),rate = Lambda_ini[j,k]) ## F^{-1}(\Phi())
      Beta_ini[j,k] <- rnorm(1,mean = 0, sd = sqrt(Psi_ini[j]*Tau_ini[j,k])) # note that the sd requires the standard deviation
    }
  }
  resultlist <- list('Z' = Z_ini,
                     'Psi'=Psi_ini,
                     'Beta'=Beta_ini,
                     'Tau'=Tau_ini,
                     'Lambda'=Lambda_ini,
                     'Sigma'=Sigma_ini)
  return(resultlist)
}}

##################

{
# Gaussian copula density under gamma marginals
library(extras)
## Z转换 检验没问题

Z_trans <- function(lambda,shape=alpha,rate=eta){ # transform \lambda_j to z_j
  z <- c()
  for(k in 1:K){
    if(pgamma(lambda[k],shape,rate)==1){ 
      p <- round(1 - 1e-16,digits = 16)}else{
        p <- pgamma(lambda[k],shape,rate)
      }
    z[k] <- qnorm(p,0,1)
  }
  return(z)
}

Lam_likelihood <- function(lambda,Sigma,shape=alpha,rate=eta){
  z<- Z_trans(lambda)
  expinner <- t(z) %*% (diag(rep(1,K))-solve(Sigma)) %*% z / 2
  GClikelihood <- pow(det(Sigma),-1/2) * exp(expinner)%>%c()
  lambdalikelihood<-lapply(lambda,
                           FUN = function(x){return(dgamma(x,shape,rate))})%>%unlist()%>%prod()
  return(GClikelihood*lambdalikelihood)
}


### f(\lambda_j | -) \propto (\tau_j | \lambda_j)(\lambda_j | \Sigma)
Lam_posterior <- function(lambda,tau,Sigma,shape=alpha,rate=eta){
  lam_likeli <- Lam_likelihood(lambda,Sigma)
  tau_likevector <- c()
  for(k in 1:K){
    tau_likevector[k] <- dexp(tau[k],rate =lambda[k])  
  }
  return(lam_likeli * prod(tau_likevector))
}

Q_tcomputation <- function(lambda,z,shape=alpha,rate = eta){
  Q_t <- matrix(0,K,K) # Q^\top
  for(k in 1:K){
    Q_t[k,k] <- (1/dnorm(z[k],0,1)) * dgamma(lambda[k],shape,rate)
  }
  return(Q_t)
}

## Gradient is obtained under log posterior !!!
Gradient_lambda <- function(lambda,Sigma,tau,shape=alpha,rate=eta){
  z <- Z_trans(lambda)
  Q_t <- Q_tcomputation(lambda,z)
  part1 <- Q_t %*% (diag(1,K)-solve(Sigma)) %*% z
  part2 <- (shape-1)/lambda-rate
  part3 <- 1/lambda - tau
  return(part1+part2+part3)
}

Hessian_lambda <- function(lambda,Sigma,shape=alpha,rate=eta){
  z <- Z_trans(lambda)
  Q_t <- Q_tcomputation(lambda,z)
  q <- c()
  for(k in 1:K){
    q[k] <- (z[k]/(dnorm(z[k],0,1))^2)*(dgamma(lambda[k],shape,rate))^2 + 
      (1/dnorm(z[k],0,1))*((rate^shape/gamma(shape))* (lambda[k]^(shape-2)) * exp(-rate*lambda[k])*
                             (shape-1+lambda[k]*(-rate)))
  }
  part1<-t(Q_t %*% (diag(1,K)-solve(Sigma))%*% t(Q_t))%>%round(digits = 8) +
    diag(q)%*%diag( (diag(1,K)-solve(Sigma))%*%z %>%c() )
  # isSymmetric(t(Q_t %*% (diag(1,K)-solve(Sigma))%*% t(Q_t))%>%round(digits = 6)): True
  part2 <- diag((1-shape)/(lambda)^2)
  part3 <- diag(-1/(lambda)^2)
  return(t(part1+part2+part3)) ## in principle, H(f) = J(gradient f)^\top
}

###########
####### posterior functions
Z_posterior_binary=function(X, Y, Beta){
  Z <- list()
  for(k in 1:K){
    z_k <- rep(0,length(Y[[k]]))
    for(i in 1:length(Y[[k]])){
      
      if(Y[[k]][i] == 0){
        
        z_k[i] <- rtruncnorm(n=1, a=-Inf, b=0, mean = X[[k]][i,] %*% Beta[,k], sd = 1)
        
      }else{
        
        z_k[i] <- rtruncnorm(n=1, a=1e-08, b=Inf, mean = X[[k]][i,] %*% Beta[,k], sd = 1)
        
      }
    }
    Z[[k]] <- z_k %>% as.matrix()
  }
  return(Z)
}

Beta_posterior_binary_probit=function(X,Z,Psi,Tau){
  smooth_invA <- function(ori_A){
    eigen_A <- eigen(ori_A)
    eigenvalue <- eigen_A$values
    eigenvec <- eigen_A$vectors
    ## 一般导致cov非正定的 并非因为eigen=0 而是因为下面的round 但round是为了确保symmteric
    #eigenvalue[which(eigenvalue==0)] <- 1e-160
    ## 用最简单的 abs(eigen) 但容易出现很小的eigen 被系统判定为0
    #newG_j <- eigen_G$vectors %*% diag(abs(eigenvalue)) %*% t(eigen_G$vectors)
    
    ## 如果用coth的方法 还能方便对很小的值进行soft regularize
    ## 
    smooth_inveigen <- diag(1 / (eigenvalue * sapply(eigenvalue,FUN=coth)))
    newinvA <- eigenvec %*% smooth_inveigen %*% t(eigenvec)
    return(newinvA)}
  
  Beta <- matrix(0,J,K)
  for(k in 1:K){## solve(diag(Psi) %*% diag(Tau[,k]))
    A <- (diag(1/Psi) %*% diag(1/Tau[,k])) + (t(X[[k]]) %*% X[[k]])
    
    invA <- smooth_invA(A)
    #Cov_beta <- as.matrix(forceSymmetric(Errvar[k] * invA))%>% round(digits = digi) #故而只能在最终结果上限制小数点位数, 还需要适应性的改变digit
    Cov_beta <- (invA %>% forceSymmetric() %>% as.matrix())
    
    mu_beta <- invA %*% t(X[[k]]) %*% Z[[k]] %>% c()
    Beta[,k] <- rmvn(n=1, mu=mu_beta, Sigma=Cov_beta)
  }
  return(Beta)
}

##### polya gamma augmentation
pg_sampling <- function(a,d,nsum = 500){
  
  g <- rgamma(nsum,a,1)
  
  gg <- c(rep(0,nsum))
  for ( k in 1:nsum){
    gg[k] <- g[k] / ((k-1/2)^2 + d^2/(4*pi^2))
  }
  return (sum(gg) / (2*pi^2) )
}

PG_posterior_logistic=function(X,Beta,n_logit = 1){
  PGlist <- list()
  for(k in 1:K){
    PGvec_k <- c()
    for(i in 1:nrow(X[[k]])){
      ## apply the PG function from source code
      PGvec_k[i] <- rpg.devroye(1, h=n_logit, z=X[[k]][i,] %*% Beta[,k])
    }
    PGlist[[k]] <- PGvec_k %>% as.matrix()
  }
  return(PGlist)
}

Beta_posterior_binary_logistic <- function(X,PG,Y,Psi,Tau,N_trial = rep(1,K)){
  Beta <- matrix(0,J,K)
  
  smooth_invA <- function(ori_A){
    eigen_A <- eigen(ori_A)
    eigenvalue <- eigen_A$values
    eigenvec <- eigen_A$vectors
    ## 一般导致cov非正定的 并非因为eigen=0 而是因为下面的round 但round是为了确保symmteric
    #eigenvalue[which(eigenvalue==0)] <- 1e-160
    ## 用最简单的 abs(eigen) 但容易出现很小的eigen 被系统判定为0
    #newG_j <- eigen_G$vectors %*% diag(abs(eigenvalue)) %*% t(eigen_G$vectors)
    
    ## 如果用coth的方法 还能方便对很小的值进行soft regularize
    ## 
    smooth_inveigen <- diag(1 / (eigenvalue * sapply(eigenvalue,FUN=coth) ))
    newinvA <- eigenvec %*% smooth_inveigen %*% t(eigenvec)
    return(newinvA)
  }
  
  for(k in 1:K){## solve(diag(Psi) %*% diag(Tau[,k]))
    Omega_k <- PG[[k]] %>% c() %>% diag()
    kappa_k <- Y[[k]] - (N_trial[k]/2)
    
    A <- (diag(1/Psi) %*% diag(1/Tau[,k])) + (t(X[[k]]) %*% Omega_k %*% X[[k]])
    #invA <- solve(A) ## A如果数值过大 则solve(A)会很小 进而有些数无法被读取
    # 但若令数值过小 则优惠出现不对称的情况 甚至有些地方就是1e-12
    invA <- A %>% smooth_invA()
    
    Cov_beta <- invA %>% forceSymmetric() %>% as.matrix() #故而只能在最终结果上限制小数点位数, 还需要适应性的改变digit
    mu_beta <- (invA %*% t(X[[k]]) %*% kappa_k) %>% c()
    Beta[,k] <- rmvn(n=1, mu=mu_beta, Sigma=Cov_beta)
  }
  return(Beta)
}


Errvar_posterior=function(X,Y,Beta,a_err = 1, b_err = 1,K=K){
  errvar <- c()
  for(k in 1:K){
    nsam <- length(Y[[k]])
    resid <- Y[[k]] - X[[k]] %*% Beta[,k]
    errvar[k] <- rinvgamma(n=1, shape = (nsam+a_err) / 2, 
                           scale = (t(resid)%*%resid+b_err) / 2)
  }
  return(errvar)
}


library(matrixcalc)
library(Matrix)


Beta_posterior=function(X,Y,Errvar,Psi,Tau){
  Beta <- matrix(0,J,K)
  
  smooth_invA <- function(ori_A){
    eigen_A <- eigen(ori_A)
    eigenvalue <- eigen_A$values
    eigenvec <- eigen_A$vectors
    ## 一般导致cov非正定的 并非因为eigen=0 而是因为下面的round 但round是为了确保symmteric
    #eigenvalue[which(eigenvalue==0)] <- 1e-160
    ## 用最简单的 abs(eigen) 但容易出现很小的eigen 被系统判定为0
    #newG_j <- eigen_G$vectors %*% diag(abs(eigenvalue)) %*% t(eigen_G$vectors)
    
    ## 如果用coth的方法 还能方便对很小的值进行soft regularize
    ## 
    smooth_inveigen <- diag(1 / (eigenvalue * sapply(eigenvalue,FUN=coth)))
    newinvA <- eigenvec %*% smooth_inveigen %*% t(eigenvec)
    return(newinvA)
  }
  
  for(k in 1:K){## solve(diag(Psi) %*% diag(Tau[,k]))
    A <- Errvar[k] * (diag(1/Psi) %*% diag(1/Tau[,k])) + (t(X[[k]]) %*% X[[k]])
    #invA <- solve(A) ## A如果数值过大 则solve(A)会很小 进而有些数无法被读取
    # 但若令数值过小 则优惠出现不对称的情况 甚至有些地方就是1e-12
    #    digi <-5 # digit怎么选？ 
    ## 必须让cov matrix对称 由于solve的计算问题 故会导致not symmetric错误
    ## 且也会出现eigen()有负值情况 是否需要warning and stop????
    #    while(sum(as.matrix(forceSymmetric(Errvar[k] * invA))%>% round(digits = digi) %>% diag()==0)>0 |
    #          !matrixcalc::is.positive.definite( as.matrix(forceSymmetric(Errvar[k] * invA))%>% round(digits = digi),tol = 1e-20)
    #          ){ # 需要cov是positive definite 且 diag上必须非0
    #      digi <- digi + 1
    
    #    }
    
    invA <- smooth_invA(A)
    #Cov_beta <- as.matrix(forceSymmetric(Errvar[k] * invA))%>% round(digits = digi) #故而只能在最终结果上限制小数点位数, 还需要适应性的改变digit
    Cov_beta <- Errvar[k] * (invA %>% forceSymmetric() %>% as.matrix())
    
    mu_beta <- invA %*% t(X[[k]]) %*% Y[[k]] %>% c()
    Beta[,k] <- rmvn(n=1, mu=mu_beta, Sigma=Cov_beta)
  }
  return(Beta)
}

Psi_posterior=function(Tau, Beta,a_psi=1,b_psi=1){
  Psi <- c()
  
  for(j in 1:J){
    v_1 = K+a_psi
    if(K == 1){
      eta_1 = t(Beta[j,]) %*% ((1 / Tau[j,])) %*% Beta[j,] + b_psi
    }else{
      eta_1 = t(Beta[j,]) %*% (diag(1 / Tau[j,])) %*% Beta[j,] + b_psi
    }
    
    Psi[j] <- rinvgamma(1,v_1/2, eta_1/2)
  }
  return(Psi)
}

Tau_posteriorGDP=function(Lambda, Beta, Psi){
  rownames(Lambda) <- NULL
  Tau <- matrix(0,nrow = J,ncol = K)
  for(j in 1:J){
    for(k in 1:K){
      mu_tau <- sqrt((Lambda[j,k]^2 * Psi[j]) / (Beta[j,k])^2 )
      shape_tau <- Lambda[j,k]^2
      ## tau很不稳定 而且存在rinvgaussian sample出0的情况
      invg <- rinvgaussian(n=1,mu_tau,shape_tau)
      if(invg == 0){
        invg <- 1e-320
      }
      Tau[j,k] <- 1 / invg
    }
  }
  return(Tau)
}

Sigma_PX_posterior_GDP<-function(Lambda,Sigma,shape,rate = eta,lowbar = 1e-07){
  Z <- matrix(0,J,K)
  for(j in 1:J){
    for(k in 1:K){
      if(pgamma(Lambda[j,k],shape[j,k],rate)==1){
        pp <- round(1 - 1e-16,digits = 16)
      }else{
        pp <- pgamma(Lambda[j,k],shape[j,k],rate)
      }
      Z[j,k] <- qnorm(pp,mean=0,sd=1)    
    }
  }
  
  D <- matrix(0,K,K)
  invSigma <- solve(Sigma)
  for(k in 1:K){
    D[k,k] <- sqrt(rinvgamma(n=1,shape=(K+1)/2, scale=diag(invSigma)[k]/2))
  }
  
  W <- Z %*% D
  
  ## 当nu < dim(S) 即会出现degenerate的情况 故此用generalized inverse Wishart来解决
  if(1+J <= K){
    Omega <- CholWishart::rGenInvWishart(n=1,df = 2+J,diag(1,K)+t(W)%*%W)[,,1]
    Omega[lower.tri(Omega)] = t(Omega)[lower.tri(Omega)]
    
    if(!Omega %>% is.positive.definite()){
      Omega_eigen <- eigen(Omega)
      Oeigenvalue <- Omega_eigen$values
      Oeigenvec <- Omega_eigen$vectors
      if(sum(Oeigenvalue<=1e-08)!=0){
        Oeigenvalue[which(Oeigenvalue<=1e-08)] <- lowbar
      }
      Omega <- Oeigenvec %*% diag(abs(Oeigenvalue)) %*% t(Oeigenvec)
      Omega[lower.tri(Omega)] = t(Omega)[lower.tri(Omega)]
    }
    
  }else{
    Omega  <- rinvwishart(nu=2+J, S=diag(1,K)+t(W)%*%W)
  }
  
  # same as newD<-diag(sqrt(diag(Omega)),K);Sigma <- solve(newD) %*% Omega %*% solve(newD)
  Sigma<-psych::cor.smooth(cov2cor(Omega)) %>% round(digits = 8)
  # 在极端状态下会出现non-diag上的+ -1 需要被替换
  
  Sigma[which(matrixcalc::upper.triangle(Sigma) - diag(diag(Sigma),K) == 1, arr.ind = T)] <- 1 - 1e-07
  Sigma[which(matrixcalc::upper.triangle(Sigma) - diag(diag(Sigma),K) == -1, arr.ind = T)] <- -(1 - 1e-07)
  
  return(matrixcalc::upper.triangle(Sigma)+t(matrixcalc::upper.triangle(Sigma))-diag(diag(Sigma),K))
}


#########################

coth <- function(x,a = 1e+06){
  ### 当x是负数时，这个函数 = -1 相当于取绝对值
  ### x是正数时 则返回1
  ### 但如果a*x的数值太大，就会报NaN, 
  ### 所以当|x| > 1/a时 直接返回+-1
  if(abs(x) >= 1/a){
    return(sign(x))
  }
  
  if(abs(x) < 1/a){
    ### 但对于close to 0 (0 <x< 1/a)的数 他会给一个正则项 使得x * coth(ax) >= 1/a 
    return((exp(a*x) + exp(-a*x)) / (exp(a*x) - exp(-a*x)))
  }
}

derivative_G_j<- function(lambda,shape,rate = eta){
  
  derivative_q <- function(lambda,z,shape,rate = eta){
    gradvec_q <- c()
    for(k in 1:K){
      p1 <- (1/dnorm(z[k],0,1))^3 * dgamma(lambda[k],shape[k],rate)^3
      p2 <- 2 * z[k]^2 * (1/dnorm(z[k],0,1))^3 * dgamma(lambda[k],shape[k],rate)^3
      p3 <- z[k] * (1/dnorm(z[k],0,1))^2 * (rate^shape[k] / gamma(shape[k]))^2 * (lambda[k]^(2*shape[k]-3) * exp(-2*rate*lambda[k]) *
                                                                                    (2*shape[k] - 2 + lambda[k] * (-2*rate)))
      p4 <- (1/dnorm(z[k],0,1))^2 * z[k] * dgamma(lambda[k],shape[k],rate) * 
        ((rate^shape[k] / gamma(shape[k])) * (lambda[k]^(shape[k]-2)) * exp(-rate*lambda[k])*
           (shape[k]-1+lambda[k]*(-rate)))
      p5 <- (1/dnorm(z[k],0,1)) * (rate^shape[k] / gamma(shape[k])) * exp(-rate*lambda[k]) * (
        (shape[k]-1)*(shape[k]-2)*lambda[k]^(shape[k]-3)
        + (shape[k]-1)*lambda[k]^(shape[k]-2)*(-rate)
        + (shape[k]-1)*lambda[k]^(shape[k]-2)*(-rate)
        + (lambda[k])^(shape[k]-1)*rate^2 
      )
      gradvec_q[k] <- p1 + p2 + p3 + p4 + p5  
    }
    return(gradvec_q)
  }
  
  z <- lam_to_z(lambda,shape)
  Q_t <- Qt_lam_GDP(lambda,z,shape)
  q <- c()
  for(k in 1:K){
    q[k] <- (z[k]/(dnorm(z[k],0,1))^2) * (dgamma(lambda[k],shape[k],rate))^2 + 
      (1/dnorm(z[k],0,1)) * ((rate^shape[k] / gamma(shape[k])) * (lambda[k]^(shape[k]-2)) * exp(-rate*lambda[k]) *
                               (shape[k]-1+lambda[k]*(-rate)))
  }
  
  der_q <- derivative_q(lambda,z,shape)
  
  target <- c()
  
  for(k in 1:K){
    p1 <- (t(z) %*% (diag(1,K)-solve(Sigma)))[k] * der_q[k] %>% c()
    p2 <- q[k] * ((diag(1,K)-solve(Sigma)) %*% c(Q_t[k,]))[k] %>% c()
    p3 <- 2 * (Q_t[k,] %*% (diag(1,K)-solve(Sigma)))[k] * q[k]
    p4 <- 2 * shape[k] * lambda[k]^(-3)
    target[k] <- p1 + p2 + p3 + p4
  }
  target <- -target %>% diag()
  
  for(k in 1:K-1){
    for(l in (k+1):K){
      target[l,k] <- -(Q_t[l,] %*% (diag(1,K)-solve(Sigma)))[k] * q[k]
    }
  }
  return(t(target) + target - diag(target)%>%diag())
}

### when G = {-H}, 即softabs 作用在negative Hessian上
Delta_H_smooth_lambda <- function(lambda,beta,Sigma,psi,shape,rate = eta,p,a = 10^6){
  
  oriG_j <- -Hessian_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)
  eigen_G_j <- eigen(oriG_j)
  
  Q_j <- eigen_G_j$vectors
  ori_eigen <- eigen_G_j$values
  R_j <- diag( 1 / (ori_eigen * sapply(ori_eigen,FUN=coth)))
  ### p要注意得是vector
  D_j <- diag( c(t(Q_j)%*%c(p)) / c(ori_eigen * sapply(ori_eigen,FUN=coth)) )
  
  J_j <- matrix(0,K,K)
  for(i in 1:(K-1)){
    for(l in (i+1):K){
      J_j[i,l] <- (ori_eigen[i]*coth(ori_eigen[i],a) - ori_eigen[l]*coth(ori_eigen[l],a)) / 
        (ori_eigen[i] - ori_eigen[l])
    }
  }
  ### const / inf = 0, 但得确保分子得是不太大的(通常来讲也不会太大) 否则就是inf / inf -> NaN
  J_j<-J_j + t(J_j) + diag(sapply(ori_eigen,FUN=coth) + ((-4 * ori_eigen) / (exp(a*ori_eigen)-exp(-a*ori_eigen))^2))
  
  M_j1 <- Q_j %*% (R_j * J_j) %*% t(Q_j)
  M_j2 <- Q_j %*% D_j %*% J_j %*% D_j %*% t(Q_j)
  ### derivative of G 只与lambda 和 alpha有关，但是使得p过大的原因
  ### 主要在part2的M_j2, 且主要体现在D_j上（因为M_j2 和 M_j1只有D_j不共享）
  ### 可是当p不大的时候 D_j也会大过头
  der_G_j <- derivative_G_j(lambda,shape)
  part1 <- c()
  part2 <- c()
  for(k in 1:K){
    
    der_G_j_k <- diag(0,K)
    der_G_j_k[k,] <- der_G_j[k,]
    der_G_j_k[,k] <- der_G_j[,k]
    
    part1[k] <- (1/2) * tr(M_j1 %*% der_G_j_k)
    part2[k] <- -(1/2) * tr(M_j2 %*% der_G_j_k)
  }
  part3 <- -Gradient_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta) %>% c()
  return(part1 + part2 + part3)
}


### 如果基于delta {H}的式子 写出非Tr的形式
Delta_H_smooth_lambda <- function(lambda,beta,Sigma,psi,shape,rate = eta,p,a = 10^6){
  
  oriG_j <- -Hessian_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)
  eigen_G_j <- eigen(oriG_j)
  
  Q_j <- eigen_G_j$vectors
  ori_eigen <- eigen_G_j$values
  R_j <- diag( 1 / (ori_eigen * sapply(ori_eigen,FUN=coth)))
  ### p要注意得是vector
  D_j <- diag( c(t(Q_j)%*%c(p)) / c(ori_eigen * sapply(ori_eigen,FUN=coth)) )
  
  J_j <- matrix(0,K,K)
  for(i in 1:(K-1)){
    for(l in (i+1):K){
      J_j[i,l] <- (ori_eigen[i]*coth(ori_eigen[i],a) - ori_eigen[l]*coth(ori_eigen[l],a)) / 
        (ori_eigen[i] - ori_eigen[l])
    }
  }
  ### const / inf = 0, 但得确保分子得是不太大的(通常来讲也不会太大) 否则就是inf / inf -> NaN
  J_j<-J_j + t(J_j) + diag(sapply(ori_eigen,FUN=coth) + ((-4 * ori_eigen) / (exp(a*ori_eigen)-exp(-a*ori_eigen))^2))
  
  ### derivative of G 只与lambda 和 alpha有关，但是使得p过大的原因
  ### 主要在part2的M_j2, 且主要体现在D_j上（因为M_j2 和 M_j1只有D_j不共享）
  ### 可是当p不大的时候 D_j也会大过头
  der_G_j <- derivative_G_j(lambda,shape)
  part1 <- c()
  part2 <- c()
  
  for(k in 1:K){
    
    der_G_j_k <- diag(0,K)
    der_G_j_k[k,] <- der_G_j[k,]
    der_G_j_k[,k] <- der_G_j[,k]
    
    shared_k <- J_j * (t(Q_j) %*% der_G_j_k %*% Q_j)
    ### 这和上面的 Tr()形式是一样的
    part1[k] <- (1/2) * tr(Q_j %*% R_j %*% shared_k %*% t(Q_j))
    part2[k] <- -(1/2) * c(t(p) %*% Q_j %*% R_j %*% shared_k %*% R_j %*% t(Q_j) %*% p) ## 但这和Tr形式不同
  }
  part3 <- -Gradient_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta) %>% c()
  return(part1 + part2 + part3)
}


### delta H_p
Delta_H_smooth_p <- function(lambda,beta,Sigma,psi,shape,rate = eta,p){
  oriG_j <- -Hessian_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)
  eigen_G_j <- eigen(oriG_j)
  
  Q_j <- eigen_G_j$vectors
  ori_eigen <- eigen_G_j$values
  R_j <- diag( 1 / (ori_eigen * sapply(ori_eigen,FUN=coth)))
  
  return( c((Q_j %*% R_j %*% t(Q_j)) %*% c(p)) ) 
}


H_smooth_func <- function(lambda,beta,Sigma,psi,shape,rate = eta,p){
  oriG_j <- -Hessian_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)
  eigen_G_j <- eigen(oriG_j)
  
  Q_j <- eigen_G_j$vectors
  ori_eigen <- eigen_G_j$values
  R_j <- diag( 1 / (ori_eigen * sapply(ori_eigen,FUN=coth)))
  smooth_eigen <- diag((ori_eigen * sapply(ori_eigen,FUN=coth)))
  newG_j <- Q_j %*% smooth_eigen %*% t(Q_j)
  newG_j_inv <- Q_j %*% R_j %*% t(Q_j)
  
  p1 <- -log(Lam_posterior_GDP(lambda,beta,Sigma,psi,shape,rate = eta))
  p2 <- (1/2) * log(det(newG_j)) # ignore the const
  p3 <- c((1/2) * t(p) %*% newG_j_inv %*% p)
  return(p1+p2+p3)
}


lam_to_z <- function(lambda,shape,rate=eta){ # transform \lambda_j to z_j
  z <- c()
  for(k in 1:K){
    if(pgamma(lambda[k],shape[k],rate)==1){ 
      p <- round(1 - 1e-16,digits = 16)}else{
        p <- pgamma(lambda[k],shape[k],rate)
      }
    z[k] <- qnorm(p,0,1)
  }
  return(z)
}

Lam_likelihood_GDP <- function(lambda,Sigma,shape,rate=eta){
  z<- lam_to_z(lambda,shape)
  expinner <- t(z) %*% (diag(rep(1,K))-solve(Sigma)) %*% z / 2
  GClikelihood <- pow(det(Sigma),-1/2) * exp(expinner)%>%c()
  lambdalikelihood<-c()
  for(k in 1:K){
    lambdalikelihood[k] <- dgamma(lambda[k],shape[k],rate)
  }
  
  return(GClikelihood*prod(lambdalikelihood))
}

Lam_posterior_GDP <- function(lambda,beta,Sigma,psi,shape,rate=eta){
  lam_likeli <- Lam_likelihood_GDP(lambda,Sigma,shape)
  beta_likevector <- c()
  for(k in 1:K){
    beta_likevector[k] <- LaplacesDemon::dlaplace( beta[k],location=0, scale = 1 / (lambda[k]/sqrt(psi)) )
  }
  return(lam_likeli * prod(beta_likevector))
}

Qt_lam_GDP <- function(lambda,z,shape,rate = eta){
  Q_t <- matrix(0,K,K) # Q^\top
  for(k in 1:K){
    Q_t[k,k] <- (1/dnorm(z[k],0,1)) * dgamma(lambda[k],shape[k],rate)
  }
  return(Q_t)
}

## Gradient is obtained under log posterior !!!

Gradient_lambda_GDP <- function(lambda,beta,Sigma,psi,shape,rate=eta){
  z <- lam_to_z(lambda,shape)
  Q_t <- Qt_lam_GDP(lambda,z,shape)
  part1 <- Q_t %*% (diag(1,K)-solve(Sigma)) %*% z %>% c()
  part2 <- (shape-1)/lambda-rate
  part3 <- 1/lambda - abs(beta) / sqrt(psi)
  return(part1+part2+part3)
}


Hessian_lambda_GDP <- function(lambda,beta,Sigma,psi,shape, rate = eta){
  z <- lam_to_z(lambda,shape)
  Q_t <- Qt_lam_GDP(lambda,z,shape)
  q <- c()
  for(k in 1:K){
    q[k] <- (z[k]/(dnorm(z[k],0,1))^2)*(dgamma(lambda[k],shape[k],rate))^2 + 
      (1/dnorm(z[k],0,1)) * ((rate^shape[k] / gamma(shape[k])) * (lambda[k]^(shape[k]-2)) * exp(-rate*lambda[k]) *
                               (shape[k]-1+lambda[k]*(-rate)))
  }
  part1 <- t(Q_t %*% (diag(1,K)-solve(Sigma)) %*% t(Q_t)) +
    diag(q) %*% diag( (diag(1,K)-solve(Sigma)) %*% z %>% c() )
  # isSymmetric(t(Q_t %*% (diag(1,K)-solve(Sigma))%*% t(Q_t))%>%round(digits = 6)): True
  part2 <- diag((-shape[k])/(lambda)^2)
  return(part1 + part2)
}

#d<-3
#-Hessian_lambda_GDP(Lambda[d,],Beta[d,],Sigma=Sigma,Psi[d],alpha[d,])%>%is.positive.definite()
#-Hessian_lambda_GDP(Lambda[d,],Beta[d,],Sigma=Sigma_block3,Psi[d],alpha[d,])%>%is.positive.definite()

Lam_sample_GDP_joint <- function(Sigma, shape, rate = eta){
  z <- rmvn(n=1, mu = rep(0,K), Sigma)%>% c()
  lambda <- c()
  for(k in 1:K){
    if(pnorm(z[k],0,1)==1){
      qq <- round(1 - 1e-16,digits = 16)}else{
        qq <- pnorm(z[k],0,1)
      }
    lambda[k] <- qgamma(qq,shape = shape[k],rate = rate)
  }
  return(lambda)
}

Lam_sample_GDP <- function(shape, rate = eta){
  lambda <- c()
  for(k in 1:K){
    lambda[k] <- rgamma(1,shape = shape[k],rate = rate)
  }
  return(lambda)
}


adaptive_simplified_M_MALA_lambda_GDP=function(beta,Sigma,psi,shape,rate = eta,
         epsilon = K^(-1/3)/10, L = 1000, ini_lambda,
         loweps = K^(-1/3)/50, uppereps = 1,
         opt_rate = 0.7, const_0 = 10, epsconst = 1){
  
  smooth_invG <- function(lambda,beta,Sigma,psi,shape,rate = eta){
    ori_G <- -Hessian_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)
    eigen_G <- eigen(ori_G)
    eigenvalue <- eigen_G$values
    eigenvec <- eigen_G$vectors
    ## 一般导致cov非正定的 并非因为eigen=0 而是因为下面的round 但round是为了确保symmteric
    #eigenvalue[which(eigenvalue==0)] <- 1e-160
    ## 用最简单的 abs(eigen) 但容易出现很小的eigen 被系统判定为0
    #newG_j <- eigen_G$vectors %*% diag(abs(eigenvalue)) %*% t(eigen_G$vectors)
    
    ## 如果用coth的方法 还能方便对很小的值进行soft regularize
    ## 
    smooth_inveigen <- diag(1 / (eigenvalue * sapply(eigenvalue,FUN=coth)))
    newinvG <- eigenvec %*% smooth_inveigen %*% t(eigenvec)
    return(newinvG)
  }
  
  truncation_eps <- function(eps,low = loweps, upper = uppereps){
    return(min(max(eps,low),upper))
  }
  
  proposal_para <- function(beta,Sigma,psi,shape,rate = eta,
                            epsilon,lambda){
    
    newinvG <- smooth_invG(lambda,beta,Sigma,psi,shape,rate = eta)
    proposal_mean = lambda + 
      (epsilon^2/2) * ((newinvG) %*% Gradient_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta))
    
    ## 如果newonvG本来就很小 结果因为eps^2，反而shrink to <1e-08
    proposal_cov = (epsilon^2 * newinvG)  ## round之后 会有很小的eigenvalue -> 0
    #proposal_cov %>% is.positive.definite()
    ## force the cov to be symmteric
    proposal_cov[lower.tri(proposal_cov)] = t(proposal_cov)[lower.tri(proposal_cov)]
    
    ## 会存在diag小于1e-08 会被系统判定为0 进而在dmvn中 整个cov是非正定的 （但是半正定）
    ## 所以可以 (<1e-08) := 1.001e-08
    if(sum(diag(proposal_cov)<1e-08)!=0){
      ind <- which(diag(proposal_cov)<1e-08)
      proposal_cov[ind,ind] <- 1.001e-08 ## 由于系统默认的 正定的判定下限是 1e-08, 故需稍大一点的数字
    }
    return(list(c(proposal_mean),proposal_cov))
  }
  
  #  eigen(proposal_cov)
  #  proposal_cov%>%is.positive.definite()
  ## 很小的值 都会被误判呈0
  #  proposal_cov%>%is.positive.semi.definite()
  #epsilon = 0.001
  old_lambda <- ini_lambda
  acc_num <- 0
  eps_set <- c()
  #tt <- 1
  for(i in 1:L){
    ### 若基于此old_lambda, 以及Sigma生成的new mean 是负的
    ### 若跳过则会一致生成无效的new mean, Sigma 使得old_lambda无法被更新
    ### 不如随机生成rgamma(5) 或者干脆将其取abs
    eps_set[i] <- epsilon
    ## 如果eps过小 则eps^2 * invG的元素就会过小 即<1e-08, 所以要防止eps过小，truncation的下界要守住
    new_para <- proposal_para(beta,Sigma,psi,shape,rate = eta,
                              epsilon = epsilon,old_lambda)
    ### 会有些落入到negative supp中 
    if(!new_para[[2]]%>%is.positive.definite()){
      #tt <- tt + 1
      next
    }
    new_lambda <- rmvn(n=1,mu=new_para[[1]],Sigma=new_para[[2]])%>%c()
    
    if(sum(new_lambda<=0)!= 0){
      next
    }
    
    new_posterior <- Lam_posterior_GDP(new_lambda,beta,Sigma,psi,shape)
    
    new_prop_density <- dmvn(new_lambda, mu=new_para[[1]], Sigma=new_para[[2]])
    
    old_posterior <- Lam_posterior_GDP(old_lambda,beta,Sigma,psi,shape)
    
    old_para <- proposal_para(beta,Sigma,psi,shape,rate = eta,
                              epsilon = epsilon,new_lambda)
    ## 会出现eigen有< 1e - 08的情况
    ## 虽然上面程序已经赋值给 < 1e-08的eigen 以 1.001e-08
    ## 但还是免不了报错
    if(!old_para[[2]]%>%is.positive.definite()){
      #tt <- tt + 1
      next
    }
    old_prop_density <- dmvn(old_lambda, mu=old_para[[1]], Sigma=old_para[[2]])
    ### proposal density太抢戏了 posterior density普遍都很低，但决定acc的 居然是proposal density
    #accrate = (new_posterior ) / (old_posterior )
    #accrate = (new_posterior * old_prop_density) / (old_posterior * new_prop_density)
    log_accrate = log(new_posterior)+log(old_prop_density) - log(old_posterior) - log(new_prop_density)
    if(log_accrate%>%is.na()){
      log_accrate <- log(10^(-7))
    }
    accrate = exp(log_accrate)
    
    if(i %% epsconst == 0){
      epsilon =  truncation_eps(eps = epsilon + (const_0/i)*(min(accrate,1) - opt_rate))
    }
    
    if(log_accrate>0 | log_accrate > log(runif(1,0,1))){
      old_lambda <- new_lambda
      
      acc_num <- acc_num + 1
    }else{old_lambda <- old_lambda}
    
    
    ## 也有建议说是 log(eps^2) = log(eps^2) + const/i * [a - opt_a]
    #if(i %% epsconst == 0){
    #  neweps =  truncation_eps(eps = 2*log(epsilon) + (const_0/i)*(min(accrate,1) - opt_rate))
    #  epsilon = exp(neweps/2)
    #  }
    
    #if(accrate >= runif(1,0,1)){
    
    #  old_lambda <- new_lambda
    
    #  acc_num <- acc_num + 1
    
    #}else{old_lambda <- old_lambda}
  }
  #tt
  return(list(old_lambda, acc_num / L, mean(eps_set)))
}


adaptive_simplified_M_MALA_Lambda_GDP=function(Beta,Sigma,Psi,shape,rate = eta,#grad_U = -Gradient_lambda, U=-Lam_posterior,
         epsilon=K^(-1/3)/10, L=300,
         loweps = K^(-1/3)/20, uppereps = 1,
         opt_rate = 0.7, const_0 = 10,
         epsconst = 5,inisd=0.25){
  
  Lam_sample_GDP_joint = function(Sigma, shape, rate = eta){
    z <- rmvn(n=1, mu = rep(0,K), Sigma)%>% c()
    lambda <- c()
    for(k in 1:K){
      if(pnorm(z[k],0,1)==1){
        qq <- round(1 - 1e-16,digits = 16)}else{
          qq <- pnorm(z[k],0,1)
        }
      lambda[k] <- qgamma(qq,shape = shape[k],rate = rate)
    }
    return(lambda)
  }
  
  Lambda <- matrix(0,nrow = J,ncol = K)
  Acc_prob <- c()
  eps_set <- c()
  for(j in 1:J){
    ini_lambda <- shape[j,] / rate + rnorm(K,0,inisd)
    #ini_lambda <- Lam_sample_GDP_joint(Sigma = Sigma,shape[j,])
    #ini_lambda <- (shape[j,] / rate )
    
    #ini_lambda <- (shape[j,] / rate ) + rnorm(K,0,0.2)
    #ini_lambda <- rgamma(K,shape=shape, rate = rate)
    samplings <- adaptive_simplified_M_MALA_lambda_GDP (beta=Beta[j,],Sigma=Sigma,psi=Psi[j],shape=shape[j,],
                                                        epsilon = epsilon,
                                                        L=L, ini_lambda = ini_lambda,
                                                        opt_rate = opt_rate, const_0 = const_0,epsconst = epsconst)
    Lambda[j,] <- samplings[[1]]##本来lambda的prior就应该是depends on Sigma的
    Acc_prob[j] <- samplings[[2]]
    eps_set[j] <- samplings[[3]]
    
  }
  return(list(Lambda,Acc_prob,eps_set))
}

adaptive_full_M_MALA_lambda_GDP <- function(beta,Sigma,psi,shape,rate = eta,
                                            epsilon = 0.01, L = 100, ini_lambda, a = 10^6,
                                            loweps = 0.005, uppereps = 1,
                                            opt_rate = 0.7, const_0 = 1, epsconst = 1){
  
  proposal_para <- function(beta,Sigma,psi,shape,rate = eta,
                            epsilon,lambda){
    
    oriG_j <- -Hessian_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)
    eigen_G_j <- eigen(oriG_j)
    
    Q_j <- eigen_G_j$vectors
    ori_eigen <- eigen_G_j$values
    R_j <- diag( 1 / (ori_eigen * sapply(ori_eigen,FUN=coth)))
    
    new_invG <- Q_j %*% R_j %*% t(Q_j)
    
    J_j <- matrix(0,K,K)
    for(i in 1:(K-1)){
      for(l in (i+1):K){
        ## 若有在off-diag的位置 出现\eigen_j = \eigen_k
        if(ori_eigen[i] ==  ori_eigen[l]){
          J_j[i,l] <- coth(ori_eigen[i],a) + ((-4 * ori_eigen[i]) / (exp(a*ori_eigen[i])-exp(-a*ori_eigen[i]))^2)
        }else{
          J_j[i,l] <- (ori_eigen[i]*coth(ori_eigen[i],a) - ori_eigen[l]*coth(ori_eigen[l],a)) / 
            (ori_eigen[i] - ori_eigen[l])
        }
      }
    }
    ### const / inf = 0, 但得确保分子得是不太大的(通常来讲也不会太大) 否则就是inf / inf -> NaN
    J_j <- J_j + t(J_j) + diag(sapply(ori_eigen,FUN=coth) + ((-4 * ori_eigen) / (exp(a*ori_eigen)-exp(-a*ori_eigen))^2))
    ### 由于oriG是neg Hess, 且derivative也是 delta (-H)
    der_neg_H_j <- derivative_G_j(lambda,shape)
    
    der_omega_j <- rep(0,K)
    for(k in 1:K){
      
      der_neg_H_j_k <- diag(0,K)
      der_neg_H_j_k[k,] <- der_neg_H_j[k,]
      der_neg_H_j_k[,k] <- der_neg_H_j[,k]
      
      #### delta {G} = delta {-H}
      der_G_j_k <- Q_j %*% (J_j * (t(Q_j) %*% der_neg_H_j_k %*% Q_j)) %*% t(Q_j)
      der_omega_j <- der_omega_j + (new_invG %*% der_G_j_k %*% new_invG)[,k]
    }
    der_omega_j <- -der_omega_j
    ## 若epsilon过大 还会导致prop mean的元素为负
    
    proposal_mean = lambda + 
      (epsilon^2/2) * c(new_invG %*% Gradient_lambda_GDP(lambda,beta,Sigma,psi,shape,rate = eta)) +
      (epsilon^2/2) * der_omega_j
    
    proposal_cov = (epsilon^2 * new_invG)  ## round之后 会有很小的eigenvalue -> 0
    #proposal_cov %>% is.positive.definite()
    ## force the cov to be symmteric
    proposal_cov[lower.tri(proposal_cov)] = t(proposal_cov)[lower.tri(proposal_cov)]
    
    ## 会存在diag小于1e-08 会被系统判定为0 进而在dmvn中 整个cov是非正定的 （但是半正定）
    ## 所以可以 (<1e-08) := 1.001e-08
    if(sum(diag(proposal_cov)<1e-08)!=0){
      ind <- which(diag(proposal_cov)<1e-08)
      proposal_cov[ind,ind] <- 1.1e-08 ## 由于系统默认的 正定的判定下限是 1e-08, 故需稍大一点的数字
      ## 但系统检验positive df 是从eigen value的角度开始的
      ## 所以最保险的方法是 U^\top |D| U
    }
    return(list(c(proposal_mean),proposal_cov))
  }
  
  truncation_eps <- function(eps,low = loweps, upper = uppereps){
    return(min(max(eps,low),upper))
  }
  
  
  #  eigen(proposal_cov)
  #  proposal_cov%>%is.positive.definite()
  ## 很小的值 都会被误判呈0
  #  proposal_cov%>%is.positive.semi.definite()
  
  # new_para[[2]]%>%is.positive.semi.definite()
  # epsilon <- 0.0002
  old_lambda <- ini_lambda
  acc_num <- 0
  eps_set <- c()
  #tt <- 1
  for(i in 1:L){
    eps_set[i] <- epsilon
    
    new_para <- proposal_para(beta,Sigma,psi,shape,rate = eta,
                              epsilon = epsilon,old_lambda)
    # new_para[[1]]
    #(new_para[[2]]) %>% eigen()
    
    if(!new_para[[2]]%>%is.positive.definite()){
      #tt <- tt + 1
      #### instead of next
      ## 若只给<1e-08的diag位置赋以1.001e-08 照样会出现semi-definite
      next
      #new_para[[2]] <- diag(1e-06,K)
    }
    new_lambda <- rmvn(n=1,mu=new_para[[1]],Sigma=new_para[[2]])%>%c()
    
    if(sum(new_lambda<0)!=0){
      #### instead of next
      #new_lambda <- Lam_sample_GDP_joint(Sigma,shape)
      new_lambda <- shape / rate + rnorm(K,0,0.1)
      #new_lambda <- rgamma(K,shape,rate=rate)
    }
    
    new_posterior <- Lam_posterior_GDP(new_lambda,beta,Sigma,psi,shape)
    
    new_prop_density <- dmvn(new_lambda, mu=new_para[[1]], Sigma=new_para[[2]])
    
    old_posterior <- Lam_posterior_GDP(old_lambda,beta,Sigma,psi,shape)
    
    old_para <- proposal_para(beta,Sigma,psi,shape,rate = eta,
                              epsilon = epsilon,new_lambda)
    ## 会出现eigen有< 1e - 08的情况
    ## 虽然上面程序已经赋值给 < 1e-08的eigen 以 1.001e-08
    ## 但还是免不了报错
    if(!old_para[[2]]%>%is.positive.definite()){
      #tt <- tt + 1
      next
      #old_para[[2]] <- diag(1e-06,K)
    }
    old_prop_density <- dmvn(old_lambda, mu=old_para[[1]], Sigma=old_para[[2]])
    
    accrate = (new_posterior * old_prop_density) / (old_posterior * new_prop_density)
    if(accrate%>%is.na()){
      accrate <- 10^(-7)
    }
    
    if(i %% epsconst == 0){
      epsilon =  truncation_eps(eps = epsilon + (const_0/i)*(min(accrate,1) - opt_rate))
    }
    
    if(accrate >= runif(1,0,1)){
      old_lambda <- new_lambda
      acc_num <- acc_num + 1
    }else{old_lambda <- old_lambda}
  }
  #tt
  return(list(old_lambda,acc_num / L,mean(eps_set)))
}


adaptive_full_M_MALA_Lambda_GDP <- function(Beta,Sigma,Psi,shape,rate = eta,#grad_U = -Gradient_lambda, U=-Lam_posterior,
                                            epsilon=1e-01, L=600,a=10^6,
                                            loweps = K^(-1/3)/20, uppereps = 1,
                                            opt_rate = 0.8, const_0 = 10,
                                            epsconst = 1,inisd=0.1){
  Lambda <- matrix(0,nrow = J,ncol = K)
  Acc_prob <- c()
  eps_set <- c()
  trial_number = c()
  for(j in 1:J){
    #ini_lambda <- shape[j,] / rate + rnorm(K,0,inisd)
    ini_lambda <- Lam_sample_GDP_joint(Sigma,shape=shape[j,]) + rnorm(K,0,inisd)
    
    Lambda_and_Accprob <- adaptive_full_M_MALA_lambda_GDP(beta=Beta[j,],Sigma=Sigma,psi=Psi[j],shape=shape[j,],rate = eta,
                                                          epsilon = epsilon,
                                                          L=L, ini_lambda = ini_lambda,
                                                          opt_rate = opt_rate, const_0 = const_0,epsconst = epsconst)##本来lambda的prior就应该是depends on Sigma的
    
    trial = 1
    while(Lambda_and_Accprob[[2]] <= opt_rate / 2 & trial <= 10){
      ini_lambda <- Lam_sample_GDP_joint(Sigma,shape=shape[j,]) + rnorm(K,0,inisd)
      
      Lambda_and_Accprob <- adaptive_full_M_MALA_lambda_GDP(beta=Beta[j,],Sigma=Sigma,psi=Psi[j],shape=shape[j,],rate = eta,
                                                            epsilon = epsilon,
                                                            L=L, ini_lambda = ini_lambda,
                                                            opt_rate = opt_rate, const_0 = const_0,epsconst = epsconst)##本来lambda的prior就应该是depends on Sigma的
      trial = trial + 1
    }
    
    Lambda[j,] <- Lambda_and_Accprob[[1]]
    Acc_prob[j] <- Lambda_and_Accprob[[2]]
    eps_set[j] <- Lambda_and_Accprob[[3]]
    trial_number[j] <- trial
  }
  return(list(Lambda,Acc_prob,eps_set,trial_number))
}


Sigma_posterior_GDP <- function(Lambda,v_0 = K+1, V_0 = diag(1,K),shape,rate = eta){
  Z <- matrix(0,J,K)
  for(j in 1:J){
    for(k in 1:K){
      if(pgamma(Lambda[j,k],shape[j,k],rate)==1){
        pp <- round(1 - 1e-16,digits = 16)
      }else{
        pp <- pgamma(Lambda[j,k],shape[j,k],rate)
      }
      Z[j,k] <- qnorm(pp,mean=0,sd=1)    
    }
  }
  v_1 = v_0 + J
  ## or v_0 * V_0
  V_1 = v_0 * V_0 + t(Z) %*% Z
  V <- rinvwishart(nu=v_1, S=V_1) ## 为什么可以生成的这么不稳定
  Sigma <- psych::cor.smooth(cov2cor(V)) %>% round(digits = 8)
  return(Sigma)
}



Sigma_PX_posterior_GDP=function(Lambda,Sigma,shape,rate = eta,lowbar = 1e-07){
  Z <- matrix(0,J,K)
  for(j in 1:J){
    for(k in 1:K){
      if(pgamma(Lambda[j,k],shape[j,k],rate)==1){
        pp <- round(1 - 1e-16,digits = 16)
      }else{
        pp <- pgamma(Lambda[j,k],shape[j,k],rate)
      }
      Z[j,k] <- qnorm(pp,mean=0,sd=1)    
    }
  }
  
  D <- matrix(0,K,K)
  invSigma <- solve(Sigma)
  for(k in 1:K){
    D[k,k] <- sqrt(rinvgamma(n=1,shape=(K+1)/2, scale=diag(invSigma)[k]/2))
  }
  
  W <- Z %*% D
  
  ## 当nu < dim(S) 即会出现degenerate的情况 故此用generalized inverse Wishart来解决
  if(1+J <= K){
    Omega <- CholWishart::rGenInvWishart(n=1,df = 2+J,diag(1,K)+t(W)%*%W)[,,1]
    Omega[lower.tri(Omega)] = t(Omega)[lower.tri(Omega)]
    
    if(!Omega %>% is.positive.definite()){
      Omega_eigen <- eigen(Omega)
      Oeigenvalue <- Omega_eigen$values
      Oeigenvec <- Omega_eigen$vectors
      if(sum(Oeigenvalue<=1e-08)!=0){
        Oeigenvalue[which(Oeigenvalue<=1e-08)] <- lowbar
      }
      Omega <- Oeigenvec %*% diag(abs(Oeigenvalue)) %*% t(Oeigenvec)
      Omega[lower.tri(Omega)] = t(Omega)[lower.tri(Omega)]
    }
    
  }else{
    Omega  <- rinvwishart(nu=2+J, S=diag(1,K)+t(W)%*%W)
  }
  
  # same as newD<-diag(sqrt(diag(Omega)),K);Sigma <- solve(newD) %*% Omega %*% solve(newD)
  Sigma<-psych::cor.smooth(cov2cor(Omega)) %>% round(digits = 8)
  # 在极端状态下会出现non-diag上的+ -1 需要被替换
  
  Sigma[which(matrixcalc::upper.triangle(Sigma) - diag(diag(Sigma),K) == 1, arr.ind = T)] <- 1 - 1e-07
  Sigma[which(matrixcalc::upper.triangle(Sigma) - diag(diag(Sigma),K) == -1, arr.ind = T)] <- -(1 - 1e-07)
  
  return(matrixcalc::upper.triangle(Sigma)+t(matrixcalc::upper.triangle(Sigma))-diag(diag(Sigma),K))
}

Lambda_posteriornew=function(Beta,Psi,shape = alpha,
         rate = eta){
  Lambda <- matrix(0, nrow = J, ncol = K)
  for(j in 1:J){
    for(k in 1:K){
      Lambda[j,k] <- rgamma( n=1,shape + 1,
                             rate + abs(Beta[j,k])/sqrt(Psi[j]) )
    }
  }
  return(Lambda)
}

{
  Y_to_binary_probit_list <- function(Z){
    Y_pred_list <- list()
    for(k in 1:K){
      prob_k <- c()
      for(i in 1:length(Z[[k]])){
        prob_k[i] <- pnorm(Z[[k]][i],0,1)
      }
      Y_pred_list[[k]] <- rbern(n=length(prob_k),prob = prob_k) %>% as.matrix()
    }
    return(Y_pred_list )
  }
  
  ## Y_to_binary: change a continuous vector to 0-1 vector
  Y_to_binary_compare0_list <- function(X,Beta){
    
    Y_pred_list <- list()
    for(k in 1:K){
      Z = X[[k]] %*% Beta[,k]
      Z_latent = rnorm(length(Z), Z, 1)
      Y_pred_list[[k]] = Y_to_binary(Z_latent)
    }
    return(Y_pred_list)
  }
  
  Y_to_binary_logistic_list <- function(X,Beta){
    
    Y_pred_list <- list()
    for(k in 1:K){
      Y_pred_list[[k]] <- Y_to_binary_logistic(X[[k]] %*% Beta[,k])
    }
    return(Y_pred_list)
  }
  
  Loss_binary_func <- function(trueY_list, predY_list){
    loss_vec <- c()
    for(k in 1:K){
      loss_vec[k] <- sum(predY_list[[k]] != trueY_list[[k]]) / length(predY_list[[k]])
    }
    return(loss_vec)
  }
  
  ###
  ## log_loss = \sum log(1 + exp(- y * (x\top \beta)))
  
  Logistic_loss<- function(Y,X,Beta){
    loss_vec <- c()
    
    log_func <- function(x){
      return(log(1+exp(-x)))
    }
    
    for(k in 1:K){
      pred_y_k <- X[[k]] %*% Beta[,k]
      true_y_k <- Y[[k]]
      loss_vec[k] = sum((pred_y_k * true_y_k) %>% log_func())
      
    }
    return(loss_vec)
  }
}


}
