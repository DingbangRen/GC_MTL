
library(dplyr)
{
#######

#sp80_1 <- c(rep(1,5),rep(0,5),rep(c(0,1,1,0,0,0,0),70/7))
## only a small subset is relevant
sp80_1 <- c(rep(1,5),rep(0,5),rep(1,10),rep(0,60))
sum(sp80_1==1)

sp80_1_1 = sp80_1 
sp80_1_2 = sp80_1
sp80_1_3 = sp80_1
## sp 2应该与sp 1有更多的重合
sp80_2 <- c(rep(1,5),rep(0,10),rep(1,10),rep(0,55))

##### 如果让task 5,6也成为majority呢？？
#sp80_2_1 = sp80_2 = sp80_1
sp80_3 <- c(rep(0,30),rep(1,10),rep(0,40))

sp80_4 <- c(rep(0,30),rep(0,40),rep(1,10))
}

{
## 讲\beta设置的更随机，例如shared_mu = 0, 且variance=2
set.seed(123)#
#shared_mu <- 5
#shared_mu <- 1
shared_mu = 0
L_80 = 80
shared_sigma <- 2

## another setting
#shared_mu = 1
#shared_sigma = 0.2

beta1_sim3 <- sp80_1  *rnorm(L_80,shared_mu,shared_sigma)
beta2_sim3 <- sp80_1_1*rnorm(L_80,shared_mu,shared_sigma)
beta3_sim3 <- sp80_1_2*rnorm(L_80,shared_mu,shared_sigma)
beta4_sim3 <- sp80_1_3*rnorm(L_80,shared_mu,shared_sigma)

## 把beta_4 设置为negative value
beta4_sim3 <- -sp80_1_3*rnorm(L_80,shared_mu,shared_sigma)

######
beta4_sim3=beta3_sim3=beta2_sim3=beta1_sim3
beta4_sim3 = - beta4_sim3

beta5_sim3 <- sp80_2  *rnorm(L_80,shared_mu,shared_sigma)
beta6_sim3 <- sp80_2 *rnorm(L_80,shared_mu,shared_sigma)

beta6_sim3 = beta5_sim3

beta7_sim3 <- sp80_3  *rnorm(L_80,shared_mu,shared_sigma)
beta8_sim3 <- sp80_4  *rnorm(L_80,shared_mu,shared_sigma)


Beta_alltasks_J_80 <- data.frame( 'task 1'=beta1_sim3,
                                  'task 2'=beta2_sim3,
                                  'task 3'=beta3_sim3,
                                  'task 4'=beta4_sim3,
                                  'task 5'=beta5_sim3,
                                  'task 6'=beta6_sim3,
                                  'task 7'=beta7_sim3,
                                  'task 8'=beta8_sim3)
Beta_alltasks_J_80%>%round(2)
}

######J = 80, K = 8, n =50
{J <- 80 ## 暂不考虑intercept
  K <- 8}
######
###################

{
library('dplyr')
library(LaplacesDemon)

#set.seed(15)
#(sample(8:18,K,replace = F) -> nsamples)

#set.seed(10)
#(sample(30:50,K,replace = F) -> nsamples)


rep(80,8) -> nsamples
## another setting 
#rep(100, 8) -> nsamples

nsamples[4] <- 20

z_score_X <- function(X){
  X_zscore <- matrix(0,nrow(X),ncol(X))
  for(i in 1:ncol(X)){
    X_zscore[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])
  }
  return(X_zscore)
}

Y_to_binary <- function(y){
  y_binary <- c()
  for(i in 1:length(y)){
    if(y[i] <= 0){
      y_binary[i] <- 0
    }else if(y[i] > 0){
      y_binary[i] <- 1
    }
  }
  return(y_binary %>% as.matrix())
}

## suggested by (y_i |\beta, x_i) \sim Bern(p_i), logit(p_i) = x^\top_i \beta
Y_to_binary_logistic <- function(y){
  
  probfunc <- function(z){
    return(1 / (1 + exp(-z)))
  }
  
  y_binary <- c()
  y_binary <- rbern(n=length(y), prob = probfunc(y))
  
  # for(i in 1:length(y)){
  #prob_i <- probfunc(y[i])
  #if(prob_i > 0.5){
  #  y_binary[i] <- 1
  #}else{
  #  y_binary[i] <- 0
  #}
  #}
  return(y_binary %>% as.matrix())
}

}

###################
{
z_score_X <- function(X){
  X_zscore <- matrix(0,nrow(X),ncol(X))
  for(i in 1:ncol(X)){
    X_zscore[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])
  }
  return(X_zscore)
}

Y_to_binary <- function(y){
  y_binary <- c()
  for(i in 1:length(y)){
    if(y[i] <= 0){
      y_binary[i] <- 0
    }else if(y[i] > 0){
      y_binary[i] <- 1
    }
  }
  return(y_binary %>% as.matrix())
}

Y_to_binary_logistic <- function(y){
  
  probfunc <- function(z){
    return(1 / (1 + exp(-z)))
  }
  
  y_binary <- c()
  y_binary <- rbern(n=length(y), prob = probfunc(y))
  
  # for(i in 1:length(y)){
  #prob_i <- probfunc(y[i])
  #if(prob_i > 0.5){
  #  y_binary[i] <- 1
  #}else{
  #  y_binary[i] <- 0
  #}
  #}
  return(y_binary %>% as.matrix())
}
}


###################
{
set.seed(1)
rmvn(J,rep(0,nsamples[1]), 10*diag(rep(1,nsamples[1]))) %>% as.matrix()%>%t() -> X1
X1 <- z_score_X(X1)
eps1 <- rmvn(1,rep(0,nsamples[1]), 1*diag(rep(1,nsamples[1])))%>%t()
colnames(eps1) <- NULL

X1 %*% beta1_sim3 + eps1-> Y1
Y_to_binary(X1 %*% beta1_sim3) -> Y1_binary
Y_to_binary_logistic(X1 %*% beta1_sim3) -> Y1_binary

#cbind(rep(1,nrow(X1)),X1) -> X1
#X1 %*% beta1_sim + eps1-> Y1
#X1 %*% beta1_sim2 + eps1-> Y1_new

#####
set.seed(2)
rmvn(J,rep(nsamples[2]), 10*diag(rep(1,nsamples[2]))) %>% as.matrix()%>%t() -> X2
X2 <- z_score_X(X2)
eps2 <- rmvn(1,rep(0,nsamples[2]), 1*diag(rep(1,nsamples[2])))%>%t()
colnames(eps2) <- NULL

X2 %*% beta2_sim3 + eps2-> Y2
Y_to_binary(X2 %*% beta2_sim3) -> Y2_binary
Y_to_binary_logistic(X2 %*% beta2_sim3) -> Y2_binary


#cbind(rep(1,nrow(X2)),X2) -> X2
#X2 %*% beta2_sim+ eps2 -> Y2
#X2 %*% beta2_sim2+ eps2 -> Y2_new

#####
set.seed(3)
rmvn(J,rep(0,nsamples[3]), 10*diag(rep(1,nsamples[3]))) %>% as.matrix()%>%t() -> X3
X3 <- z_score_X(X3)
eps3 <- rmvn(1,rep(0,nsamples[3]), 1*diag(rep(1,nsamples[3])))%>%t()
colnames(eps3) <- NULL


X3 %*% beta3_sim3 + eps3-> Y3
Y_to_binary(X3 %*% beta3_sim3) -> Y3_binary
Y_to_binary_logistic(X3 %*% beta3_sim3) -> Y3_binary


#cbind(rep(1,nrow(X3)),X3) -> X3
#X3 %*% beta3_sim +eps3  -> Y3
#X3 %*% beta3_sim2 +eps3  -> Y3_new

#####
set.seed(4)
rmvn(J,rep(0,nsamples[4]), 10*diag(rep(1,nsamples[4]))) %>% as.matrix()%>%t() -> X4
X4 <- z_score_X(X4)
eps4 <- rmvn(1,rep(0,nsamples[4]), 1*diag(rep(1,nsamples[4])))%>%t()
colnames(eps4) <- NULL


X4 %*% beta4_sim3 + eps4-> Y4
Y_to_binary(X4 %*% beta4_sim3) -> Y4_binary
Y_to_binary_logistic(X4 %*% beta4_sim3) -> Y4_binary

#cbind(rep(1,nrow(X4)),X4) -> X4
#X4 %*% beta4_sim+ eps4 -> Y4
#X4 %*% beta4_sim2+ eps4 -> Y4_new

#####
set.seed(5)
rmvn(J,rep(0,nsamples[5]), 10*diag(rep(1,nsamples[5]))) %>% as.matrix()%>%t() -> X5
X5 <- z_score_X(X5)
eps5 <- rmvn(1,rep(0,nsamples[5]), 1*diag(rep(1,nsamples[5])))%>%t()
colnames(eps5) <- NULL


X5 %*% beta5_sim3 + eps5-> Y5
Y_to_binary(X5 %*% beta5_sim3) -> Y5_binary
Y_to_binary_logistic(X5 %*% beta5_sim3) -> Y5_binary

#cbind(rep(1,nrow(X5)),X5) -> X5
#X5 %*% beta5_sim+eps5 -> Y5
#X5 %*% beta5_sim2+eps5 -> Y5_new

######
set.seed(6)
rmvn(J,rep(0,nsamples[6]), 10*diag(rep(1,nsamples[6]))) %>% as.matrix()%>%t() -> X6
X6 <- z_score_X(X6)
eps6 <- rmvn(1,rep(0,nsamples[6]), 1*diag(rep(1,nsamples[6])))%>%t()
colnames(eps6) <- NULL


X6 %*% beta6_sim3 + eps6-> Y6
Y_to_binary(X6 %*% beta6_sim3) -> Y6_binary
Y_to_binary_logistic(X6 %*% beta6_sim3) -> Y6_binary

#cbind(rep(1,nrow(X6)),X6) -> X6
#X6 %*% beta6_sim+eps6 -> Y6
#X6 %*% beta6_sim2+eps6 -> Y6_new

#####
set.seed(7)
rmvn(J,rep(0,nsamples[7]), 10*diag(rep(1,nsamples[7]))) %>% as.matrix()%>%t() -> X7
X7 <- z_score_X(X7)
eps7 <- rmvn(1,rep(0,nsamples[7]), 1*diag(rep(1,nsamples[7])))%>%t()
colnames(eps7) <- NULL


X7 %*% beta7_sim3 + eps7-> Y7
Y_to_binary(X7 %*% beta7_sim3) -> Y7_binary
Y_to_binary_logistic(X7 %*% beta7_sim3) -> Y7_binary

#cbind(rep(1,nrow(X7)),X7) -> X7
# X7 %*% beta7_sim+eps7 -> Y7
#X7 %*% beta7_sim2+eps7 -> Y7_new

#####
set.seed(8)
rmvn(J,rep(0,nsamples[8]), 10*diag(rep(1,nsamples[8]))) %>% as.matrix()%>%t() -> X8
X8 <- z_score_X(X8)
eps8 <- rmvn(1,rep(0,nsamples[8]), 1*diag(rep(1,nsamples[8])))%>%t()
colnames(eps8) <- NULL


X8 %*% beta8_sim3 + eps8-> Y8
Y_to_binary(X8 %*% beta8_sim3) -> Y8_binary
Y_to_binary_logistic(X8 %*% beta8_sim3) -> Y8_binary

#cbind(rep(1,nrow(X8)),X8) -> X8
# X7 %*% beta7_sim+eps7 -> Y7
#X8 %*% beta8_sim2+eps8 -> Y8_new
}

#################
###### J = 80, n = 50, K = 8
{
Xlist_sim <- list(  'X1' = X1,
                    'X2' = X2,
                    'X3' = X3,
                    'X4' = X4,
                    'X5' = X5,
                    'X6' = X6,
                    'X7' = X7,
                    'X8' = X8)
Ylist_sim <- list('Y1' = Y1,
                   'Y2' = Y2,
                   'Y3' = Y3,
                   'Y4' = Y4,
                   'Y5' = Y5,
                   'Y6' = Y6,
                   'Y7' = Y7,
                   'Y8' = Y8)

Ylist_sim_binary <- list('Y1' = Y1_binary,
                          'Y2' = Y2_binary,
                          'Y3' = Y3_binary,
                          'Y4' = Y4_binary,
                          'Y5' = Y5_binary,
                          'Y6' = Y6_binary,
                          'Y7' = Y7_binary,
                          'Y8' = Y8_binary)
}

{
  #J <- 80
  #J <- 30
  #J  = 40
  #K <- 8  
  
  X_train = list()
  Y_train = list()
  
  X_original = Xlist_sim
  Y_original = Ylist_sim
  
  #Y_original = Ylist_sim_binary
  
  #X_original = Xlist_J_50
  #Y_original = Ylist_J_50
  
  nsam_original <- nsamples[1]
  
  train_prop = 0.8
  
  set.seed(1)
  nsample_train <- floor(nsam_original*train_prop)
  sample_train <- sample(x=seq(1,nsam_original,1), size=nsample_train)
  
  for(k in 1:K){
    if(k == 4){
      next
    }
    X_train[[k]] <- X_original[[k]][sample_train,] %>% as.matrix()
    Y_train[[k]] <- Y_original[[k]][sample_train] %>% as.matrix()
  }
  
  set.seed(1)
  nsample_train_4 <- floor(nsamples[4]*train_prop)
  sample_train_4 <- sample(x=seq(1,nsamples[4],1), size=nsample_train_4)
  
  X_train[[4]] <- X_original[[4]][sample_train_4,] %>% as.matrix()
  Y_train[[4]] <- Y_original[[4]][sample_train_4] %>% as.matrix
  
  #######
  
  Y_test = list()
  X_test = list()
  
  sample_test <- seq(1,nsam_original,1)[!(seq(1,nsam_original,1) %in% sample_train)]
  
  ###
  
  for(k in 1:K){
    if(k == 4){
      next
    }
    X_test[[k]] <- X_original[[k]][sample_test,]%>% as.matrix()
    Y_test[[k]] <- Y_original[[k]][sample_test]%>% as.matrix()
  }
  
  sample_test_4 <- seq(1,nsamples[4],1)[!(seq(1,nsamples[4],1) %in% sample_train_4)]
  X_test[[4]] <- X_original[[4]][sample_test_4,]%>% as.matrix()
  Y_test[[4]] <- Y_original[[4]][sample_test_4]%>% as.matrix()
}



