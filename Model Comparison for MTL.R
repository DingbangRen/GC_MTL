####### Simulation
Xlist_test = X_test
Ylist_test = Y_test
simul_if = T

if(rmtl_type == 'Classification'){
  for(k in 1:K){
    Ylist_test[[k]][Ylist_test[[k]] == 0] = -1
  }
}

####### Sacros
#Xlist_test = Xlist_sarcos_test
#Ylist_test = Ylist_sarcos_test
####### School
Xlist_test = X_test
Ylist_test = Y_test
#### isolet
Xlist_test = X_test
Ylist_test = Y_test
########
{
Beta_ours_all_mean_sim = matrix(0,nrow = J, ncol = K)
Beta_linear_all_sim = matrix(0,nrow = J, ncol = K)

library(coda)

mcmcbeta<-function(n,b,task,intv = interval){
  
  vari <-  matrix(0,ncol = J, nrow=(n/intv)-b/intv+1)
  for(i in 1:(n/intv-b/intv+1)){
    vari[i,] <- resultlist[[i]]$Beta[,task]%>%t()
  } 
  #### how can we diagnose the matrix
  return(mcmc(vari))
}

for(d in 1:K){
  # 需要把原来的intercept项取出， 因为linear regression自带intercept考虑
  X_d <- Xlist_test[[d]][,-1] %>% data.frame()
  datalr_d <- cbind('Y_d' = Ylist_test[[d]],
                    X_d)
  
  #lm(Y_d ~. ,data = datalr_d) -> lm_simulation_task
  ## 使用posterior mean or median??
  Beta_mean_simulation <- apply(mcmcbeta(n=niter,b=burnin,task = d), MARGIN = 2, function(x){return(mean(x))})
  
  
  #betalinear_simulation <- lm_simulation_task$coefficients
  #betalinear_simulation[is.na(betalinear_simulation)] <- 0
  ####
  ##若满足条件  此处需要把intercept去掉
  #if(ncol(Xlist_test[[d]]) != length(betalinear_simulation)){
  #  betalinear_simulation <- betalinear_simulation[-1] %>% as.vector()
  #}
  
  #Beta_linear_all_sim[,d] = betalinear_simulation
  Beta_ours_all_mean_sim[,d] = Beta_mean_simulation
  
}
}
###############

nMSE <- function(X,Y,beta){
  estimate <- X %*% c(beta)
  nmse <- (Y-estimate)^2 %>% sum() / 
    (nrow(X) * var(Y) )
  return(nmse)
}

rMSE <- function(X,Y,beta){
  estimate <- X %*% c(beta)
  n = length(Y)
  rmse <- sqrt(( sum((Y-estimate)^2) ) / n)
  return(rmse)
}

#### to avoid the negative R square case, and we should decompose the SS total into two parts.
Rsquare = function(X,Y,beta){
  estimate <- X %*% c(beta)
  
  SS_explained = (estimate - mean(Y))^2 %>% sum()
  SS_res = (Y - estimate)^2 %>% sum()
  
  return(SS_explained / (SS_explained+SS_res))
}
###############
## simulation
rsquare_allmodels_sim = data.frame(L21 = 0,
                                      Lasso = 0,
                                      Trace = 0,
                                      Mean_regular = 0,
                                      #K_means = 0,
                                      DMFS = 0,
                                      ARWUL = 0, # 在isolet数据下 该算法会爆掉
                                      Ours = 0
                                      #Single_Linear = 0
                                      )

## 
for(d in 1:K){
  i = 1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  
  i = i+1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_mean_reg_simulation$W[,d])
  
  ## K specification is needed before
  #i = i+1
  #rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_clustering$W[,d])
  
  i = i+1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  i = i+1
  rsquare_allmodels_sim[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rsquare_allmodels_sim%>% round(4)
rsquare_allmodels_sim %>% colMeans() %>% round(4)
########
rmse_allmodels_sim = data.frame(L21 = 0,
                                   Lasso = 0,
                                   Trace = 0,
                                   Mean_regular = 0,
                                   #K_means = 0,
                                   DMFS = 0,
                                   ARWUL = 0, # 在isolet数据下 该算法会爆掉
                                   Ours = 0
                                   #Single_Linear = 0
)

for(d in 1:K){
  i = 1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  
  i = i+1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_mean_reg_simulation$W[,d])
  
  ## K specification is needed before
  #i = i+1
  #rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_clustering$W[,d])
  
  i = i+1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  i = i+1
  rmse_allmodels_sim[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rmse_allmodels_sim%>% round(4)
rmse_allmodels_sim %>% colMeans() %>% round(4)



###########
## sacros

rsquare_allmodels_sacros = data.frame(L21 = 0,
                                   Lasso = 0,
                                   Trace = 0,
                                   Mean_regular = 0,
                                   #K_means = 0,
                                   DMFS = 0,
                                   ARWUL = 0, # 在isolet数据下 该算法会爆掉
                                   Ours = 0
                                   #Single_Linear = 0
)

## 
for(d in 1:K){
  i = 1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  i = i+1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_mean_reg_simulation$W[,d])
  
  ## K specification is needed before
  #i = i+1
  #rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_clustering$W[,d])
  
  i = i+1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  i = i+1
  rsquare_allmodels_sacros[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rsquare_allmodels_sacros%>% round(4)
rsquare_allmodels_sacros %>% colMeans() %>% round(4)

########
rmse_allmodels_sacros = data.frame(L21 = 0,
                                Lasso = 0,
                                Trace = 0,
                                Mean_regular = 0,
                                #K_means = 0,
                                DMFS = 0,
                                ARWUL = 0, # 在isolet数据下 该算法会爆掉
                                Ours = 0
                                #Single_Linear = 0
)

## Sacros
for(d in 1:K){
  i = 1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  
  i = i+1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_mean_reg_simulation$W[,d])
  
  ## K specification is needed before
  #i = i+1
  #rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_clustering$W[,d])
  
  i = i+1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  i = i+1
  rmse_allmodels_sacros[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rmse_allmodels_sacros %>% round(4)
rmse_allmodels_sacros %>% colMeans() %>% round(4)



######## 
##for 'school' performance
{
rMSE_allmodels = data.frame(L21 = 0,
                               Lasso = 0,
                               Trace = 0,
                               DMFS = 0,
                               ARWUL = 0,
                               Ours = 0
                               #Single_Linear = 0
)
for(d in 1:K){
  i = 1
  rMSE_allmodels[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  i = i+1
  rMSE_allmodels[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rMSE_allmodels[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rMSE_allmodels[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rMSE_allmodels[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  
  i = i+1
  rMSE_allmodels[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rMSE_allmodels %>% colMeans()

####
rsquare_allmodels = data.frame(L21 = 0,
                            Lasso = 0,
                            Trace = 0,
                            DMFS = 0,
                            ARWUL = 0,
                            Ours = 0
                            #Single_Linear = 0
)
for(d in 1:K){
  i = 1
  rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  i = i+1
  rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  
  i = i+1
  rsquare_allmodels[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rsquare_allmodels
rsquare_allmodels %>% colMeans()

}
#########


#########
## isolet
rMSE_allmodels_isolet = data.frame(L21 = 0,
                                       Lasso = 0,
                                       Trace = 0,
                                       Mean_regular = 0,
                                       # K_means = 0,
                                       DMFS = 0,
                                       ARWUL = 0, # 在isolet数据下 该算法会爆掉 select a good step size will help
                                       Ours = 0
                                       #Single_Linear = 0)
                            )

for(d in 1:K){
  i = 1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  i = i+1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  i = i+1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_mean_reg_simulation$W[,d])
  
  i = i+1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  
  i = i+1
  rMSE_allmodels_isolet[d,i] <- rMSE(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}

rMSE_allmodels_isolet %>% colMeans() %>% round(4)


rsquare_allmodels_isolet = data.frame(L21 = 0,
                               Lasso = 0,
                               Trace = 0,
                               Mean_regular = 0,
                              # K_means = 0,
                               DMFS = 0,
                               ARWUL = 0, # 在isolet数据下 该算法会爆掉 select a good step size will help
                               Ours = 0
                               #Single_Linear = 0
)

## 
for(d in 1:K){
  i = 1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_L21_simulation$W[,d])
  i = i+1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Lasso_simulation$W[,d])
  i = i+1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_Trace_simulation$W[,d])
  
  i = i+1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_mean_reg_simulation$W[,d])
  
  i = i+1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = dmfs_results_simulation$W[,d])
  i = i+1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = ARWUL_simulation_vanilla_coef_matrix[,d])
  i = i+1
  rsquare_allmodels_isolet[d,i] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
  #rsquare_allmodels[d,7] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_linear_all_sim[,d])                               
}
rsquare_allmodels_isolet %>% round(5)
rsquare_allmodels_isolet %>% colMeans() %>% round(4)

########


if(simul_if == T){
  priornet_fit = c()
  clustering_fit = c()
  
  for(d in 1:K){
    
    priornet_fit[d] = Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_network_priornet_simulation$W[,d])
    
    clustering_fit[d] = Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = modelRMTL_clustering$W[,d])
  }
  
  rsquare_allmodels = data.frame(rsquare_allmodels,
                                 Priornet = priornet_fit,
                                 Clustering = clustering_fit)
  
}

rsquare_allmodels %>% round(4)
apply(rsquare_allmodels,2,mean)
#rsquare_schoolfor_ours = c()
#for(d in 1:K){
#  rsquare_schoolfor_ours [d] <- Rsquare(X=Xlist_test[[d]],Y=Ylist_test[[d]],beta = Beta_ours_all_mean_sim[,d])
#}


library(ggplot2)
library(tidyr)
rsquare_allmodels <- data.frame(xindex = seq(1,K),
                                rsquare_allmodels) %>% round(4)


df <- rsquare_allmodels %>%
  gather(key = "methods", value = "rsquare", -xindex)

max_yaxis <- max(df$rsquare) * 1.2

ggplot(df, aes(x = xindex, y = rsquare)) + 
  geom_line(aes(color = methods, linetype = methods)) + 
  scale_color_manual(values = 
                       c("darkred", "navyblue",'seagreen2','chocolate',
                         'deeppink','slateblue','cyan1','orangered','grey2'))+
  ylim(0,max_yaxis) +labs(x='task index')


