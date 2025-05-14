#data preparation for real data

#######Sacros
library(dplyr)
{library('R.matlab')
  
  R.matlab::readMat("C:\\Users\\scott\\Downloads\\sarcos_inv.mat") -> sarcos_data
  sarcos_df <- data.frame(sarcos_data)
  
  Y_sarcos <- sarcos_df[,22:28]
  str(Y_sarcos)
  
  X_sarcos <- sarcos_df[,1:21]
  str(X_sarcos)
  
  names(Y_sarcos) <- NULL
  names(X_sarcos) <- NULL
  
  Y_sarcos %>%as.matrix() -> Y_sarcos
  X_sarcos %>%as.matrix() -> X_sarcos
}
nrow(Y_sarcos)

# sacros data
## keep the original Y, Z score to the X[,j]
## random 500 samples for training dataset
{set.seed(123)
  randomrows <- sample(nrow(X_sarcos),size=500,replace=F)
  
  X_sarcos_semifinal <- X_sarcos[randomrows,]
  
  X_sarcos_semifinal %>% str()
  
  X_sarcos_zscore <- matrix(0,nrow(X_sarcos_semifinal),ncol(X_sarcos_semifinal))
  for(i in 1:ncol(X_sarcos_semifinal)){
    X_sarcos_zscore[,i] <- (X_sarcos_semifinal[,i] - mean(X_sarcos_semifinal[,i])) / sd(X_sarcos_semifinal[,i])
  }
  
  X_sarcos_final <- cbind(rep(1,nrow(X_sarcos_zscore))%>%c(),X_sarcos_zscore)
  ###考虑是否要带intercept
  
  #X_sarcos_final <- X_sarcos_zscore
  X_sarcos_final <- X_sarcos_final
  
  
  X_sarcos_final %>% glimpse()
  
  Y_sarcos_final <- Y_sarcos[randomrows,]
  
  Y_sarcos_final %>% str()
  Y_sarcos_final %>% head()
  ####
  # 构造training data
  
  Xlist_sarcos <- list()
  for(i in 1:ncol(Y_sarcos_final)){
    Xlist_sarcos[[i]] <- X_sarcos_final
  }
  Xlist_sarcos%>%str()
  
  Ylist_sarcos <- list()
  for(i in 1:ncol(Y_sarcos_final)){
    Ylist_sarcos[[i]] <- Y_sarcos_final[,i]%>% as.matrix()
  }
  Ylist_sarcos%>%str()
  ##########
}
# 构造testing data

{ 
  allsam_vec <- 1:nrow(X_sarcos)
  allsam_vec_rest =  allsam_vec[!(allsam_vec %in% randomrows)]
  set.seed(123)
  
  randomrows_test <- sample(allsam_vec_rest,size = 1000,replace = F)
  sum(randomrows_test %in% randomrows)
  
  X_sarcos_semifinal_test <- X_sarcos[randomrows_test,]
  
  X_sarcos_semifinal_test %>% str()
  
  X_sarcos_zscore_test <- matrix(0,nrow(X_sarcos_semifinal_test),ncol(X_sarcos_semifinal_test))
  for(i in 1:ncol(X_sarcos_semifinal_test)){
    X_sarcos_zscore_test[,i] <- (X_sarcos_semifinal_test[,i] - 
                                   mean(X_sarcos_semifinal_test[,i])) / sd(X_sarcos_semifinal_test[,i])
  }
  
  X_sarcos_final_test <- cbind(rep(1,nrow(X_sarcos_zscore_test))%>%c(),X_sarcos_zscore_test)
  ###考虑是否要带intercept
  
  #X_sarcos_final_test <- X_sarcos_zscore
  X_sarcos_final_test <- X_sarcos_final_test
  
  
  X_sarcos_final_test %>% glimpse()
  
  Y_sarcos_final_test <- Y_sarcos[randomrows_test,]
  
  Y_sarcos_final_test %>% str()
  
  ####
  # 构造testing data
  
  Xlist_sarcos_test <- list()
  for(i in 1:ncol(Y_sarcos_final_test)){
    Xlist_sarcos_test[[i]] <- X_sarcos_final_test
  }
  Xlist_sarcos_test%>%str()
  
  Ylist_sarcos_test <- list()
  for(i in 1:ncol(Y_sarcos_final_test)){
    Ylist_sarcos_test[[i]] <- Y_sarcos_final_test[,i]%>% as.matrix()
  }
  Ylist_sarcos_test%>%str()
  
}

#save(Xlist_sarcos, file = 'C:\\Users\\scott\\Xlist_sarcos.RData')
#save(Ylist_sarcos, file = 'C:\\Users\\scott\\Ylist_sarcos.RData')

load(file = 'C:\\Users\\scott\\Xlist_sarcos.RData')
load(file = 'C:\\Users\\scott\\Ylist_sarcos.RData')

#save(Xlist_sarcos_test, file = 'C:\\Users\\scott\\Xlist_sarcos_test.RData')
#save(Ylist_sarcos_test, file = 'C:\\Users\\scott\\Ylist_sarcos_test.RData')

load(file = 'C:\\Users\\scott\\Xlist_sarcos_test.RData')
load(file = 'C:\\Users\\scott\\Ylist_sarcos_test.RData')

X_train = Xlist_sarcos
Y_train = Ylist_sarcos

Xlist_test = Xlist_sarcos_test
Ylist_test = Ylist_sarcos_test

J = ncol(Xlist_sarcos[[1]])
K = length(Ylist_sarcos)

#-------------------------------------------------------------------------------

## schools data experiments
{
{schoolsdata <- read.delim('C:\\Users\\scott\\Downloads\\schools original data.txt')
str(schoolsdata)
df <- schoolsdata %>% as.matrix()
## 第一行数据读不进来 需要补充
rbind('1  117241803 111',df) -> df
colnames(df) <- NULL
glimpse(df)


i=1
school_df1 <- data.frame()
while(i <= nrow(df)){
  data.frame(Year = df[i,]%>%substr(start=1,stop=1)%>%as.numeric(),
             School = df[i,]%>%substr(start=2,stop=4)%>%as.numeric(),
             Score = df[i,]%>%substr(start=5,stop=6)%>%as.numeric(),
             FSM = df[i,]%>%substr(start=7,stop=8)%>%as.numeric(),
             VR1band = df[i]%>%substr(start=9,stop=10)%>%as.numeric(),
             Gender = df[i,]%>%substr(start=11,stop=11)%>%as.numeric(),
             VRstudent = df[i]%>%substr(start=12,stop=12)%>%as.numeric(),
             Ethnics = df[i,]%>%substr(start=13,stop=14)%>%as.numeric(),
             SchoolGender = df[i]%>%substr(start=15,stop=15)%>%as.numeric(),
             Denomination = df[i,]%>%substr(start=16,stop=16)%>%as.numeric() ) -> newrow_school
  school_df1 <- rbind(school_df1,newrow_school)
  i = i+1
}

school_df1 %>% head()
school_df1 %>% glimpse()

##########  将categorical data 变成binary
library(nnet)

## year
school_df1$Year %>% n_distinct()
school_df1$Year %>% class.ind() -> Year_dummy
Year_dummy %>% head()

colnames(Year_dummy) <- sapply(as.numeric(colnames(Year_dummy)),FUN = function(x){
  return(paste('Year',x,sep='_'))
})
Year_dummy %>% head()

## VR band
## 在介绍中 说只有3个值，但实际有4个，还有=0的value
school_df1$VRstudent %>% n_distinct()
school_df1$VRstudent%>% class.ind() -> VRstudent_dummy
VRstudent_dummy %>% head()

which(school_df1$VRstudent==0)
school_df1$VRstudent[9978]


VRstudent_dummy[,1] %>% n_distinct()

colnames(VRstudent_dummy) <- sapply(as.numeric(colnames(VRstudent_dummy)),FUN = function(x){
  return(paste('VRstudent',x,sep='_'))
})
VRstudent_dummy %>% head()

## Ethnics
school_df1$Ethnics%>% n_distinct()
school_df1$Ethnics%>% class.ind() -> Ethnics_dummy
Ethnics_dummy %>% head()

colnames(Ethnics_dummy) <- sapply(as.numeric(colnames(Ethnics_dummy)),FUN = function(x){
  return(paste('Ethnics',x,sep='_'))
})
Ethnics_dummy%>% head()

## schoolgender
school_df1$SchoolGender%>%n_distinct()
school_df1$SchoolGender%>% class.ind() -> SchoolGender_dummy
SchoolGender_dummy %>% head()

colnames(SchoolGender_dummy) <- sapply(as.numeric(colnames(SchoolGender_dummy)),FUN = function(x){
  return(paste('SchoolGender',x,sep='_'))
})
SchoolGender_dummy %>% head()

##school denomination
school_df1$Denomination%>%n_distinct()
school_df1$Denomination%>% class.ind() -> Denomination_dummy
Denomination_dummy %>% head()

colnames(Denomination_dummy) <- sapply(as.numeric(colnames(Denomination_dummy)),FUN = function(x){
  return(paste('Denomination',x,sep='_'))
})
Denomination_dummy %>% head()

## 合并
cbind(Year_dummy,School=school_df1$School,Score=school_df1$Score,
      FSM=school_df1$FSM,VR1band=school_df1$VR1band,Gender=school_df1$Gender,
      VRstudent_dummy,Ethnics_dummy,SchoolGender_dummy,Denomination_dummy) %>% data.frame() -> schools_withdummy
schools_withdummy %>% head()
schools_withdummy %>% glimpse()


#### 根据school name将data分类成list中的不同元素
split(schools_withdummy,factor(schools_withdummy$School)) -> schools_list
schools_list[[1]] %>% glimpse()
schools_list[[ floor(n_distinct(schools_withdummy$School)/2) ]] %>% glimpse()
schools_list[[unique(schools_withdummy$School)%>%tail(n=1)]] %>% glimpse()

}
#########
Ylist_schools <- list()
for(i in 1:n_distinct(schools_withdummy$School)){
  Ylist_schools[[i]] <- schools_list[[i]]$Score %>% as.matrix()
}
str(Ylist_schools)

Xlist_schools <- list()
for(j in 1:n_distinct(schools_withdummy$School)){
  prep_X <- subset(schools_list[[j]],select = -c(School,Score)) %>% as.matrix()
  colnames(prep_X) <- NULL
  rownames(prep_X) <- NULL
  Xlist_schools[[j]] <- cbind(rep(1,nrow(prep_X)), prep_X)
}

Xlist_schools %>% glimpse()
varnames_schooldata <- colnames(subset(schools_list[[1]],select = -c(School,Score)) %>% as.matrix())

#save(Xlist_schools,file = 'C:\\Users\\scott\\Xlist_schools.RData')
#save(Ylist_schools,file = 'C:\\Users\\scott\\Ylist_schools.RData')

load('C:\\Users\\scott\\Xlist_schools.RData')
load('C:\\Users\\scott\\Ylist_schools.RData')

#######################

{
  J=ncol(Xlist_schools[[1]]) 
  #K=n_distinct(schools_withdummy$School)
  K= 139
  
  X_train = list()
  Y_train = list()
  X_test = list()
  Y_test = list()
  
  #### 要不就随机选20个tasks 看看效果
  set.seed(123)
  K=20
  
  random_schools = sample(length(Xlist_schools), size=K,replace = F)
  
  ### 尝试normalize the reponse, 不然过于scattered了
  z_norm_vec = function(vec){
    var_vec = var(vec)
    mean_vec = mean(vec)
    normed_vec = (vec - mean_vec) / sqrt(var_vec)
    return(normed_vec)
  }
  ### 但是z-score y的效果也不好
  ### 所以要不z-score predictors
  ### 对于predictors来说 5, 6的数值最诡异 需要特殊处理下？
  ### 但是同样也会报错
  z_norm_to_X_56= function(X_mat){
    X_mat[,5] =  X_mat[,5] %>% z_norm_vec()
    X_mat[,6] =  X_mat[,6] %>% z_norm_vec()
    return(X_mat)
  }
  ###
  ### 又有论文说 可以只用27个feature不需要intercept
  #J = J-1
  ###
  
  train_perc = 0.7
  
  for(i in 1:K){
    set.seed(i)
    prep_rows = 1:nrow(Xlist_schools[[ random_schools[i] ]])
    ran_train_rows = sample(prep_rows ,size = floor(train_perc * length(prep_rows)),replace = F)
    ran_test_rows = prep_rows[!(prep_rows%in%ran_train_rows)]
    
    
    X_train_i_prep = Xlist_schools[[ random_schools[i] ]][ran_train_rows,]
    #X_train_i_prep = X_train_i_prep %>% z_norm_to_X_56()
    #X_train_i_prep = X_train_i_prep[,-1]
    X_train[[i]] = X_train_i_prep
    
    Y_train_i_prep = Ylist_schools[[ random_schools[i] ]][ran_train_rows] 
    #Y_train_i_prep  = Y_train_i_prep %>% z_norm_y_school()
    Y_train[[i]] = Y_train_i_prep %>% as.matrix()
    
    ###
    
    
    X_test_i_prep = Xlist_schools[[ random_schools[i] ]][ran_test_rows,]
    #X_test_i_prep = X_test_i_prep %>% z_norm_to_X_56()
    #X_test_i_prep = X_test_i_prep[,-1]
    X_test[[i]] = X_test_i_prep
    
    Y_test_i_prep = Ylist_schools[[ random_schools[i] ]][ran_test_rows]
    #Y_test_i_prep = Y_test_i_prep %>% z_norm_y_school()
    Y_test[[i]] = Y_test_i_prep %>% as.matrix()
  }
}

str(X_train)
str(X_test)

}

## isolet
#-------------------------------------------------------------
isolet_original_1234 <- read.table(file='C:\\Users\\scott\\Downloads\\isolet1+2+3+4.txt')
ncol(isolet_original_1234)
(nrow(isolet_original_1234)+2) / 4

isolet_original_5 <- read.table(file='C:\\Users\\scott\\Downloads\\isolet5.txt')
nrow(isolet_original_5)


(nrow(isolet_original_1234)+2) / (52*4) # =120
## 相当于30 * 4, 即每组都有30个 说了52遍单词的 speakers

isolet_numerical_1234 <- matrix(0,nrow=nrow(isolet_original_1234),ncol=ncol(isolet_original_1234))
for(i in 1:nrow(isolet_numerical_1234)){
  isolet_numerical_1234[i,] <- 
    gsub(',','',isolet_original_1234[i,]  %>% droplevels()  %>% ## factor以及data.frame的形式 都会产生问题
           as.matrix() %>% as.character()) %>% as.numeric() 
}


gsub(',','',isolet_original_1234[i,]  %>% droplevels()  %>%
       as.matrix() %>% as.character()) %>% as.numeric() 

isolet_numerical_1234 %>% str()

isolet_numerical_5 <- matrix(0,nrow=nrow(isolet_original_5),ncol=ncol(isolet_original_5))
for(i in 1:nrow(isolet_numerical_5)){
  isolet_numerical_5[i,] <- 
    gsub(',','',isolet_original_5[i,]  %>% droplevels()  %>% ## factor以及data.frame的形式 都会产生问题
           as.matrix() %>% as.character()) %>% as.numeric() 
}

isolet_numerical_5 %>% str()
## 先在最后添加两行 方便循环读取数据最后再删掉

##############
## PCA???

isolet_preready_1234<-  rbind(isolet_numerical_1234,matrix(0,nrow=2,ncol = ncol(isolet_numerical_1234)))
isolet_ready_all <- rbind(isolet_preready_1234,isolet_numerical_5)
isolet_ready_all %>% str()

isolet_pca_ready <- prcomp(isolet_ready_all[,-ncol(isolet_ready_all)],
                           center = TRUE,scale = TRUE) #scale = T?
names(isolet_pca_ready) 

## percent of var has been explaint
summary(isolet_pca_ready)$importance[3,] %>% plot() ## 100的variables大概能解释0.8¬0.9
#The prcomp object stores the eigenvectors as a matrix in the rotation attribute.
summary(isolet_pca_ready)$importance[3,100]

## eigen vec of PC
pca_eigenvec_iso <- isolet_pca_ready$rotation
pca_eigenvec_iso %>% str()

## PC scores which is the projection of the original data
pca_scores_iso <- isolet_pca_ready$x
pca_scores_iso %>% str()


pcnum <- 100
isolet_pca_all <- cbind(pca_scores_iso[,1:pcnum],isolet_ready_all[,ncol(isolet_ready_all)])
isolet_pca_all %>% str()
#####

isolet_ready_1234 <- isolet_pca_all[1:nrow(isolet_preready_1234),]
isolet_ready_5 <- isolet_pca_all[(nrow(isolet_preready_1234)+1):nrow(isolet_pca_all),]

isolet_ready_1234[(1*3+51*2  + 52*4*2):(1*3+51*3 + 52*4*2),ncol(isolet_ready_1234)]
(nrow(isolet_ready_1234) - (1*4+51*4)) / (52*4) # 29

iso1 <- rep(0,ncol(isolet_ready_1234)) 
iso2 <- rep(0,ncol(isolet_ready_1234))
iso3 <- rep(0,ncol(isolet_ready_1234))
iso4 <- rep(0,ncol(isolet_ready_1234)) 


for(i in 1:(1+(nrow(isolet_ready_1234) - (1*4+51*4)) / (52*4))){
  iso1 <- rbind(iso1,isolet_ready_1234[(1   +        52*4*(i-1)): (1   + 51   +  52*4*(i-1)),])
  iso2 <- rbind(iso2,isolet_ready_1234[(1*2 + 51   + 52*4*(i-1)): (1*2 + 51*2 +  52*4*(i-1)),])
  iso3 <- rbind(iso3,isolet_ready_1234[(1*3 + 51*2 + 52*4*(i-1)): (1*3 + 51*3 +  52*4*(i-1)),])
  iso4 <- rbind(iso4,isolet_ready_1234[(1*4 + 51*3 + 52*4*(i-1)): (1*4 + 51*4 +  52*4*(i-1)),])
}


iso1 <- iso1[-1,]
iso2 <- iso2[-1,]
iso3 <- iso3[-1,]
## task 4末尾的两个rows是missing的 要去掉
iso4 <- iso4[-1,]
iso4 <- iso4[1:(nrow(iso4)-2),]

iso5 <-  isolet_ready_5
iso5 %>% str()

#############

uni_rows_iso = function(dat){
  apply(dat, 1, sum) %>% unique() -> uni_sumrow
  apply(dat, 1, sum) -> sumrow
  dat.mat = matrix(0,ncol=ncol(dat),nrow=length(uni_sumrow))
  for(i in 1:length(uni_sumrow)){
    dat.mat[i,] = dat[which(sumrow==uni_sumrow[i])[1],]
  }
  return(dat.mat)
}

iso1_uni = uni_rows_iso(iso1)
iso2_uni = uni_rows_iso(iso2)
iso3_uni = uni_rows_iso(iso3)
iso4_uni = uni_rows_iso(iso4)
iso5_uni = uni_rows_iso(iso5)


############
## training data
Xlist_isolet_train = list()
Xlist_isolet_test = list()

Ylist_isolet_train = list()
Ylist_isolet_test = list()

set.seed(113)
train_perc = 0.5
ran_rows_train_1 = sample(1:nrow(iso1_uni), size = floor(train_perc * nrow(iso1_uni)), replace = F)
ran_rows_train_2 = sample(1:nrow(iso2_uni), size = floor(train_perc * nrow(iso2_uni)), replace = F)
ran_rows_train_3 = sample(1:nrow(iso3_uni), size = floor(train_perc * nrow(iso3_uni)), replace = F)
ran_rows_train_4 = sample(1:nrow(iso4_uni), size = floor(train_perc * nrow(iso4_uni)), replace = F)
ran_rows_train_5 = sample(1:nrow(iso5_uni), size = floor(train_perc * nrow(iso5_uni)), replace = F)

ran_rows_test_1 <- (1:nrow(iso1_uni))[!( (1:nrow(iso1_uni)) %in% ran_rows_train_1)]
ran_rows_test_2 <- (1:nrow(iso2_uni))[!( (1:nrow(iso2_uni)) %in% ran_rows_train_2)]
ran_rows_test_3 <- (1:nrow(iso3_uni))[!( (1:nrow(iso3_uni)) %in% ran_rows_train_3)]
ran_rows_test_4 <- (1:nrow(iso4_uni))[!( (1:nrow(iso4_uni)) %in% ran_rows_train_4)]
ran_rows_test_5 <- (1:nrow(iso5_uni))[!( (1:nrow(iso5_uni)) %in% ran_rows_train_5)]

############
Xlist_isolet_train[[1]] <- matrix(cbind(rep(1,length(ran_rows_train_1)) , iso1_uni[ran_rows_train_1,-(pcnum+1)]), 
                                  nrow = length(ran_rows_train_1), ncol = pcnum+1)

Xlist_isolet_train[[2]] <- matrix(cbind(rep(1,length(ran_rows_train_2)) , iso2_uni[ran_rows_train_2,-(pcnum+1)]), 
                                  nrow = length(ran_rows_train_2), ncol = pcnum+1)

Xlist_isolet_train[[3]] <- matrix(cbind(rep(1,length(ran_rows_train_3)) , iso3_uni[ran_rows_train_3,-(pcnum+1)]), 
                                  nrow = length(ran_rows_train_3), ncol = pcnum+1)

Xlist_isolet_train[[4]] <- matrix(cbind(rep(1,length(ran_rows_train_4)) , iso4_uni[ran_rows_train_4,-(pcnum+1)]), 
                                  nrow = length(ran_rows_train_4), ncol = pcnum+1)

Xlist_isolet_train[[5]] <- matrix(cbind(rep(1,length(ran_rows_train_5)) , iso5_uni[ran_rows_train_5,-(pcnum+1)]), 
                                  nrow = length(ran_rows_train_5), ncol = pcnum+1)
#############
Xlist_isolet_test[[1]] <- matrix(cbind(rep(1,length(ran_rows_test_1)) , iso1_uni[ran_rows_test_1,-(pcnum+1)]), 
                                 nrow = length(ran_rows_test_1), ncol = pcnum+1)

Xlist_isolet_test[[2]] <- matrix(cbind(rep(1,length(ran_rows_test_2)) , iso2_uni[ran_rows_test_2,-(pcnum+1)]), 
                                 nrow = length(ran_rows_test_2), ncol = pcnum+1)

Xlist_isolet_test[[3]] <- matrix(cbind(rep(1,length(ran_rows_test_3)) , iso3_uni[ran_rows_test_3,-(pcnum+1)]), 
                                 nrow = length(ran_rows_test_3), ncol = pcnum+1)

Xlist_isolet_test[[4]] <- matrix(cbind(rep(1,length(ran_rows_test_4)) , iso4_uni[ran_rows_test_4,-(pcnum+1)]), 
                                 nrow = length(ran_rows_test_4), ncol = pcnum+1)

Xlist_isolet_test[[5]] <- matrix(cbind(rep(1,length(ran_rows_test_5)) , iso5_uni[ran_rows_test_5,-(pcnum+1)]), 
                                 nrow = length(ran_rows_test_5), ncol = pcnum+1)

#############
Ylist_isolet_train[[1]] <- iso1_uni[ran_rows_train_1,pcnum+1] %>% as.matrix()
Ylist_isolet_train[[2]] <- iso2_uni[ran_rows_train_2,pcnum+1] %>% as.matrix()
Ylist_isolet_train[[3]] <- iso3_uni[ran_rows_train_3,pcnum+1] %>% as.matrix()
Ylist_isolet_train[[4]] <- iso4_uni[ran_rows_train_4,pcnum+1] %>% as.matrix()
Ylist_isolet_train[[5]] <- iso5_uni[ran_rows_train_5,pcnum+1] %>% as.matrix()

Ylist_isolet_test[[1]] <- iso1_uni[ran_rows_test_1,pcnum+1] %>% as.matrix()
Ylist_isolet_test[[2]] <- iso2_uni[ran_rows_test_2,pcnum+1] %>% as.matrix()
Ylist_isolet_test[[3]] <- iso3_uni[ran_rows_test_3,pcnum+1] %>% as.matrix()
Ylist_isolet_test[[4]] <- iso4_uni[ran_rows_test_4,pcnum+1] %>% as.matrix()
Ylist_isolet_test[[5]] <- iso5_uni[ran_rows_test_5,pcnum+1] %>% as.matrix()

#save(Xlist_isolet_train,file='C:\\Users\\scott\\Xlist_isolet_train.RData')
#save(Xlist_isolet_test,file='C:\\Users\\scott\\Xlist_isolet_test.RData')
#save(Ylist_isolet_train,file='C:\\Users\\scott\\Ylist_isolet_train.RData')
#save(Ylist_isolet_test,file='C:\\Users\\scott\\Ylist_isolet_test.RData')

load(file='C:\\Users\\scott\\Xlist_isolet_train.RData')
load(file='C:\\Users\\scott\\Xlist_isolet_test.RData')
load(file='C:\\Users\\scott\\Ylist_isolet_train.RData')
load(file='C:\\Users\\scott\\Ylist_isolet_test.RData')
###### 为什么X中有这么多相同的rows?????????
##### 这会导致ARWUL 生成NA

X_train = Xlist_isolet_train
X_test = Xlist_isolet_test

Y_train = Ylist_isolet_train
Y_test = Ylist_isolet_test

Xlist_test = X_test
Ylist_test = Y_test
#######################
## MCMC for isolet

J = ncol(X_train[[1]])
K = 5
