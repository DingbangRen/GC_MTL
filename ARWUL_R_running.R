
### 为了让ARWU读取 需要将数据转为json

##
#install.packages("https://cran.r-project.org/src/contrib/Archive/reticulate/reticulate_1.4.tar.gz", repos=NULL, type="source")
##
{
  library(jsonlite)
  
  write_json(X_train,  'C:\\Users\\scott\\Downloads\\GCVS_simulatedata_forpy\\X_train.json')
  write_json(Y_train,  'C:\\Users\\scott\\Downloads\\GCVS_simulatedata_forpy\\Y_train.json')
  write_json(X_test,  'C:\\Users\\scott\\Downloads\\GCVS_simulatedata_forpy\\X_test.json')
  write_json(Y_test,  'C:\\Users\\scott\\Downloads\\GCVS_simulatedata_forpy\\Y_test.json')
}
##需要在ipython上跑
## http://localhost:8888/notebooks/Downloads/ARMUL-main/ARWUL_model.ipynb

########
### 如果ARWUL 数值爆掉了 或者趋于inf, 
### 就会报错NA

######
##ARWUL 跑出的结果
## running on python
ARWUL_simulation_vanilla_coef = read_json('C:/Users/scott/Downloads/ARMUL-main/vanilla_simulation_coef_matrix.json')


ARWUL_simulation_vanilla_coef_matrix = matrix(0,nrow = J, ncol = K)
for(k in 1:K){
  ARWUL_simulation_vanilla_coef_matrix[,k] = ARWUL_simulation_vanilla_coef[k]%>% unlist()
}
ARWUL_simulation_vanilla_coef_matrix %>% head(5)

