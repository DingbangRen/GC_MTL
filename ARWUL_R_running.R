
### 为了让ARWU读取 需要将数据转为json
script_file <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else "ARWUL_R_running.R"
repo_root <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
arwul_exchange_dir <- file.path(repo_root, "results", "legacy_exchange", "ARWUL")
dir.create(arwul_exchange_dir, recursive = TRUE, showWarnings = FALSE)
exchange_file <- function(...) file.path(arwul_exchange_dir, ...)
armul_root <- function(){
  path <- Sys.getenv("GCVS_ARMUL_ROOT", "")
  if (!nzchar(path)) {
    stop("Set GCVS_ARMUL_ROOT to the ARMUL repository before importing ARMUL outputs.")
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

##
#install.packages("https://cran.r-project.org/src/contrib/Archive/reticulate/reticulate_1.4.tar.gz", repos=NULL, type="source")
##
{
  library(jsonlite)
  
  write_json(X_train, exchange_file("X_train.json"))
  write_json(Y_train, exchange_file("Y_train.json"))
  write_json(X_test, exchange_file("X_test.json"))
  write_json(Y_test, exchange_file("Y_test.json"))
}
##需要在ipython上跑
## Run the Python notebook from the directory pointed to by GCVS_ARMUL_ROOT.

########
### 如果ARWUL 数值爆掉了 或者趋于inf, 
### 就会报错NA

######
##ARWUL 跑出的结果
## running on python
ARWUL_simulation_vanilla_coef = read_json(file.path(armul_root(), "vanilla_simulation_coef_matrix.json"))


ARWUL_simulation_vanilla_coef_matrix = matrix(0,nrow = J, ncol = K)
for(k in 1:K){
  ARWUL_simulation_vanilla_coef_matrix[,k] = ARWUL_simulation_vanilla_coef[k]%>% unlist()
}
ARWUL_simulation_vanilla_coef_matrix %>% head(5)
