
#draw the correlation hotmap
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#corr_matrix = covnmat_post_median

corr_matrix = covnmat_post_mean

diag(corr_matrix) <- NA
rownames(corr_matrix) <- paste0("task ", 1:8)
colnames(corr_matrix) <- paste0("task ", 1:8)

# 把对角线设为 NA（可选）
diag(corr_matrix) <- 0

# 转换为 long format
df_long <- melt(corr_matrix, varnames = c("Var1", "Var2"), value.name = "Correlation")

# 设定统一顺序
task_order <- paste0("task ", 1:8)
df_long$Var1 <- factor(df_long$Var1, levels = rev(task_order))
df_long$Var2 <- factor(df_long$Var2, levels = task_order)

# 然后照常画图

### Task Correlation Heatmap for simulation data
ggplot(df_long, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "red3", midpoint = 0,
    limits = c(-1, 1), na.value = "white", name = "Correlation"
  ) +
  theme_minimal(base_size = 12) +
  coord_fixed() +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#########

sim_Beta_summary <- function(n, b, task,iniburn = burnin){
  Beta_mat = matrix(0,nrow = 80, ncol = 8)
  
  t = 1
  for(j in ((b-iniburn)/interval + 1):((n-iniburn)/interval) ){
    Beta_mat<-  Beta_mat + resultlist[[j]]$Beta  
    t = t+1
  } 
  return((Beta_mat / (t)))
}


### true beta
Beta_alltasks_J_80


# 示例：创建一个 80 x 8 的矩阵（你可以替换为自己的矩阵）
set.seed(123)
coef_matrix_est <- sim_Beta_summary(n=niter,b=burnin,iniburn=burnin) 
coef_matrix_true = Beta_alltasks_J_80

# 给行和列命名（可选）
rownames(coef_matrix_est) <- paste0("var_", 1:80)
colnames(coef_matrix_est) <- paste0("task_", 1:8)
####
df_long_est <- melt(coef_matrix_est, varnames = c("Variable", "Task"), value.name = "Coefficient")
df_long_est$Variable <- factor(df_long_est$Variable, levels = rev(unique(df_long$Variable)))
df_long_est$AbsCoef <- abs(df_long_est$Coefficient)

#########

rownames(coef_matrix_true) <- paste0("var_", 1:80)
colnames(coef_matrix_true) <- paste0("task_", 1:8)
coef_df_true <- as.data.frame(coef_matrix_true)
coef_df_true$Variable <- rownames(coef_matrix_true)

# 然后 melt 的时候指定 ID 变量为 Variable
df_long_true <- melt(coef_df_true, id.vars = "Variable", variable.name = "Task", value.name = "Coefficient")
df_long_true$Variable <- factor(df_long_true$Variable, levels = rev(unique(df_long_true$Variable)))

# 取绝对值
df_long_true$AbsCoef <- abs(df_long_true$Coefficient)

####### 需要统一颜色的口径
library(patchwork)

library(reshape2)
library(ggplot2)
library(patchwork)  # 用于并排画图

# 假设你已有两个矩阵 coef1 和 coef2（80×8）
# 示例生成：
set.seed(123)

df_long_est$Matrix <- "Matrix 1"

df_long_true$Matrix <- "Matrix 2"

# 合并两个 data frame
df_both <- rbind(df_long_est, df_long_true)

# 保证 Variable 顺序一致
df_both$Variable <- factor(df_both$Variable, levels = rev(unique(df1$Variable)))
df_both$Task <- factor(df_both$Task, levels = unique(df1$Task))

max_abs <- max(df_both$AbsCoef)

p1 <- ggplot(df_long_est, aes(x = Task, y = factor(Variable, levels = rev(unique(Variable))), fill = AbsCoef)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(low = "black", high = "white", limits = c(0, max_abs), name = "|Coef|") +
  ggtitle("Estimated Absolute Coefficient Magnitudes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none" ) + labs(x = "", y = "")

p2 <- ggplot(df_long_true, aes(x = Task, y = factor(Variable, levels = rev(unique(Variable))), fill = AbsCoef)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(low = "black", high = "white", limits = c(0, max_abs), name = "|Coef|") +
  ggtitle("True Absolute Coefficient Magnitudes") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(x = "", y = "")
p1+p2        
