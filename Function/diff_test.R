#' @title perform_pairwise_tests
#' 
#' @description This function testing diff between groups, focus on one groups vs the rest groups, and return the results.
#' 
#' @param df data.frame for expr,mets or any other data,the rownames is the features, 
#' the colnames is the samples,and all samples must be in the group_table rownames.
#' @param group_table data.frame for group_table, the rownames is the samples,you can choose any column as the group_name.
#' @param group_name the column name in the group_table, which is used to define the groups.
#' @param diff_test_method the method to test the diff between groups, now only support "t.test","w.test","auc.test".
#' @param p.adjust.method the method to adjust the p value, default is "BH".
#' @return all_results.results[['group1']] is the result for group1 vs the rest groups.
#
#

#循环实现多组差异检验
perform_pairwise_tests <- function(df = flux, group_table = res, group_name = "Phenotype", 
                                   diff_test_method = "t.test", p.adjust.method = "BH",remove_sum0_features = TRUE) {
  # 获取所有的分组类别
  all_groups <- unique(group_table[[group_name]])
  
  # 创建一个列表保存结果
  all_results <- list()
  
  # 对每个元素与其余的分组进行差异检验
  for (group in all_groups) {
    # 其余分组作为control
    other_groups <- all_groups[all_groups != group]
    
    # 打印当前正在处理的group
    message("Performing test for case =", group, "; control = other groups")
    
    # 将group作为case，其他剩余的作为control进行差异分析
    group_table_temp <- group_table
    group_table_temp[[group_name]][group_table_temp[[group_name]] %in% other_groups] <- "Other"
    
    # 调用 perform_diff_test
    result <- perform_diff_test(df = df, 
                                group_table = group_table_temp, 
                                group_name = group_name, 
                                case = group, 
                                control = "Other", 
                                diff_test_method = diff_test_method, 
                                p.adjust.method = p.adjust.method,
                                remove_sum0_features = remove_sum0_features)
    
    # 将结果存储到列表中
    all_results[[group]] <- result
  }
  
  # 返回所有结果
  return(all_results)
}



#' @title perform_diff_test
#' @description This function testing diff between groups, focus on case and contol, and return the results.
#' @param df data.frame for expr,mets or any other data,the rownames is the features,
#' the colnames is the samples,and all samples must be in the group_table rownames.
#' @param group_table data.frame for group_table, the rownames is the samples,you can choose any column as the group_name.
#' @param group_name the column name in the group_table, which is used to define the groups.
#' @param case the case group name.
#' @param control the control group name.
#' @param diff_test_method the method to test the diff between groups, now only support "t.test","w.test","auc.test".
#' @param p.adjust.method the method to adjust the p value, default is "BH".
#' @return results.results$p.value is the p value,results$statistic is the statistic,results$mean_case is the mean of case,results$mean_control is the mean of control.

perform_diff_test <- function(df=flux, 
                              group_table = res, 
                              group_name = "Phenotype",
                              case="Hot",
                              control="Cold",
                              diff_test_method = "t.test",##目前仅支持t_test,Wilcoxon test(w.test),auc test
                              p.adjust.method = "BH",
                              remove_sum0_features = TRUE
) {
  message("case =",case,"; ","control =",control)
  # 提取分组信息
  #group <- group_table[[group_name]]
  ##提取group_table中group=case的行名，不是行编号
  case_samples <- rownames(group_table)[group_table[[group_name]] == case]
  control_samples <- rownames(group_table)[group_table[[group_name]] == control]
  ##提取df中存在的数据
  valid_case_samples <- case_samples[case_samples %in% colnames(df)]
  #print(valid_case_samples[1:5])
  
  valid_control_samples <- control_samples[control_samples %in% colnames(df)]
  
  #print(valid_control_samples[1:5])
  print(paste0("valid_case_samples:",length(valid_case_samples)))
  print(paste0("valid_control_samples:",length(valid_control_samples)))
  #删除整行都是0的缺失值
  if (remove_sum0_features == TRUE){df <- df[rowSums(df) > 0, ]}
  
  if (diff_test_method == "t.test") {
    message("Performing t-test ...")
    # 应用函数到df的每一行
    results <- apply(df, 1, function(x) {
      # 进行t检验
      test <- t.test(x[valid_case_samples],x[valid_control_samples])
      # 返回p值和其他统计量
      #direction <- ifelse(test$statistic >0, "up", "down")
      c(p.value = test$p.value, 
        statistic = test$statistic, 
        mean_case = mean(x[valid_case_samples],na.rm = TRUE), 
        mean_control = mean(x[valid_control_samples],na.rm = TRUE))
    })
    # 转置结果，使得每一列代表一个统计量
    results <- as.data.frame(t(results))
    #增加p.adjust by BH
    results$p.adjust <- p.adjust(results$p.value, method = p.adjust.method)
    results$log2FC <- log2(results$mean_case/results$mean_control)
    # 返回结果
    return(results)
  } else if (diff_test_method == "w.test") {
    message("Performing Wilcoxon test ...")
    results <- apply(df, 1, function(x) {
      #print(x[valid_case_samples])
      test <- wilcox.test(x[valid_case_samples], x[valid_control_samples])
      #direction <- ifelse(test$statistic >0, "up", "down")
      c(p.value = test$p.value, 
        statistic = test$statistic, 
        mean_case = mean(x[valid_case_samples],na.rm = TRUE), 
        mean_control = mean(x[valid_control_samples],na.rm = TRUE))
      #direction = direction
    })
    results <- as.data.frame(t(results))
    results$p.adjust <- p.adjust(results$p.value, method = p.adjust.method)
    results$log2FC <- log2(results$mean_case/results$mean_control)
    return(results)
  } else if (diff_test_method == "auc.test") {
    library(pROC)
    message("Performing AUC scoring ...")
    results <- apply(df, 1, function(x) {
      case_values <- x[valid_case_samples]
      control_values <- x[valid_control_samples]
      # 创建标签和预测值
      response <- c(rep(1, length(case_values)), rep(0, length(control_values)))
      predictor <- c(case_values, control_values)
      ##计算roc
      roc_obj <- roc(response, predictor)
      direction <- ifelse(roc_obj$direction == "<", "up", "down")
      c(auc = auc(roc_obj),
        mean_case = mean(case_values),
        mean_control = mean(control_values))
      direction = direction
    })
    results <- as.data.frame(t(results))
    return(results)
  } else {
    stop("diff_test_method must be 't.test', 'w.test' or 'auc.test'")
  }
}
##使用perform_diff_test
#perform_diff_test(df=flux, 
#                  group_table = res, 
#                  group_name = "Phenotype",
#                  case="Hot",
#                  control="Cold",
#                  diff_test_method = "t.test",
#                  p.adjust.method = "BH")



cor_test<-function(data =data, 
                   name="TCGA",
                   method = "spearman"){
  # 初始化矩阵存储相关性和 p 值
  cor_matrix <- matrix(NA, ncol = ncol(data), nrow = ncol(data))
  p_values <- matrix(NA, ncol = ncol(data), nrow = ncol(data))
  
  rownames(cor_matrix) <- colnames(data)
  colnames(cor_matrix) <- colnames(data)
  rownames(p_values) <- colnames(data)
  colnames(p_values) <- colnames(data)
  
  # 计算 Spearman 相关性和 p 值
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      if (i != j) {
        test <- cor.test(data[, i], data[, j], method = method)
        cor_matrix[i, j] <- test$estimate
        p_values[i, j] <- test$p.value
      } else {
        cor_matrix[i, j] <- NA
        p_values[i, j] <- NA
      }
    }
  }
  
  # 将相关性和显著性数据转换为长格式
  cor_data <- melt(cor_matrix, na.rm = TRUE)
  pval_data <- melt(p_values, na.rm = TRUE)
  
  # 合并相关性和显著性数据
  heatmap_data <- cbind(cor_data, pval_data$value)
  colnames(heatmap_data) <- c("Var1", "Var2", "Correlation", "P_value")
  
  # 创建显著性分级：将 p 值转换为颜色分组
  heatmap_data$Significance <- cut(
    heatmap_data$P_value,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns") # ns = not significant
  )
  
  # 绘制相关性热图（包括显示相关性值）
  p<-ggplot(heatmap_data, aes(Var1, Var2)) +
    geom_point(aes(size = abs(Correlation), color = Significance), alpha = 0.8) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4, color = "black") +
    scale_size_continuous(range = c(3, 10), name = "Correlation\n(abs)") +
    scale_color_manual(values = c("***" = "red", "**" = "orange", "*" = "yellow", "ns" = "gray"),
                       name = "Significance\n(p-value)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    labs(title = "Correlation Heatmap",
         x = "", y = "") +
    coord_fixed()
  write.csv(heatmap_data, file = paste0(ResultPath,name,"correlation_heatmap.csv"), row.names = FALSE)
  ggsave(plot = p,paste0(ResultPath,name,"correlation_heatmap.pdf"),
         bg = "white", width = 5, height = 5)
  return(p)
}
