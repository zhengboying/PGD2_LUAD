library(forestplot)
data <- read.csv("./04_sur/Univar_result_OS.csv")
head(data)
#2,3,4,5列取4位有效数字
data$HR <- round(data$HR, 4)
data$lower.95 <- round(data$lower.95, 4)
data$upper.95 <- round(data$upper.95, 4)
data$pvalue <- round(data$pvalue, 4)

data1 <- data %>% dplyr::filter(HR > 1)
data2 <- data %>% dplyr::filter(HR < 1)
hrtable <- rbind(  c("Increase in Hazard", NA, NA, NA, NA, NA),
                   data1,  c("Reduce in Hazard", NA, NA, NA, NA, NA), 
                   data2)
# 将要在图中展示的文本
tabletext <- cbind(  
  c("Gene", "Increase in Hazard", data1$gene, "Reduce in Hazard", data2$gene),  
  c("P-value", NA, data1$pvalue, NA, data2$pvalue), 
  c("Hazard Ratio", NA, data1$`HR`, NA, data2$`HR`))
which(tabletext[, 1] == "Reduction in Hazard")
nrow(tabletext) + 1
quantile(data$lower)
quantile(data$upper)
pdf("./04_sur/Univar_forestplot.pdf", width =6, height = 6) 


##所有tabletext数值取4位有效数字

p<-forestplot(  
  labeltext = tabletext, 
  # 图中的文本  
  mean = c(NA, as.numeric(hrtable$HR)), 
  # HR  
  lower = c(NA, as.numeric(hrtable$lower)), 
  upper = c(NA, as.numeric(hrtable$upper)),
  graph.pos = 4, 
  # 图在表中的列位置  graphwidth = unit(.4, "npc"), 
  # 图在表中的宽度比例  fn.ci_norm = "fpDrawDiamondCI", 
  # box类型选择钻石
  col = fpColors(box = "#CD534CFF",lines="#0073C2FF", zero = "#0073C2FF"),
  #box颜色  
  boxsize = 0.7, 
  # box大小根据样本量设置  
  lwd.ci = 3, 
  ci.vertices.height = 0.125,  
  ci.vertices = T, 
  # 置信区间用线宽、高、型  
  zero = 1, 
  # zero线横坐标  
  lwd.zero = 2,  
  xticks = seq(0, 4, 0.5), 
  # 横坐标刻度根据需要可随意设置  
  lwd.xaxis = 2, # X轴线宽  
  xlab = "Hazard Ratio",  
  hrzl_lines = list("1" = gpar(lwd = 2, col = "black"), 
                     # 第二行顶部加黑实线    
                     "2" = gpar(lwd = 1, col = "grey50", lty = 2),
                      # 第二行顶部加灰色虚线    
                     "3" = gpar(lwd = 1, col = "grey50", lty = 2), 
                    "4" = gpar(lwd = 1, col = "grey50", lty = 2),
                    # 第二行顶部加灰色虚线    
                    "5" = gpar(lwd = 1, col = "grey50", lty = 2), 
                      # 第九行顶部加灰色虚线    
                      "8" = gpar(lwd = 2, col = "black")  ),
  # 最后一行底部加黑线，""中数字为nrow(tabletext) + 1  
  lineheight = unit(0.8, "cm"),   
  txt_gp = fpTxtGp(label = gpar(cex = 1), 
                  # 各种字体大小设置    
                  ticks = gpar(cex = 1.25),    
                  xlab = gpar(cex = 1.85),    
                  title = gpar(cex = 1.5)),  
  clip = c(0, 2),  
  colgap = unit(0.15, "cm"),   
  mar = unit(rep(1.25, times = 4), "cm"), 
  #  new_page = F # 是否新页
  )
p

dev.off()