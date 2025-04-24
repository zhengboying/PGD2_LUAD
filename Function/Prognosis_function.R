# Function: getPk; deal_dup
#           get_cox; get_forest
#           get_km
#           fun_to_corr; get_chrod
# 区分肿瘤/正常样本
getPk <- function(data){
  barCode <- data
  group <- sapply(barCode,function(i){
    tmp <- unlist(strsplit(i,"-"))[4]
    z <- ""
    if (grepl("(01|02|03|04|05|06|07|08|09)(A|B|C)",tmp,ignore.case = T)){
      z <- "Tumor"
    }
    if (grepl("(11|12|13|14|15|16|17|18|19)(A|B|C)",tmp,ignore.case = T)){
      z <- "Normal"
    }
    return(z)
  })
  id <- sapply(barCode,function(i){
    tmp <- paste0(unlist(strsplit(i,"-"))[1:3],collapse = "-")
    return(tmp)
  })
  orders <- sapply(barCode,function(i){
    tmp <- unlist(strsplit(i,"-"))[4]
    tmp <- gsub("^[0-9]+","",tmp)
    return(tmp)
  })
  res <- cbind.data.frame(barCode,id,orders,group)
  res <- res[order(res$orders,decreasing = F),]
  labs <- paste(res$id,res$group,sep = "_")
  res <- res[!duplicated(labs),]
  return(res)
}

# 去除重复基因（重复基因表达取平均值）
deal_dup <- function(data,desc=c(1:2),index="Gene"){
  if(sum(duplicated(data[,names(data)==index]))!=0){
    Desc <- data[,desc]
    Desc <- Desc[!duplicated(Desc[,names(data)==index]),]
    data <- data[,c(which(names(data)==index),(desc[length(desc)]+1):ncol(data))]
    
    uni_gene <- unique(data[,names(data)==index])
    dup_gene <- data[,names(data)==index][duplicated(data[,names(data)==index])]
    dup_gene <- unique(dup_gene)
    dat1 <- data[!(data[,names(data)==index] %in% dup_gene),]
    
    dup_df <- lapply(1:length(dup_gene),function(i){
      temp_df <- data[data[,names(data)==index]==dup_gene[i],]
      temp <- apply(temp_df[,-1],2,mean,na.rm = T)
      temp <- t(temp)
      return(temp)
    })
    dat2 <- do.call(rbind.data.frame,dup_df)
    dat2 <- cbind.data.frame(dup_gene,dat2)
    names(dat2)[1] <- index
    dat <- rbind.data.frame(dat1,dat2)
    data <- merge(Desc,dat,by=index)
    return(data) 
  } else {
    return(data)
  }
}

get_km <- function(df,methods=c("mean","median","sur")[3],output, data.name = data.name,cutoffs= NULL){
  methods <- tolower(methods)
  # cox_model <- coxph(Surv(OS, Status)~., data = df)
  # cox_tab <- summary(cox_model)$coefficients
  # cox_tab <- as.data.frame(cox_tab)
  # row.names(cox_tab) <- gsub("`","",row.names(cox_tab))
  # write.xlsx(cox_tab,paste0(output,"/",data.name,"_cox_table.xlsx"),rowNames = T)
  # # df$Risk <- predict(cox_model,newdata = df,type='risk')
  
  for (i in 3:ncol(df)){
    surdf <- df
    osdf <- surdf[,c(1:2,i)]
    names(osdf)[3] <- "value"
    osdf$value[is.na(osdf$value)] <- min(osdf$value,na.rm = T)/2
    if (methods=="sur"){
      res.cut <- surv_cutpoint(osdf,
                               time = "OS",
                               event = "Status",
                               variables = "value"
      )
      cutoffs <- res.cut$cutpoint$cutpoint
    } else if (methods=="median"){
      cutoffs <- median(osdf$value)
    } else if (methods=="mean"){
      cutoffs <- mean(osdf$value)
    } else if (methods== "custom"){
      cutoffs <- cutoffs
    } else {
      stop("incorrect methods!")
    }
    message(paste0(data.name," ",i,"cutoffs: ",cutoffs))
    #print(osdf)
    osdf$group <- sapply(osdf$value,function(i){
      if (i>cutoffs){
        i <- "High"
      } else {i <- "Low"}
    })
    fit <- survfit(Surv(OS, Status)~group, data = osdf)
    print(fit)
    jpeg(filename = paste0(output,"/",data.name, '_', names(surdf)[i], '_', methods, "_OS.jpg"),width = 2500,height = 3400,res = 600)
    p <- ggsurvplot(fit,data = osdf, conf.int = T, risk.table = T,
                    palette = c("red","blue"),surv.median.line = "hv",pval = T,font.x = 14,
                    font.y = 14,font.tickslab = 12,font.legend = 10,xlab = paste0("Time (day)\n",names(surdf)[i]),ylab = "Survival probability (OS)",
                    ggtheme = theme_bw())
    print(p)
    dev.off()
    
    pdf(file = paste0(output,"/",data.name, '_', names(surdf)[i], '_', methods, "_OS.pdf"),width = 6,height = 8,pointsize = 6)
    p <- ggsurvplot(fit,data = osdf, conf.int = T, risk.table = T,
                    palette = c("red","blue"),surv.median.line = "hv",pval = T,font.x = 14,
                    font.y = 14,font.tickslab = 12,font.legend = 10,xlab = paste0("Time (day)\n",names(surdf)[i]),ylab = "Survival probability (OS)",
                    ggtheme = theme_bw())
    print(p,newpage = FALSE)
    dev.off()
  }
}

fun_to_corr <- function(df1, df2,methods = "pearson") {
  require(psych)
  require(reshape2)
  t <- corr.test(df1, df2, use = "pairwise", method = methods, adjust = "none",ci=F,alpha=0.05)
  t_cor <- data.frame(t$r, check.names = FALSE)
  t_p <- data.frame(t$p, check.names = FALSE)
  cut_sig <- function(p) {
    out <- cut(p, breaks = c(0, 0.001,0.01,0.05,1), include.lowest = T, labels = c("***", "**", "*", ""))
    return(out)
  }
  t_p_sig <- apply(t_p, 2, cut_sig)
  rownames(t_p_sig) <- rownames(t_p)
  
  get_df <- function(data,names="Corr"){
    data <- as.data.frame(data)
    data <- cbind.data.frame(row.names(data),data)
    names(data)[1] <- "node1"
    data <- melt(data,id.vars = "node1",variable.name = "node2",value.name = names)
    data$id <- paste(data$node1,data$node2,sep = ";;")
    return(data)
  }
  
  cordf <- get_df(t_cor,names = "Corr")
  pdf <- get_df(t_p,names = "P")
  sigdf <- get_df(t_p_sig,names = "Lable")
  
  merge_fun <- function(x,y){
    z <- merge.data.frame(x,y,by="id")
    return(z)
  }
  df <- Reduce(merge_fun,list(cordf,pdf,sigdf))
  df <- df[,c("node1","node2","Corr","P","Lable")]
  df$node1 <- as.character(df$node1)
  df$node2 <- as.character(df$node2)
  return(df)
}

get_chrod <- function(nets,nodes,sds=0,N = T,output){
  require(circlize)
  names(nets) <- c("node1","node2","scores")
  names(nodes) <- c("nodes","group")
  toCheck <- unique(c(nets$node1,nets$node2))
  if (all(toCheck %in% nodes$nodes)){
    print("network data paired with nodes data, start to plots.")
  } else {
    stop("your network data was not all in nodes data! please check you data.")
  }
  get_col <- function(n,seed = 123,random=T){
    c102 <- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941",
              "#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6",
              "#63FFAC","#B79762","#004D43","#8FB0FF","#997D87",
              "#5A0007","#809693","#6A3A4C","#1B4400","#4FC601",
              "#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900",
              "#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA",
              "#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299",
              "#300018","#0AA6D8","#013349","#00846F","#372101",
              "#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2",
              "#C2FF99","#001E09","#00489C","#6F0062","#0CBD66",
              "#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66",
              "#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459",
              "#456648","#0086ED","#886F4C","#34362D","#B4A8BD",
              "#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F",
              "#938A81","#575329","#00FECF","#B05B6F","#8CD0FF",
              "#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7",
              "#A77500","#6367A9","#A05837","#6B002C","#772600",
              "#D790FF","#9B9700","#549E79","#FFF69F","#201625",
              "#72418F","#BC23FF","#99ADC0","#3A2465","#922329",
              "#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98",
              "#A4E804","#324E72"
    )
    if (random){
      if (n>102){
        return(rainbow(n))
      } else {
        set.seed(seed)
        z <- sample(c102,102,replace = F)
        return(z[1:n])
      }
    } else {
      if (n>102){
        return(rainbow(n))
      } else {
        return(c102[1:n])
      }
    }
  }
  
  mem <- unique(nodes$nodes)
  col_var <- get_col(length(mem),seed = sds,random = F)
  names(col_var) <- mem
  
  mem1 <- unique(nodes$group)
  col_grp <- get_col(length(mem1),seed = sds,random = T)
  names(col_grp) <- mem1
  
  if (N){
    # mid <- (min(nets$scores)+max(nets$scores))/2
    # col_line <- colorRamp2(c(min(nets$scores),mid, max(nets$scores)), c("darkblue","gray50", "darkred"))
    col_line <- colorRamp2(c(-1,0, 1), c("darkblue","gray50", "darkred"))
    cols <- sapply(nets$scores,col_line)
  } else {
    cols <- "gray50"
  }
  
  pdf(file=paste0(output,"/circlize.pdf"), width=12, height=12, pointsize=8)
  chordDiagram(nets, 
               col='#00000000',
               grid.col = col_var,
               directional = 1, direction.type='arrows',
               link.arr.lwd = 2,
               link.arr.width = 0.2,
               link.arr.lty = NULL,
               link.arr.type = "triangle",
               link.arr.col = cols,
               diffHeight = 0.2,
               annotationTrack = 'grid',
               annotationTrackHeight = c(0.07, 0.07),
               preAllocateTracks = list(
                 track.height = mm_h(6),
                 track.margin = c(mm_h(8), 0)
               ) 
  )
  
  cls <- unique(nodes$group)
  for (i in 1:length(cls)){
    tmpNode <- nodes$nodes[nodes$group==cls[i]]
    txt <- cls[i]
    highlight.sector(tmpNode, track.index = 1,col = col_grp[i],facing = "clockwise",
                     text.col = "gray50",border = "black",lwd = 1.5,
                     text = txt, cex = 1.5, font= 2,niceFacing = TRUE
    )
  }
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=.8)
  }, bg.border = NA)
  circos.clear()
  dev.off()
  
  jpeg(file=paste0(output,"/circlize.jpg"), width=12, height=12, res = 600,units = "in")
  chordDiagram(nets, 
               col='#00000000',
               grid.col = col_var,
               directional = 1, direction.type='arrows',
               link.arr.lwd = 2,
               link.arr.width = 0.2,
               link.arr.lty = NULL,
               link.arr.type = "triangle",
               link.arr.col = cols,
               diffHeight = 0.2,
               annotationTrack = 'grid',
               annotationTrackHeight = c(0.07, 0.07),
               preAllocateTracks = list(
                 track.height = mm_h(6),
                 track.margin = c(mm_h(8), 0)
               ) 
  )
  
  cls <- unique(nodes$group)
  for (i in 1:length(cls)){
    tmpNode <- nodes$nodes[nodes$group==cls[i]]
    txt <- cls[i]
    highlight.sector(tmpNode, track.index = 1,col = col_grp[i],facing = "clockwise",
                     text.col = "gray50",border = "black",lwd = 1.5,
                     text = txt, cex = 1.5, font= 2,niceFacing = TRUE
    )
  }
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=.8)
  }, bg.border = NA)
  circos.clear()
  dev.off()
}

get_cox <- function(df,starts = 3){
  require(survival)
  names(df)[1:2] <- c("OS","Status")
  covariates <- names(df)[starts:ncol(df)]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(OS, Status)~', "`",x,"`",sep = "")))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = df)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value <- signif(x$wald["pvalue"], digits=4)
                           wald.test <- signif(x$wald["test"], digits=4)
                           beta <- signif(x$coef[1], digits=4);#coeficient beta
                           HR <- signif(x$coef[2], digits=4);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                           # HR <- paste0(HR, " (", 
                           #              HR.confint.lower, "-", HR.confint.upper, ")")
                           # names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                           #               "p.value")
                           res <- c(beta, HR,HR.confint.lower,HR.confint.upper, wald.test, p.value)
                           names(res) <- c("Beta","HR","Lower","Upper","Wald.test","Pvalue")
                           return(res)
                         })
  res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  res$Pvalue <- as.numeric(res$Pvalue)
  res <- cbind.data.frame(row.names(res),res)
  names(res)[1] <- "Var"
  return(res)
}

get_forest <- function(data,output){
  require(forestplot)
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  df <- data
  df$Pvalue <- signif(df$Pvalue,4)
  p <- forestplot(df[,c("Var","Pvalue")], 
                  mean=df$HR,   
                  lower=df$Lower,  
                  upper=df$Upper,
                  zero = 1,
                  fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                  boxsize = 0.3, ##误差条中的圆心点大小
                  col=fpColors(line = "gray50", #误差条的线的颜色
                               box="darkred"), #误差条的圆心点的颜色
                  lty.ci = 7,   # 误差条的线的线型
                  lwd.ci = 3,   # 误差条的线的宽度
                  ci.vertices.height = 0.15, # # 误差条末端的长度
                  txt_gp = fpTxtGp(ticks = gpar(cex = 1.2), xlab = gpar(cex = 0.8), cex = 1.2), #文本大小设置
                  lineheight = "auto") |>
    fp_add_header(Var = c("Genes"),
                  Pvalue = c("P value"))
  heights <- nrow(df)*.2+2
  widths <- heights+.3
  jpeg(paste0(
    output, '/Ca_gene_cor_forests.jpg'),
    width = widths, height = heights, units = 'in', res = 600, pointsize = 8
  )
  print(p)
  dev.off()
  
  pdf(paste0(
    output, '/Ca_gene_cor_forests.pdf'),
    width = widths, height = heights, pointsize = 8
  )
  print(p)
  dev.off()
}