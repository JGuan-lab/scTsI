library(Rtsne)
library(tsne)
library(ggplot2)
library(Rmisc)
library(ggpubr)
library(lemon)
suppressMessages(library('monocle'))
library(TSCAN)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggrepel)
library(parallel)
library(igraph)
library(Seurat)
library(mclust)
library(devtools)
library(fpc)
library(ngram)
library(FNN)
library(Matrix)
library(Metrics)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Rcpp::sourceCpp("G:/Zhanghongyu/data/julei.cpp")

cluster_evalu<-function(H,realcluster){
  cluster <- kmeans(H,3)[["cluster"]]
  reallist<-list()
  cluslist<-list()
  for (i in 1:length(unique(realcluster))) {
    reallist[[i]]<-which(realcluster==unique(realcluster)[i])
    reallist[[i]]<-as.vector(reallist[[i]])
  }
  for (i in 1:length(unique(cluster))) {
    cluslist[[i]]<-which(cluster==unique(cluster)[i])
    cluslist[[i]]<-as.vector(cluslist[[i]])
  }
  nmi <- NMI(reallist,cluslist)
  ari <- adjustedRandIndex(realcluster,cluster)
  stats <- cluster.stats(dist(H)^2, cluster)
  sc <- stats$avg.silwidth
  clu_eva<-list(ari,sc,nmi)
  return(clu_eva)
}

# ARI、NMI、SC
evalu_mean <- function(impute,realcluster){
  H = t(impute)
  ari<-vector(length=100)
  nmi<-vector(length=100)
  for(j in c(1:5)){
    test=cluster_evalu(H,realcluster)
    ari[j]<-test[[1]]
    nmi[j]<-test[[2]]
  }
  ari_ave=mean(ari)
  nmi_ave=mean(nmi)
  result <- list(ari_ave,nmi_ave)
  names(result)=c('ARI','NMI')
  return(result)
}


########simulation data matrix pcc/RMSE,cell pss.
get_cor_result <- function(data_true, impute_data){
  options( warn = -1 ) 
  
  # cell-cell correlation
  
  data_true_cell = cor(as.matrix((data_true)))
  
  data_true_cell[is.na(data_true_cell)] = 0
  
  # gene-gene correlation
  
  data_true_gene = cor(t((data_true)), method = "pearson")
  
  data_true_gene[is.na(data_true_gene)] = 0
  
  
  
  # cell-cell correlation
  
  impute_data_cell = cor((impute_data), method = "pearson")
  
  impute_data_cell[is.na(impute_data_cell)] = 0
  
  # gene-gene correlation
  
  impute_data_gene = cor(t((impute_data)), method = "pearson")
  
  impute_data_gene[is.na(impute_data_gene)] = 0
  
  
  vector_true<- as.vector(unlist(data_true))
  vector_impute_data <- as.vector(unlist(impute_data))
  
  c_true<- as.vector(unlist(data_true_cell))
  c_impute_data <- as.vector(unlist(impute_data_cell))
  
  g_true<- as.vector(unlist(data_true_gene))
  g_impute_data <- as.vector(unlist(impute_data_gene))
  
  cell_pcc<-cor(c_true,c_impute_data)
  
  gene_pcc<-cor(g_true,g_impute_data)
  
  matrix_pcc <- cor(vector_true,vector_impute_data)
  
  RMSE <- rmse(data_true,impute_data)
  
  
  
  result <- vector()
  
  result[[1]] = RMSE
  
  result[[2]] = matrix_pcc
  
  result[[3]] = cell_pcc
  
  result[[4]] = gene_pcc
  
  
  return(result)
  
}


##################### trajectory ####################################
ctlevel <- data.frame(ct=c('0','1','2','3','4','5','6'),stringsAsFactors = F)
row.names(ctlevel) <- ctlevel[,1]

correctorder <- wrongorder <- NULL
evct <- ctlevel[,1]
pair <- expand.grid(evct,evct)
pair[,1] <- as.character(pair[,1])
pair[,2] <- as.character(pair[,2])
pair <- pair[pair[,1]!=pair[,2],]
corid <- which(ctlevel[pair[,1],1] < ctlevel[pair[,2],1])
wroid <- which(ctlevel[pair[,1],1] > ctlevel[pair[,2],1])
correctorder <- c(correctorder,sapply(corid,function(si) paste0(pair[si,],collapse = '_')))
wrongorder <- c(wrongorder,sapply(wroid,function(si) paste0(pair[si,],collapse = '_')))

prep_traj_order <- function(cds,group){
  cds$H2228_to_H1975 = NA
  cds$H2228_to_H1975[group=="0 9 0"] = 0
  cds$H2228_to_H1975[group=="1 7 1"] = 1
  cds$H2228_to_H1975[group=="2 5 2"] = 2
  cds$H2228_to_H1975[group=="3 3 3"] = 3
  cds$H2228_to_H1975[group=="5 2 2"] = 4
  cds$H2228_to_H1975[group=="7 1 1"] = 5
  cds$H2228_to_H1975[group=="9 0 0"] = 6
  
  cds$H2228_to_HCC827 = NA
  cds$H2228_to_HCC827[group=="0 9 0"] = 0
  cds$H2228_to_HCC827[group=="1 7 1"] = 1
  cds$H2228_to_HCC827[group=="2 5 2"] = 2
  cds$H2228_to_HCC827[group=="3 3 3"] = 3
  cds$H2228_to_HCC827[group=="2 2 5"] = 4
  cds$H2228_to_HCC827[group=="1 1 7"] = 5
  cds$H2228_to_HCC827[group=="0 0 9"] = 6
  
  cds$H1975_to_HCC827 = NA
  cds$H1975_to_HCC827[group=="9 0 0"] = 0
  cds$H1975_to_HCC827[group=="7 1 1"] = 1
  cds$H1975_to_HCC827[group=="5 2 2"] = 2
  cds$H1975_to_HCC827[group=="3 3 3"] = 3
  cds$H1975_to_HCC827[group=="2 2 5"] = 4
  cds$H1975_to_HCC827[group=="1 1 7"] = 5
  cds$H1975_to_HCC827[group=="0 0 9"] = 6
  return(cds)
}

prep_RNA_traj_order <- function(cds,group){
  cds$H2228_to_H1975 = NA
  cds$H2228_to_H1975[cds$group=="1 0 0"] = 0
  cds$H2228_to_H1975[cds$group=="0.68 0.16 0.16"] = 1
  cds$H2228_to_H1975[cds$group=="0.33 0.33 0.33"] = 2
  cds$H2228_to_H1975[cds$group=="0.16 0.68 0.16"] = 3
  cds$H2228_to_H1975[cds$group=="0 1 0"] = 4
  
  cds$H2228_to_HCC827 = NA
  cds$H2228_to_HCC827[cds$group=="1 0 0"] = 0
  cds$H2228_to_HCC827[cds$group=="0.68 0.16 0.16"] = 1
  cds$H2228_to_HCC827[cds$group=="0.33 0.33 0.33"] = 2
  cds$H2228_to_HCC827[cds$group=="0.16 0.16 0.68"] = 3
  cds$H2228_to_HCC827[cds$group=="0 0 1"] = 4
  
  cds$H1975_to_HCC827 = NA
  cds$H1975_to_HCC827[cds$group=="0 1 0"] = 0
  cds$H1975_to_HCC827[cds$group=="0.16 0.68 0.16"] = 1
  cds$H1975_to_HCC827[cds$group=="0.33 0.33 0.33"] = 2
  cds$H1975_to_HCC827[cds$group=="0.16 0.16 0.68"] = 3
  cds$H1975_to_HCC827[cds$group=="0 0 1"] = 4
  return(cds)
}

prep_traj_wrapper <- function(cds,group){
  if ('0 0 9' %in% cds$group){
    return(prep_traj_order(cds,group))
  } else {
    return(prep_RNA_traj_order(cds,group))
  }
}


get_RNAmix_cds <- function(expression_matrix){
  set.seed(12345)
  # preproccess
  
  group <- sub('\\..*\\.','',colnames(expression_matrix))
  # Change the cellmix colnames
  cn <- gsub('\\.','\\:',colnames(expression_matrix))
  colnames(expression_matrix) = cn
  
  # Change the RNAmix colnames
  #cn <- gsub('\\.',':',colnames(expression_matrix))
  #cn <- gsub('0:16','0.16',cn)
  #cn <- gsub('0:68','0.68',cn)
  #cn <- gsub('0:33','0.33',cn)
  #colnames(expression_matrix) = cn
  
  prefix = unique(sub(':.*','',colnames(expression_matrix)))  
  if (FALSE){
    raw = readRDS(paste0('./data/processed/cellbench/',sub('.rds','',f),'/genebycell.rds'))
    colnames(expression_matrix) = colnames(raw)
    # cn <- sapply(cn,function(i) {
    #   paste0(round(as.numeric(strsplit(i,'_')[[1]]) * 9),collapse = '_')
    # },USE.NAMES = F)
    cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
  } else {
    cn = sapply(colnames(expression_matrix), function(i) strsplit(i,':')[[1]][2])
    cn <- sapply(cn,function(i) {
      j <- as.numeric(strsplit(i,'_')[[1]])
      paste0(j,collapse='_')
    },USE.NAMES = F)
  } 
  
  group <- sub('.*:','',colnames(expression_matrix))
  colnames(expression_matrix) <- paste0(colnames(expression_matrix),'_',1:ncol(expression_matrix))
  print(str(group))
  prop <- t(sapply(group,function(i) {
    as.numeric(strsplit(sub('.*:','',i),'_')[[1]])
  }))
  cell_metadata <- data.frame(cell=colnames(expression_matrix),p1=prop[,1],p2=prop[,2],p3=prop[,3])
  row.names(cell_metadata) <- colnames(expression_matrix)
  gene_annotation <- data.frame(gene_short_name=row.names(expression_matrix))
  row.names(gene_annotation) <- row.names(expression_matrix)
  
  
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  cds <- newCellDataSet(as.matrix(expression_matrix),phenoData = pd, featureData = fd,expressionFamily=uninormal())
  
  cds$group = gsub('_',' ',group)
  cds <- prep_traj_wrapper(cds,cds$group)
  flag <- 0
  tryCatch({cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0);flag <- 1},warning=function(w){},error=function(e){})
  if (flag == 0) {
    cds = reduceDimension(cds, method = "DDRTree",norm_method="none",pseudo_expr=0,auto_param_selection=F)
  }
  
  cds = orderCells(cds)
  tmp = table(pData(cds)$State[pData(cds)$group=='0_1_0'])  # specify H2228 as root state.
  # tmp = table(pData(cds)$State[pData(cds)$group=='0 9 0'])  # specify H2228 as root state.

  tmp = tmp[order(tmp,decreasing = T)]
  cds = orderCells(cds,root_state=names(tmp)[1],num_paths=3) # reorder the cells
  return(cds)
}


get_max_score = function(col_anno, states){
  col_st = col_anno[col_anno$State %in% states,]
  # cor1 = cor(col_st$Pseudotime, col_st$H2228_to_H1975, use = "pairwise.complete.obs",method='spearman')
  # cor2 = cor(col_st$Pseudotime, col_st$H2228_to_HCC827, use = "pairwise.complete.obs",method='spearman')
  # cor3 = cor(col_st$Pseudotime, col_st$H1975_to_HCC827, use = "pairwise.complete.obs",method='spearman')
  ov1 = sum(col_anno[!is.na(col_anno$H2228_to_H1975),"State"] %in% states)/sum(!is.na(col_anno$H2228_to_H1975))
  ov2 = sum(col_anno[!is.na(col_anno$H2228_to_HCC827),"State"] %in% states)/sum(!is.na(col_anno$H2228_to_HCC827))
  ov3 = sum(col_anno[!is.na(col_anno$H1975_to_HCC827),"State"] %in% states)/sum(!is.na(col_anno$H1975_to_HCC827))
  #sp1 = 1-sum(col_anno[is.na(col_anno$H2228_to_H1975),"State"] %in% states)/sum(is.na(col_anno$H2228_to_H1975))
  #sp2 = 1-sum(col_anno[is.na(col_anno$H2228_to_HCC827),"State"] %in% states)/sum(is.na(col_anno$H2228_to_HCC827))
  #sp3 = 1-sum(col_anno[is.na(col_anno$H1975_to_HCC827),"State"] %in% states)/sum(is.na(col_anno$H1975_to_HCC827))
  #res_df = data.frame(corr=abs(c(cor1,cor2,cor3)),overlap=c(ov1,ov2,ov3),sp=c(sp1,sp2,sp3))
  # res_df = data.frame(corr=abs(c(cor1,cor2,cor3)),overlap=c(ov1,ov2,ov3))
  res_df = data.frame(overlap=c(ov1,ov2,ov3))
  
  res_df = res_df[order(res_df,decreasing = T),]
  return(res_df[1])
}
# calculate overlap
get_RNAmix_ov <- function(cds){
  max_score_df = list()
  if (length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points) > 0) {
    for(i in 1:length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)){
      print(i)
      # rm(cds_reduced)
      tryCatch({cds_reduced <- buildBranchCellDataSet(cds, branch_point=i)},error=function(e) {}) # get branch info
      if (exists('cds_reduced')) {
        col_data_df = as.data.frame(pData(cds))
        max_score_df[[i]] = lapply(1:length(unique(pData(cds_reduced)$Branch)),function(j){
          state1 = table(pData(cds_reduced)$State[pData(cds_reduced)$Branch == unique(pData(cds_reduced)$Branch)[j]])
          state1 = state1[state1>=1] ## Aug30,19. change > to >=
          print(state1)
          max_score_df1 = get_max_score(col_data_df, names(state1))
        })
        max_score_df[[i]] = Reduce(rbind,max_score_df[[i]])      
      }
    }      
    if (length(max_score_df) ==0) {
      c(NA)
    } else {
      max_score_df <- max_score_df[!sapply(max_score_df,is.null)]
      tmp = unlist(lapply(max_score_df,function(x){colMeans(x)[1]}))         ### autoimpute
      finalres = max_score_df[[which(tmp==max(tmp))]]
      return(colMeans(finalres))
    }
  } else {
    c(NA)
  }
  
}


kendall_func <- function(x,y,tot=NULL){
  if(length(x)>1&length(y)>1){
    P   <- 0
    Q   <- 0
    for(i in 1:(length(x)-1)){
      for(j in (i+1):length(x)){
        if(sign(x[i]-x[j])*sign(y[i]-y[j])<0){
          Q = Q + 1
        }
        
        if(sign(x[i]-x[j])*sign(y[i]-y[j])>0){
          P = P + 1
        }
        
      }
    }
    if(is.null(tot)){tot=length(x)}
    out <- (P-Q)/choose(tot,2) # option 1, slingshot, max
  }else{
    out <- 0
  }
  
  return(out)
}

# calculate kendall
score_kendall <- function(cds) {
  if (length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points) > 0) {
    sl <- NULL
    s1 <- NULL
    for(i in 1:length(cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points)){
      tryCatch({cds_reduced <- buildBranchCellDataSet(cds,branch_point=i)},error=function(e) {})
      df = data.frame(pData(cds_reduced),stringsAsFactors = F)
      df <- df[order(df$Pseudotime),]
      sl <- rbind(sl,sapply(unique(df$Branch),function(ub) {
        
        so <- as.character(df[df[,'Branch']==ub,1])
        # soct <- ct[match(so,ct[,1]),3]
        x1 = as.numeric(na.omit(cds$H2228_to_H1975[match(so,cds$cell)]))
        x2 = as.numeric(na.omit(cds$H2228_to_HCC827[match(so,cds$cell)]))
        x3 = as.numeric(na.omit(cds$H1975_to_HCC827[match(so,cds$cell)]))
        z<-c()
        # The first column loops through each element in sequence, 
        # The second column loops through each element by number of lengths
        for (pid in list(x1,x2,x3)){
          y = matrix(1:length(pid))
          z <- rbind( kendall_func(pid,y), z )
        } 
        z <- max(z)
        # eid <- expand.grid(1:length(soct),1:length(soct))
        # eid <- eid[eid[,1]<eid[,2],]
        # eid <- sprintf('%s_%s',soct[eid[,1]],soct[eid[,2]])
        # c(sum(eid %in% correctorder),sum(eid %in% wrongorder))
        
      }))
    }
    s1 = mean(sl[which.max(rowMeans(sl)),])
  } else {
    NA
  }
}



###############plot########################
scaleFUN <- function(x) sprintf("%.2f", x)

pcc_plot <- function(mydata){
  mydata <- mydata[-2,]
  mydata <- melt(mydata,id="Method")
  colnames(mydata) <- c("Method","Dropout_rate","value")
  pcc<- ggplot(data = mydata,aes(x=Method,y=value,group =Dropout_rate,color=Dropout_rate))+
    geom_point()+
    geom_line()+
    xlab('800×1000')+
    ylab("PCC")+
    theme_bw() +
    theme(
      panel.grid.major=element_line(colour=NA),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.minor = element_blank(),
      axis.title =  element_text(size=16,face = "bold"),
      axis.text = element_text(size=10), 
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"), 
      #text = element_text(family = "STXihei"),
      #legend.position = c(.075,.075),
      #legend.box.background = element_rect(color="black")
    )+
      scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061","#440154FF"))+
      # scale_y_continous(limits=c(min(y),max(y)),breaks=round(seq(min(y),max(y),length.out=3),2))
    scale_y_continuous(labels=scaleFUN)
  return(pcc)
}

rmse_plot <- function(mydata){
  
  mydata <- mydata[-2,]
  mydata <- melt(mydata,id="Method")
  colnames(mydata) <- c("Method","Dropout_rate","value")
  rmse<- ggplot(data = mydata,aes(x=Method,y=value,group =Dropout_rate,color=Dropout_rate))+
    geom_point()+
    geom_line()+
    xlab('800×1000')+
    ylab("RMSE")+
    theme_bw() +
    theme(
      panel.grid.major=element_line(colour=NA),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.minor = element_blank(),
      axis.title =  element_text(size=16,face = "bold"),
      axis.text = element_text(size=10), 
      legend.text = element_text(size = 14), 
      legend.title = element_text(size = 16, face = "bold"), 
      #text = element_text(family = "STXihei"),
      #legend.position = c(.075,.075),
      #legend.box.background = element_rect(color="black")
    )+
    scale_color_manual(values = c("#FDE5A7","#FDE725FF","#FFC200","#FFA500","#FF8C00","#FF6600"))
    # scale_y_continuous(labels=scaleFUN)
  return(rmse)
}

cell_cell_plot <- function(mydata){
  cc <- ggplot(mydata, aes(x=Dropout, y=value, fill=Method)) +
        geom_bar(stat="identity", position=position_dodge(),
           color="white", width=.9) +
        scale_y_continuous(expand = c(0, 0),limits = c(0,1.1))+
        labs(x = '800×1000', y = 'cell-cell')+
        theme_bw()+
        theme(
            panel.grid = element_blank(),
            axis.title =  element_text(size=10,face = "bold")
        )  +   
        guides(fill = guide_legend(title = NULL))+
        scale_fill_manual(values = c("#7F2704", "#963003", "#B13A02", "#D14501", "#E45709", "#F16B16", "#F9812F", "#FD974A", "#FDAB67", "#FDC088", "#FDD3A8", "#FDE0C3", "#FEEBD9"))
  return(cc)

}
time_plot <- function(mydata){
  mydata <- melt(mydata,id="Method")
  colnames(mydata) <- c("Method","dimension","value")#更改列名
  time<-ggplot(data=mydata, aes(x=Method,y=log(value+1))) +     
  geom_boxplot(aes(fill = Method),outlier.shape = NA)+# 移除异常值点
  theme_bw()+
  geom_jitter(width = 0.2, alpha = 0.7)+
  scale_fill_manual(values = colors) +  # 手动设置填充颜色
  labs(x = "",  # 设置x轴标题
       y = "Time/log(s+1)")+  theme(
         panel.grid = element_blank(),
         axis.title =  element_text(size=16,face = "bold"),#设置标题字体大小
         axis.text = element_text(size=10), 
         legend.text = element_text(size = 14), # 设置图例标签字体大小
         legend.title = element_text(size = 16, face = "bold"), # 设置图例标题字体大小
         axis.text.x = element_text(angle = 25, hjust = 1),  # 控制 x 轴标签角度和对齐,
         # axis.text.x = element_blank(),  # 隐藏 x 轴的标签)
       )
  return(time)
}

simulated_data_tsne_plot <- function(data,realcluster,name){
  options( warn = -1 ) 
  impute_data <- data
  tsne.coords <- Rtsne(t(impute_data),pca=FALSE,perplexity = 30,theta=0.5,check_duplicates = F)$Y
  #rownames(tsne.coords[[k]]) <- rownames(t(data_true))
  colnames(tsne.coords) <- c("tSNE1","tSNE2")
  tsne.coords <- as.data.frame(tsne.coords)
  p <- ggplot(tsne.coords,aes(tSNE1,tSNE2,color=factor(realcluster))) +
    geom_point(size=0.2) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major=element_line(colour=NA),
          panel.grid.minor = element_blank(),
          axis.title =  element_text(size=10,face = "bold"),
          axis.text.x = element_blank(),  # 隐藏 x 轴的标签
          axis.text.y = element_blank(),  # 隐藏 y 轴的标签
          axis.ticks.x = element_blank(),  # 隐藏 x 轴的刻度
          axis.ticks.y = element_blank()   # 隐藏 y 轴的刻度
    ) +
    labs(title = name,color="Species")
  
  
  return(p)
  
}

truedata_tsne_plot <- function(data,realcluster,name){
  options( warn = -1 ) 
  impute_data <- data
  tsne.coords <- Rtsne(t(impute_data),pca=FALSE,perplexity = 30,theta=0.5,check_duplicates = F)$Y
  #rownames(tsne.coords[[k]]) <- rownames(t(data_true))
  colnames(tsne.coords) <- c("tSNE1","tSNE2")
  tsne.coords <- as.data.frame(tsne.coords)
  p <- ggplot(tsne.coords,aes(tSNE1,tSNE2,color=factor(realcluster))) +
    geom_point(size=0.2) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major=element_line(colour=NA),
          panel.grid.minor = element_blank(),
          axis.title =  element_text(size=10,face = "bold"),
    ) +
    labs(title = name,color=NULL)
  
  
  return(p)
}

a <- ggplot(data = mydata,aes(x=Method,y=value,group =Dropout,color=Dropout))+
  geom_point()+
  geom_line()+
  xlab('800×1000')+
  ylab("cel-cell")+
  theme_bw() +
  theme(
    panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    axis.title =  element_text(size=16,face = "bold"),
    axis.text = element_text(size=10), 
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16, face = "bold"), 
    #text = element_text(family = "STXihei"),
    #legend.position = c(.075,.075),
    #legend.box.background = element_rect(color="black")
  )+
  scale_color_manual(values = c("#84C7E1","#1A95C8","#0064B2","#254F82","#053061","#440154FF"))+
  # scale_y_continous(limits=c(min(y),max(y)),breaks=round(seq(min(y),max(y),length.out=3),2))
  scale_y_continuous(labels=scaleFUN)

plot__boxplot <- function(mydata){
  bp<-ggplot(data=mydata, aes(x=Method,y=value)) +     
      geom_boxplot(aes(fill = Method),outlier.shape = NA)+
      theme_bw()+
      geom_jitter(width = 0.2, alpha = 0.7)+
      scale_fill_manual(values = colors) + 
      labs(x = "",  y = "NMI")+  theme(
         panel.grid = element_blank(),)
  return(bp)
}


plot__trajectory <- function(cds,lab){
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  colorvec = getPalette(8)
  x = colorvec[1]
  colorvec[1] = colorvec[8]
  colorvec[8] = x
  cds$celltype = cds$H2228_to_H1975
  cds$celltype[is.na(cds$celltype)] = 'Na'
  p1 <- plot_cell_trajectory(cds, cell_size = 0.5, color_by = 'celltype')  + theme(legend.position = 'bottom',axis.text.x = element_blank(),  
                                                                                   axis.text.y = element_blank(),  
  )  + 
    scale_color_manual(values=colorvec ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,override.aes = list(size=2,alpha=1))) +
    theme(legend.spacing.x = unit(0, 'cm'),legend.spacing.y = unit(-0.5, 'cm'))+
    labs(color='') + xlab(NULL) + ylab(NULL)+ggtitle(lab)
  return(p1)
}


