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









##########################################plot#########################################


tsne_plot <- function(data,realcluster,name){
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

plot__trajectory <- function(cds,lab){
  
  cds$celltype = cds$H2228_to_H1975
  cds$celltype[is.na(cds$celltype)] = 'Na'
  p1 <- plot_cell_trajectory(cds, cell_size = 0.5, color_by = 'celltype')  + theme(legend.position = 'bottom',axis.text.x = element_blank(),  # 隐藏 x 轴的标签
                                                                                   axis.text.y = element_blank(),  # 隐藏 y 轴的标签
  )  + 
    scale_color_manual(values=colorvec ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,override.aes = list(size=2,alpha=1))) +
    theme(legend.spacing.x = unit(0, 'cm'),legend.spacing.y = unit(-0.5, 'cm'))+
    labs(color='') + xlab(NULL) + ylab(NULL)+ggtitle(lab)
  return(p1)
}


# 创建一个空列表，用于保存每个图像
p_list <- list()
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colorvec = getPalette(8)
x = colorvec[1]
colorvec[1] = colorvec[8]
colorvec[8] = x



