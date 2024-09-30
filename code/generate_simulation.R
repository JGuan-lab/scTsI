library(splatter)
library(scater)
rm(list=ls())


params = newSplatParams()
params = setParams(params, list(batchCells = 1000, nGenes = 5000,group.prob = c(0.20, 0.35, 0.45),de.prob = c(0.045, 0.045, 0.045),de.facLoc = 0.1,de.facScale = 0.4))
sim = splatSimulateGroups(params,dropout.type = "experiment",dropout.shape = -1,dropout.mid =3,seed = 6)


# genereate the cpm levels of the true simulation data5
data_true = data_true = edgeR::cpm(sim@assays@data$TrueCounts)

count = sim@assays@data$counts
DEFacGroup1 = sim@rowRanges@elementMetadata$DEFacGroup1
DEFacGroup2 = sim@rowRanges@elementMetadata$DEFacGroup2
DEFacGroup3 = sim@rowRanges@elementMetadata$DEFacGroup3


data_dropout = data_true
data_dropout[sim@assays@data@listData$Dropout] = 0
data_dropout=as.matrix(data_dropout)

# generate the bulk RNAseq data
data_bulk = data.frame(val = rowMeans(data_true))

# define the data list for the simulation data
# indcluding: data_true: true data
#          data_dropout: dropout data
#             data_bluk: bulk data
#      percentage_zeros: dropout rate
#                 group: the group label

group =colData(sim)@listData$Group





#stratified sampling

jwcol <-function(n,m,data_sc,data_true,group){
  # n is the dimension of change
  # m is the original dimension
  subset_list <- list()
  
  for (category in unique(group)) {

    indices <- which(group == category)
    
    subset_list[[as.character(category)]] <- data_sc[,indices, drop = FALSE]
  }
  
  # Set the sampling ratio
  sampling_ratio <- n/m
  sampled_list <- list()
  for (category in names(subset_list)) {
    subset_data <- subset_list[[category]]
    
    num_samples <- round(ncol(subset_data) * sampling_ratio)
    
    sampled_data <- subset_data[,sample(1:ncol(subset_data), size = num_samples, replace = FALSE) , drop = FALSE]
    
    sampled_list[[category]] <- sampled_data
  }
  
  data_new <- do.call(cbind, sampled_list)
  data_dro <- data_sc[, colnames(data_sc) %in% colnames(data_new), drop = FALSE]
  data_true_new <- data_true[, colnames(data_true) %in% colnames(data_new), drop = FALSE]
  
  data = list(data_new = data_dro,data_true_new = data_true_new)
  return(data)
}

jwrow <-function(n,data_sc,data_true,DEFacGroup1,DEFacGroup2,DEFacGroup3){
  #  n is the dimension of change
  
  Group1 <- which(DEFacGroup1 != 1)
  Group2 <- which(DEFacGroup2 != 1)
  Group3 <- which(DEFacGroup3 != 1)

  merged_vector <- union(Group1,Group2)
  merged_vector <- union(merged_vector,Group3)
  total_nums <- length(Group1)+length(Group2)+length(Group3)
  ratio_1 <- length(Group1) / total_nums
  ratio_2 <- length(Group2) / total_nums
  ratio_3 <- length(Group3) / total_nums
  
  intersect_positions1_2 <- intersect(Group1,Group2)
  intersect_positions1_3 <- intersect(Group1,Group3)
  intersect_positions2_3 <- intersect(Group2,Group3)
  intersect_positions1_2_3 <- intersect(intersect_positions1_2,Group3)
  # if(length(intersect_positions1_2_3){
  #   
  # }else{
  total_intersect <- union(intersect_positions1_2,intersect_positions1_3)
  total_intersect <- union(total_intersect,intersect_positions2_3)
  if(n<length(total_intersect)){
    
  }else{
    
    processed_vector_1 <- Group1[!Group1 %in% (union(intersect_positions1_2,intersect_positions1_3))]
    processed_vector_2 <- Group2[!Group2 %in% (union(intersect_positions1_2,intersect_positions2_3))]
    processed_vector_3 <- Group3[!Group3 %in% (union(intersect_positions1_3,intersect_positions2_3))]
    
    if(n<length(merged_vector)){
      m=n-length(total_intersect)
      sample_1 <- sample(processed_vector_1, round(ratio_1 * m),replace = FALSE)
      sample_2 <- sample(processed_vector_2, round(ratio_2 * m),replace = FALSE)
      sample_3 <- sample(processed_vector_3, round(ratio_3 * m),replace = FALSE)
      
      combined_positions <- sort(c(sample_1,sample_2,sample_3,total_intersect))
      data_new <- data_sc[combined_positions, , drop = FALSE]
      data_true_new <- data_true[combined_positions, , drop = FALSE]
      
      data = list(data_new = data_new,data_true_new = data_true_new)
      return(data)
     
    }else{

      m=n-length(merged_vector)
      
      selected_rows_indices <- setdiff(1:nrow(data_sc), merged_vector)
      # selected_rows <- data_sc[selected_rows_indices, , drop = FALSE]
      
      selected_indices <- sample(selected_rows_indices, m, replace = FALSE) 
      
      combined_all <- c(merged_vector,selected_indices)
      data_new <- data_sc[combined_all, , drop = FALSE]  
      data_true_new <- data_true[combined_all, , drop = FALSE]
      
      data = list(data_new = data_new,data_true_new = data_true_new)
      return(data)
    }
    
  }
  # }
}
