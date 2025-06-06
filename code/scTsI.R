library(mclust)
library(devtools)
library(fpc)
library(ngram)
library(FNN)
library(Matrix)
library(Metrics)
library(glmnet)

scTsI <- function(data_sc,threshold=0,data_bulk,k1=25,k2=25){
  
  data_sc[data_sc < threshold] <- 0
 
  dimensions <- dim(data_sc)
  m = dimensions[1]
  n = dimensions[2]
  ones_n <- rep(1, times = n)
  ones_m <- rep(1, times = m)
  diag_mn <- Matrix(diag(ones_m)/n, sparse = TRUE)
  
  # Construct the B matrix(Refer to the formula introduced)
  B = kronecker(t(ones_n),diag_mn)
  
  # twice knn
  data_before <- t(data_sc)
  data_before_tt <- data_sc
  
  # cell
  idx <-get.knn(data_before,  k = k1) 
  
  for (i in 1:n) {
    data <- colMeans(data_before[idx$nn.index[i,],])
    zero_positions <- which(data_before[i,] == 0)
    data_before[i,zero_positions] = data[zero_positions]
  }
  
  # gene
  idx <-get.knn(data_before_tt,  k = k2) 
  
  for (i in 1:m) {
    data <- colMeans(data_before_tt[idx$nn.index[i,],])
    zero_positions <- which(data_before_tt[i,] == 0)
    data_before_tt[i,zero_positions] = data[zero_positions]
  }
  # first stage
  data_knn <- (t(data_before)+data_before_tt)/2
  
  # Vectorize
  combined_vector <- c(data_sc)
  
  n = length(combined_vector)
  
  # Construct the M matrix
  M <- sparseMatrix(i = 1:n, j = 1:n, x = 1)
  nonzero_indices <- which(combined_vector != 0)
  zero_indices <- which(combined_vector == 0)
  x = c(zero_indices,nonzero_indices)
  M1 = M[x,]
  M_ = M[,x]
  B_result = B %*%M_
  vector_length <- length(zero_indices)
  x <- B_result[,1:vector_length]
  y <- as.matrix(data_bulk - rowMeans(data_knn))
  
  #  lower bound for output
  combined_knn <- c(data_knn)  
  combined_knn = M1 %*% combined_knn
  combined_knn = as(combined_knn,'matrix')
  limit_vector = -combined_knn
  
  # output y2
  fit = glmnet(x=x, y=y, alpha = 0,family = "gaussian",lower.limits = limit_vector,intercept = FALSE,nlambda = 10)
  
  # The median result is chosen by default
  # The beta result is in order
  B <- as.matrix(fit$beta[,5])
  dimensions <- dim(data_knn)
  m = dimensions[1]
  n = dimensions[2]
  u = numeric((m*n-vector_length))
  Y = c(B,u)
  Y_result <- as(M_%*%Y, "matrix")

  # second stage
  result <- matrix(Y_result, nrow = m, ncol = n)

  
  data_result <- result+data_knn
  # Output the final result
  data_result[data_result < 0] <- 0
  return(data_result)
}


find_threshold <-function(data_sc){

  vec <- as.vector(data)
  unique_sorted <- sort(unique(vec))
  q_vals <- quantile(unique_sorted, probs = c(0.001,0.002,0.005,0.01,0.02,0.05, 0.10, 0.15, 0.20))
  vec_sorted <- sort(vec)
  n <- length(vec_sorted)

  data.frame(
    quantile = names(q_vals),
    value = round(as.numeric(q_vals), 4),
  )
}
