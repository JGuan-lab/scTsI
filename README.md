# scTsI
We propose a two-stage scTsI algorithm for single cell gene expression imputation.  In the first stage, scTsI finds neighboring cells and genes for initial imputation through KNN, and in the second stage, scTsI leverages bulk RNA-seq data as a constraint and uses ridge regression for adjusting the initial imputed values.

## 1.Introduction
scTsI is a two-stage algorithm for single cell gene expression imputation.  In the first stage, scTsI finds neighboring cells and genes for initial imputation through KNN, and in the second stage, scTsI leverages bulk RNA-seq data as a constraint and uses ridge regression for adjusting the initial imputed values.

## 2.Requirements:
  library(glmnet)
  source('scTsI.R')

## 3.Quick start
### 3.1 Prepare data
    The inputs include scRNA-seq data and bulk RNA-seq data.
    path = ''
    setwd(path)
    data_sc<as.matrix(readRDS(paste0(path,'data_sc.rds'))
    data_bulk<as.matrix(readRDS(paste0(path,'data_bulk.rds'))
### 3.2 Run scTsI
scTsI selects the result of the first λ of glmnet by default. 
The better results can be obtained by calculating other λ results, adjusting the number of λ, or changing the range of λ.

    result <- scTsI(data_sc,data_bulk)
