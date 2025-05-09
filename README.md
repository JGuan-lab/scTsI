# scTsI
## 1. Introduction
scTsI is a two-stage algorithm for single cell gene expression imputation. In the first stage, scTsI finds neighboring cells and genes for initial imputation through KNN, and in the second stage, scTsI leverages bulk RNA-seq data as a constraint and uses ridge regression for adjusting the initial imputed values.

The simulated datasets and the analyzed experimental single-cell datasets can be accessed at: https://zenodo.org/uploads/13859610

## 2. Requirements:
    library(glmnet)
    source('scTsI.R')
## 3. Quick start
### 3.1 Prepare data
The inputs include scRNA-seq data and bulk RNA-seq data. 

    data_sc <- as.matrix(readRDS('data_sc.rds'))
    data_bulk <- as.matrix(readRDS('data_bulk.rds'))
### 3.2 Run scTsI
Threshold represents the threshold of high expression genes (the expression values greater than which would be unchanged by scTsI).
k1 and k2 represent the number of neighbor cells and neighbor genes in the first stage, respectively.

    # run demo.R
    # The annotations for the two-stage algorithm are in scTsI.R
    
    result <- scTsI(data_sc,threshold=0,data_bulk,k1=25,k2=25)
    write.csv(result, file="scTsI_result.csv")
