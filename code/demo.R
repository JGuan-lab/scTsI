#######Read the data#############
path="" 
setwd(path)  
data_sc <- as.matrix(readRDS(paste0(path,'data_sc.rds'))
data_bulk <- as.matrix(readRDS(paste0(path,'data_bulk.rds'))
###########Run scTsI#######
source('scTsI.R')
result <- scTsI(data_sc,threshold=0,data_bulk,k1=25,k2=25)
