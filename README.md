# scTsI
We propose a two-stage scTsI algorithm for single cell gene expression imputation.  In the first stage, scTsI finds neighboring cells and genes for initial imputation through KNN, and in the second stage, scTsI leverages bulk RNA-seq data as a constraint and uses ridge regression for adjusting the initial imputed values.
