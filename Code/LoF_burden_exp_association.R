# Load and filter expression matrix --------------------------------------------

load("/path/to/Expression_matrix_NAremoved_NaNreplaced")
delete<-vector()
for (u in 1:ncol(Expression_matrix_NAremoved_NaNreplaced)){
  if (var(Expression_matrix_NAremoved_NaNreplaced[,u])==0) delete<-c(delete,u)
}
Expression_matrix_NAremoved_NaNreplaced<-Expression_matrix_NAremoved_NaNreplaced[,-delete]

# Load and filter LoF matrix ---------------------------------------------------

LoF_matrix_collapsed <- read.csv("/path/to/S1_LoF_matrix_collapsed.csv")
LoF_matrix_collapsed <- LoF_matrix_collapsed[,-1]
Allele_count<-colSums(LoF_matrix_collapsed)
LoF_matrix_collapsed<-LoF_matrix_collapsed[,-which(Allele_count<3)]
LoF_matrix_collapsed_filter<-LoF_matrix_collapsed[rownames(Expression_matrix_NAremoved_NaNreplaced),]

# Load and filter kinship matrix -----------------------------------------------

library("rhdf5")
kinship <- h5read(file = "/path/to/kinship_ibs_mac5.hdf5", name = "kinship")
accessions <- h5read(file = "/path/to/kinship_ibs_mac5.hdf5", name = "accessions")

colnames(kinship) <- accessions
rownames(kinship) <- accessions

kinship_sub <- kinship[rownames(LoF_matrix_collapsed_filter),rownames(LoF_matrix_collapsed_filter)]

# The function to perform INT if applied ---------------------------------------

INT <- function(x) {
  valid_x <- x[!is.na(x)]  # Remove NA values
  ranked <- rank(valid_x, ties.method = "average")
  int_transformed <- qnorm((ranked - 0.5) / length(valid_x))
  
  # Create a result vector matching the original input size
  result <- rep(NA, length(x))
  result[!is.na(x)] <- int_transformed
  
  return(result)
}

# Perform LoF burden tests with expression using EMMAX -------------------------

pvalue_matrix_Exp <- matrix(0,ncol(Expression_matrix_NAremoved_NaNreplaced),ncol(LoF_matrix_collapsed_filter))
beta_matrix_Exp <- matrix(0,ncol(Expression_matrix_NAremoved_NaNreplaced),ncol(LoF_matrix_collapsed_filter))
marker_variance_Exp <- vector(length = ncol(Expression_matrix_NAremoved_NaNreplaced))
residual_variance_Exp <- vector(length = ncol(Expression_matrix_NAremoved_NaNreplaced))
library(cpgen)
# This was run using R version 3.6.1. cGWAS.emmax may fail with R version 4.2 or later.
for (i in 1:ncol(Expression_matrix_NAremoved_NaNreplaced)){
  mod <- cGWAS.emmax(y=Expression_matrix_NAremoved_NaNreplaced[,i], M=as.matrix(LoF_matrix_collapsed_filter), A=as.matrix(kinship_sub)) # if applying INT, y=INT(Expression_matrix_NAremoved_NaNreplaced[,i])
  pvalue_matrix_Exp[i,]<-mod$p_value
  beta_matrix_Exp[i,]<-mod$beta
  marker_variance_Exp[i]<-mod$marker_variance
  residual_variance_Exp[i]<-mod$residual_variance
}
rownames(pvalue_matrix_Exp)<-colnames(Expression_matrix_NAremoved_NaNreplaced)
colnames(pvalue_matrix_Exp)<-colnames(LoF_matrix_collapsed_filter)
rownames(beta_matrix_Exp)<-colnames(Expression_matrix_NAremoved_NaNreplaced)
colnames(beta_matrix_Exp)<-colnames(LoF_matrix_collapsed_filter)

save(pvalue_matrix_Exp, file = "pvalue_matrix_Exp")
save(beta_matrix_Exp, file = "beta_matrix_Exp")
save(marker_variance_Exp, file = "marker_variance_Exp")
save(residual_variance_Exp, file = "residual_variance_Exp")