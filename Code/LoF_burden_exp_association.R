# LoF_burden_exp_association.R
# Tests for associations between per-gene LoF burden and genome-wide gene
# expression using the EMMAX mixed model (via cpgen::cGWAS.emmax), which
# accounts for population structure through a kinship matrix.
#
# NOTE: cGWAS.emmax was tested under R 3.6.1. It may fail under R >= 4.2.

library(rhdf5)
library(cpgen)

# Load and filter expression matrix --------------------------------------------
# Remove genes with zero variance (constant expression across accessions),
# as these cannot be tested.

load("/path/to/Expression_matrix_NAremoved_NaNreplaced")

zero_var_cols <- which(apply(Expression_matrix_NAremoved_NaNreplaced, 2, var) == 0)
if (length(zero_var_cols) > 0) {
  Expression_matrix_NAremoved_NaNreplaced <- Expression_matrix_NAremoved_NaNreplaced[, -zero_var_cols]
}

# Load and filter LoF matrix ---------------------------------------------------
# Retain only genes with LoF allele count >= 3 across accessions.

LoF_matrix_collapsed <- read.csv("/path/to/Supplementary Table S4.csv", row.names = 1)

Allele_count <- colSums(LoF_matrix_collapsed)
LoF_matrix_collapsed <- LoF_matrix_collapsed[, Allele_count >= 3]

# Subset to accessions present in the expression matrix
LoF_matrix_collapsed_filter <- LoF_matrix_collapsed[
  rownames(Expression_matrix_NAremoved_NaNreplaced), ,
  drop = FALSE
]

# Load and subset kinship matrix -----------------------------------------------

kinship    <- h5read(file = "/path/to/kinship_ibs_mac5.hdf5", name = "kinship")
accessions <- h5read(file = "/path/to/kinship_ibs_mac5.hdf5", name = "accessions")

colnames(kinship) <- accessions
rownames(kinship) <- accessions

kinship_sub <- kinship[
  rownames(LoF_matrix_collapsed_filter),
  rownames(LoF_matrix_collapsed_filter)
]

# Inverse normal transformation (INT) helper -----------------------------------
# Applied optionally; see comment in the EMMAX loop below.

INT <- function(x) {
  valid_x <- x[!is.na(x)]
  ranked  <- rank(valid_x, ties.method = "average")
  int_transformed <- qnorm((ranked - 0.5) / length(valid_x))

  result <- rep(NA, length(x))
  result[!is.na(x)] <- int_transformed
  return(result)
}

# EMMAX association testing ----------------------------------------------------
# Rows = expression traits; columns = LoF genes.
# To apply INT to expression phenotypes, replace y= argument with:
#   y = INT(Expression_matrix_NAremoved_NaNreplaced[, i])

n_exp <- ncol(Expression_matrix_NAremoved_NaNreplaced)
n_lof <- ncol(LoF_matrix_collapsed_filter)

pvalue_matrix_Exp      <- matrix(0, n_exp, n_lof)
beta_matrix_Exp        <- matrix(0, n_exp, n_lof)
marker_variance_Exp    <- vector(length = n_exp)
residual_variance_Exp  <- vector(length = n_exp)

for (i in seq_len(n_exp)) {
  mod <- cGWAS.emmax(
    y = Expression_matrix_NAremoved_NaNreplaced[, i],
    M = as.matrix(LoF_matrix_collapsed_filter),
    A = as.matrix(kinship_sub)
  )
  pvalue_matrix_Exp[i, ]   <- mod$p_value
  beta_matrix_Exp[i, ]     <- mod$beta
  marker_variance_Exp[i]   <- mod$marker_variance
  residual_variance_Exp[i] <- mod$residual_variance
}

rownames(pvalue_matrix_Exp) <- colnames(Expression_matrix_NAremoved_NaNreplaced)
colnames(pvalue_matrix_Exp) <- colnames(LoF_matrix_collapsed_filter)
rownames(beta_matrix_Exp)   <- colnames(Expression_matrix_NAremoved_NaNreplaced)
colnames(beta_matrix_Exp)   <- colnames(LoF_matrix_collapsed_filter)

# Save results -----------------------------------------------------------------

save(pvalue_matrix_Exp,     file = "pvalue_matrix_Exp")
save(beta_matrix_Exp,       file = "beta_matrix_Exp")
save(marker_variance_Exp,   file = "marker_variance_Exp")
save(residual_variance_Exp, file = "residual_variance_Exp")
