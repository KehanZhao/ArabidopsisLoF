# Expression_matrix.R
# Constructs the expression matrix aligned to the LoF matrix accessions and genes.
# NAs (missing expression data) are removed row/column-wise, and NaNs are
# replaced with the column minimum before saving.

# Load transcriptome data -------------------------------------------------------

load("/path/to/TG_data_20180606.Rdata")

Deletion_gene_name       <- colnames(LoF_matrix_collapsed)
Deletion_accession_name  <- rownames(LoF_matrix_collapsed)

# Build expression matrix -------------------------------------------------------
# Iterates over every (accession, gene) pair and looks up the batch-corrected
# log2 expression value from the TG dataset.

n_acc  <- length(Deletion_accession_name)
n_gene <- length(Deletion_gene_name)

Expression_matrix <- matrix(NA, n_acc, n_gene)

for (u in seq_len(n_acc * n_gene)) {
  Gene_series      <- ceiling(u / n_acc)
  Accession_series <- u - n_acc * (Gene_series - 1)

  Gene_name      <- Deletion_gene_name[Gene_series]
  Accession_name <- Deletion_accession_name[Accession_series]

  TG_Gene_series      <- which(TG.genes$name == Gene_name)
  TG_Accession_series <- which(TG.meta$index == Accession_name)

  Gene_expression_level <- TG.genes$d_log2_batch[TG_Gene_series, TG_Accession_series]

  if (length(Gene_expression_level) < 1) Gene_expression_level <- NA
  if (length(Gene_expression_level) > 1) Gene_expression_level <- mean(Gene_expression_level)

  Expression_matrix[Accession_series, Gene_series] <- Gene_expression_level
}

colnames(Expression_matrix) <- Deletion_gene_name
rownames(Expression_matrix) <- Deletion_accession_name

save(Expression_matrix, file = "Expression_matrix")

# Remove all-NA rows and columns ------------------------------------------------

Expression_matrix_NAremoved <- Expression_matrix[
  apply(Expression_matrix, 1, function(y) any(!is.na(y))),
]
Expression_matrix_NAremoved <- Expression_matrix_NAremoved[,
  apply(Expression_matrix_NAremoved, 2, function(y) any(!is.na(y)))
]

# Replace NaN with column minimum -----------------------------------------------
# NaN can arise from log2(0) or 0/0 in upstream batch correction.
# Each NaN is replaced with the minimum non-NaN value in its column.

n_row <- nrow(Expression_matrix_NAremoved)
n_col <- ncol(Expression_matrix_NAremoved)

Expression_matrix_NAremoved_NaNreplaced <- matrix(0, n_row, n_col)

for (u in seq_len(n_col)) {
  col_vals <- Expression_matrix_NAremoved[, u]
  valid    <- col_vals[!is.nan(col_vals)]
  col_vals[is.nan(col_vals)] <- min(valid)
  Expression_matrix_NAremoved_NaNreplaced[, u] <- col_vals
}

colnames(Expression_matrix_NAremoved_NaNreplaced) <- colnames(Expression_matrix_NAremoved)
rownames(Expression_matrix_NAremoved_NaNreplaced) <- rownames(Expression_matrix_NAremoved)

save(Expression_matrix_NAremoved_NaNreplaced, file = "Expression_matrix_NAremoved_NaNreplaced")
