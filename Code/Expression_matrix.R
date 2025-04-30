# This script generates the expression matrix used for LoF-expression association testing

# Load transcriptome data
load("/path/to/TG_data_20180606.Rdata")

Deletion_gene_name<-colnames(LoF_matrix_collapsed)
Deletion_accession_name<-rownames(LoF_matrix_collapsed)

Expression_matrix <- matrix(NA,length(Deletion_accession_name),length(Deletion_gene_name))
for (u in 1:(length(Deletion_accession_name)*length(Deletion_gene_name))){
  Gene_series<-ceiling(u/length(Deletion_accession_name))
  Accession_series<-u-length(Deletion_accession_name)*(Gene_series-1)
  Gene_name<-Deletion_gene_name[Gene_series]
  Accession_name<-Deletion_accession_name[Accession_series]
  TG_Gene_series<-which(TG.genes$name==Gene_name)
  TG_Accession_series<-which(TG.meta$index==Accession_name)
  Gene_expression_level<-TG.genes$d_log2_batch[TG_Gene_series,TG_Accession_series]
  if (length(Gene_expression_level)<1) Gene_expression_level<-NA
  if (length(Gene_expression_level)>1) Gene_expression_level<-mean(Gene_expression_level)
  Expression_matrix[Accession_series,Gene_series]<-Gene_expression_level
}

colnames(Expression_matrix)<-colnames(LoF_matrix_collapsed)
rownames(Expression_matrix)<-rownames(LoF_matrix_collapsed)

save(Expression_matrix, file="Expression_matrix")

# Remove NA values 
Expression_matrix_NAremoved <- Expression_matrix[apply(Expression_matrix,1,function(y) any(!is.na(y))),]
Expression_matrix_NAremoved <- Expression_matrix_NAremoved[,apply(Expression_matrix_NAremoved,2,function(y) any(!is.na(y)))]

# Replace NaN in Expression matrix with lowest value of each column

Expression_matrix_NAremoved_NaNreplaced <- matrix(0,nrow(Expression_matrix_NAremoved),ncol(Expression_matrix_NAremoved))
for (u in 1:ncol(Expression_matrix_NAremoved)){
  c<-Expression_matrix_NAremoved[,u]
  No_NaN <- c[!is.nan(c)]
  Lowest_value <- min(No_NaN)
  c[which(is.nan(c))] <- Lowest_value
  Expression_matrix_NAremoved_NaNreplaced[,u]<-c
}

colnames(Expression_matrix_NAremoved_NaNreplaced)<-colnames(Expression_matrix_NAremoved)
rownames(Expression_matrix_NAremoved_NaNreplaced)<-rownames(Expression_matrix_NAremoved)

save(Expression_matrix_NAremoved_NaNreplaced, file = "Expression_matrix_NAremoved_NaNreplaced")