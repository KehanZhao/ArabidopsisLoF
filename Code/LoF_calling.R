# Load necessary libraries -----------------------------------------------------
library(vcfR)
library(data.table)
library(stringr)

# Calling LoF variants from structural variance (deletion) ---------------------

# Load gene annotation file 
# Read in gff annotation file and give it useful column names 
gff<-fread("/path/to/TAIR10_GFF3_genes.gff")
colnames(gff)<-c("chr","source","type","start","stop","X","direction","X2","info")

# Subest to only coding regions (for now we only care about deletions that overlap with coding sequences)
CDS<-gff[type=="CDS"]
# Remove Chr from chromosome names to match the vcf
CDS$chr<-gsub("Chr","",CDS$chr)

# Use some grep type commands to extract the gene model and gene name from the info column provided in the annotaiton gff format
CDS$model<-gsub("Parent=(.+),.+","\\1",CDS$info)
CDS$gene<-gsub("\\..+","",CDS$model)
CDS_split<-split(CDS, by="chr")

# Subset to gene regions
gene<-gff[type=="gene"]
gene$chr<-gsub("Chr","",gene$chr)
gene$model<-gsub("ID=(.+);.+","\\1",gene$info)

# Load deletion data 
vcf<-read.vcfR("/path/to/Deletions.vcf.gz")
variants<-data.table(vcf@fix)

# Clean genotype calls
gt<-data.table(vcf@gt[,-1])
gt_processed<-apply(gt,2,function(col){
  clean <- gsub("/.+","",col)
  col<-gsub("\\.","0",clean)
  return(col)
})

# Looping through variants dataframe to:
# 1. Identify which genes/models are impacted by a variant (frameshifts, deletions).
# 2. Record information about these impacts (gene, gene model, frameshift?, large deletion?).
# 3. Summarize how many isoforms/models are affected per gene.
# 4. Calculate the percentage of the CDS affected by the variant.
# 5. Output the final annotated dataframe as variants.

# Pre-allocate result vectors
lof_new_vector <- rep(NA, nrow(variants))
frameshift_vector <- rep(NA, nrow(variants))
large_deletion_vector <- rep(NA, nrow(variants))

pb <- txtProgressBar(min = 0, max = nrow(variants), style = 3)
for (i in seq_len(nrow(variants))) {
  row <- unlist(variants[i, ])
  chr <- row[1]
  pos <- as.numeric(row[2])
  end <- as.numeric(gsub("END=(\\d+);.+", "\\1", row[8]))
  
  # Extract CDS entries for this chromosome
  CDS_chr <- CDS_split[[chr]]
  
  # Find overlaps with CDS models
  overlaps <- which(
    (CDS_chr$start <= pos & CDS_chr$stop >= pos) |
      (CDS_chr$start <= end & CDS_chr$stop >= end) |
      (pos <= CDS_chr$start & end >= CDS_chr$stop)
  )
  
  overlapping_models <- unique(CDS_chr$model[overlaps])
  
  real_lof_new <- frameshift <- large_deletion <- character(0)
  
  for (candidate in overlapping_models) {
    gene_CDS <- unlist(apply(CDS[model == candidate], 1, function(x) x[4]:x[5]))
    gene_CDS <- sort(gene_CDS)
    middle_region <- gene_CDS[-c(1:floor(0.05 * length(gene_CDS)), floor(0.95 * length(gene_CDS)):length(gene_CDS))]
    
    affected <- gene_CDS %in% pos:end
    in_middle <- middle_region %in% pos:end
    
    if ((sum(affected) %% 3 != 0 || mean(affected) >= 0.1) && any(in_middle)) {
      real_lof_new <- c(real_lof_new, candidate)
    }
    if ((sum(affected) %% 3 != 0) && any(in_middle)) {
      frameshift <- c(frameshift, candidate)
    }
    if (mean(affected) >= 0.1 && any(in_middle)) {
      large_deletion <- c(large_deletion, candidate)
    }
  }
  
  # Assign to vectors
  lof_new_vector[i] <- ifelse(length(real_lof_new) > 0, paste(real_lof_new, collapse = ","), NA)
  frameshift_vector[i] <- ifelse(length(frameshift) > 0, paste(frameshift, collapse = ","), NA)
  large_deletion_vector[i] <- ifelse(length(large_deletion) > 0, paste(large_deletion, collapse = ","), NA)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Assign results to variants
variants <- as.data.frame(variants)
variants$lof_model_5_95 <- lof_new_vector
variants$frameshift <- frameshift_vector
variants$large_deletion <- large_deletion_vector

# Extract unique gene names (first 9 characters of model names)
variants$gene <- sapply(lof_new_vector, function(x) {
  if (!is.na(x)) {
    genes <- unique(substr(unlist(strsplit(x, ",")), 1, 9))
    return(paste(genes, collapse = ","))
  } else {
    return(NA)
  }
})

# Add model counts and percent affected
All_gene_model<-read.delim("/path/to/TAIR10_all_gene_models.txt")
variants$model_number <- NA
variants$affected_model_number <- NA
variants$affected_model_percentage <- NA

pb <- txtProgressBar(min = 0, max = nrow(variants), style = 3)
for (i in seq_len(nrow(variants))) {
  gene <- variants$gene[i]
  if (!is.na(gene)) {
    affected_models <- variants$lof_model_5_95[i]
    gene_list <- unlist(strsplit(gene, ","))
    
    model_counts <- sapply(gene_list, function(g) {
      sum(grepl(g, All_gene_model$This.file.lists.ALL.TAIR10.gene.models..both.representative_gene.models.and.all.other.splice.variants..a.total.of.41.671.models.))
    })
    
    affected_counts <- sapply(gene_list, function(g) {
      str_count(affected_models, g)
    })
    
    percents <- round(affected_counts / model_counts, 3)
    
    variants$model_number[i] <- paste(model_counts, collapse = ",")
    variants$affected_model_number[i] <- paste(affected_counts, collapse = ",")
    variants$affected_model_percentage[i] <- paste(percents, collapse = ",")
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Compute % of CDS affected
variants$affected_CDS_percentage <- NA
affected_rows <- which(!is.na(variants$lof_model_5_95))

for (row_n in affected_rows) {
  chr <- variants[row_n, 1]
  pos <- as.numeric(variants[row_n, 2])
  end <- as.numeric(gsub("END=(\\d+);.+", "\\1", variants[row_n, 8]))
  models <- strsplit(variants[row_n, 9], ",")[[1]]
  
  length_percents <- sapply(models, function(m) {
    gene_CDS <- unlist(apply(CDS[model == m], 1, function(x) x[4]:x[5]))
    sum(gene_CDS %in% pos:end) / length(gene_CDS)
  })
  
  variants$affected_CDS_percentage[row_n] <- paste(round(length_percents, 3), collapse = ",")
}

# Calling LoF variants from the 1001 Genome Project ----------------------------

# Load high impact VCF data
# High impact variants were extracted from the 1001 Genome VCF using command line: zgrep "HIGH" 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz > high_impact.txt
high_impact <- read.delim("/path/to/high_impact.txt", header = FALSE)
header <- read.delim("/path/to/1001GenomeVCF_header.txt", header = FALSE) # This file contains the header row from 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz 
colnames(high_impact) <- header$V1

# Split data into info and genotype parts
high_impact_info <- high_impact[, 1:9]
high_impact_gt <- high_impact[, -c(1:9)]

# Clean genotype calls: extract first character and replace '.' with '0'
high_impact_gt_processed <- apply(high_impact_gt, 2, function(col) {
  clean <- substr(col, 1, 1)
  gsub("\\.", "0", clean)
})
high_impact_gt_processed <- as.data.frame(high_impact_gt_processed)

# Handle multiple alternate alleles
multiple_changes <- grep(",", high_impact_info$ALT)

for (i in seq_along(multiple_changes)) {
  row <- multiple_changes[i]
  changes <- strsplit(as.character(high_impact_info[row, 5]), ",")[[1]]
  
  # Add new rows for each alternate allele (starting from the 2nd)
  for (w in 2:length(changes)) {
    new_row_vcf <- high_impact_info[row, ]
    new_row_vcf$ALT <- changes[w]
    high_impact_info <- rbind(high_impact_info, new_row_vcf)
    
    new_row_gt <- rep(0, 1135)
    new_row_gt[which(high_impact_gt_processed[row, ] == w)] <- w
    high_impact_gt_processed <- rbind(high_impact_gt_processed, new_row_gt)
  }
  
  # Update original row for ALT and genotype
  high_impact_info[row, 5] <- changes[1]
  original_gt <- rep(0, 1135)
  original_gt[which(high_impact_gt_processed[row, ] == 1)] <- 1
  high_impact_gt_processed[row, ] <- original_gt
  
  message(sprintf("Processed multi-ALT row %d of %d", i, length(multiple_changes)))
}

# Assign genotype for each row
for (i in seq_len(nrow(high_impact_gt_processed))) {
  gt <- unique(as.character(high_impact_gt_processed[i, ]))
  gt <- setdiff(gt, "0")
  if (length(gt) > 0) {
    high_impact_info[i, 10] <- gt
  }
  if (i %% 1000 == 0) message(sprintf("Annotated GT row %d", i))
}
colnames(high_impact_info)[10] <- "gt"

# Parse gene models from annotation
genes_model <- apply(high_impact_info, 1, function(row) {
  variant_vec <- strsplit(row[8], ",")[[1]]
  LoF <- variant_vec[grepl("HIGH", variant_vec)]
  LoF_real_p <- gsub(".*CODING|", "", LoF)
  LoF_real_p <- gsub(").*", "", LoF_real_p)
  gt_info <- substr(LoF_real_p, nchar(LoF_real_p), nchar(LoF_real_p))
  gt <- row[10]
  LoF_model <- substr(LoF_real_p, 2, 12)
  match <- which(gt_info == gt)
  if (length(match) != 0) {
    LoF_model <- LoF_model[match]
  } else if (is.na(gt)) {
    LoF_model <- "No_LoF"
  } else {
    LoF_model <- NA
  }
  paste(LoF_model, collapse = ",")
})
high_impact_info[, 11] <- genes_model
colnames(high_impact_info)[11] <- "gene_model"

# Assign affected genes
for (i in seq_along(genes_model)) {
  if (!is.na(genes_model[i]) && genes_model[i] != "NA") {
    models <- strsplit(genes_model[i], ",")[[1]]
    unique_genes <- unique(substr(models, 1, 9))
    high_impact_info[i, 12] <- paste(unique_genes, collapse = ",")
  }
}
colnames(high_impact_info)[12] <- "gene"

# Calculate affected percentages
for (i in seq_len(nrow(high_impact_info))) {
  gene <- high_impact_info$gene[i]
  if (!is.na(gene)) {
    affected_model <- high_impact_info$gene_model[i]
    if (nchar(gene) == 9) {
      model_number <- sum(grepl(gene, All_gene_model[[1]]))
      affected_number <- str_count(affected_model, gene)
      high_impact_info[i, 13:15] <- c(model_number, affected_number, affected_number / model_number)
    }
    if (nchar(gene) > 9) {
      genes <- strsplit(gene, ",")[[1]]
      result <- sapply(genes, function(g) {
        total <- sum(grepl(g, All_gene_model[[1]]))
        affected <- str_count(affected_model, g)
        percent <- affected / total
        c(total, affected, percent)
      })
      high_impact_info[i, 13:15] <- lapply(1:3, function(j) paste(result[j, ], collapse = ","))
    }
  }
  if (i %% 1000 == 0) message(sprintf("Calculated percentages row %d", i))
}
colnames(high_impact_info)[13:15] <- c("model_number", "affected_model_number", "affected_percentage")

# Determine LoF type
LoF_type <- apply(high_impact_info, 1, function(row) {
  variant_vec <- strsplit(row[8], ",")[[1]]
  LoF <- variant_vec[grepl("HIGH", variant_vec)]
  type <- gsub(".*=", "", gsub("\\(.*", "", LoF))
  LoF_real_p <- gsub(".*CODING|", "", LoF)
  LoF_real_p <- gsub(").*", "", LoF_real_p)
  gt_info <- substr(LoF_real_p, nchar(LoF_real_p), nchar(LoF_real_p))
  gt <- row[10]
  match <- which(gt_info == gt)
  if (length(match) != 0) {
    type[match]
  } else if (is.na(gt)) {
    "No_LoF"
  } else {
    NA
  }
})

# Make sure each element is a character string, even if there are multiple types
LoF_type_collapsed <- sapply(LoF_type, function(x) paste(x, collapse = ";"))

high_impact_info[, 16] <- LoF_type_collapsed
colnames(high_impact_info)[16] <- "LoF_type"

# Determine gene model from 5–95% CDS region
genes_model_5_95 <- apply(high_impact_info, 1, function(row) {
  chr <- as.numeric(row[1])
  pos <- as.numeric(row[2])
  original <- as.character(row[4])
  new <- strsplit(as.character(row[5]), ",")[[1]]
  old_end <- pos + nchar(original) - 1
  new_end <- pos + nchar(new) - 1
  CDS_chr <- CDS_split[[chr]]
  
  overlaps <- which(
    (CDS_chr$start <= pos & CDS_chr$stop >= pos) |
      (CDS_chr$start <= old_end & CDS_chr$stop >= old_end) |
      (CDS_chr$start <= new_end & CDS_chr$stop >= new_end)
  )
  candidates <- unique(CDS_chr$model[overlaps])
  real_lof <- vector()
  for (candidate in candidates) {
    gene_CDS <- unlist(apply(CDS[model == candidate], 1, function(x) x[4]:x[5]))
    gene_CDS_middle <- sort(gene_CDS)[floor(0.05 * length(gene_CDS)):floor(0.95 * length(gene_CDS))]
    if (pos %in% gene_CDS_middle || old_end %in% gene_CDS_middle || new_end %in% gene_CDS_middle) {
      real_lof <- c(real_lof, candidate)
    }
  }
  if (length(real_lof) == 0) real_lof <- NA
  paste(real_lof, collapse = ",")
})
high_impact_info[, 17] <- genes_model_5_95
colnames(high_impact_info)[17] <- "gene_model_5_95"

# Create LoF matrix (structural variants - deletion) ---------------------------

LoF_matrix_SV <- lapply(unique(CDS$gene), function(g) {
  cat(paste0(g, " "))
  
  lof_var_index <- which(grepl(g, variants$gene))
  
  if (length(lof_var_index) > 0) {
    info_subset <- variants[lof_var_index, ]
    
    # Resolve multi-gene annotations
    multiple_genes <- which(grepl(",", info_subset$gene))
    if (length(multiple_genes) > 0) {
      for (i in multiple_genes) {
        genes <- strsplit(info_subset$gene[i], ",")[[1]]
        affected <- strsplit(info_subset$affected_model_percentage[i], ",")[[1]]
        n <- which(genes == g)
        
        info_subset$gene[i] <- genes[n]
        info_subset$affected_model_percentage[i] <- affected[n]
      }
    }
    
    # Filter variants not fully affecting gene models
    delete <- which(as.numeric(info_subset$affected_model_percentage) != 1)
    lof_var_index_new <- if (length(delete) > 0) lof_var_index[-delete] else lof_var_index
    
    # Assign genotype call
    if (length(lof_var_index_new) > 1) {
      lof_gt <- gt_processed[lof_var_index_new, ]
      lof_call <- unlist(apply(lof_gt, 2, function(x) max(as.numeric(x))))
    } else if (length(lof_var_index_new) == 1) {
      lof_call <- as.numeric(gt_processed[lof_var_index_new, ])
    } else {
      lof_call <- rep(0, 1301)
    }
  } else {
    lof_call <- rep(0, 1301)
  }
  
  lof_call <- data.table(lof_call)
  colnames(lof_call) <- g
  return(lof_call)
})

# Combine SV matrix
LoF_matrix_SV <- do.call(cbind, LoF_matrix_SV)
LoF_matrix_SV <- as.matrix(LoF_matrix_SV)
rownames(LoF_matrix_SV) <- colnames(gt_processed)

# Create LoF matrix (the 1001 Genome Project high impact variants) -------------

# Create a binary LoF matrix (genes × accessions)
LoF_matrix_1001 <- lapply(unique(CDS$gene), function(g) {
  cat(paste0(g, " "))
  
  # Get variant rows associated with the gene
  lof_var_index <- which(grepl(g, high_impact_info$gene))
  
  if (length(lof_var_index) > 0) {
    info_subset <- high_impact_info[lof_var_index, ]
    
    # If variants affect multiple genes, extract info specific to this gene
    multiple_genes <- which(grepl(",", info_subset$gene))
    if (length(multiple_genes) > 0) {
      for (i in multiple_genes) {
        genes <- strsplit(info_subset$gene[i], ",")[[1]]
        affected <- strsplit(info_subset$affected_percentage[i], ",")[[1]]
        type <- strsplit(info_subset$LoF_type[i], ",")[[1]]
        n <- which(genes == g)
        
        info_subset$gene[i] <- genes[n]
        info_subset$affected_percentage[i] <- affected[n]
        info_subset$LoF_type[i] <- type[n]  
      }
    }
    
    # Filter out variants that do not meet LoF criteria
    delete1 <- which(as.numeric(info_subset$affected_percentage) != 1)
    
    delete2 <- vector()
    for (o in 1:nrow(info_subset)) {
      lof_type <- paste(unique(strsplit(info_subset$LoF_type[o], ",")[[1]]), collapse = ",")
      if ((lof_type %in% c("frameshift_variant", "stop_gained")) && info_subset$gene_model_5_95[o] == "NA") {
        delete2 <- c(delete2, o)
      }
    }
    
    delete3 <- which(info_subset$LoF_type == "NA")
    
    delete <- union(union(delete1, delete2), delete3)
    lof_var_index_new <- if (length(delete) > 0) lof_var_index[-delete] else lof_var_index
    
    # Assign genotype call
    if (length(lof_var_index_new) > 1) {
      lof_gt <- high_impact_gt_processed[lof_var_index_new, ]
      lof_call <- unlist(apply(lof_gt, 2, function(x) max(as.numeric(x))))
    } else if (length(lof_var_index_new) == 1) {
      lof_call <- as.numeric(high_impact_gt_processed[lof_var_index_new, ])
    } else {
      lof_call <- rep(0, 1135)
    }
  } else {
    lof_call <- rep(0, 1135)
  }
  
  lof_call <- data.table(lof_call)
  colnames(lof_call) <- g
  return(lof_call)
})

LoF_matrix_1001 <- do.call(cbind, LoF_matrix_1001)
LoF_matrix_1001 <- as.matrix(LoF_matrix_1001)
rownames(LoF_matrix_1001) <- colnames(high_impact_gt_processed)

# Collapse two matrices (1001 high impact + SV) --------------------------------

# Sanity check: compare genes and accessions
Gene_list_1 <- colnames(LoF_matrix_1001)
Gene_list_2 <- colnames(LoF_matrix_SV)
identical(Gene_list_1, Gene_list_2)

Accession_list_1 <- rownames(LoF_matrix_1001)
Accession_list_2 <- rownames(LoF_matrix_SV)

# Intersect of accessions
Accession_list_intersect <- intersect(Accession_list_1, Accession_list_2)

# Subset matrices to intersected accessions
Matrix_1_intersect <- LoF_matrix_1001[Accession_list_intersect, ]
Matrix_2_intersect <- LoF_matrix_SV[Accession_list_intersect, ]

# Confirm same order
identical(colnames(Matrix_1_intersect), colnames(Matrix_2_intersect))
identical(rownames(Matrix_1_intersect), rownames(Matrix_2_intersect))

# Combine matrices (collapsed LoF = 1001 high impact + SV)
LoF_matrix_collapsed <- Matrix_1_intersect + Matrix_2_intersect

# Cap values at 1 (LoF presence/absence)
LoF_matrix_collapsed[LoF_matrix_collapsed > 1] <- 1
write.csv(LoF_matrix_collapsed, file = "/path/to/LoF_matrix_collapsed.csv")


