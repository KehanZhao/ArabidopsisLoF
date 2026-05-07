# ArabidopsisLoF

Loss-of-function burden–gene expression association analyses in *Arabidopsis thaliana*

Paper: Kehan Zhao, Mariele Lensink, J Grey Monroe, Functional insights into dispensable genes using genome-wide loss-of-function burden tests in Arabidopsis, The Plant Cell, Volume 38, Issue 4, April 2026, koag087, https://doi.org/10.1093/plcell/koag087

---

## Overview

This repository contains R code for testing associations between loss-of-function (LoF) variant burden and genome-wide gene expression across a large natural population of *Arabidopsis thaliana* accessions. Analyses integrate multi-omic data including SNPs, structural variants, RNA-seq expression, kinship, and phenotypic data from publicly available resources.

The pipeline proceeds in four stages: (1) LoF variant calling from two independent variant sources and construction of the LoF genotype matrix, (2) expression matrix construction aligned to the LoF matrix, (3) EMMAX-based LoF burden–expression association testing, and (4) network visualization of significant associations.

---

## Scripts

### `LoF_calling.R`

Identifies loss-of-function variants in *A. thaliana* from two sources and merges them into a single binary LoF matrix (accessions × genes).

**Source 1 — Structural variants (deletions) from the European Variation Archive:** For each deletion, overlapping CDS models are identified using the TAIR10 GFF3 annotation. Variants are restricted to the 5–95% middle region of the CDS, and flagged as LoF if they cause a frameshift (deletion length not divisible by 3) or delete ≥10% of the coding sequence.

**Source 2 — High-impact SNPs/indels from the 1001 Genomes Project:** High-impact variants (pre-filtered with SnpEff's `HIGH` impact flag) are parsed for frameshift and stop-gain annotations. The same 5–95% CDS filter is applied. Multi-allelic sites are expanded to one row per alternate allele before genotype calling.

The two resulting matrices are intersected on shared accessions and combined into a single binary LoF matrix (`LoF_matrix_collapsed.csv`).

> **Note:** `LoF_calling.R` is specifically tailored for the Arabidopsis TAIR10 annotation and the datasets described in the paper. For a universal LoF caller, see [LoFMatrixBuilder](https://github.com/KehanZhao/LoFMatrixBuilder).

**Key inputs:**
- `TAIR10_GFF3_genes.gff` — gene/CDS annotation
- `TAIR10_all_gene_models.txt` — all TAIR10 gene model IDs
- `Deletions.vcf.gz` — SV/deletion calls (European Variation Archive, PRJEB38975)
- `high_impact.txt` — SnpEff HIGH-impact variants extracted from the 1001 Genomes VCF
- `1001GenomeVCF_header.txt` — header row from the 1001 Genomes VCF

**Key output:** `LoF_matrix_collapsed.csv`

---

### `Expression_matrix.R`

Constructs a gene expression matrix aligned to the accessions and genes present in the LoF matrix. Expression values are batch-corrected log2 values (`d_log2_batch`) from the Kawakatsu *et al*. 2016 transcriptome dataset.

For each (accession, gene) pair, the expression value is looked up from the TG dataset. If an accession has multiple expression measurements for the same gene, they are averaged. Rows and columns that are entirely `NA` are removed, and remaining `NaN` values, which can arise from log2 transformations of zero in upstream batch correction, are replaced with the column minimum.

**Key inputs:** `TG_data_20180606.Rdata`, `LoF_matrix_collapsed` (for accession/gene ordering)

**Key outputs:** `Expression_matrix`, `Expression_matrix_NAremoved_NaNreplaced`

---

### `LoF_burden_exp_association.R`

Tests for associations between per-gene LoF burden and genome-wide gene expression using the EMMAX mixed model, which accounts for population structure via a kinship matrix. For each expression trait, a single EMMAX model is fit with all LoF genes as markers simultaneously.

Prior to testing, expression genes with zero variance are dropped, LoF genes with fewer than 3 alleles across accessions are excluded, and both matrices are subsetted to their shared accessions. An inverse normal transformation (INT) function is included but not applied by default.

> **Note:** Developed and tested under R 3.6.1 using the `cpgen` package. `cGWAS.emmax` may fail under R ≥ 4.2.

**Key inputs:** `Expression_matrix_NAremoved_NaNreplaced`, `S1_LoF_matrix_collapsed.csv`, `kinship_ibs_mac5.hdf5`

**Key outputs:** `pvalue_matrix_Exp`, `beta_matrix_Exp`, `marker_variance_Exp`, `residual_variance_Exp`

---

### `LoF_exp_gene_network.R`

Builds and plots a directed gene network from significant LoF–expression associations. Edges point from LoF gene (source, circle node) to expression gene (target, square node) and are colored by direction of effect: pink for positive β and blue for negative β.

The network is filtered to LoF genes with ≥10 significant expression associations. Highlighted genes (all expression targets of *FRI* (AT4G00650) plus *FRI* itself) are rendered in orange with labels.

**Key input:** `S5_LoF_burden_expression_significant_associations.csv`

**Key output:** `LoF_exp_network.pdf`

---

## Public Data Sources

| Data | Accessions | Source |
|------|-----------|--------|
| SNP & small indel calls | 1,135 | [1001 Genomes Project v3.1](https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snpeff_v3.1/) |
| Structural variants & indels | 1,301 | [European Variation Archive — PRJEB38975](https://www.ebi.ac.uk/ena/browser/view/PRJEB38975) |
| Gene expression (RNA-seq) | — | [Kawakatsu *et al*. 2016, *Cell* 166(2):492–505](http://dx.doi.org/10.1016/j.cell.2016.06.044) |
| Kinship matrix | — | [1001 Genomes — SNP matrix (imputed HDF5)](https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/) |
| Flowering time (FT10, FT16) | — | [1001 Genomes phenotype tables](http://1001genomes.org/tables/1001genomes-FT10-FT16_and_1001genomes-accessions.html) |

---

## Dependencies

All scripts are written in R. Key packages:

| Package | Used in |
|---------|---------|
| `vcfR` | `LoF_calling.R` |
| `data.table` | `LoF_calling.R`, `LoF_exp_gene_network.R` |
| `stringr` | `LoF_calling.R` |
| `rhdf5` | `LoF_burden_exp_association.R` |
| `cpgen` | `LoF_burden_exp_association.R` |
| `igraph` | `LoF_exp_gene_network.R` |

`cpgen::cGWAS.emmax` requires R 3.6.1 and may not work on R ≥ 4.2.

---

## Citation

If you use this code, please cite:

> Kehan Zhao, Mariele Lensink, J Grey Monroe, Functional insights into dispensable genes using genome-wide loss-of-function burden tests in Arabidopsis, The Plant Cell, Volume 38, Issue 4, April 2026, koag087, https://doi.org/10.1093/plcell/koag087

---

## Related

Universal LoF caller for any species: [KehanZhao/LoFMatrixBuilder](https://github.com/KehanZhao/LoFMatrixBuilder)
