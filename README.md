This repository contains code for analyses of LoF burden - gene expression association testing in Arabidopsis.
https://doi.org/10.1093/plcell/koag087

Public data used for this project:

SNP variants and small indel calls for 1,135 Arabidopsis accessions from the 1001 Genomes Project: https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snpeff_v3.1/.
Structural variants and indel calls for 1,301 accessions were downloaded from European Variation Archive (PRJEB38975).
Expression data were obtained from Kawakatsu et al., Epigenomic Diversity in a Global Collection of Arabidopsis thaliana Accessions, Cell 166 (2), p492-505, 14 July 2016, http://dx.doi.org/10.1016/j.cell.2016.06.044.
Kinship: https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/.
Flowering time data: http://1001genomes.org/tables/1001genomes-FT10-FT16_and_1001genomes-accessions.html.

LoF_calling.R in this repository is specific to the Arabidopsis dataset. A universal LoF caller that can be used for any species can be found at: https://github.com/KehanZhao/LoFMatrixBuilder.
