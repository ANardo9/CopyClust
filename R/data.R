#' CopyClust Example TCGA Integrative Cluster Labels
#'
#' Gold standard integrative cluster labels of 1,075 example TCGA samples from 'test_data' determined using combined copy number and gene expression classifier from iC10 package.
#'
#' @format ## `test_labels`
#' A 1075 by 2 matrix
"test_labels"

#' CopyClust Example Data from TCGA Pre-Formatted
#'
#' Example data from 1,075 TCGA samples used to test the package pre-formatted for the CopyClust function.
#'
#' @format ## `test_data`
#' A 1075 by 478 matrix
"test_data"

#' CopyClust Example Data from TCGA Raw
#'
#' Example raw copy number values in DNACopy format from 100 TCGA samples used to test the CC_format function.
#'
#' @format ## `test_data_raw`
#' A 20413 by 6 matrix
"test_data_raw"

#' hg18 Ranges
#'
#' Genomic ranges used as features for the CopyClust function linked to the hg18 reference genome. Utilized for formatting with CC_format.
#' Columns: "chrom", "start", "end", "width", "range".
#'
#' @format ## `hg18_ranges`
#' A 478 by 5 matrix
"hg18_ranges"

#' hg19 Ranges
#'
#' Genomic ranges used as features for the CopyClust function linked to the hg19 reference genome. Utilized for formatting with CC_format.
#' Columns: "chrom", "start", "end", "width", "range".
#'
#' @format ## `hg19_ranges`
#' A 478 by 5 matrix
"hg19_ranges"

#' hg38 Ranges
#'
#' Genomic ranges used as features for the CopyClust function linked to the hg38 reference genome. Utilized for formatting with CC_format.
#' Columns: "chrom", "start", "end", "width", "range".
#'
#' @format ## `hg38_ranges`
#' A 478 by 5 matrix
"hg38_ranges"


