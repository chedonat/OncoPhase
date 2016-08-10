
#' chr22_11152 : chromosome 22 of Patient 11152.
#'
#' A dataset containing allele counts, haplotype phasing and copy number information on chromosome 22 
#' of Patient 11152 of the SOX2 Study.
#'
#' Also available on the same patient : chr10, chr15, chr18 and chr22
#'
#'
#' @format Contains the following data :  :
#' \describe{
#'   \item{tumoursamples}{The list of tumor samples of the study}
#'   \item{SNP_allelecount_df}{Data frame containing the count of allele supporting the variant of each mutations
#'      \describe{
#'        \item{}{Chrom : Chromosomes }
#'         \item{}{Start : Starting position}
#'         \item{}{End : End  position}
#'         \item{}{Vartype : variant Type}
#'         \item{}{IsGermline : is the mutation a Germline SNP or a Somatic mutation}
#'         \item{}{Ref : Reference sequence}
#'         \item{}{All : Variant sequence}
#'         \item{}{One entry per tumor samples}
#'      }
#'   }
#'   \item{ref_allelecount_df}{Data frame containing the count of allele supporting the reference at each mutation}
#'   \item{phasing_association_df}{A data frame containing for each somatic mutations, a colon separated  list of Germline mutations phased to it. }
#'   \item{major_copynumber_df}{A data frame containing the major copy number of the mutations at each tumor samples }
#'   \item{minor_copynumber_df}{A data frame containing the minor copy number of the mutations at each tumor samples }
#'   \item{minor_copynumber_df}{A data frame containing the normal cell contamination rate for each mutations  at each tumor samples }
#' }
"chr22_11152" 

