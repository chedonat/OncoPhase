#' chr15_OP1019 : chromosome 15 of Patient OP1019.
#'
#' A dataset containing Allele counts, haplotype phaisng and copy number information on chromosome 10 
#' of Patient OP1019 of the SOX2 Study.
#'
#' @format Contains the following data :  :
#' \describe{
#'   \item{tumoursamples}{The list of tumours/samples of the study}
#'   \item{SNP_allelecount_df}{Data frame containing the count of allele supporting the variant of each mutations
#'      \describe{
#'        \item{}{Chrom : Chromosomes}
#'         \item{}{Start : Starting position}
#'         \item{}{End : End  position}
#'         \item{}{Vartype : variant Type}
#'         \item{}{IsGermline : is the mutation a Germline SNP or a Somatic mutation}
#'         \item{}{Ref : Reference sequence}
#'         \item{}{All : Variant sequence}
#'         \item{}{One entry per sample/tumour}
#'      }
#'   }
#'   \item{ref_allelecount_df}{Data frame containing the count of allele supporting the reference at each mutation}
#'   \item{phasing_association_df}{A data frame containing for each somatic mutations, a colon separated  list of Germline mutations phased to it. }
#'   \item{major_copynumber_df}{A data frame containing the major copy number of the mutations at each tumours/samples}
#'   \item{minor_copynumber_df}{A data frame containing the minor copy number of the mutations at each tumours/samples}
#'   \item{minor_copynumber_df}{A data frame containing the normal cell contamination rate for each mutations  at each tumours/samples}
#' }
"chr15_OP1019"

#' chr22_11152 : chromosome 22 of Patient 11152.
#'
#' A dataset containing Allele counts, haplotype phaisng and copy number information on chromosome 22 
#' of Patient 11152 of the SOX2 Study.
#'
#' @format Contains the following data :  :
#' \describe{
#'   \item{tumoursamples}{The list of tumours/samples of the study}
#'   \item{SNP_allelecount_df}{Data frame containing the count of allele supporting the variant of each mutations
#'      \describe{
#'        \item{}{Chrom : Chromosomes }
#'         \item{}{Start : Starting position}
#'         \item{}{End : End  position}
#'         \item{}{Vartype : variant Type}
#'         \item{}{IsGermline : is the mutation a Germline SNP or a Somatic mutation}
#'         \item{}{Ref : Reference sequence}
#'         \item{}{All : Variant sequence}
#'         \item{}{One entry per sample/tumour}
#'      }
#'   }
#'   \item{ref_allelecount_df}{Data frame containing the count of allele supporting the reference at each mutation}
#'   \item{phasing_association_df}{A data frame containing for each somatic mutations, a colon separated  list of Germline mutations phased to it. }
#'   \item{major_copynumber_df}{A data frame containing the major copy number of the mutations at each tumours/samples}
#'   \item{minor_copynumber_df}{A data frame containing the minor copy number of the mutations at each tumours/samples}
#'   \item{minor_copynumber_df}{A data frame containing the normal cell contamination rate for each mutations  at each tumours/samples}
#' }
"chr22_11152" 
