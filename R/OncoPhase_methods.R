


#' chr22_XYZ101 : Patient XYZ101 data from chromosome 22 .
#'
#' A generated dataset containing allele counts, haplotype phasing and copy number information on chromosome 22 
#' for a patient XYZ101 (Data created) .
#'
#' Also available on the same patient : chr10, chr15 and chr18
#'
#'
#' @format Contains the following data :  :
#' \describe{
#'  \item{tumoursamples}{The list of tumor samples of the study}
#'  \item{SNP_allelecount_df}{Data frame containing the count of allele supporting the variant of each mutations
#'     \describe{
#'       \item{}{Chrom : Chromosomes }
#'        \item{}{Start : Starting position}
#'        \item{}{End : Position of the mutation}
#'        \item{}{Vartype : variant Type}
#'        \item{}{IsGermline : is the mutation a Germline SNP or a Somatic mutation}
#'        \item{}{Ref : Reference sequence}
#'        \item{}{All : Variant sequence}
#'        \item{}{One entry per tumor samples}
#'     }
#'  }
#'  \item{ref_allelecount_df}{Data frame containing the count of allele supporting the reference at each mutation}
#'  \item{phasing_association_df}{A data frame containing for each somatic mutations, a colon separated  list of Germline mutations phased to it. }
#'  \item{major_copynumber_df}{A data frame containing the major copy number of the mutations at each tumor samples }
#'  \item{minor_copynumber_df}{A data frame containing the minor copy number of the mutations at each tumor samples }
#'  \item{minor_copynumber_df}{A data frame containing the normal cell contamination rate for each mutations  at each tumor samples }
#' }
"chr22_XYZ101"










#' 
#' 
#' Computes cellular prevalence at a single mutation point
#' 
#' This is a generic function to compute the cellular prevalence of a somatic mutation point using OncoPhase method. The method computes the prevalence of the somatic mutation relatively to phased nearby SNPs whose prevalence are known to be 1.  \code{\link{getPrevalence}} requires the allelic-information of the somatic mutation and the aggregated information of its Phased SNP but the function can also be run in the absence of phasing information (FlankingSNP mode) or nearby SNP (SNVOnly mode).
#' 
#'     The method particularly exhibits an increase in the accuracy when the locus of the SNP is also affected by a somatic copy number alteration (SCNA). The method detect the temporal relationship between the two alterations (C1: SNV occurred after the SCNA; C2: SNV occurred before the SCNA) and computes the detailed prevalence of each of the following group of cells (if detail is set to TRUE) : 
#'       \describe{
#'        \item{Germ}{ Cells having a germline genotype  at the locus of the SNV. That is No SNV, no SCNA}
#'        \item{Alt}{ Cells having one alternative of  the two somatic alteration. That is either the SCNA, either the SNV not both.}
#'        \item{Both}{ Cells having both somatic alterations. That is the SNV and the SCNA}
#'     }
#' 
#' OncoPhase can be run under three modes:
#' 
#'  \describe{
#'        \item{PhasedSNP}{ The prevalence is computed relatively to a Phased SNP}
#'        \item{FlankingSNP}{ In the absence of phasing information, the prevalence is computed relatively to a neighbor SNP located on the same locus with the somatic SNV. NA is returned if the prevalence cannot be resolved without knowing the phasing information between the SNP and the SNV.}
#'        \item{SNVOnly}{The prevalence is computed using only the SNV information without the usage of any nearby SNP}
#'     }
#' 
#' @param varcounts_snv A count (or a vector of counts  if multiple samples ) of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param refcounts_snv A count (or a vector of counts  if multiple samples ) of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn  major copy number (or a vector  if multiple samples ) 
#' at the locus of the mutation
#' @param minor_cn minor copy number (or a vector  if multiple samples)  
#'  at the locus of the mutation 
#' @param varcounts_snp  A count (or a vector of counts  if multiple samples) of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param refcounts_snp  A count (or a vector of counts  if multiple samples) of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param mode The mode under which the prevalence is computed  (default : PhasedSNP , alternatives methods  are FlankingSNP and SNVOnly). Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP
#' @param detail when set to 0, the function simply output the cellular prevalence of the somatic mutation. if set to 1,  a detailed output is generated containing:
#'  \describe{
#'        \item{Context}{ The associated context (C=1 or C=2) }
#'        \item{Prevalence}{The computed cellular prevalence}
#'        \item{DetailedPrevvalence}{the detailed prevalence for each group of cells (germline cells (Germ), cells affected by one of the two genomic alterations (Alt), cells affected  by both genomic alterations (Both)}
#'        \item{ResidualNorm}{Norm of the residuals when an approximation is performed}
#'        \item{CondensedPrevalence}{A colon separated list of the above fields (Context, Prevalence, Detailedprevalence and ResidualNorm). The detailed prevalence are separated by "|"}
#' }
#'   if detail is set to 2, the function outputs the CondensedPrevalence field above.     
#' @param Trace if set to TRUE, print the trace of the computation.    
#'  
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage. 
#' @param SomaticCountAdjust when set to 1, varcounts_snv and refcounts_snv might be adjusted if necessary so that they meet the rules varcounts_snv <= varcounts_snp, refcounts_snv >= refcounts_snp and varcounts_snv + refcounts_snv ~ Poiss(varcounts_snp + refcounts_snp). Not used if mode=SNVOnly,  
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment.  
#' 
#'
#' @return   The cellular prevalence if detail =0, a detailed output if detail = 1, and a condensed output if detail =2. See the usage of the parameter detail above.
#' 
#'     
#'      
#' @examples
#' 
#' #Example 1
#' 
#' prevalence=getPrevalence(varcounts_snv=14,refcounts_snv=10,major_cn=3,minor_cn=1,
#' varcounts_snp=16, refcounts_snp=8  )
#' 
#' print(prevalence)
#' 
#' #Sample_1 
#' #0.75   
#' #The above example  gives the same  prevalence with mode FlankingSNP but not with mode SNVOnly
#' prevalence=getPrevalence(varcounts_snv=14,refcounts_snv=10,major_cn=3,minor_cn=1,
#' varcounts_snp=16, refcounts_snp=8 ,mode="FlankingSNP")
#' print(prevalence)
#' #Sample_1 
#' #0.75 
#' prevalence=getPrevalence(varcounts_snv=14,refcounts_snv=10,major_cn=3,minor_cn=1,
#' varcounts_snp=16, refcounts_snp=8 ,mode="SNVOnly")
#' print(prevalence)
#' #Sample_1 
#' #0.79
#'  
#' #Example 2 Case Study A (see paper)
#' prevalence = getPrevalence(varcounts_snv=6,refcounts_snv=8,major_cn=2,minor_cn=1,
#' varcounts_snp=8, refcounts_snp=6, detail=1)
#' print(prevalence)
#' #$Sample_1
#' #$Sample_1$Context
#' #[1] "C2"
#' #
#' #$Sample_1$Prevalence
#' #[1] 0.66
#' #
#' #$Sample_1$DetailedPrevalence
#' #Germ  Alt Both 
#' #0.33 0.33 0.33 
#' #
#' #$Sample_1$ResidualNorm
#' #residual 
#' #2e-31 
#' #
#' #$Sample_1$CondensedPrevalence
#' #[1] "C2:0.66:0.33|0.33|0.33:2e-31"
#'   
#' ## The above example is not resolvable without phasing information
#'  prevalence =  getPrevalence(varcounts_snv=6,refcounts_snv=8,major_cn=2,minor_cn=1,
#'  varcounts_snp=8, refcounts_snp=6, mode="FlankingSNP",detail=TRUE)
#'  print(prevalence)
#' #Warning message:
#' #  In getFlankingSNPPrevalence(varcounts_snv, refcounts_snv, major_cn, minor_cn,  :
#' #                                The prevalence is not resolved without the knowledge of the Phased Germline
#' #NA
#' #'## it gives an inaccurate  prevalence under mode "SNVOnly"    
#' prevalence= getPrevalence(varcounts_snv=6,refcounts_snv=8,major_cn=2,minor_cn=1,
#' varcounts_snp=8, refcounts_snp=6, mode="SNVOnly",detail=0)
#' print(prevalence)
#' #Sample_1 
#' #1 
#' #Example 3 Case Study B (see paper) Not resolvable without phasing information
#' prevalence=getPrevalence(varcounts_snv=4,refcounts_snv=8,major_cn=2,minor_cn=0,
#' varcounts_snp=8, refcounts_snp=4, detail=2)
#' print(prevalence)
#' #Sample_1 
#' # "C2:0.33:0.67|0|0.Prevalence(varcounts_snv=6,refcounts_snv=8,major_cn=2,minor_cn=1,
#' varcounts_snp=8, refcounts_snp=6, detail=TRUE,Trace=TRUE )33:5.2e-32" 
#' 
#' #Example 4 Case Study A (see paper) Not resolvable without phasing information
#' prevalence = getPrevalence(varcounts_snv=6,refcounts_snv=8,major_cn=2,minor_cn=1,
#' varcounts_snp=8, refcounts_snp=6 )
#' print(prevalence)      
#' #Sample_1 
#' #0.66  
#' 
#' #Example 5 We group case study A, B and C above to form a multi-sample  case
#' prevalence=getPrevalence(varcounts_snv=c(6,4,6),refcounts_snv=c(8,8,14),major_cn=c(2,2,2),
#' minor_cn=c(1,0,1),varcounts_snp=c(8,8,8), refcounts_snp=c(6,4,12) )                   
#' print(prevalence)
#' #Sample_1 Sample_2 Sample_3 
#' #0.66     0.33     0.75   
#' #'
#' 
#' @seealso \code{\link{getPrevalence}},  \code{\link{getSamplePrevalence}},   \code{\link{getSinglePrevalence}}, \code{\link{getPrevalenceSNVOnly}}                                                 
#' @export
getPrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, varcounts_snp=NULL, refcounts_snp=NULL,  detail=0, mode="PhasedSNP",Trace=FALSE,LocusCoverage=TRUE,SomaticCountAdjust=TRUE,NormalCellContamination=NULL,Optimal=TRUE)
{
  
  
  
  
  N=length(varcounts_snv) # Number of samples
  
  if((length(refcounts_snv)!=N) || 
     (!is.null(varcounts_snp) && (mode !="SNVOnly") &&((length(varcounts_snp) !=N) || (length(refcounts_snp) !=N))) ||
     (length(major_cn) !=N) || (length(minor_cn)!=N) )
    stop("\n\nThe vectors passed as input should have the same size\n\n")
  
  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=FlankingSNP
  numeric_mode=c("SNVOnly", "PhasedSNP","FlankingSNP")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1 or 2")
    }
  }
  
  
  
  
  
  if(mode=="PhasedSNP"){
    prev_somatic=getPhasedSNPPrevalence( varcounts_snv,refcounts_snv , major_cn,minor_cn, varcounts_snp , refcounts_snp,detail,Trace=Trace,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination= getSinglePrevalence,Optimal=Optimal)
  }else if(mode=="FlankingSNP"){
    prev_somatic=getFlankingSNPPrevalence( varcounts_snv,refcounts_snv ,  major_cn,minor_cn,varcounts_snp , refcounts_snp, detail,Trace=Trace,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination=NormalCellContamination,Optimal=Optimal)
  }else if(mode=="SNVOnly"){
    prev_somatic=getSNVOnlyPrevalence(varcounts_snv,refcounts_snv ,major_cn,minor_cn, detail, Trace=Trace,NormalCellContamination=NormalCellContamination)
  }else {
    stop("parameter mode should be either FlankingSNP, PhasedSNP or SNVOnly")
  }
  
  
  
  
  
  prev_somatic
}









#' Somatic mutations cellular prevalence on a sample.
#' 
#' This function computes the cellular prevalence of a list of somatic mutations in
#'  cancer.  The function applies the model to a range of mutations located at a given genomic region or at the whole genome scale. 
#' The function invokes \code{\link{getPrevalence}}  to compute the cellular prevalence for each mutation of the list.
#' The model computes the prevalence of a somatic
#'   mutation relatively to close and eventually phased germline SNP as specified in \code{\link{getPrevalence}}. 
#' 
#' @param input_df A data frame containing for each mutations :
#' \describe{
#'        \item{varcounts_snv}{Alelle counts supporting the SNV}
#'        \item{refcounts_snv}{Alelle counts supporting the reference at the SNV locus}
#'        \item{major_cn}{Major copy number at the SNV locus}
#'        \item{minor_cn}{Minor copy number at the SNV locus}
#'        \item{varcounts_snp}{Alelle counts supporting the SNP}
#'        \item{refcounts_snp}{Alelle counts supporting the reference at the SNP}
#'      }
#' 
#' @param nbFirstColumns Number of first columns in input_df to reproduce in the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation.
#' 
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end   e.g "chr22:179800-98767. 
#' @param mode The mode under which the prevalence is computed  (Default : PhasedSNP , alternatives methods  are FlankingSNP, OptimalSNP,and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP and 3 = OptimalSNP
#' @param detail when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations (SNV or CNV) but not both, cells affected by  both copy number alteration and SNV ). 
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage.
#' @param SomaticCountAdjust when set to TRUE, varcounts_snv and refcounts_snv might be adjusted if necessary so that they meet the rules varcounts_snv <= varcounts_snp, refcounts_snv >= refcounts_snp and varcounts_snv + refcounts_snv ~ Poiss(varcounts_snp + refcounts_snp). Not used if mode=SNVOnly,  
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment.  
#' 
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame input_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence dataframe (e.g REF, ALL, ...) }
#'         \item{}{and the following information}
#'   \describe{
#'        \item{Prev}{The Cellular Prevalence of the mutation}
#'        \item{Germ}{The proportion of cells with a normal genotype}
#'        \item{Alt}{The proportion of cells with only the CNA if the context C=C1 or with only the SNV if the context C=C2}
#'        \item{Both}{The proportion of cells with both the SNV and the SCNA}
#'        \item{Context}{Context at the mutation. If C1 then the SNV occurred after the SCNA, if C=c2 then the SNV occurred before the SCNA}
#'        \item{residual}{Residual after limSolve approximation.}
#'      }
#'      }
#'      
#' @examples
#' 
#' #Example 1: 
#' 
#' input_file=system.file("extdata","phylogeny1_d300_n80.tsv", package = "OncoPhase")
#' input_df<-read.table(input_file,header=TRUE)
#' rownames(input_df) = input_df$mutation_id
#' print(input_df)
#' #  mut_id varcounts_snv refcounts_snv major_cn minor_cn varcounts_snp refcounts_snp
#' #a      a      151  152        1        1      151  135
#' #b      b      123  176        1        1      161  150
#' #c      c       94  209        2        1      176  134
#' #d      d       23  283        1        1      155  144
#' #e      e       60  228        2        0      174  125
#' 
#' prevalence_df=getSamplePrevalence(input_df,nbFirstColumns = 1)
#' 
#' print(prevalence_df)
#' #  mut_id   Prev   Germ    Alt   Both Residual Context
#' #  a      a 0.9967 0.0017 0.0017 0.9967  3.1e-03      C1
#' #  b      b 0.8230 0.0890 0.0890 0.8230  1.3e-03      C1
#' #  c      c 0.4010 0.6000 0.0910 0.3100  3.9e-33      C2
#' #  d      d 0.1500 0.4200 0.4200 0.1500  1.4e-03      C1
#' #  e      e 0.2490 0.7500 0.0890 0.1600  5.1e-31      C2
#' 
#' 
#' #'@seealso \code{\link{getPrevalence}}
#' @export
#' 
getSamplePrevalence<-function(input_df,mode="PhasedSNP",  nbFirstColumns=0, region=NULL,detail=TRUE,  LocusCoverage=TRUE,SomaticCountAdjust=TRUE,NormalCellContamination=NULL,Optimal=TRUE)
{
  
  #Check the compulsory columns
  compulsory_columns=c("varcounts_snv","refcounts_snv","major_cn","minor_cn")
  
  if (length(setdiff(compulsory_columns,colnames(input_df)))>0){
    stop(" The allele count master matrices should have at least the following headers
         columns : ", compulsory_columns)
    # print(compulsory_columns)
  }
  
  
  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=FlankingSNP, 3 = OptimalSNP
  numeric_mode=c("SNVOnly", "PhasedSNP","FlankingSNP","OptimalSNP")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2,3))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1,  2 or 3")
    }
  }
  
  
  #Prepare the master matrices for the prevalence.
  masterprevalence=as.data.frame(matrix(nrow=nrow(input_df), ncol=nbFirstColumns+6))
  rownames(masterprevalence) = rownames(input_df)
  if(nbFirstColumns>0){
    masterprevalence[1:nbFirstColumns] = input_df[1:nbFirstColumns]
    colnames(masterprevalence) = c(colnames(input_df[1:nbFirstColumns]),"Prevalence", "Germ","Alt","Both","Residual","Context")
  }else{
    colnames(masterprevalence) = c("Prevalence", "Germ","Alt","Both","Residual","Context")
  }
  
  
  for(imut in 1: nrow(input_df))
  {
    varcounts_snv = input_df[imut,"varcounts_snv"]
    refcounts_snv = input_df[imut,"refcounts_snv"]
    minor_cn=input_df[imut,"minor_cn"]
    major_cn=input_df[imut,"major_cn"]
    varcounts_snp=input_df[imut,"varcounts_snp"]
    refcounts_snp=input_df[imut,"refcounts_snp"]
    
    InputValues=paste(varcounts_snv,refcounts_snv,minor_cn,major_cn,varcounts_snp,refcounts_snp,sep=":")
    prevalence= getPrevalence(varcounts_snv, refcounts_snv, major_cn, minor_cn, varcounts_snp, refcounts_snp, detail=T, mode=mode ,LocusCoverage=LocusCoverage, SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination=NormalCellContamination,Optimal=Optimal)
    
    if(is.null(prevalence) ||  length(prevalence)==0 || is.na(prevalence[[1]]) )
      next
    
    masterprevalence[imut,"Prevalence"] = prevalence[[1]]$Prevalence
    masterprevalence[imut,"Germ"] = prevalence[[1]]$DetailedPrevalence["Germ"]
    masterprevalence[imut,"Alt"] = prevalence[[1]]$DetailedPrevalence["Alt"]
    masterprevalence[imut,"Both"] = prevalence[[1]]$DetailedPrevalence["Both"]
    masterprevalence[imut,"Context"] = prevalence[[1]]$Context
    masterprevalence[imut,"Residual"] = prevalence[[1]]$ResidualNorm
    masterprevalence[imut,"InputValues"] = InputValues
  }
  
  
  
  masterprevalence
  
}







#' Somatic mutations cellular prevalence computation using haplotype phasing on a multiple sample study.
#' 
#' This is a generic function to compute the cellular prevalence of somatic mutations in
#'  cancer using haplotype phasing.  The function applies the model to a range of mutations located at a given genomic region or at the whole genome scale. The model computes the prevalence of a somatic
#'   mutation relatively to close and eventually phased germline mutations. It uses three sources
#'    of information as input : The allelic counts, the phasing information and the 
#'    copy number alteration.  Multiple tumor samples can be provided for the prevalence computation.
#' 
#' @param snp_allelecount_df A data frame containing for each mutation the  allelic 
#' Counts of the variant at each tumor samples. The data frame should contain at least the following three columns among its firsts columns: Chrom (The mutation
#'  chromosome) , Pos or End (The mutation position) and IsGermline (is the mutation a germline
#'   or somatic mutation).
#' @param ref_allelecount_df A data frame containing for each mutation the allelic count
#'  of the reference at each tumor sample. The data frame should contain at least the following three columns among its firsts columns:  Chrom (The mutation
#' chromosome) , Pos or End  and IsGermline (is the mutation a Germline
#'    or Somatic mutation)  
#' @param major_copynumber_df A data frame containing for each mutation, its  major
#' chromosomal copy number at each tumor samples. Should contain at least the following three columns among its firsts columns:  Chrom (The mutation
#' chromosome) , Pos or End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic mutation)  
#' 
#' @param minor_copynumber_df A data frame containing for each mutation the minor
#'  chromosomal copy number at each tumor samples. Should contain at least the following three columns among its firsts columns:  Chrom (The mutation
#' chromosome) , Pos or End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic mutation)  
#' 
#' @param mode The mode under which the prevalence is computed  (default : PhasedSNP , alternatives modes  are FlankingSNP, OptimalSNP and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP and 3 = OptimalSNP
#' @param phasing_association_df A data frame containing for each somatic mutation, 
#' a colon separated list of germline SNP phased to it.
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment.  
#' @param nbFirstColumns Number of first columns in snp_allelecount_df to reproduce in
#'  the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation at each tumor sample
#' @param tumoursamples The list of tumor samples to consider for the prevalence
#' computation.  These samples should be present as column headers in the data frame
#'   snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df. If not provided, the headers from nbFirstColumns + 1 to 
#'   the last column of snp_allelecount_df are retrieved and their  intersection with the
#' other inputted data frames headers is considered
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end 
#' e.g "chr22:179800-98767 
#' @param detail when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations (SNV or CNV) but not both, cells affected by  both copy number alteration and SNV ). 
#' @param LocusRadius Only phased SNPs located within LocusRadius bp from the somatic mutation will be considered.
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage. 
#' @param SomaticCountAdjust when set to TRUE, varcounts_snv and refcounts_snv might be adjusted if necessary so that they meet the rules varcounts_snv <= varcounts_snp, refcounts_snv >= refcounts_snp and varcounts_snv + refcounts_snv ~ Poiss(varcounts_snp + refcounts_snp). Not used if mode=SNVOnly.
#' 
#' 
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence data frame (e.g REF, ALL, ...) }
#'         \item{}{One column per tumor sample reporting the prevalence of the mutation 
#'         at each samples}
#'      }
#'      
#' @examples
#' 
#' #Example 1: Loading a simple example data set with two somatic mutations, 5 germlines SNP
#' # and 3 tumor samples
#' data(simpleExample2)
#' se=simpleExample2
#' prevalence_df=getMultiSamplesPrevalence(se$snp_allelecount_df, se$ref_allelecount_df,
#'  se$major_copynumber_df,se$minor_copynumber_df,phasing_association_df=se$phasing_association_df, )
#' print(prevalence_df)
#' 
#' #Chrom     End IsGermline  Tumour1        Tumour2        Tumour3
#' #mutation2  chr2 3003000          0 C2:0|0|1 C2:0.15|0|0.85 C2:0.12|0|0.88
#' #mutation6  chr2 4008000          0 C1:1|0|0       C1:1|0|0 C2:0|0.24|0.76
#' 
#' #Example 2 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient XYZ101 
#' # (data created)
#' 
#' data("chr22_XYZ101")
#' ds=chr22_XYZ101
#' masterprevalence_df=getMultiSamplesPrevalence(ds$snp_allelecount_df, ds$ref_allelecount_df,
#'   ds$major_copynumber_df,ds$minor_copynumber_df,phasing_association_df = ds$phasing_association_df,
#'    nbFirstColumns=6,detail=FALSE)
#' print(head(masterprevalence_df))
#' 
#' data("chr18_XYZ101")
#' df=chr18_XYZ101
#' masterprevalence_df=getMultiSamplesPrevalence(df$snp_allelecount_df, df$ref_allelecount_df,
#'  df$major_copynumber_df,df$minor_copynumber_df,phasing_association_df=df$phasing_association_df, 
#'  nbFirstColumns=6, region="chr18:10000000-80000000")
#' print(head(masterprevalence_df))
#' 
#' 
#' 
#' #'@seealso \code{\link{getPrevalence}}
#' @export
getMultiSamplesPrevalence<-function(snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode="PhasedSNP", phasing_association_df=NULL,NormalCellContamination=NULL,nbFirstColumns=3, tumoursamples=NULL,   region=NULL,detail=TRUE,  LocusRadius = 10000,LocusCoverage=TRUE,SomaticCountAdjust=TRUE,Optimal=TRUE)
{
  
  
  if("End" %in% colnames(snp_allelecount_df)) 
    colnames(snp_allelecount_df)[which(colnames(snp_allelecount_df)=="End")] = "Pos"
  if("End" %in% colnames(ref_allelecount_df)) 
    colnames(ref_allelecount_df)[which(colnames(ref_allelecount_df)=="End")] = "Pos"
  if("End" %in% colnames(minor_copynumber_df)) 
    colnames(minor_copynumber_df)[which(colnames(minor_copynumber_df)=="End")] = "Pos"
  if("End" %in% colnames(major_copynumber_df)) 
    colnames(major_copynumber_df)[which(colnames(major_copynumber_df)=="End")] = "Pos"
  

  # Extract the somatic mutations 
  
  cat("\n\n Cellular prevalence computation using OncoPhase")
  
  
  compulsory_columns=c("Chrom","Pos","IsGermline")
  
  print(head(snp_allelecount_df))
  if (length(setdiff(compulsory_columns,colnames(snp_allelecount_df)))>0){
    stop(paste(" The allele count master matrices should have at least the following headers columns : ",paste(compulsory_columns,collapse=" ")))
  }
  
  
  if (length(setdiff(compulsory_columns,colnames(ref_allelecount_df)))>0){
    stop(paste(" The allele count master matrices should have at least the following headers columns : ",paste(compulsory_columns,collapse=" ")))
  }
  
  #Restriction to the region
  
  
  # If a region is provided, a restriction is performed on the given region
  if(!is.null(region)){
    
    cat("\n\n Restriction of the matrices within the region ", region,"...\n\n")
    region_parts= unlist(strsplit(region,":"))
    
    chrom = region_parts[1]
    startPosition = 1
    endPosition = hg19_dfsize[chrom]

    snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$Chrom == chrom , ]
    ref_allelecount_df = ref_allelecount_df[ref_allelecount_df$Chrom == chrom , ]
    major_copynumber_df = major_copynumber_df[major_copynumber_df$Chrom == chrom , ]
    minor_copynumber_df = minor_copynumber_df[minor_copynumber_df$Chrom == chrom , ]
    if(!is.null(phasing_association_df))
      phasing_association_df = phasing_association_df[phasing_association_df$Chrom == chrom , ]
   # if(!is.null(NormalCellContamination_df))
   #   NormalCellContamination_df = NormalCellContamination_df[NormalCellContamination_df$Chrom == chrom , ]
    
    if(length(region_parts)>1){
      coordinates = unlist(strsplit(region_parts[2],"-"))
      startPosition = as.numeric(coordinates[1])
      endPosition = as.numeric(coordinates[2])
      
      snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$Pos >= startPosition & snp_allelecount_df$Pos <= endPosition,]
      ref_allelecount_df = ref_allelecount_df[ref_allelecount_df$Pos >= startPosition & ref_allelecount_df$Pos <= endPosition,]
      major_copynumber_df = major_copynumber_df[major_copynumber_df$Pos >= startPosition & major_copynumber_df$Pos <= endPosition,]
      minor_copynumber_df = minor_copynumber_df[minor_copynumber_df$Pos >= startPosition & minor_copynumber_df$Pos <= endPosition,]

      if(!is.null(phasing_association_df))
        phasing_association_df= phasing_association_df[phasing_association_df$Pos >= startPosition & phasing_association_df$Pos <= endPosition,]
     # if(!is.null(NormalCellContamination_df))
     #   NormalCellContamination_df = NormalCellContamination_df[NormalCellContamination_df$Pos >= startPosition & NormalCellContamination_df$Pos <= endPosition,]
      
    }
  }
  
  
  
  
  
  
  
  if (is.null(tumoursamples)){
    tumoursamples = colnames(snp_allelecount_df[(nbFirstColumns+1):ncol(snp_allelecount_df)])
  }
  
  #print(colnames(snp_allelecount_df))
  
  
  tumoursamples = Reduce(intersect,list(tumoursamples,colnames(snp_allelecount_df),
                                        colnames(ref_allelecount_df),
                                        colnames(major_copynumber_df), 
                                        colnames(minor_copynumber_df)
  ))
  

  if(length(tumoursamples) ==0)
  {
    stop(" None of the tumour samples provided is present in the five  matrices :
         snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df")
  }
  
  snp_allelecount_df=numeric_column(snp_allelecount_df,tumoursamples)
  ref_allelecount_df=numeric_column(ref_allelecount_df,tumoursamples)  
  major_copynumber_df=numeric_column(major_copynumber_df,tumoursamples)
  minor_copynumber_df=numeric_column(minor_copynumber_df,tumoursamples)

  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=FlankingSNP, 3 = OptimalSNP
  numeric_mode=c("SNVOnly", "PhasedSNP","FlankingSNP","OptimalSNP")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2,3))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1,  2 or 3")
    }
  }
  
 
  
    
    #Extraction of the list of somatic mutations the cellular prevalence will be computed.
    somatic_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==0, ]

    
    
    #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=FlankingSNP, 3 = OptimalSNP
    numeric_mode=c("SNVOnly", "PhasedSNP","FlankingSNP","OptimalSNP")
    if(is.numeric(mode))
    {
      if(mode %in% c(0,1,2,3))
      {
        mode = numeric_mode[mode +1 ]
      }else{
        stop("\n\n Mode parameter, if numeric,  should be either 0, 1,  2 or 3")
      }
    }
    
    
    
    #We then  retrieve for each somatic mutation the list of germline mutations to consider for the prevalence computation
    # a) if PhasedSNP mode, then the considered germline are the germline mutations phased to the somatic mutation and located within the same locus
    # b) if FlankingSNP mode then the considered germline are the close germlines located within LocusRadius distance from the somatic mutation.
    #c) if OptimalSNP mode tho columns are provided, the first for the phased germline, and if only there is not phasing information for this mutation, then the second column contains the close germlines located within LocusRadius
    
    
    
    
    #Initialising the masterprevalence matrice, setting the headers
    prevalence_columns=tumoursamples
    if(detail ==1)
    {
      prevalence_columns=c()
      for(sample in tumoursamples){
        prevalence_columns=c(prevalence_columns,paste(sample,c("Prevalence","Context","Germ","Alt","Both","ResidualNorm","InputValues"),sep="_"))
      }
    }
    masterprevalence<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=nbFirstColumns + length(prevalence_columns))
    masterprevalence<-as.data.frame(masterprevalence)
    colnames(masterprevalence) <- c(colnames(somatic_snp_allelecount_df[1:nbFirstColumns]),prevalence_columns)
    rownames(masterprevalence) <- rownames(somatic_snp_allelecount_df)
    masterprevalence[1:nbFirstColumns] = somatic_snp_allelecount_df[1:nbFirstColumns]
    
    
    
    
    # masterprevalence=masterprevalence[listover_estimated,]
    #  mut=""
    Nbmutations = nrow(masterprevalence)
    
    
    cat("\n\n Number of mutations : ", Nbmutations)
    
      TraceProgress=T
      cat(" \n A progression message giving information about the mutation under processing will be displayed each  (N/100)th mutation. Set ProgressOutputs to FALSE to turn off this. ")
      trace_step= ceiling(Nbmutations/100)

    
    for (imut in 1:nrow(masterprevalence))
      #for (imut in 1:5)
    {
      
      
      #Mutation name and mutation position
      mut <- rownames(masterprevalence[imut,]); 
      mut_pos=as.numeric(masterprevalence[imut,"Pos"])
      
      if(imut%%trace_step==0){
        cat("\n Mutation ", mut, " ", imut, "/",Nbmutations,sep="" )
      }
      
      
      #For each mutation, we need to extract one value or one vector (if multiple samples)  of :
      # - varcounts_snp and refcounts_snp : Respectively Variant and reference coverage/count of the phased/nearby Germline Mutations
      # - varcounts_snv and refcounts_snv : Respectively Variant and reference coverage/count of the somatic mutations
      # - major_cn and minor_cn : major and minor copy number at the locus of the somatic mutations 
      # - CNV_fraction : If provided, fraction of cells affected by the CNV
      # - NormalCell_fraction : If provided, fraction of normal cell contamination.
      
      
      
      mode_locus=mode  # Mode to consider while computing the mutation prevalence and retrievning the germline on the same locus with the mutation
      if(mode=="OptimalSNP") # If Mode = Optimal, then we  willtry first PhasedSNP, latter if no germline found, we  will try FlankingSNP
        mode_locus="PhasedSNP"
      
      
      ############################
      #####Source of information 1: The list of linked germline mutations 
      ##############
      
      if(mode!="SNVOnly")
      {
        
        
        
        linked_germline_df=getLocusGermlineMutations(somatic_snp_allelecount_df[mut,], snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,phasing_association_df,  tumoursamples,mode_locus,  LocusRadius)
        linked_germline=as.character(linked_germline_df[mut,"LinkedGermlines"])
        
        #Id mode=Optimal then in case no germline is found with PhasedSNP the search is relaunched with FlankingSNP
        
        if(is.null(linked_germline)|| is.na(linked_germline))
          if((mode_locus=="PhasedSNP")&&(mode=="OptimalSNP")){
            mode_locus="FlankingSNP"
            
            linked_germline_df=getLocusGermlineMutations(somatic_snp_allelecount_df[mut,], snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,phasing_association_df,  tumoursamples,mode_locus,  LocusRadius)
            linked_germline=as.character(linked_germline_df[mut,"LinkedGermlines"])
            
            
          }
        
        
        #       cat("\n\n Linked for mutation ", mut,"\n")
        #       print(linked_germline)
        #       cat("\n\n mode locus ", mut,"\n")
        #       print(mode_locus)
        
        if(is.null(linked_germline)|| is.na(linked_germline))
          next
        
        
        
        
        #List of germline mutations linked to the considered somatic mutation
        linkedGermlines_list<-unlist(strsplit(linked_germline,":"))
        if (length(unlist(linkedGermlines_list))==0)
          next
        
        
      }
      
      
      ############################
      #####Source of information 2: The Copy Number 
      ##############
      
      #For the somatic mutation
      major_cn_sample_somatic = major_copynumber_df[mut,tumoursamples, drop=F]
      minor_cn_sample_somatic = minor_copynumber_df[mut,tumoursamples,drop=F]
   #   NormalContamination_sample_somatic = minor_copynumber_df[mut,tumoursamples,drop=F]
      
      ############################
      #####Source of information 3: The Allele Count 
      
      #Somatic
      #wellfraction_somatic =snp_allelecount_df[mut,cifs:nbcolumns_wellfraction]
      snpwellcount_somatic =  data.matrix(snp_allelecount_df[mut,tumoursamples,drop=F])
      refwellcount_somatic = data.matrix(ref_allelecount_df[mut,tumoursamples,drop=F])
      
      #Germline
      # a) if mode_locus=PhasedSNP, according to the MLE (EM) estimation, at each sample, we consider the average counts of the linked germline.
      # b) if mode_locus=FlankingSNP, at each sample, we consider the closest germline mutation having all the required information (Alleles Count and Copy number Information)
      

      if(mode_locus=="PhasedSNP")
      {
        snpwellcount_germlines=  colMeans(data.matrix(snp_allelecount_df[linkedGermlines_list,tumoursamples, drop=F]), 
                                          na.rm=T)
        refwellcount_germlines =colMeans( data.matrix(ref_allelecount_df[linkedGermlines_list,tumoursamples, drop=F]), 
                                          na.rm=T)
        
      }else   if(mode_locus=="FlankingSNP")
      {
        snpwellcount_germlines=  vector("numeric",length=length(tumoursamples))
        refwellcount_germlines = vector("numeric",length=length(tumoursamples))
        
        names(snpwellcount_germlines) = tumoursamples
        names(refwellcount_germlines) = tumoursamples
        
        
        #For the Linked germline mutations
        major_cn_sample_LinkedGermline_df = major_copynumber_df[linkedGermlines_list,tumoursamples,drop=F ]
        snp_allelecount_LinkedGermline_df= data.matrix(snp_allelecount_df[linkedGermlines_list,tumoursamples, drop=F])
        ref_allelecount_LinkedGermline_df= data.matrix(ref_allelecount_df[linkedGermlines_list,tumoursamples, drop=F])
        #To fix a bug when only one germline, 
        snp_allelecount_LinkedGermline_df=as.data.frame(snp_allelecount_LinkedGermline_df)
        ref_allelecount_LinkedGermline_df=as.data.frame(ref_allelecount_LinkedGermline_df)
        
        
        for(sample in tumoursamples)
        {
          
          selectedlinkedGermlines_list = rownames(major_cn_sample_LinkedGermline_df[!is.na(major_cn_sample_LinkedGermline_df[,sample]) ,sample,drop=F])
          selectedlinkedGermlines_list = intersect(selectedlinkedGermlines_list,rownames(snp_allelecount_LinkedGermline_df[!is.na(snp_allelecount_LinkedGermline_df[,sample]),sample,drop=F] ))
          selectedlinkedGermlines_list =  intersect(selectedlinkedGermlines_list,rownames(ref_allelecount_LinkedGermline_df[!is.na(ref_allelecount_LinkedGermline_df[,sample]),sample,drop=F] ))
          
          #Compute the distance to the somatic mutation
          distance_germlines= snp_allelecount_df[selectedlinkedGermlines_list,"Pos", drop=F]
          distance_germlines["distance"] = abs(as.numeric(distance_germlines$Pos) - mut_pos)
          #select the closest germline
          closest_germline = rownames(distance_germlines[distance_germlines$distance== min(distance_germlines$distance, na.rm=T), ])
          
          if(length(closest_germline)==0)
            next
          
          #retrieve its allele counts
          snpwellcount_germlines[sample] = snp_allelecount_LinkedGermline_df[closest_germline, sample]
          refwellcount_germlines[sample] = ref_allelecount_LinkedGermline_df[closest_germline, sample]
        }
        
      }else  if(mode=="SNVOnly"){ 
        snpwellcount_germlines=NA
        refwellcount_germlines=NA
      }else{
        stop("\n\n mode shpu;d be one of SNVOnly, PhasedSNP, FlankingSNP or OptimalSNP")
      }
      
      
      # wellfraction_germlines=snp_allelecount_df[linkedGermlines_list,cifs:nbcolumns_wellfraction]
      
      ##### Preparing the input for the formulw
      
      ##Allele Counts
      varcounts_snvomatic=snpwellcount_somatic[mut, tumoursamples,drop=F] # Somatic variant counts 
      refcounts_snvomatic=refwellcount_somatic[mut, tumoursamples,drop=F] # Somatic reference counts 
      lambda_LinkedGermline<-snpwellcount_germlines #  Germline variant counts 
      mu_LinkedGermline<- refwellcount_germlines # Germline reference counts 
      
      ###Copy number profile (simply the one of the somatic mutation locus)

      major_cn= unlist(major_cn_sample_somatic)
      minor_cn = unlist(minor_cn_sample_somatic)
      
      #Summarising the inputs
      # stop(30)
      Trace=F
      if(Trace ){
        cat("\n\n The inputs are : ")
        cat("\n\t varcounts_snvomatic :\n");print( varcounts_snvomatic)
        cat("\n\t refcounts_snvomatic  :\n");print(  refcounts_snvomatic )
        cat("\n\t  lambda_LinkedGermline :\n");print(lambda_LinkedGermline )
        #stop(20)
        cat("\n\t  mu_LinkedGermline :\n");print(mu_LinkedGermline  )
        cat("\n\t  major_cn :\n");print( major_cn )
        cat("\n\t minor_cn  :\n");print( minor_cn  )
      }
      
      # stop(10)
      ###Calling the prevalence quantification
     # if(detail)
       # detail=2
      
      
      Trace=F
      prev_somatic=getPrevalence(varcounts_snvomatic,refcounts_snvomatic,major_cn,minor_cn, lambda_LinkedGermline , mu_LinkedGermline,  detail ,mode_locus,Trace,LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination=NormalCellContamination,Optimal=Optimal)
      
      
      
      if (detail!=1)
        {
        masterprevalence[mut,names(prev_somatic)] = prev_somatic
        }else{
        
        
        for(sample in tumoursamples)
        {
          prevalence=prev_somatic[[sample]]  
          if(is.null(prevalence) || is.na(prevalence)){
            masterprevalence[mut,paste(sample,"Prevalence",sep="_")] = NA
            next()
          }

          
          masterprevalence[mut,paste(sample,"Prevalence",sep="_")] = prevalence$Prevalence
          masterprevalence[mut,paste(sample,"Context",sep="_")] = prevalence$Context
          masterprevalence[mut,paste(sample,"Germ",sep="_")] = prevalence$DetailedPrevalence["Germ"]
          masterprevalence[mut,paste(sample,"Alt",sep="_")] = prevalence$DetailedPrevalence["Alt"]
          masterprevalence[mut,paste(sample,"Both",sep="_")] = prevalence$DetailedPrevalence["Both"]
          masterprevalence[mut,paste(sample,"ResidualNorm",sep="_")] = prevalence$ResidualNorm
          masterprevalence[mut,paste(sample,"InputValues",sep="_")] = prevalence$InputValues
          
          
        }
        
        
      }
      
      
      
    }
    
 
  masterprevalence
  
  }




#' @export
getPhasedSNPPrevalence<-function( varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp, detail=0,Trace=FALSE,LocusCoverage=TRUE,SomaticCountAdjust=TRUE,NormalCellContamination=NULL, Optimal=TRUE)
  {
  
  
  
    #if(length(varcounts_snp)>1)
    {
      tumoursamples= names(varcounts_snp)
      if (is.null(tumoursamples)){
        tumoursamples = paste("Sample",c(1:length(varcounts_snp)),sep="_")
      }
      
      #To avoid some side error, we transform them into vector
      varcounts_snp=as.vector(varcounts_snp)
      refcounts_snp=as.vector(refcounts_snp)
      varcounts_snv=as.vector(varcounts_snv)
      refcounts_snv = as.vector(refcounts_snv)
      major_cn = as.vector(major_cn)
      minor_cn=as.vector(minor_cn)
      
      names(varcounts_snp) =  tumoursamples
      names(refcounts_snp) =  tumoursamples
      names(varcounts_snv) =  tumoursamples
      names(refcounts_snv) =  tumoursamples
      names(major_cn) =  tumoursamples
      names(minor_cn) =  tumoursamples
      
    }
  

  #Initialisation of prev_S
  if(detail==1)
  { prev_S = list()
  }else{
    prev_S =vector("numeric", length=length(varcounts_snp))
    names(prev_S) = names(varcounts_snp)
    prev_S[prev_S==0]<-NA  
    }

    
    for(sample in tumoursamples)
    {
      
    #  Trace=(sample=="ABpre")
     # if(Trace) cat("\n\n\n Computing the prevalence on sample :", sample)
      args_list=list(
        varcounts_snv=varcounts_snv[sample],
        refcounts_snv=refcounts_snv[sample],
        major_cn=major_cn[sample],
        minor_cn=minor_cn[sample],
        varcounts_snp=varcounts_snp[sample],
        refcounts_snp=refcounts_snp[sample],#/ omega_G[sample] - varcounts_snp[sample],
        detail=1,
        Trace=Trace,
        LocusCoverage=LocusCoverage,
        SomaticCountAdjust=SomaticCountAdjust)
      
      if(anyNA(args_list))
        next
      if(varcounts_snp[sample]+refcounts_snp[sample] ==0)
        next
      if(varcounts_snv[sample]+ refcounts_snv[sample] ==0)
        next
      
      if(Optimal){
        args_list[["LocusCoverage"]] = F
        args_list[["SomaticCountAdjust"]] = F
        prevalence_00 = do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
        
        args_list[["LocusCoverage"]] = F
        args_list[["SomaticCountAdjust"]] = T
        prevalence_01 = do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
        
        args_list[["LocusCoverage"]] = T
        args_list[["SomaticCountAdjust"]] = F
        prevalence_10 = do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
        
        #We Normalise the count first.
        locus_snp=varcounts_snp[sample] + refcounts_snp[sample]
        locus_snv = varcounts_snv[sample] + refcounts_snv[sample]
        newvarcount_snv= varcounts_snv[sample] * locus_snp/locus_snv
        newrefcount_snv= refcounts_snv[sample] * locus_snp/locus_snv
        args_list[["varcounts_snv"]] = newvarcount_snv
        args_list[["refcounts_snv"]] = newrefcount_snv
        prevalence_100 = do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
        
        Bestprevalence= prevalence_00
        if(prevalence_01$Residual < Bestprevalence$Residual )
          Bestprevalence= prevalence_01
        if(prevalence_10$Residual < Bestprevalence$Residual )
          Bestprevalence= prevalence_10
        if(prevalence_100$Residual < Bestprevalence$Residual )
          Bestprevalence= prevalence_100
        
          
        if (Trace){
          cat("\n\n Prevalence 00 :\n"); print(prevalence_00) 
          cat("\n\n Prevalence 01 :\n"); print(prevalence_01) 
          cat("\n\n Prevalence 10 :\n"); print(prevalence_10) 
          cat("\n\n Prevalence 100 :\n"); print(prevalence_100)
          
          cat("\n\n Best Prevalence  :\n"); print(Bestprevalence) 
        }
        
        prevalence = Bestprevalence
        
      }else{
        prevalence=do.call(getPhasedSNPPrevalence_on_singlemutation, args_list)
      }
      

      
     # print(prevalence)
      
      #  if detail, the context and  tree type of prevalence are collapsed else only the prevalence is outputed
      
      prev_S[sample] =NA
      
     # if(length(prevalence) > 0 && !is.na(prevalence))
     #   {
        if(detail==2){
          if(length(prevalence) > 0)
          prev_S[sample] = prevalence$CondensedPrevalence
        } else if(detail==1){
          if(length(prevalence) > 0)
          prev_S[sample] =list(prevalence)
        }else{
          if(length(prevalence) > 0)
          prev_S[sample] =prevalence$Prevalence
        }
     # }else{
     #   prev_S[sample] =NA
     # }

    }
    
    #if(Trace) {cat("\n\n\n Computed prevalences are  :\n")
   # print(prev_S)}
    
    
    prev_S
  }
  
  

#' @export
getFlankingSNPPrevalence<-function( varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp, detail=FALSE,Trace=FALSE,LocusCoverage=1,SomaticCountAdjust=TRUE,NormalCellContamination=NULL)
  {
  
  #For each case, we compute the prevalence twice. By considering the somatic to be phased to the germline SNP or phased with the alternative chromosome harboring the reference of the Germline. The latter is achieved just by 
 
    prevalence_phasedSNP = getPhasedSNPPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp, detail=1,Trace=0,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination=NormalCellContamination)
     prevalence_phasedREF= getPhasedSNPPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn, refcounts_snp,varcounts_snp, detail=1, Trace=0,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination=NormalCellContamination)
     
     if(Trace){
       cat("\n\n Prevalence case phased with SNP \n ")
       print( prevalence_phasedSNP)
       cat("\n Prevalence case phased with REF \n ")
       print( prevalence_phasedREF)
     }
     
     
     
     #Initialisation
     prevalence_flankingSNP=prevalence_phasedSNP
     
     
     # we now have choose the prevalence having the highest residual, if Same tumour we sum the residual across all the sample and we compare, else We go sample by sample.
     #if both have the same residual, we check into sample prevalence if they have the same prevalence, in that case the value of that prevalence will be assigned. Otherwise NA will be assigned.
     
     
     #Previously an argument to the function..
     SameTumour=T
     
     if(SameTumour)
     {
       #We sum the residuals and we choose the prevalence with the lower residuals.
       SNP_residual=sum(as.numeric(unlist(lapply(prevalence_phasedSNP, function(x)  if(!is.na(x[["ResidualNorm"]])) as.numeric(x[["ResidualNorm"]]) else NA ))) ,na.rm=T)
       REF_residual=sum(as.numeric(unlist(lapply( prevalence_phasedREF, function(x)  if(!is.na(prevalence_phasedREF[1])) as.numeric(x[["ResidualNorm"]]) else NA ))),na.rm=T)
       
       if(SNP_residual<1e-8)
         SNP_residual=0
       if(REF_residual<1e-8)
         REF_residual=0
       
       if(SNP_residual<REF_residual){
         prevalence_flankingSNP=prevalence_phasedSNP
       }else if(SNP_residual>REF_residual) {
         prevalence_flankingSNP=prevalence_phasedREF
       } else{
         for(sample in names(prevalence_flankingSNP))
           if (prevalence_phasedSNP[[sample]]$Prevalence != prevalence_phasedREF[[sample]]$Prevalence)
             prevalence_flankingSNP[[sample]] =NA 
       }
       
       
     }else{
       for(sample in names(prevalence_flankingSNP))
       {
         #First any residual less than 1e-8 is considered to be 0.
         
         if(prevalence_phasedSNP[[sample]]$ResidualNorm <1e-8)
           prevalence_phasedSNP[[sample]]$ResidualNorm=0
         if(prevalence_phasedREF[[sample]]$ResidualNorm <1e-8)
           prevalence_phasedREF[[sample]]$ResidualNorm=0  
         
         if(prevalence_phasedSNP[[sample]]$ResidualNorm < prevalence_phasedREF[[sample]]$ResidualNorm){
           prevalence_flankingSNP[[sample]] = prevalence_phasedSNP[[sample]]  
         } else if (prevalence_phasedSNP[[sample]]$ResidualNorm > prevalence_phasedREF[[sample]]$ResidualNorm) {
           prevalence_flankingSNP[[sample]] = prevalence_phasedREF[[sample]]  
         }else if (prevalence_phasedSNP[[sample]]$Prevalence != prevalence_phasedREF[[sample]]$Prevalence){
           prevalence_flankingSNP[[sample]] =NA
           
           warning("The prevalence is not resolved without the knowledge of the Phased Germline at this sample")
         }
       }
     }
     
     
     
     
     
     
     
       if(detail==1){
         prev_S= prevalence_flankingSNP
       }else{
         prev_S =vector("numeric", length=length(prevalence_flankingSNP))
         names(prev_S) = names(prevalence_flankingSNP)
         prev_S[prev_S==0]<-NA 
         for(sample in names(prevalence_flankingSNP))
         {
           if(length(prevalence_flankingSNP[[sample]])==1)
             if(is.na(prevalence_flankingSNP[[sample]]))
             next
           
           if(detail==2){
             prev_S[sample] = prevalence_flankingSNP[[sample]]$CondensedPrevalence             
           }else{
             prev_S[sample] = prevalence_flankingSNP[[sample]]$Prevalence               
           }
         }
       }
       
  if(Trace){
        cat("\n Final prevalence \n ")
    print( prev_S)
  }
  
 prev_S
}



#' @export
getPhasedSNPPrevalence_on_singlemutation<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, varcounts_snp, refcounts_snp, detail=FALSE,Trace=FALSE,LocusCoverage=TRUE,SomaticCountAdjust=TRUE,NormalCellContamination=NULL)
  {
  
  
  varcounts_snv=as.numeric(varcounts_snv)
  refcounts_snv=as.numeric(refcounts_snv)
  major_cn = as.numeric(major_cn)
  minor_cn = as.numeric(minor_cn)
  varcounts_snp = as.numeric(varcounts_snp)
  refcounts_snp = as.numeric(refcounts_snp)
  

  Prevalence=NA
  DetailedPrevvalence=NA

  if(SomaticCountAdjust){
      if(varcounts_snv > varcounts_snp)
        varcounts_snv=varcounts_snp
      if(refcounts_snv < refcounts_snp)
        refcounts_snv = refcounts_snp 
      
      if(refcounts_snv + varcounts_snv < qpois(0.05,max(0,refcounts_snp + varcounts_snp ))||
         refcounts_snv + varcounts_snv > qpois(0.95,max(0,refcounts_snp + varcounts_snp ))){
        
        locus_snp=varcounts_snp + refcounts_snp
        locus_snv = varcounts_snv + refcounts_snv
        
        varcounts_snv = varcounts_snv * locus_snp/locus_snv
        refcounts_snv = refcounts_snv * locus_snp/locus_snv
      }

    }
  

    
    
    
    #We compute the prevalence for the two contexts and we choose the one with the less residual
    if(Trace) cat("\n\n\n Context : C1 (C=0)  SNV after CNA \n **********")
    PrevalenceCond_C1 = getSinglePrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,"C1",Trace,LocusCoverage,NormalCellContamination)
    if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA \n **********")
    PrevalenceCond_C2 = getSinglePrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,"C2",Trace,LocusCoverage,NormalCellContamination)
    
    
    
    PrevalenceCond = PrevalenceCond_C1
    context="C1"
    AlternateResidualNorm = PrevalenceCond_C2["residual"]
    
    if(as.numeric(PrevalenceCond_C2["residual"]) < as.numeric(PrevalenceCond_C1["residual"])){
      PrevalenceCond = PrevalenceCond_C2
      context="C2"
      AlternateResidualNorm = PrevalenceCond_C1["residual"]
    }
    
    
    PrevalenceCond=as.numeric(format(PrevalenceCond,digits=2))
    names(PrevalenceCond) = names(PrevalenceCond_C2)
    
    
    
    AllPrevalences=PrevalenceCond[1:3]
    if(context=="C2"){
      Prevalence=sum(PrevalenceCond["Alt"],PrevalenceCond["Both"],na.rm=T)
    }
    if(context=="C1"){
      Prevalence=PrevalenceCond["Both"]
    }
    
    residualNorm = PrevalenceCond["residual"]
    
    condensedPrevalence=paste( context,Prevalence,paste(AllPrevalences,collapse="|"),residualNorm, sep=":")
    varcounts_snp=as.numeric(format(round(varcounts_snp, 2), nsmall = 2))
    refcounts_snp=as.numeric(format(round(refcounts_snp, 2), nsmall = 2))
    input_values = paste(varcounts_snv,refcounts_snv,major_cn,minor_cn, varcounts_snp, refcounts_snp,sep=":")
    
    if(detail){
      Prevalence_output = list(Context=context,Prevalence=Prevalence,DetailedPrevalence=AllPrevalences,ResidualNorm=residualNorm, AlternateResidualNorm=AlternateResidualNorm, CondensedPrevalence = condensedPrevalence,InputValues=input_values)
    }else{
      Prevalence_output = Prevalence
    }
    
    if(Trace) 
    {
      cat("\n\n\n\t\t ***** Final Prevalence is \n")
      print(Prevalence_output)
    }

  
  
  
  
  
  Prevalence_output
  
}



#' @export
getSNVOnlyPrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, detail=FALSE,Trace=FALSE,NormalCellContamination=NULL)
{
  
    tumoursamples= colnames(varcounts_snv)
  if(is.null(tumoursamples))
    tumoursamples= names(varcounts_snv)   


  if (is.null(tumoursamples))
    tumoursamples = paste("Sample",c(1:length(varcounts_snv)),sep="_")
 
  
  
  
  #To avoid some side error, we transform them into vector
  varcounts_snv=as.vector(varcounts_snv)
  refcounts_snv = as.vector(refcounts_snv)
  major_cn = as.vector(major_cn)
  minor_cn=as.vector(minor_cn)
  
  names(varcounts_snv) =  tumoursamples
  names(refcounts_snv) =  tumoursamples
  names(major_cn) =  tumoursamples
  names(minor_cn) =  tumoursamples

  
  
  
  if(detail==1)
  { prev_S = list()
  }else{
    prev_S =vector("numeric", length=length(varcounts_snv))
    names(prev_S) = names(varcounts_snv)
    prev_S[prev_S==0]<-NA  
  }
  
  
  
  for(sample in tumoursamples)
  {
    if(Trace) cat("\n\n\n Computing the prevalence on sample :", sample)
    args_list=list(
      varcounts_snv=varcounts_snv[sample],
      refcounts_snv=refcounts_snv[sample],
      major_cn=major_cn[sample],
      minor_cn=minor_cn[sample],
      detail=1,
      Trace=Trace,
      NormalCellContamination=NormalCellContamination)
    
    if(anyNA(args_list))
      next
    if(varcounts_snv[sample]+ refcounts_snv[sample] ==0)
      next
    
    prevalence=do.call(getSNVOnlyPrevalence_on_singlemutation, args_list)
 
    # print(prevalence)
    
    #  if detail, the context and  tree type of prevalence are collapsed else only the prevalence is outputed
    if(detail==2){
      prev_S[sample] = prevalence$CondensedPrevalence
    } else if(detail==1){
      prev_S[sample] =list(prevalence)
    }else{
      prev_S[sample] =prevalence$Prevalence
    }
  }
  
  if(Trace) {cat("\n\n\n Computed prevalences are  :\n")
    print(prev_S)}
  
  
  prev_S
}


#' @export
getSNVOnlyPrevalence_on_singlemutation<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, detail=FALSE,Trace=FALSE,NormalCellContamination=NULL)
  {
  

  Prevalence=NA
  DetailedPrevvalence=NA
  sigma=1
  if(Trace) cat("\n The parameters are : ", c(varcounts_snv,refcounts_snv,major_cn,minor_cn))
  
  
  #We compute the prevalence for the two contexts and we choose the one with the less residual
  if(Trace) cat("\n\n\n Context : C1 (C=0)  SNV after CNA \n **********")
  PrevalenceCond_C1 = getPrevalenceSNVOnly(varcounts_snv,refcounts_snv,major_cn,minor_cn,"C1")
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA and sigma=major copy number \n **********")
  PrevalenceCond_C2_major = getPrevalenceSNVOnly(varcounts_snv,refcounts_snv,major_cn,minor_cn,"C2",major_cn)
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA and sigma=minor copy number \n **********")
  PrevalenceCond_C2_minor = getPrevalenceSNVOnly(varcounts_snv,refcounts_snv,major_cn,minor_cn,"C2",minor_cn)
  
  PrevalenceCond = PrevalenceCond_C1
  context="C1"
  if(min(as.numeric(PrevalenceCond_C2_major["residual"]),as.numeric(PrevalenceCond_C2_minor["residual"]) ) < as.numeric(PrevalenceCond_C1["residual"])){

    context="C2"
    PrevalenceCond = PrevalenceCond_C2_major
    sigma=major_cn
    if(as.numeric(PrevalenceCond_C2_minor["residual"]) < as.numeric(PrevalenceCond_C2_major["residual"])){
      PrevalenceCond = PrevalenceCond_C2_minor
      sigma=minor_cn
    }
    
  }
  
  
  PrevalenceCond=as.numeric(format(PrevalenceCond,digits=2))
  names(PrevalenceCond) = names(PrevalenceCond_C1)

  if(Trace){
    cat("\n Prevalence: \n")
    print(PrevalenceCond)
  }
  
  
  AllPrevalences=PrevalenceCond[1:3]
  
  
  if(context=="C2"){
    Prevalence=sum(PrevalenceCond["Alt"],PrevalenceCond["Both"],na.rm=T)
  }
  if(context=="C1"){
    Prevalence=PrevalenceCond["Both"]
  }
  
  residualNorm = PrevalenceCond["residual"]
  
  
  condensedPrevalence=paste( context,Prevalence,paste(AllPrevalences,collapse="|"),residualNorm, sep=":")
  
  
  input_values = paste(varcounts_snv,refcounts_snv,major_cn,minor_cn,sep=":")
  
  if(detail){
    Prevalence_output = list(Context=context,Prevalence=Prevalence,DetailedPrevalence=AllPrevalences,ResidualNorm=residualNorm, CondensedPrevalence = condensedPrevalence,InputValues=input_values)
  }else{
    Prevalence_output = Prevalence
  }
  
  
  if(Trace) {cat("\n\n\n\t\t ***** Final Prevalence is \n")
  print(Prevalence_output)
  }
  Prevalence_output
  
  
  
  
}





#' Generate the matrices C, W and M from a set of parameters.
#' 
#' This is a generic function to generate the matrices of the linear system (see the paper) from the allele counts and the copy number information.
#' 
#' 
#' @param varcounts_snv  A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param refcounts_snv  A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn major copy number  at the locus of the mutation
#' @param minor_cn  minor copy number   at the locus of the mutation 
#' @param varcounts_snp  A count of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param refcounts_snp  A count of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage. 
#' 
#' 
#' @return   the matrices W, C and M for the linear system of prevalence computation.
#'     
#'    
#'     
#' @examples
#' 
#' Matrices = getMatrices(3, 10,2,1,8,5,"C1")
#' 
#' print(Matrices)
#' #$context
#' #[1] "C1"
#' #
#' #$W
#' #          SNP       SNV
#' #SNP 0.6153846 0.0000000
#' #SNV 0.0000000 0.2307692
#' #
#' #$M
#' #    Germ Alt Both
#' #SNP    1   2    2
#' #SNV    0   0    1
#' #
#' #$C
#' #    Germ Alt Both
#' #SNP    2   3    3
#' #SNV    2   3    3
#' 
#' @export
getMatrices<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,LocusCoverage=FALSE){
  total_cn = major_cn + minor_cn
  
#   if((LocusCoverage>0)){
#     if(LocusCoverage==1){
#       Locus_Coverage= refcounts_snp+varcounts_snp
#     }else if (LocusCoverage==2){
#       Locus_Coverage=(refcounts_snp+varcounts_snp + refcounts_snv +varcounts_snv) /2
#     }
# 
#     omega_G = min(1, varcounts_snp/Locus_Coverage)
#     omega_S= min(1, varcounts_snv/Locus_Coverage)
#   }else {
#     if (LocusCoverage!=0)
#       warnings("\n\n LocusCoverage should be one of 0, 1 or 2. 0 will be considered")
#     omega_G = varcounts_snp/(refcounts_snp+varcounts_snp)
#     omega_S= varcounts_snv/(refcounts_snv +varcounts_snv) 
#   }
  
  

    if(LocusCoverage){
      Locus_Coverage= refcounts_snp+varcounts_snp
    omega_G = min(1, varcounts_snp/Locus_Coverage)
    omega_S= min(1, varcounts_snv/Locus_Coverage)
  }else {
    omega_G = varcounts_snp/(refcounts_snp+varcounts_snp)
    omega_S= varcounts_snv/(refcounts_snv +varcounts_snv) 
  }
  
  W=matrix(c(omega_G,0,0,omega_S),ncol=2,nrow=2)
  colnames(W)= c("SNP","SNV")
  rownames(W) = c("SNP","SNV")
  
  C=matrix(nrow=2,ncol=3)
  colnames(C) = c("Germ","Alt","Both")
  rownames(C) = c("SNP","SNV")
  M=C
  
  if(major_cn<minor_cn)
  {
    stop("\n The major copy number ", major_cn, " can not be less than the minor copy number", minor_cn)
  }
  #if the germline VAF is < 0.5 then sigma = major_Cn else minor_CN
  if(omega_G > 0.5) sigma =major_cn
  if(omega_G <= 0.5) sigma =minor_cn
  
  if(context=="C1"){
    M["SNP",] = c(1,sigma,sigma)
    M["SNV",] = c(0,0,1)
    C["SNP",] = c(2,total_cn,total_cn)
    C["SNV",] = c(2,total_cn,total_cn)
  }else{
    M["SNP",] = c(1,1,sigma)
    M["SNV",] = c(0,1,sigma)
    C["SNP",] = c(2,2,total_cn)
    C["SNV",] = c(2,2,total_cn)
  }
  
  list(context=context,W=W,M=M,C=C)
}




#' Generate the matrices C, W and M from a set of parameters under the mode "SNVOnly.
#' 
#' This is a generic function to generate the matrices of the linear system (see the paper) from the allele counts and the copy number information under the SNVOnly mode.
#' 
#' 
#' @param varcounts_snv A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param refcounts_snv A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn major Copy number  at the locus of the mutation
#' @param minor_cn minor copy number   at the locus of the mutation 
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#' @param sigma Copy number of the parental chromosome harboring the mutation. 
#' 
#' 
#' @return   the matrices W, C and M for the linear system of prevalence computation.
#'     
#'    
#'     
#' @examples
#' 
#' Matrices = getMatricesSNVOnly(3,10,2,1,"C1")
#' 
#' print(Matrices)
#' #$context
#' #[1] "C1"
#' #
#' #$W
#' #SNV
#' #SNV 0.2307692
#' #
#' #$M
#' #Germ Alt Both
#' #SNV    0   0    1
#' #
#' #$C
#' #Germ Alt Both
#' #SNV    2   3    3
#' 
#' @export
getMatricesSNVOnly<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,context, sigma=NULL){
  total_cn = major_cn + minor_cn
  #omega_G = varcounts_snp/(refcounts_snp+varcounts_snp)
  omega_S= varcounts_snv/(refcounts_snv +varcounts_snv) 
  W=matrix(c(omega_S),ncol=1,nrow=1)
  colnames(W)= c("SNV")
  rownames(W) = c("SNV")
  
  C=matrix(nrow=1,ncol=3)
  colnames(C) = c("Germ","Alt","Both")
  rownames(C) = c("SNV")
  M=C
  
  if(major_cn<minor_cn)
  {
    stop("\n The major copy number ", major_cn, " can not be less than the minor copy number", minor_cn)
  }
  
  
  #if the germline VAF is < 0.5 then sigma = major_Cn else minor_CN
 # if(omega_G > 0.5) sigma =major_cn
 # if(omega_G <= 0.5) sigma =minor_cn
  
  if(context=="C1"){
   # M["SNP",] = c(1,sigma,sigma)
    M["SNV",] = c(0,0,1)
    #C["SNP",] = c(2,total_cn,total_cn)
    C["SNV",] = c(2,total_cn,total_cn)
  }else{
    if (is.null(sigma)){
      stop("\n\n If the context is C_2, the value of sigma is required. Please provide it or set the context to C1\n\n")
    }else{
      if((sigma!= major_cn) && (sigma != minor_cn))
      {
        stop("\n\n sigma should be either the major copy number either the minor copy number. Please check your parameters and try again\n\n")
      }
    }

    
    #M["SNP",] = c(1,1,sigma)
    M["SNV",] = c(0,1,sigma)
   # C["SNP",] = c(2,2,total_cn)
    C["SNV",] = c(2,2,total_cn)
  }
  
  list(context=context,W=W,M=M,C=C)
}




#' Compute the cellular prevalence of each group of cells
#' 
#' This is a generic function to compute the detailed prevalence of a single mutation using the linear system of  the model.
#' 
#' 
#' @param varcounts_snp A count of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param refcounts_snp A count of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param varcounts_snv A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param refcounts_snv A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn major copy number at the locus of the mutation
#' @param minor_cn minor copy number (or a vector of copy number if multiple tumor samples)
#’ at the locus of the mutation 
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage. 
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment.  
#' @param Trace Print a trace of the eecution.
#' 
#' @return   A list of the three cellular prevalence of each of the three groups of cells
#' 
#' @examples
#' 
#' Prevalences = getSinglePrevalence(3, 10,2,1,8,5,"C1")
#' 
#'  print(Prevalences)
#' # Germ  Alt Both 
#' # 0.4  0.0  0.6 
#' 
#' @seealso \code{\link{getPrevalence}},   \code{\link{getMatrices}}
#' @export
getSinglePrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,Trace=FALSE,LocusCoverage=TRUE,NormalCellContamination=NULL)
  {
  

  
  if(is.na(refcounts_snv) || is.na(refcounts_snp)){
    Prevalence_output=NA
  }else{
    
    
    
  }
  
  
  
  if(Trace){
    cat("\n\n The input :\n ")
    cat(" varcounts_snv :", varcounts_snv," refcounts_snv :", refcounts_snv," major_cn :", major_cn," minor_cn :", minor_cn," varcounts_snp :", varcounts_snp, " refcounts_snp :", refcounts_snp," context :", context)
  }

   
  matrix=getMatrices(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,LocusCoverage=LocusCoverage)
  if(Trace){
    cat("\n\n The matrices :\n ")
    print(matrix)
  }
  
  #We solve a system A*X =B
 
  
  
  #matrix A
  A=matrix(ncol=3,nrow=3)
  A[1,] = c(1,1,1)
  A[c(2,3),] = matrix$W %*% matrix$C-matrix$M
  B=c(1,0,0)
  

  
  # equality constraints
  e = matrix( c(1, 1, 1), nrow=1, byrow=TRUE) # theta_G + theta_A + theta_B = 1
  f = 1
  
  # bounds (first 3 rows theta_G, theta_A, theta_B > 0, second 3 rows -theta_G, -theta_A, -theta_B > -1 (equiv. theta_G, theta_A, theta_b < 1))
  g = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, -1), nrow=6, byrow=TRUE)
  h = matrix(c(0, 0, 0, -1, -1, -1), nrow=6, byrow=TRUE)
  
  #   cat("\n A \n") 
  #   print( A)
  #   cat("\n B \n") 
  #   print( B)
  #   cat("\n E \n") 
  #   print( e)
  #   cat("\n F \n") 
  #   print( f)
  #   cat("\n G \n") 
  #   print( g)
  #   cat("\n H \n") 
  #   print( h)
  # use constrained linear solver
  
  if(Trace){
    
    
    cat("\n \n matrices in  system A*X = B\n")
    cat("\n A \n")
    print(A)
    cat("\n B\n")
    print(B)
   # cat("\n ** Result obtained with solve() of  A*X=B\n")
    #P=solve(A,B)
    #names(P) = c("Germ","Alt","Both")
    #print(P)
    
    
    #cat("\n\n ** Result obtained from lsei without bounds\n") 
   # print(  lsei( A = A, B = B, E = e, F = f))
    
    cat("\n\n ** Result obtained from lsei with bounds\n") 
    print(  lsei( A = A, B = B, E = e, F = f, G = g, H = h))
    
  }

  
  #coutput_with_bounds 
  Prevalence= lsei( A = A, B = B, E = e, F = f, G = g, H = h)
  #print(as.numeric(format(Prevalence$X,digits=2)))
  
  Prevalence=c(as.numeric(format(Prevalence$X,digits=2)),Prevalence$solutionNorm)
  #print(Prevalence)
  Prevalence=as.numeric(format(Prevalence,digits=2))
  names(Prevalence) = c("Germ","Alt","Both","residual")
  
  if(Trace){
    cat("\n\n The prevalences (from lsei() of  limSolve with bounds) :\n ")
    print(Prevalence)
  }
  
 if(!is.null(NormalCellContamination)){
   normalc=as.numeric(NormalCellContamination)
   if(0<=normalc && normalc<=1){
     prevalence=max(1,prevalence/(1-normalc))
   }else{
     stop(" The normal cell contamination rate should take value between 0 and 1.0")
   }

 }
  
  Prevalence
  
  
  
}





#' Compute the cellular prevalence of each group of cells in case of SNVOnly mode
#' 
#' This is a generic function to compute the detailed prevalence of a single mutation using the linear system making the model.
#' 
#' 
#' @param varcounts_snv A count of 
#' alleles supporting the variant sequence of  the somatic mutation 
#' @param refcounts_snv A count of 
#' alleles supporting the reference sequence  of  the somatic mutation
#' @param major_cn  major copy number at the locus of the mutation
#' @param minor_cn minor copy number (or a vector of copy number if multiple tumor samples)
#’ at the locus of the mutation 
#' @param context represents either the situation of a mutation which occurred after the CNV ("C1") or the context of a mutation which occurred before the CNV ("C2"). If not provided, the right context will be estimated from the input
#' @param sigma The parental copy number of  the chromosome harboring the mutation locus. Only needed if the context = C2. Should be either the major copy number either minor copy number
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment. 
#' @param Trace If TRUE, a trace of the execution will be printed
#' 
#' 
#' @return   A list of the three cellular prevalence of each of the three groups of cells
#'     
#' @examples
#' 
#' Prevalences = getPrevalenceSNVOnly(3,10,2,1,"C2",2)
#' 
#' print(Prevalences)
#' #Germ      Alt     Both residual 
#' #0.60     0.31     0.09     0.00 
#' 
#' @seealso \code{\link{getPrevalence}},   \code{\link{getMatricesSNVOnly}}
#' @export
getPrevalenceSNVOnly<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,context,sigma=NULL,Trace=FALSE,NormalCellContamination=NULL)
  {
  
   if(Trace){
     cat("\n\n The input :\n ")
     cat(" varcounts_snv :", varcounts_snv," refcounts_snv :", refcounts_snv," major_cn :", major_cn," minor_cn :", minor_cn,"sigma", sigma, " context :", context)
   }
   
   
  matrix=getMatricesSNVOnly(varcounts_snv,refcounts_snv,major_cn,minor_cn,context,sigma)
  
   if(Trace){
     cat("\n\n The matrices :\n ")
     print(matrix)
   }
   
  #We solve a system A*X =B
  Trace=FALSE
  
  
  if(Trace){
    print(  matrix)
  }
  #matrix A
  A=matrix(ncol=3,nrow=2)
  A[1,] = c(1,1,1)
  A[2,] = matrix$W %*% matrix$C-matrix$M
  B=c(1,0)
  
  # equality constraints
  e = matrix( c(1, 1, 1), nrow=1, byrow=TRUE) # theta_G + theta_A + theta_B = 1
  f = 1
  
  # bounds (first 3 rows theta_G, theta_A, theta_B > 0, second 3 rows -theta_G, -theta_A, -theta_B > -1 (equiv. theta_G, theta_A, theta_b < 1))
  g = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, -1), nrow=6, byrow=TRUE)
  h = matrix(c(0, 0, 0, -1, -1, -1), nrow=6, byrow=TRUE)
  
#   cat("\n A \n") 
#   print( A)
#   cat("\n B \n") 
#   print( B)
#   cat("\n E \n") 
#   print( e)
#   cat("\n F \n") 
#   print( f)
#   cat("\n G \n") 
#   print( g)
#   cat("\n H \n") 
#   print( h)
  # use constrained linear solver
  #coutput_with_bounds 
  Prevalence= lsei( A = A, B = B, E = e, F = f, G = g, H = h)
  Prevalence=c(as.numeric(format(Prevalence$X,digits=2)),Prevalence$solutionNorm)
  names(Prevalence) = c("Germ","Alt","Both","residual")
  
  
  
  if(!is.null(NormalCellContamination)){
    normalc=as.numeric(NormalCellContamination)
    if(0<=normalc && normalc<=1){
      prevalence=max(1,prevalence/(1-normalc))
    }else{
      stop(" The normal cell contamination rate should take value between 0 and 1.0")
    }
    
  }
  
  
  Prevalence
  
  #cat("\n\n\n")
 # print(Prevalence)
}











#' @export
getLocusGermlineMutations<-function(somatic_snp_allelecount_df, snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,phasing_association_df,  tumoursamples,mode="PhasedSNP",  LocusRadius)
{
  
  #We then  retrieve for each somatic mutation the list of germline mutations to consider for the prevalence computation
  # a) if PhasedSNP mode, then the considered germline are the germline mutations phased to the somatic mutation and located within the same locus
  # b) if FlankingSNP mode then the considered germline are the close germlines located within LocusRadius distance from the somatic mutation.
  #c) if OptimalSNP mode tho columns are provided, the first for the phased germline, and if only there is not phasing information for this mutation, then the second column contains the close germlines located within LocusRadius
  
  
  # LinkedGermlineMutation=getLocusGermlineMutations(somatic_snp_allelecount_df, snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,cnv_fraction,phasing_association_df,  samples_to_pool,mode,  LocusRadius)
  
  
  
  Association=("Phased_List" %in% colnames(phasing_association_df)) 
  if(!Association)
  PhasingCode=TRUE
  
  
  #Preparing the data frame to contains the Germlines linked to each somatic mutation.
  LinkedGermlines<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=4)
  LinkedGermlines<-as.data.frame(LinkedGermlines)
  colnames(LinkedGermlines) <- c(colnames(somatic_snp_allelecount_df[1:3]),"LinkedGermlines")
  rownames(LinkedGermlines) <- rownames(somatic_snp_allelecount_df)
  LinkedGermlines[1:3] = somatic_snp_allelecount_df[1:3]
  
  imutation=0 #Simply an iterator
  
  germline_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==1, ]
  
  for (imut in 1:nrow(LinkedGermlines))
  {
    #Mutation name and mutation position
    mut <- rownames(LinkedGermlines[imut,]); 
    mut_pos=as.numeric(LinkedGermlines[imut,"Pos"])
    
    
    
    if(mode=="PhasedSNP"){
      
      if(Association){
        
        if(is.null(phasing_association_df))
          stop("\n\n if mode=PhasedSNP the phasing association matrice should be provided")
        CandidateGermlines = as.character(phasing_association_df[mut,"Phased_List"])
        if(is.null(CandidateGermlines))
          next
        CandidateGermlines<-unlist(strsplit(CandidateGermlines,":"));  
        if (length(unlist(CandidateGermlines))==0)
          next
        
      }else if(PhasingCode){
        
        PhasedSamples=intersect(colnames(snp_allelecount_df[4:ncol(snp_allelecount_df)]), colnames(phasing_association_df))
        if(length(PhasedSamples)==0){
          stop("\n\n\t No samples in common between the columns of the Allele Count and the phasing Association.")
        }
        CandidateGermlines=c()
        for(sample in PhasedSamples)
        {
          mut_phasingcode=as.character(phasing_association_df[mut,sample])
          if(length(mut_phasingcode)==0 || is.na(mut_phasingcode))
            next
          
          CandidateGermlines = unique(c(CandidateGermlines,rownames(phasing_association_df[!is.na(as.character(unlist(phasing_association_df[sample]))) & as.character(unlist(phasing_association_df[sample]))==mut_phasingcode,])))
          
        }
        
      }else{
        stop("\n'n Unknown Phasing Information format. Format should be an Association or a matrice of phasing codes. Check your phasing information matrice and try again")
      }
      
    }else if(mode=="FlankingSNP"){
      CandidateGermlines =rownames(germline_snp_allelecount_df[germline_snp_allelecount_df$Pos >= mut_pos - LocusRadius & germline_snp_allelecount_df$Pos <= mut_pos + LocusRadius, ])
    }else if(mode=="SNVOnly"){
      CandidateGermlines =NA
    }else{
      stop("\n\n The mode parameters shuld be either PhasedSNP either FlankingSNP")
    }
    
    CandidateGermlines=setdiff(CandidateGermlines, mut)
    
    #Now, Among the candidate germline, we want to keep only those present on the same locus with the somatic mutation
    
    submatrix_CandidateGermlines = snp_allelecount_df[CandidateGermlines, 1:3 , drop=F] #retrieve a submatrix of the CandidateGermlines with the chromosome and the position.
    
    submatrix_CandidateGermlines_leftside = submatrix_CandidateGermlines[submatrix_CandidateGermlines$Pos <mut_pos,1:3 , drop=F ] #germlines at the left of the mutation 
    submatrix_CandidateGermlines_rightside = submatrix_CandidateGermlines[submatrix_CandidateGermlines$Pos >mut_pos,1:3 , drop=F ] # germline at the right
    
    #order the leftside from highest to smallest position, leave the rightside from smallest to highest
    orders_phase=order(submatrix_CandidateGermlines_leftside["Pos"],decreasing=T)
    submatrix_CandidateGermlines_leftside=submatrix_CandidateGermlines_leftside[orders_phase,]
    
    
    at_least_one_good_germline=F # is there atleast one good germline kept?
    
    at_least_one_good_germline=F # is there atleast one good germline kept?
    germline_to_exclude=c()
    for(sample in tumoursamples)
    {
      samelocus_germline=c()
      major_som=major_copynumber_df[mut,sample]
      minor_som=minor_copynumber_df[mut,sample]
      
      absence_copynumberprofile=c()
      count_lower_than_somatic<-c() # For more accuracy in case of abundance of germline, someone can choose to exclude germline having an allele count less than the somatic allele count.
      
      for (germ in rownames(submatrix_CandidateGermlines_leftside))
      {
        major_germ=major_copynumber_df[germ , sample]
        minor_germ=minor_copynumber_df[germ, sample]
       
        if (  (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) &&
              (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          break
        }
      }
      
      for (germ in rownames(submatrix_CandidateGermlines_rightside))
      {
        major_germ=major_copynumber_df[germ , sample]
        minor_germ=minor_copynumber_df[germ, sample]
        
        
        if ( (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) && 
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          break
        }
      }
      
      notsamelocus_germline = setdiff(CandidateGermlines, samelocus_germline)
      germline_to_exclude=unique(c(germline_to_exclude,notsamelocus_germline))
    }
    samelocus_germlines= setdiff(CandidateGermlines,germline_to_exclude)
    LinkedGermlines[mut,"LinkedGermlines"] = paste(samelocus_germlines, collapse=":")
    
  }
  
  LinkedGermlines 
  
}







