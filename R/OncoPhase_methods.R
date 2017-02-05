

#' Somatic mutations cellular prevalence on a sample set of mutations
#'
#' This function computes the cellular prevalence of a list of somatic mutations of a tumor.  The function applies OncoPhase linear model to a range of mutations located at a given genomic region or at the whole genome scale.
#' It invokes the function  \code{\link{getPrevalence}}  to compute the cellular prevalence for each mutation of the set.
#' When phasing information are available, the method can computes the prevalence of a somatic
#'   mutation relatively to phased germline SNP under the mode “PhasedSNP”. If the phasing information are not available the mode “SNVOnly” will be used to derive the cellular prevalence. as specified in \code{\link{getPrevalence}}.
#'
#' @param input_df A data frame containing for each mutation the following information (columns or fields) :
#' \describe{
#'        \item{varcounts_snv}{Allele counts supporting the SNV}
#'        \item{refcounts_snv}{Allele counts supporting the reference at the SNV locus}
#'        \item{major_cn}{Major copy number at the SNV locus}
#'        \item{minor_cn}{Minor copy number at the SNV locus}
#'        \item{varcounts_snp}{(Optional) Allele counts supporting the nearby phased SNP. Required if mode= PhasedSNP}
#'        \item{refcounts_snp}{(Optional) Allele counts supporting the reference at the nearby phased SNP. Required if mode= PhasedSNP }
#'      }
#'
#' @param nbFirstColumns Number of first columns in input_df to reproduce in the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation.
#'
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end   e.g "chr22:179800-98767.
#' @param mode The mode under which the prevalence is computed  (Default : Ultimate , alternatives methods  are PhasedSNP and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=Ultimate.
#' @param detail when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations SNV or copy number alteration and cells affected by  both copy number alteration and SNV ). The residual and the linear models inputs and parameters are also reported.
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage.
#' @param SomaticCountAdjust when set to TRUE, varcounts_snv and refcounts_snv might be adjusted if necessary so that they meet the reqirements varcounts_snv <= varcounts_snp, refcounts_snv >= refcounts_snp and varcounts_snv + refcounts_snv ~ Poiss(varcounts_snp + refcounts_snp). Not used if mode=SNVOnly.
#' @param Optimal The model will be run under different configurations  of the parameters LocusCoverage and SomaticCountAdjust. The configuration yielding the optimal residual is then selected and returned.
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment.
#' @param c2_max_residual_treshold Maximum residual threshold under which the context C2 can be inferred.
#' @param c1_ultimate_c2_replacing_treshold Context C1 is inferred if its linear model residual is less than the specified threshold.
#' @param snvonly_max_treshold Maximum threshold the linear model under SNVOnly is considered valid. Is the residual is greater than the value, then PhasedSNP is considered in case the phasing information are available.
#'
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame input_df.
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence dataframe (e.g REF and ALL sequences, ...) }
#'         \item{}{and the following information}
#'   \describe{
#'        \item{Prevalence}{The Cellular Prevalence of the mutation}
#'        \item{Germ}{The proportion of cells with a normal genotype}
#'        \item{Alt}{The proportion of cells with only the CNA if the context C=C1 or with only the SNV if the context C=C2}
#'        \item{Both}{The proportion of cells with both the SNV and the SCNA}
#'        \item{Context}{Context at the mutation. If C1 then the SNV occurred after the SCNA, if C=c2 then the SNV occurred before the SCNA}
#'        \item{solutionNorm}{Residual of the linear model.}
#'        \item{residualNorm}{Constraints residual  representing the sum of absolute values of solutionNorms of equalities and violated inequalities.}
#'        \item{Quality}{Quality of the prevalence calling. H if residual < 1e-05, F if residual < 1e-03 and L if residual > 1e-03}
#'        \item{Alt_Prevalence}{Prevalence estimated by the model if the context were to be the alternative context}
#'        \item{Alt_solutionNorm}{Residual of the linear modelunder the alternative context}
#'        \item{Alt_residualNorm}{Constraints residual  representing  the sum of absolute values of solutionNorms of equalities and violated inequalities under the alternative context}
#'        \item{Mode}{The mode considered for the cellular prevalence computation (either SNVOnly of PhasedSNP)}
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
#' #mutation_id Prevalence   Germ    Alt   Both Context solutionNorm residualNorm Quality Alt_Prevalence
#' #a           a     0.9967 0.0017 0.0017 0.9967      C1 7.718183e-32 2.220446e-16       H         0.9966
#' #b           b     0.8230 0.0890 0.0890 0.8230      C1 1.925930e-32 0.000000e+00       H         0.8200
#' #c           c     0.9000 0.1000 0.0000 0.9000      C1 5.238529e-32 2.220446e-16       H         0.7300
#' #d           d     0.1500 0.4200 0.4200 0.1500      C1 1.972152e-31 4.440892e-16       H         0.1500
#' #e           e     0.4200 0.2900 0.2900 0.4200      C1 5.007418e-32 2.220446e-16       H         0.7100
#' #Alt_solutionNorm Alt_residualNorm         InputValues    Mode      lm_inputs lm_params
#' #a     1.222984e-32     0.000000e+00 151:152:1:1:151:135 SNVOnly 151:152:1:1:C1  0.5:NA:2
#' #b     5.623715e-32     2.220446e-16 123:176:1:1:161:150 SNVOnly 123:176:1:1:C1 0.41:NA:2
#' #c     6.933348e-33     0.000000e+00  94:209:2:1:176:134 SNVOnly  94:209:2:1:C1 0.31:NA:3
#' #d     3.081488e-33     0.000000e+00  23:283:1:1:155:144 SNVOnly  23:283:1:1:C1 0.08:NA:2
#' #e     5.007418e-32     2.220446e-16  60:228:2:0:174:125 SNVOnly  60:228:2:0:C1 0.21:NA:2
#'
#'
#' #'@seealso \code{\link{getPrevalence}}
#' @export
#'
getSamplePrevalence<-function(input_df,mode="Ultimate",  nbFirstColumns=0, region=NULL,detail=TRUE,  LocusCoverage=FALSE,SomaticCountAdjust=FALSE,Optimal=TRUE, c2_max_residual_treshold=Inf, c1_ultimate_c2_replacing_treshold=0.0, snvonly_max_treshold=0.01,verbose=TRUE)
{
  
  #Check the compulsory columns
  compulsory_columns=c("varcounts_snv","refcounts_snv","major_cn","minor_cn")
  
  if (length(setdiff(compulsory_columns,colnames(input_df)))>0){
    stop(" The allele count master matrices should have at least the following headers
         columns : ", compulsory_columns)
    # print(compulsory_columns)
  }
  
  
  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=Ultimate, 3 = OptimalSNP
  numeric_mode=c("SNVOnly", "PhasedSNP","Ultimate")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1, or 2")
    }
  }
  
  if(!is.null(snvonly_max_treshold) && mode != "Ultimate")
    warnings(" The prevalences will be computed under the mode : ", mode, " Your provided parameter snvonly_max_treshold will be ignored")
  
  #Prepare the master matrices for the prevalence.
  masterprevalence=as.data.frame(matrix(nrow=nrow(input_df), ncol=nbFirstColumns+5))
  rownames(masterprevalence) = rownames(input_df)
  if(nbFirstColumns>0){
    masterprevalence[1:nbFirstColumns] = input_df[1:nbFirstColumns]
    colnames(masterprevalence) = c(colnames(input_df[1:nbFirstColumns]),"Prevalence", "Germ","Alt","Both","Context")
  }else{
    colnames(masterprevalence) = c("Prevalence", "Germ","Alt","Both","Context")
  }
  
  if(verbose)
  cat("\n Computing the cellular prevalence of ", nrow(input_df), " somatic mutations (mode = ", mode,")\n")
  for(imut in 1: nrow(input_df))
  {
    varcounts_snv = input_df[imut,"varcounts_snv"]
    refcounts_snv = input_df[imut,"refcounts_snv"]
    minor_cn=input_df[imut,"minor_cn"]
    major_cn=input_df[imut,"major_cn"]
    varcounts_snp=input_df[imut,"varcounts_snp"]
    refcounts_snp=input_df[imut,"refcounts_snp"]
    NormalCellContamination =input_df[imut,"normal_fraction"]
    
    
    InputValues=paste(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp,refcounts_snp,sep=":")

    
    prevalence= getPrevalence(varcounts_snv, refcounts_snv, major_cn, minor_cn, varcounts_snp, refcounts_snp, detail=T, mode=mode ,LocusCoverage=LocusCoverage, SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination=NormalCellContamination,Optimal=Optimal,c2_max_residual_treshold=c2_max_residual_treshold,  c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold,snvonly_max_treshold=snvonly_max_treshold)
    

    if(is.null(prevalence) ||  length(prevalence)==0 || is.na(prevalence) )
      next
    
   
    
    masterprevalence[imut,"Prevalence"] = prevalence$Prevalence
    masterprevalence[imut,"Germ"] = prevalence$DetailedPrevalence["Germ"]
    masterprevalence[imut,"Alt"] = prevalence$DetailedPrevalence["Alt"]
    masterprevalence[imut,"Both"] = prevalence$DetailedPrevalence["Both"]
    masterprevalence[imut,"Context"] = prevalence$Context
    masterprevalence[imut,"solutionNorm"] = prevalence$solutionNorm
    masterprevalence[imut,"residualNorm"] = prevalence$residualNorm
    masterprevalence[imut,"Quality"] = prevalence$Quality
    masterprevalence[imut,"Alt_Prevalence"] = prevalence$Alt_Prevalence
    masterprevalence[imut,"Alt_solutionNorm"] = prevalence$Alt_solutionNorm
    masterprevalence[imut,"Alt_residualNorm"] = prevalence$Alt_residualNorm
   masterprevalence[imut,"TumourPrevalence"] = prevalence$TumourPrevalence
      
    masterprevalence[imut,"InputValues"] = InputValues
    masterprevalence[imut,"Mode"] = prevalence$Mode
    masterprevalence[imut,"lm_inputs"] = prevalence$lm_inputs
    masterprevalence[imut,"lm_params"] = prevalence$lm_params

    
    
  }
  
  if(verbose)
  cat("\n End of somatic mutations cellular prevalence computation\n")
  
  masterprevalence
  
}



#'
#'
#' Computes cellular prevalence at a single mutation point
#'
#' This is a generic function to compute the cellular prevalence of a somatic mutation point using OncoPhase method. The method computes the prevalence of the somatic mutation relatively to phased nearby SNPs whose prevalence are known to be 1.  \code{\link{getPrevalence}} requires the allelic-information of the somatic mutation and the aggregated information of its Phased SNP but the function can also be run in the absence of phasing information (Ultimate mode) or nearby SNP (SNVOnly mode).
#'
#’ The method particularly exhibits an increase in the accuracy when the locus of the SNP is also affected by a somatic copy number alteration (SCNA). The method infer the temporal relationship between the two alterations (C1: SNV occurred after the SCNA; C2: SNV occurred before the SCNA) and computes the detailed prevalence of each of the following group of cells (if detail is set to TRUE) :
#'       \describe{
#'        \item{Germ}{ Cells having a germline genotype  at the locus of the SNV. That is No SNV, no SCNA}
#'        \item{Alt}{ Cells having one alternative of  the two somatic alteration. That is either the SCNA, either the SNV not both.}
#'        \item{Both}{ Cells having both somatic alterations. That is the SNV and the SCNA}
#'     }
#'
#' OncoPhase can be run under three modes:
#'
#'  \describe{
#'        \item{PhasedSNP}{ Phasing information is required. The prevalence is computed relatively to a nearby Phased SNP whose allelic counts should be provided}
#'        \item{SNVOnly}{The prevalence is computed using only the SNV information without the usage of any nearby SNP}
#'        \item{Ultimate}{ This is the default mode. For a given mutation, the method checks if the phasing information is required to compute an accurate cellular prevalence. If it is not, the SNVOnly mode  is used. If instead the phasing information is required the mode is then set to PhasedSNP if allelic counts of a phased nearby SNP are provided.  This is done by first computing the prevalence under the SNVOnly mode. If the data do not fit into this mode (hiogh residual of the linear model), then the prevalence is computed using PhasedSNP mode.
#'     }
#'     }
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
#' @param mode The mode under which the prevalence is computed  (default : Ultimate , alternatives modes are PhasedSNP and SNVOnly). Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=Ultimate
#' @param detail when set to FALSE, the function simply output the cellular prevalence of the somatic mutation. if set to TRUE,  a detailed output is generated containing:
#'  \describe{
#'        \item{Context}{ The inferred associated context : C1 if the SNV occurred after the copy number alteration  or C2 if the SNV occurred after the CNA) }
#'        \item{Prevalence}{The computed somatic mutation cellular prevalence}
#'        \item{DetailedPrevvalence}{the detailed prevalence for each subpopulation of cells  (germline cells (Germ), cells affected by one of the two genomic alterations (Alt), cells affected  by both genomic alterations (Both)}
#'        \item{solutionNorm}{The residual of the linear model representing   the value of the minimized quadratic function at the solution, i.e. ||Ax-b||^2.}
#'        \item{residualNorm}{Residuals from the constraints of the linear model(The sum of absolute values of solutionNorms of equalities and violated inequalities.)}
#'        \item{Quality}{Quality of the prevalence calling. H if residual < 1e-05, F if residual < 1e-03 and L if residual > 1e-03}
#'        \item{Alt_Prevalence}{Prevalence estimated by the model if the context were to be the alternative context}
#'        \item{Alt_solutionNorm}{Residual of the linear modelunder the alternative context}
#'        \item{Alt_residualNorm}{Constraints residual representing  the sum of absolute values of solutionNorms of equalities and violated inequalities under the alternative context}
#'        \item{CondensedPrevalence}{A colon separated list of the above fields (Context, Prevalence, Detailedprevalence and solutionNorm). The detailed prevalence are separated by "|"}
#'        \item{lm_inputs}{Inputs to the linear models separated by “|” and containing the allele count supporting the variant at the SNV, Allele count supporting the reference at the SNV, major and minor copy number, Alleles counts supporting respectively the variant and the reference at the phased SNP if mode=”PhasedSNP” and the context associated to the mutation.  }
#'        \item{lm_params}{Parameters of the linear model separated by “|” and containing the SNV allele fraction, the SNP allele fraction if mode= PhasedSNP, the copy number of the allele harboring the mutation (sigma)}
#'        \item{lm_params}{Parameters of the linear model separated by “|” and containing the SNV allele fraction, the SNP allele fraction if mode= PhasedSNP, the copy number of the allele harboring the mutation (sigma)}
#' }
#' @param Trace if set to TRUE, print the trace of the computation.
#' @param LocusCoverage when set to TRUE, the SNV locus coverage is estimated to the average coverage of the phased SNP and the variant allele fraction is the ratio of the variant allele count over the estimated locus coverage.
#' @param SomaticCountAdjust when set to 1, varcounts_snv and refcounts_snv might be adjusted if necessary so that they meet the rules varcounts_snv <= varcounts_snp, refcounts_snv >= refcounts_snp and varcounts_snv + refcounts_snv ~ Poiss(varcounts_snp + refcounts_snp). Not used if mode=SNVOnly.
#' @param NormalCellContamination If provided, represents the rate of normal cells contaminations in the experiment.
#' @param Optimal If TRUE, the prevalence is computed under all combination of the options SomaticCountAdjust, LocusCoverage and NormalisedCount, the value with the lower residual  is returned as the best prevalence
#' @param Context if provided, the prevalence will be computed strictly under the given context, if not the prevalence is computed under both context and the one yielding the smallest solutionNorm is retained. Default : NULL
#' @param c2_max_residual_treshold Maximum residual threshold under which the context C2 can be inferred. Default: INF.
#' @param c1_ultimate_c2_replacing_treshold Context C1 is inferred if its linear model residual is less than the specified threshold. Default: 0.1.
#' @param snvonly_max_treshold Maximum threshold the linear model under SNVOnly is considered valid. Is the residual is greater than the value, and PhasedSNP is considered in case the phasing information are available.
#' @param SearchContext  When set to true, an optimal  search of the context is done in ta  region of values around the SNV allele fraction and SNP Allele fraction if mode= PhasedSNP. 
#'
#' @return   The cellular prevalence if detail =0, a detailed output if detail = 1, and a condensed output if detail =2. See the usage of the parameter detail above.
#'
#'
#'
#' @examples
#' #Example 1
#' prevalence=getPrevalence(5,10,3,1,16,8)
#' print(prevalence)
#' # 0.86
#' #The above example under mode the mode ultimate compute the prevalence under SNVOnly. 
#' #We can set the mode to PhasedSNP and force the usage of the phasing information.
#' prevalence=getPrevalence(5,10,3,1,16,8,mode="PhasedSNP")
#' print(prevalence)
#' # 0.56
#'
#' #Example 2
#' prevalence = getPrevalence(varcounts_snv=2,refcounts_snv=8,major_cn=2,minor_cn=1,
#' varcounts_snp=8, refcounts_snp=6, detail=TRUE)
#' print(prevalence)
#' # Context
#' # [1] "C1"
#' # 
#' # $Prevalence
#' # Both 
#' # 0.55 
#' # 
#' # $DetailedPrevalence
#' # Germ  Alt Both 
#' # 0.26 0.19 0.55 
#' # 
#' # $solutionNorm
#' # [1] 6.933348e-33
#' # 
#' # $residualNorm
#' # [1] 0
#' # 
#' # $Quality
#' # [1] "H"
#' # 
#' # $Alt_Prevalence
#' # [1] 0.45
#' # 
#' # $Alt_solutionNorm
#' # [1] 4.506676e-31
#' # 
#' # $Alt_residualNorm
#' # [1] 6.661338e-16
#' # 
#' # $CondensedPrevalence
#' # [1] "C1:0.55:0.26|0.19|0.55:6.93334779979405e-33"
#' # 
#' # $lm_inputs
#' # [1] "2:8:2:1:C1"
#' # 
#' # $lm_params
#' # [1] "0.2:NA:3"
#' # 
#' # $Mode
#' # [1] "SNVOnly"
#'
#' #Example: 3
#' prevalence =  getPrevalence(13,5,2,0,47,3,detail=TRUE)
#' print(prevalence)
#' 
#' 
#' # $Context
#' # [1] "C1"
#' # 
#' # $Prevalence
#' # Both 
#' # 0.52 
#' # 
#' # $DetailedPrevalence
#' # Germ  Alt Both 
#' # 0.12 0.36 0.52 
#' # 
#' # $solutionNorm
#' # [1] 5.4e-32
#' # 
#' # $residualNorm
#' # [1] 2.2e-16
#' # 
#' # $Quality
#' # [1] "H"
#' # 
#' # $Alt_Prevalence
#' # [1] 0.38
#' # 
#' # $Alt_solutionNorm
#' # [1] 0.31
#' # 
#' # $Alt_residualNorm
#' # [1] 6.7e-16
#' # 
#' # $CondensedPrevalence
#' # [1] "C1:0.52:0.12|0.36|0.52:5.4e-32"
#' # 
#' # $lm_inputs
#' # [1] "13:5:2:0:47:3:C1"
#' # 
#' # $lm_params
#' # [1] "0.26:0.94:2:2"
#' # 
#' # $Mode
#' # [1] "PhasedSNP"
#' # 
#' #' # Example 4:
#' prevalence= getPrevalence(varcounts_snv=c(6,4,6),refcounts_snv=c(8,8,14),major_cn=c(2,2,2),
#' minor_cn=c(1,0,1),varcounts_snp=c(8,8,8), refcounts_snp=c(6,4,12))
#' print(prevalence)
#' #Sample_1 Sample_2 Sample_3 
#' #1.00     0.67     0.86 
#' #
#' # Example 5:
#' prevalence= getPrevalence(c(6,4,6),c(8,8,14),c(2,2,2),c(1,0,1),c(8,8,8), c(6,4,12), mode="PhasedSNP")
#' print(prevalence)
#' # Sample_1 Sample_2 Sample_3 
#' #0.66     0.67     0.90 
#'
#' @seealso \code{\link{getPrevalence}},  \code{\link{getSamplePrevalence}},   \code{\link{getSinglePhasedSNPPrevalence}}, \code{\link{getSingleSNVOnlyPrevalence}}
#' @export
getPrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, varcounts_snp=NULL, refcounts_snp=NULL,  detail=FALSE, mode="Ultimate",Trace=FALSE,LocusCoverage=FALSE,SomaticCountAdjust=FALSE,Optimal=TRUE, NormalCellContamination=NULL, Context=NULL, SearchContext=TRUE,  c2_max_residual_treshold=Inf, c1_ultimate_c2_replacing_treshold=0.0, snvonly_max_treshold=0.01)
{

  N=length(varcounts_snv) # Number of samples
  
  used_mode=mode
  
  if((length(refcounts_snv)!=N) ||
     (!is.null(varcounts_snp) && (mode !="SNVOnly") &&((length(varcounts_snp) !=N) || (length(refcounts_snp) !=N))) ||
     (length(major_cn) !=N) || (length(minor_cn)!=N) )
    stop("\n\nThe vectors passed as input should have the same size\n\n")

  #set the mode if numeric, 0=SNVOnly, 1 = PhasedSNP, 2=Ultimate
  numeric_mode=c("SNVOnly", "PhasedSNP","Ultimate")
  if(is.numeric(mode))
  {
    if(mode %in% c(0,1,2))
    {
      mode = numeric_mode[mode +1 ]
    }else{
      stop("\n\n Mode parameter, if numeric,  should be either 0, 1 or 2")
    }
  }
  
  if(is.null(varcounts_snp) & mode !="SNVOnly"){
    warning("\n Mode will be set to SNVOnly since no phased SNP information provided")
    mode="SNVOnly"
  }


  if(mode=="PhasedSNP"){
    prev_somatic=getPhasedSNPPrevalence( varcounts_snv,refcounts_snv , major_cn,minor_cn, varcounts_snp , refcounts_snp,detail=TRUE,Trace=Trace,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination= NormalCellContamination,Optimal=Optimal,Context=Context,SearchContext=SearchContext,c2_max_residual_treshold=c2_max_residual_treshold,  c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold)
  }else if(mode=="SNVOnly"){
    prev_somatic=getSNVOnlyPrevalence(varcounts_snv,refcounts_snv ,major_cn,minor_cn, detail=TRUE, Trace=Trace,NormalCellContamination=NormalCellContamination,Context=Context,SearchContext=SearchContext,c2_max_residual_treshold=c2_max_residual_treshold,  c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold)
     }else if(mode=="Ultimate") {
       # We run SNV Only first
       prev_somatic_SNVOnly=getSNVOnlyPrevalence(varcounts_snv,refcounts_snv ,major_cn,minor_cn, detail=TRUE, Trace=Trace,NormalCellContamination=NormalCellContamination,Context=Context,SearchContext=SearchContext,c2_max_residual_treshold=c2_max_residual_treshold,  c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold)
       
       prev_somatic = prev_somatic_SNVOnly
       used_mode="SNVOnly"
       residual_snvonly=Inf
       

       if(length(prev_somatic_SNVOnly)>0){
         if(length(prev_somatic_SNVOnly)==1){
           residual_snvonly=prev_somatic_SNVOnly[[1]]$solutionNorm 
         }else{
           print(prev_somatic_SNVOnly)
           residual_snvonly=0
           for(iprev in 1:length(prev_somatic_SNVOnly))
             residual_snvonly= residual_snvonly+ prev_somatic_SNVOnly[[iprev]]$solutionNorm 
         }
       }
      
      
       


       #If the residual is greater  than the max SNVOnly residual then PhasedSNP is run and the mode with less residual is chosed
       
       if(residual_snvonly > snvonly_max_treshold ){
         prev_somatic_PhasedSNP=getPhasedSNPPrevalence( varcounts_snv,refcounts_snv , major_cn,minor_cn, varcounts_snp , refcounts_snp,detail=TRUE,Trace=Trace,LocusCoverage=LocusCoverage,SomaticCountAdjust=SomaticCountAdjust,NormalCellContamination= NormalCellContamination,Optimal=Optimal,Context=Context,SearchContext=SearchContext,c2_max_residual_treshold=c2_max_residual_treshold,  c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold)
         
         residual_phasedsnp=Inf
         if(length(prev_somatic_PhasedSNP)>0){
           if(length(prev_somatic_PhasedSNP)==1){
             residual_phasedsnp=prev_somatic_PhasedSNP[[1]]$solutionNorm 
           }else{
             residual_phasedsnp=0
             for(iprev in 1:length(prev_somatic_PhasedSNP))
               residual_phasedsnp= residual_phasedsnp+ prev_somatic_PhasedSNP[[iprev]]$solutionNorm 
           }
           
         }


  
         
         
         if( residual_snvonly > residual_phasedsnp ){
           prev_somatic = prev_somatic_PhasedSNP
           used_mode="PhasedSNP"
         }

         
       }
       
     }else{
       stop("parameter mode should be either Ultimate, PhasedSNP or SNVOnly")  
     }
  
  

  #Assign quality
  if(length(prev_somatic)>0){
    for(iprev in 1:length(prev_somatic))
    prev_somatic[[iprev]]["Mode"] = used_mode
 
  if (length(varcounts_snv)==1){
    prev_somatic= prev_somatic[[1]]
    if(detail==FALSE){
      prev_somatic =prev_somatic$Prevalence
    }
  }else{ #MultiSample
    if(detail==FALSE){
      prev_list=c()
      for(iprev in 1:length(prev_somatic))
        prev_list =c(prev_list,prev_somatic[[iprev]]$Prevalence)
      names(prev_list) = names(prev_somatic)
      prev_somatic=prev_list
    }
  }
    
  }
 
  prev_somatic
}












#' @export
getPhasedSNPPrevalence<-function( varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp, detail=0,Trace=FALSE,LocusCoverage=FALSE,SomaticCountAdjust=FALSE,NormalCellContamination=NULL, Optimal=TRUE,Context=NULL,  SearchContext=T, c2_max_residual_treshold=Inf, c1_ultimate_c2_replacing_treshold=0.1)
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
      NormalCellContamination = as.vector(NormalCellContamination)
      names(varcounts_snp) =  tumoursamples
      names(refcounts_snp) =  tumoursamples
      names(varcounts_snv) =  tumoursamples
      names(refcounts_snv) =  tumoursamples
      names(major_cn) =  tumoursamples
      names(minor_cn) =  tumoursamples
      if(!is.null(NormalCellContamination)) names(NormalCellContamination) =  tumoursamples

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
        NormalCellContamination=NormalCellContamination[sample],
        Trace=Trace,
        LocusCoverage=LocusCoverage,
        SomaticCountAdjust=SomaticCountAdjust,
        Context=Context,
        SearchContext=SearchContext,
        c2_max_residual_treshold=c2_max_residual_treshold, 
        c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold)


      if(anyNA(args_list))
        next
      if(varcounts_snp[sample]+refcounts_snp[sample] ==0)
        next
      if(varcounts_snv[sample]+ refcounts_snv[sample] ==0)
        next

      OptimalString=""
      
      if(Optimal){
        
    
        
       # NewBestprevalence= prevalence_1000
       # if(prevalence_1000$solutionNorm >0.01)
      #    {
          args_list[["LocusCoverage"]] = F
          args_list[["SomaticCountAdjust"]] = F
          prevalence_00 = do.call(getPhasedSNPPrevalence_singlesample, args_list)
          OptimalString=paste(OptimalString,"|(B:",prevalence_00$Prevalence," ", prevalence_00$solutionNorm," ", prevalence_00$residualNorm,")|",sep="")
          
          args_list[["LocusCoverage"]] = F
          args_list[["SomaticCountAdjust"]] = T
          prevalence_01 = do.call(getPhasedSNPPrevalence_singlesample, args_list)
          OptimalString=paste(OptimalString,"|(C:",prevalence_01$Prevalence," ", prevalence_01$solutionNorm," ", prevalence_01$residualNorm,")|",sep="")
          
          
          args_list[["LocusCoverage"]] = T
          args_list[["SomaticCountAdjust"]] = F
          prevalence_10 = do.call(getPhasedSNPPrevalence_singlesample, args_list)
          OptimalString=paste(OptimalString,"|(D: ",prevalence_10$Prevalence," ", prevalence_10$solutionNorm," ", prevalence_10$residualNorm,")|",sep="")
          
          
          #We Normalise the count first.
#           locus_snp=varcounts_snp[sample] + refcounts_snp[sample]
#           locus_snv = varcounts_snv[sample] + refcounts_snv[sample]
#           newvarcount_snv= round(varcounts_snv[sample] * locus_snp/locus_snv)
#           newrefcount_snv= round(refcounts_snv[sample] * locus_snp/locus_snv)
#           args_list[["varcounts_snv"]] = newvarcount_snv
#           args_list[["refcounts_snv"]] = newrefcount_snv
#           args_list[["LocusCoverage"]] = F
#           args_list[["SomaticCountAdjust"]] = F
          
          # prevalence_100 = do.call(getPhasedSNPPrevalence_singlesample, args_list)
        
          
          #  OptimalString=paste(OptimalString,"|(D: ",prevalence_100$Prevalence," ", prevalence_100$solutionNorm," ", prevalence_100$residualNorm,")|",sep="")
          
          
          Bestprevalence= prevalence_00
          
          
          if(prevalence_01$solutionNorm < Bestprevalence$solutionNorm )
            Bestprevalence= prevalence_01
          if(prevalence_10$solutionNorm < Bestprevalence$solutionNorm )
            Bestprevalence= prevalence_10
         # if(prevalence_100$solutionNorm < Bestprevalence$solutionNorm )
         #   Bestprevalence= prevalence_100
          
          
          
          NewBestprevalence= Bestprevalence
          
      #  }
        
    


        if (Trace){
          cat("\n\n Prevalence 00 :\n"); print(prevalence_00)
          cat("\n\n Prevalence 01 :\n"); print(prevalence_01)
          cat("\n\n Prevalence 10 :\n"); print(prevalence_10)
        #  cat("\n\n Prevalence 100 :\n"); print(prevalence_100)
          cat("\n\n Best Prevalence  :\n"); print(Bestprevalence)
        }
        

       # cat("\n end 2\n")
        prevalence = NewBestprevalence

      }else{
        prevalence=do.call(getPhasedSNPPrevalence_singlesample, args_list)
      }


   #  if detail, the context and  tree type of prevalence are collapsed else only the prevalence is outputed

      prev_S[sample] =NA


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
  

    }



    prev_S
  }


#' @export
getPhasedSNPPrevalence_singlesample<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, varcounts_snp, refcounts_snp, detail=FALSE,Trace=FALSE,LocusCoverage=TRUE,SomaticCountAdjust=TRUE,NormalCellContamination=NULL, Context=NULL, SearchContext=T, c2_max_residual_treshold=Inf, c1_ultimate_c2_replacing_treshold=0.1)
  {
  varcounts_snv=as.numeric(varcounts_snv)
  refcounts_snv=as.numeric(refcounts_snv)
  major_cn = as.numeric(major_cn)
  minor_cn = as.numeric(minor_cn)
  varcounts_snp = as.numeric(varcounts_snp)
  refcounts_snp = as.numeric(refcounts_snp)
  NormalCellContamination=as.numeric(NormalCellContamination)

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

        varcounts_snv = round(varcounts_snv * locus_snp/locus_snv)
        refcounts_snv = round(refcounts_snv * locus_snp/locus_snv)
      }

    }


  
  #We compute the prevalence for the two contexts and if context=NULL, we choose the one with the less solutionNorm
  if(Trace) cat("\n\n\n Context : C1 (C=0)  SNV after CNA \n **********")
  PrevalenceCond_C1 = getSinglePhasedSNPPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,"C1",Trace,LocusCoverage,NormalCellContamination)
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA \n **********")
  PrevalenceCond_C2 = getSinglePhasedSNPPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,"C2",Trace,LocusCoverage,NormalCellContamination)
  #By default we assign context c1
  
  PrevalenceCond = PrevalenceCond_C1
  
  

  
  
  my_context=Context
  if(SearchContext && is.null(my_context))
  {
    my_context=  BruteForceContextPhasedSNP(varcounts_snv, refcounts_snv, major_cn, minor_cn, varcounts_snp, refcounts_snp,LocusCoverage=LocusCoverage)
    Context=my_context

    #If searchContext On and Context=C2, we check if C2 requirement are met
    if( Context=="C2") #Is C2 requirement met?
    {
      if((PrevalenceCond_C1[[1]]["solutionNorm"] <=  c1_ultimate_c2_replacing_treshold) || (as.numeric(PrevalenceCond_C2[[1]]["solutionNorm"])> c2_max_residual_treshold))
        Context="C1"
    }
    
    #If searchContext On and Context=C1, if C1 residual > SNVOnly treshold, we return between C1 and C2 the context with the smallest residual
    if( Context=="C1") #Is C2 requirement met?
    {
      if(PrevalenceCond_C1[[1]]["solutionNorm"] > 1)
      {
        Context = NULL  #By setting NULL at the next code block the context with the smallet residual will be returned
      }
    }
    
    
    
  }
    
    selectedcontext="C1"
    AlternativeContextsolutionNorm = PrevalenceCond_C2[[1]]["solutionNorm"]
    AlternativeContextresidualNorm = PrevalenceCond_C2[[1]]["residualNorm"]
    AlternativeContextPrevalence = PrevalenceCond_C2[[1]]["Alt"] + PrevalenceCond_C2[[1]]["Both"]

    if(is.null(Context)){
      #Here instead of simply taking the context with the smalle dsolutionNorm we opt for a brute force approach to decide of the context.
     
        
      if(((PrevalenceCond_C1[[1]]["solutionNorm"] >  c1_ultimate_c2_replacing_treshold) && (as.numeric(PrevalenceCond_C2[[1]]["solutionNorm"])<=c2_max_residual_treshold))
         || (PrevalenceCond_C1[[1]]["solutionNorm"] > 1) ){
        
        if(as.numeric(PrevalenceCond_C2[[1]]["solutionNorm"]) < as.numeric(PrevalenceCond_C1[[1]]["solutionNorm"])){
          
          PrevalenceCond = PrevalenceCond_C2
          selectedcontext="C2"
          AlternativeContextsolutionNorm = PrevalenceCond_C1[[1]]["solutionNorm"]
          AlternativeContextresidualNorm = PrevalenceCond_C1[[1]]["residualNorm"]
          AlternativeContextPrevalence = PrevalenceCond_C1[[1]]["Both"]
        }
        
      }
   
      
    }else{
      
       if(Context=="C2"){
        PrevalenceCond = PrevalenceCond_C2
        selectedcontext="C2"
        AlternativeContextsolutionNorm = PrevalenceCond_C1[[1]]["solutionNorm"]
        AlternativeContextresidualNorm = PrevalenceCond_C1[[1]]["residualNorm"]
        AlternativeContextPrevalence = PrevalenceCond_C1[[1]]["Both"]
      }


    }



    AllPrevalences=PrevalenceCond[[1]][1:3]
    if(selectedcontext=="C2"){
      Prevalence=sum(PrevalenceCond[[1]]["Alt"],PrevalenceCond[[1]]["Both"],na.rm=T)
    }
    if(selectedcontext=="C1"){
      Prevalence=PrevalenceCond[[1]]["Both"]
    }
    

    solutionNorm = PrevalenceCond[[1]]["solutionNorm"]
    residualNorm = PrevalenceCond[[1]]["residualNorm"]
    lm_inputs=PrevalenceCond[[2]]
    lm_params=PrevalenceCond[[3]]
    
    condensedPrevalence=paste( selectedcontext,Prevalence,paste(AllPrevalences,collapse="|"),solutionNorm, sep=":")

    #Quality
    if(as.numeric(solutionNorm)<1e-05){
      Qual="H"
    }else if(as.numeric(solutionNorm)<1e-03){
      Qual="F"
    }else{
      Qual="L"
    }
  
    TumourPrevalence = Prevalence
  if(!is.null(NormalCellContamination))
    TumourPrevalence = min(1, Prevalence * 1/(1-NormalCellContamination) )
    
    if(detail){
      Prevalence_output = list(Context=selectedcontext,Prevalence=Prevalence,DetailedPrevalence=AllPrevalences,solutionNorm=as.numeric(solutionNorm), residualNorm= as.numeric(residualNorm), Quality=Qual,Alt_Prevalence=as.numeric(AlternativeContextPrevalence), Alt_solutionNorm=as.numeric(AlternativeContextsolutionNorm),Alt_residualNorm=as.numeric(AlternativeContextresidualNorm),CondensedPrevalence = condensedPrevalence,TumourPrevalence =TumourPrevalence ,lm_inputs=lm_inputs, lm_params=lm_params)
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
getSNVOnlyPrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, detail=FALSE,Trace=FALSE,NormalCellContamination=NULL, Context=NULL,SearchContext=T, c2_max_residual_treshold=Inf, c1_ultimate_c2_replacing_treshold=0.1)
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
      NormalCellContamination=NormalCellContamination,
      Context=Context,
      SearchContext=SearchContext,
      c2_max_residual_treshold=c2_max_residual_treshold, 
      c1_ultimate_c2_replacing_treshold=c1_ultimate_c2_replacing_treshold)

    if(anyNA(args_list))
      next
    if(varcounts_snv[sample]+ refcounts_snv[sample] ==0)
      next

    prevalence=do.call(getSNVOnlyPrevalence_singlesample, args_list)

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
getSNVOnlyPrevalence_singlesample<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn, detail=FALSE,Trace=FALSE,NormalCellContamination=NULL, Context=NULL, SearchContext=T, c2_max_residual_treshold=Inf, c1_ultimate_c2_replacing_treshold=0.1)
  {


  Prevalence=NA
  DetailedPrevvalence=NA
 
  if(Trace) cat("\n The parameters are : ", c(varcounts_snv,refcounts_snv,major_cn,minor_cn))

  
  
  
  #We compute the prevalence for the two contexts and we choose the one with the less solutionNorm
  if(Trace) cat("\n\n\n Context : C1 (C=0)  SNV after CNA \n **********")
  PrevalenceCond_C1 = getSingleSNVOnlyPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,"C1")
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA and sigma=major copy number \n **********")
  PrevalenceCond_C2_major = getSingleSNVOnlyPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,"C2",major_cn)
  if(Trace) cat("\n\n\n Context : C2 (C=1)  SNV before CNA and sigma=minor copy number \n **********")
  PrevalenceCond_C2_minor = getSingleSNVOnlyPrevalence(varcounts_snv,refcounts_snv,major_cn,minor_cn,"C2",minor_cn)
  
  #cat("\n\nPrevalenceCond_C1 \n")
 # print(PrevalenceCond_C1)
  
 # cat("\n\nPrevalenceCond_C2_major \n")
 # print(PrevalenceCond_C2_major)
  
 # cat("\n\nPrevalenceCond_C2_major\n")
 # print(PrevalenceCond_C2_minor)
  
  PrevalenceCond_C2=PrevalenceCond_C2_minor
  if (as.numeric(PrevalenceCond_C2_major[[1]]["solutionNorm"])< as.numeric(PrevalenceCond_C2_minor[[1]]["solutionNorm"]))
    PrevalenceCond_C2=PrevalenceCond_C2_major
    
 # cat("\n\nPrevalenceCond_C2\n")
 # print(PrevalenceCond_C2)
  
  PrevalenceCond = PrevalenceCond_C1
  
  
  
  
  
  sigma=NULL
  my_context=Context
  if(SearchContext & is.null(my_context))
  {
    my_context=  BruteForceContextSNVOnly(varcounts_snv, refcounts_snv, major_cn, minor_cn)
    Context=my_context["Context"]
    sigma = my_context["sigma"]
    
    
    #If searchContext On and Context=C2, we check if C2 requirement are met
    if( Context=="C2") #Is C2 requirement met?
    {
      if((PrevalenceCond_C1[[1]]["solutionNorm"] <=  c1_ultimate_c2_replacing_treshold) || (as.numeric(PrevalenceCond_C2[[1]]["solutionNorm"])> c2_max_residual_treshold))
        Context="C1"
    }
    
    #If searchContext On and Context=C1, if C1 residual > SNVOnly treshold, we return between C1 and C2 the context with the smallest residual
    if( Context=="C1") #Is C2 requirement met?
    {
      if(PrevalenceCond_C1[[1]]["solutionNorm"] > 1)
      {
        Context = NULL  #By setting NULL at the next code block the context with the smallet residual will be returned
      }
    }
  
    
  }
  
  selectedcontext="C1"
  AlternativeContextsolutionNorm = PrevalenceCond_C2[[1]]["solutionNorm"]
  AlternativeContextresidualNorm = PrevalenceCond_C2[[1]]["residualNorm"]
  AlternativeContextPrevalence = PrevalenceCond_C2[[1]]["Alt"] + PrevalenceCond_C2[[1]]["Both"]
  
  if(is.null(Context)){
    #Here instead of simply taking the context with the smalle dsolutionNorm we opt for a brute force approach to decide of the context.

    if(((PrevalenceCond_C1[[1]]["solutionNorm"] >  c1_ultimate_c2_replacing_treshold) && (as.numeric(PrevalenceCond_C2[[1]]["solutionNorm"])<=c2_max_residual_treshold))
       || (PrevalenceCond_C1[[1]]["solutionNorm"] > 1) ){
      

      if(as.numeric(PrevalenceCond_C2[[1]]["solutionNorm"]) < as.numeric(PrevalenceCond_C1[[1]]["solutionNorm"])){
        
        PrevalenceCond = PrevalenceCond_C2
        selectedcontext="C2"
        AlternativeContextsolutionNorm = PrevalenceCond_C1[[1]]["solutionNorm"]
        AlternativeContextresidualNorm = PrevalenceCond_C1[[1]]["residualNorm"]
        AlternativeContextPrevalence = PrevalenceCond_C1[[1]]["Both"]
      }
      
    }

    
  }else{
    
    if(Context=="C2"){
      PrevalenceCond = PrevalenceCond_C2
      selectedcontext="C2"
      AlternativeContextsolutionNorm = PrevalenceCond_C1[[1]]["solutionNorm"]
      AlternativeContextresidualNorm = PrevalenceCond_C1[[1]]["residualNorm"]
      AlternativeContextPrevalence = PrevalenceCond_C1[[1]]["Both"]
    }
    
    
  }
  
  
  
  AllPrevalences=PrevalenceCond[[1]][1:3]
  if(selectedcontext=="C2"){
    Prevalence=sum(PrevalenceCond[[1]]["Alt"],PrevalenceCond[[1]]["Both"],na.rm=T)
  }
  if(selectedcontext=="C1"){
    Prevalence=PrevalenceCond[[1]]["Both"]
  }
  
  
  solutionNorm = PrevalenceCond[[1]]["solutionNorm"]
  residualNorm = PrevalenceCond[[1]]["residualNorm"]
  lm_inputs=PrevalenceCond[[2]]
  lm_params=PrevalenceCond[[3]]
  
  condensedPrevalence=paste( selectedcontext,Prevalence,paste(AllPrevalences,collapse="|"),solutionNorm, sep=":")
  
  #Quality
  if(as.numeric(solutionNorm)<1e-05){
    Qual="H"
  }else if(as.numeric(solutionNorm)<1e-03){
    Qual="F"
  }else{
    Qual="L"
  }
  
  TumourPrevalence = Prevalence
  if(!is.null(NormalCellContamination))
    TumourPrevalence = min(1, Prevalence * 1/(1-NormalCellContamination) )
  
  
  if(detail){
    Prevalence_output = list(Context=selectedcontext,Prevalence=Prevalence,DetailedPrevalence=AllPrevalences,solutionNorm=as.numeric(solutionNorm), residualNorm= as.numeric(residualNorm), Quality=Qual,Alt_Prevalence=as.numeric(AlternativeContextPrevalence), Alt_solutionNorm=as.numeric(AlternativeContextsolutionNorm),Alt_residualNorm=as.numeric(AlternativeContextresidualNorm),CondensedPrevalence = condensedPrevalence,TumourPrevalence=TumourPrevalence,lm_inputs=lm_inputs, lm_params=lm_params)
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
#' Matrices = getMatricesPhasedSNP(3, 10,2,1,8,5,"C1")
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
getMatricesPhasedSNP<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,LocusCoverage=FALSE){
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


  if(major_cn<minor_cn)
  {
    stop("\n The major copy number ", major_cn, " can not be less than the minor copy number", minor_cn)
  }
  #if the germline VAF is < 0.5 then sigma = major_Cn else minor_CN
  if(omega_G > 0.5) sigma =major_cn
  if(omega_G <= 0.5) sigma =minor_cn

  
  #Maintenance, lets fix Omega_G to sigma/(major +minor)
  

  
  W=matrix(c(omega_G,0,0,omega_S),ncol=2,nrow=2)
  colnames(W)= c("SNP","SNV")
  rownames(W) = c("SNP","SNV")
  
  C=matrix(nrow=2,ncol=3)
  colnames(C) = c("Germ","Alt","Both")
  rownames(C) = c("SNP","SNV")
  M=C


  
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
  
  lm_inputs=paste(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,sep=":")
  lm_params=paste(round(omega_S,digits=2),round(omega_G,digits=2),sigma,total_cn,sep=":")

  list(context=context,W=W,M=M,C=C,lm_inputs=lm_inputs,lm_params=lm_params )
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
  if(is.null(sigma))
  sigma=NA

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

  lm_inputs=paste(varcounts_snv,refcounts_snv,major_cn,minor_cn, context,sep=":")
  lm_params=paste(round(omega_S,digits=2),sigma,total_cn,sep=":")
  list(context=context,W=W,M=M,C=C, lm_inputs=lm_inputs,lm_params=lm_params)
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
#' Prevalences = getSinglePhasedSNPPrevalence(3, 10,2,1,8,5,"C1")
#'
#'  print(Prevalences)
#' # Germ  Alt Both
#' # 0.4  0.0  0.6
#'
#' @seealso \code{\link{getPrevalence}},   \code{\link{getMatricesPhasedSNP}}
#' @export
getSinglePhasedSNPPrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,Trace=FALSE,LocusCoverage=TRUE,NormalCellContamination=NULL)
  {


  if(is.na(refcounts_snv) || is.na(refcounts_snp)){
    Prevalence_output=NA
  }else{



  }



  if(Trace){
    cat("\n\n The input :\n ")
    cat(" varcounts_snv :", varcounts_snv," refcounts_snv :", refcounts_snv," major_cn :", major_cn," minor_cn :", minor_cn," varcounts_snp :", varcounts_snp, " refcounts_snp :", refcounts_snp," context :", context)
  }


  matrix=getMatricesPhasedSNP(varcounts_snv,refcounts_snv,major_cn,minor_cn,varcounts_snp, refcounts_snp,context,LocusCoverage=LocusCoverage)
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

  
 # print(Prevalence)
  
  Prevalence=c(as.numeric(format(Prevalence$X,digits=2)),Prevalence$solutionNorm,Prevalence$residualNorm)
#  print(Prevalence)
  Prevalence[1:5]=as.numeric(format(Prevalence[1:5],digits=2))
  names(Prevalence) = c("Germ","Alt","Both","solutionNorm","residualNorm")

  if(Trace){
    cat("\n\n The prevalences (from lsei() of  limSolve with bounds) :\n ")
    print(Prevalence)
  }

 if(!is.null(NormalCellContamination)){
   
  
   #To implement
   
   

 }

  
  list(Prevalence=Prevalence, lm_inputs = matrix$lm_inputs, lm_params=matrix$lm_params)



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
#' Prevalences = getSingleSNVOnlyPrevalence(3,10,2,1,"C2",2)
#'
#' print(Prevalences)
#' #Germ      Alt     Both solutionNorm
#' #0.60     0.31     0.09     0.00
#'
#' @seealso \code{\link{getPrevalence}},   \code{\link{getMatricesSNVOnly}}
#' @export
getSingleSNVOnlyPrevalence<-function(varcounts_snv,refcounts_snv,major_cn,minor_cn,context,sigma=NULL,Trace=FALSE,NormalCellContamination=NULL)
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
  Prevalence=c(as.numeric(format(Prevalence$X,digits=2)),Prevalence$solutionNorm,Prevalence$residualNorm)
  names(Prevalence) = c("Germ","Alt","Both","solutionNorm","residualNorm")



  if(!is.null(NormalCellContamination)){
    
    #For latter : 
#     normalc=as.numeric(NormalCellContamination)
#     if(0<=normalc && normalc<=1){
#       prevalence=max(1,prevalence/(1-normalc))
#     }else{
#       stop(" The normal cell contamination rate should take value between 0 and 1.0")
#     }

  }

  


  list(Prevalence=Prevalence, lm_inputs = matrix$lm_inputs, lm_params=matrix$lm_params)

  #cat("\n\n\n")
 # print(Prevalence)
}










#' @export
BruteForceContextPhasedSNP<-function(varcount_snv,refcount_snv,major_cn,minor_cn,varcount_snp=NULL,refcount_snp=NULL,LocusCoverage=F)
  {
  #we compute for +/- the allele fraction to catch up the majority context.
  lowerbound = max(0,round(varcount_snv/(refcount_snv+varcount_snv),digits=2) - 0.1)
  upperbound = min(1,round(varcount_snv/(refcount_snv+varcount_snv),digits=2) + 0.1)
  
  #cat("\n",varcount_snv,refcount_snv,major_cn,minor_cn,varcount_snp,refcount_snp)
  #cat("\n lowerbound", lowerbound)
  #cat("\n upperbound", upperbound)
  
  if(is.na(lowerbound) || is.na(upperbound)){
    mycontext=NULL
  }else{
    #generate the admissible SNV counts
    snv_counts= unique(round(seq(lowerbound,upperbound, 0.01 ) * (refcount_snv+varcount_snv)))
    
    sumResidC1=0
    sumResidC2=0
    numC1=0
    numC2=0
    for (varcount in snv_counts){


        p_C1=getSinglePhasedSNPPrevalence(varcount,refcount_snv,major_cn,minor_cn,varcount_snp,refcount_snp,"C1",LocusCoverage=LocusCoverage)
        p_C2=getSinglePhasedSNPPrevalence(varcount,refcount_snv,major_cn,minor_cn,varcount_snp,refcount_snp,"C2",LocusCoverage=LocusCoverage)
        
      #  print(p_C1)
        
        if(length(p_C1)==0 ||length(p_C2)==0 )
          next
        
        p = p_C1
        selectContext="C1"
        
        if(as.numeric(p_C2[[1]]["solutionNorm"]) < as.numeric(p_C1[[1]]["solutionNorm"])){
          p  = p_C2    
          selectContext="C2"
        }
        
      #  print(p)
 

      if(selectContext=="C1"){
        numC1=numC1+1
        sumResidC1 = sumResidC1 + p[[1]]["solutionNorm"]
      }else{
        numC2=numC2+1
        sumResidC2 = sumResidC2 + p[[1]]["solutionNorm"]
      }
    }
    if(numC2==0){
      mycontext = "C1"
    }else if (numC1==0){
      mycontext = "C2"
    }else {
      avgResidC1= sumResidC1/numC1
      avgResidC2 = sumResidC2/numC2
      if(avgResidC1 <= avgResidC2){
        mycontext = "C1"
      }else{
        mycontext = "C2"
      }
    }
  }
 
  
  
  mycontext
}



#' @export
BruteForceContextSNVOnly<-function(varcount_snv,refcount_snv,major_cn,minor_cn)
{
  #we compute for +/- the allele fraction to catch up the majority context.
  lowerbound = max(0,round(varcount_snv/(refcount_snv+varcount_snv),digits=2) - 0.1)
  upperbound = min(1,round(varcount_snv/(refcount_snv+varcount_snv),digits=2) + 0.1)
  
  #cat("\n lowerbound", lowerbound)
  #cat("\n upperbound", upperbound)
  
  mycontext=NULL
  sigma=NULL
  if(!is.na(lowerbound) && !is.na(upperbound)){
    #generate the admissible SNV counts
    snv_counts= unique(round(seq(lowerbound,upperbound, 0.01 ) * (refcount_snv+varcount_snv)))
    
    sumResidC1=0
    sumResidC2major=0
    sumResidC2minor=0
    numC1=0
    numC2major=0
    numC2minor=0
    
    for (varcount in snv_counts){

        p_C1=getSingleSNVOnlyPrevalence(varcount_snv,refcount_snv,major_cn,minor_cn,"C1")
        p_C2_major=getSingleSNVOnlyPrevalence(varcount_snv,refcount_snv,major_cn,minor_cn,"C2",major_cn)
        p_C2_minor=getSingleSNVOnlyPrevalence(varcount_snv,refcount_snv,major_cn,minor_cn,"C2",minor_cn)
        
        p=p_C1
        
    
        
        selectContext="C1"
        if(min(as.numeric(p_C2_major[[1]]["solutionNorm"]),as.numeric(p_C2_minor[[1]]["solutionNorm"]) ) < as.numeric(p_C1[[1]]["solutionNorm"])){
          selectContext="C2"
          if(as.numeric(p_C2_major[[1]]["solutionNorm"]) < as.numeric(p_C2_minor[[1]]["solutionNorm"] )){
            p = p_C2_major
            sigma=major_cn
          }else{
            p = p_C2_minor
            sigma=minor_cn
          }
        }
          
        
        
      
      if(selectContext=="C1"){
        numC1=numC1+1
        sumResidC1 = sumResidC1 + p[[1]]["solutionNorm"]
      }else{
        if(sigma==major_cn){
          numC2major=numC2major+1
          sumResidC2major = sumResidC2major + p[[1]]["solutionNorm"]
        }else{
          numC2minor=numC2minor+1
          sumResidC2minor = sumResidC2minor + p[[1]]["solutionNorm"]
        }

      }
    }
    
    
  #  print(c(numC1,numC2major,numC2minor))
    
    avgResidC1= sumResidC1/numC1
    avgResidC2major = sumResidC2major/numC2major
    avgResidC2minor = sumResidC2minor/numC2minor
    
    
 #   print(c(avgResidC1, avgResidC2minor, avgResidC2major))
    
    if(!is.na(avgResidC2major) && min(avgResidC1, avgResidC2minor, avgResidC2major,na.rm=T ) == avgResidC2major){
      mycontext = "C2"
      sigma=major_cn
    }else  if(!is.na(avgResidC2minor) &&  min(avgResidC1, avgResidC2minor, avgResidC2major,na.rm=T  ) == avgResidC2minor){
      mycontext = "C2"
      sigma=minor_cn
    }else if(!is.na(avgResidC1) && min(avgResidC1, avgResidC2minor, avgResidC2major,na.rm=T  ) == avgResidC1)
      mycontext = "C1"
 
  }
  
  list(Context=mycontext,sigma=as.numeric(sigma))
}

















# 
# 
# 
# 
# #' chr22_XYZ101 : Patient XYZ101 data from chromosome 22 .
# #'
# #' A generated dataset containing allele counts, haplotype phasing and copy number information on chromosome 22
# #' for a patient XYZ101 (Data created) .
# #'
# #' Also available on the same patient : chr10, chr15 and chr18
# #'
# #'
# #' @format Contains the following data :  :
# #' \describe{
# #'  \item{tumoursamples}{The list of tumor samples of the study}
# #'  \item{SNP_allelecount_df}{Data frame containing the count of allele supporting the variant of each mutations
# #'     \describe{
# #'       \item{}{Chrom : Chromosomes }
# #'        \item{}{Start : Starting position}
# #'        \item{}{End : Position of the mutation}
# #'        \item{}{Vartype : variant Type}
# #'        \item{}{IsGermline : is the mutation a Germline SNP or a Somatic mutation}
# #'        \item{}{Ref : Reference sequence}
# #'        \item{}{All : Variant sequence}
# #'        \item{}{One entry per tumor samples}
# #'     }
# #'  }
# #'  \item{ref_allelecount_df}{Data frame containing the count of allele supporting the reference at each mutation}
# #'  \item{phasing_association_df}{A data frame containing for each somatic mutations, a colon separated  list of Germline mutations phased to it. }
# #'  \item{major_copynumber_df}{A data frame containing the major copy number of the mutations at each tumor samples }
# #'  \item{minor_copynumber_df}{A data frame containing the minor copy number of the mutations at each tumor samples }
# #'  \item{minor_copynumber_df}{A data frame containing the normal cell contamination rate for each mutations  at each tumor samples }
# #' }
# "chr22_XYZ101"






