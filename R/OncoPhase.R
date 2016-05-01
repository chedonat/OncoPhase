############################################################
#
# OncoPhase : files Organization 
# data.R - file to generate the Roxygen help documents for the data
# OncoPhase.R -  this file is the main file of the package
# OncoPhase_methods.R - This files contains the main implementations of the prevalence computation
# OncoPhase_EMmodel.R - This files contains the expecattion-maximisation model




#Roxigen help for the package (top help)

#' OncoPhase package for somatic mutations cellular prevalence quantification using haplotype phasing
#' 
#' The main function for somatic mutation cellular prevalence computation are    \code{\link{getPrevalence}} and  \code{\link{getPrevalenceMultiSamples}} . \code{\link{getPrevalence}}  compute the cellular prevalence at a particular mutations under 4 distincts modes : PhasedSNP, FlankingSNP, OptimalSNP and SNVOnly.  \code{\link{getPrevalenceMultiSamples}}  computes the cellular prevalence of a list of mutations located at a given region of the genome. It can also work on a whole genome scale.
#' See the manual and examples at \code{\link{getPrevalence}} and  \code{\link{getPrevalenceMultiSamples}}  for more details.
#' 
#'  The function \code{\link{getPrevalenceLinear}} compute the prevalence of a given mutation by directly solving the linear system associated to the model.
#' 
#' Input data for simple case studies can be generated with the function \code{\link{build_casestudy}}.
#' 
#' The package include experimental data for chromosome 10, 15, 18 and 22 for two patients retrieved from a parallel clinical study.
#'  (see for example  \code{\link{chr10_OP1019}} and \code{\link{chr22_OP1019}} )
#' 
#' For more detailed information on usage, see the package vignette, by typing
#' \code{vignette("OncoPhase")}. All support questions should be emailed to the authors.
#'
#' @references
#'
#' OncoPhase reference:
#' 
#' OncoPhase: A package for computing Somatic Mutation cellular Prevalence in cancer using haplotype phasing. Bioinformatics 2016. Submitted
#'
#' OncoPhase reference:
#' 
#' ” Ovarian cancer haplotype sequencing reveals ubiquitous SOX2 overexpression in the premalignant fallopian tube epithelium”
#'
#' @author Donatien Chedom-Fotso, Ahmed Ahmed, Christopher Yau.
#' 
#' @docType package
#' @name aOncoPhase
#' @aliases OncoPhase-package
#' @keywords package
#' @import limSolve
NULL



#'  @export
hg19_dfsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, 
                  chr5=180915260, chr6=171115067, chr7=159138663, chr8=146364022, 
                  chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, 
                  chr13=115169878, chr14=107349540,chr15=102531392, chr16=90354753,
                  chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,
                  chr21=48129895,chr22=51304566, chrX=155270560, chrY=59373566,
                  chrM=16571)





# Welcome message
#################

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\n************\nWelcome to OncoPhase package\n************")
}

# Set up some custom options
############################

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.OncoPhase <- list(
    OncoPhase.path = "~/R",
    OncoPhase.install.args = "",
    OncoPhase.name = "Donatien Chedom-Fotso",
    OncoPhase.desc.author = '"Donatien Chedom-Fotso <donatien.chedom@gmail.com> 
    [aut, cre]"',
    OncoPhase.desc.license = "GPL-2",
    OncoPhase.desc.suggests = NULL,
    OncoPhase.desc = list(),
    OncoPhase.PatientName = "NoName",
    OncoPhase.WD = getwd()
  )
  toset <- !(names(op.OncoPhase) %in% names(op))
  if(any(toset)) options(op.OncoPhase[toset])
  
  invisible()
}



# Main prevalence function
##########################
##########################
##########################



setClass("SampleList" )

setClass("SampleCnvAssociation" )



Sample <- setClass(
  #Set the name of the class
  "Sample",
  
  #Define the slots
  slots = c(
    name = "character",
    ID = "character",
    type = "character",
    cnvprofile = "character",
    group = "character",
    alleletype = "character",
    variantfiletype = "character",
    variantfilename = "character",
    phasinganalysis = "logical",
    prevalenceanalysis = "logical"
  ),
  
  #Set the defaults values for the slots
  prototype = list(
    name = "",
    type = "tumour",
    variantfiletype = "varfile",
    alleletype = "wells",
    phasinganalysis = TRUE,
    prevalenceanalysis = TRUE
  ),
  
  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  
  validity = function(object)
  {
    if(! object@type %in% c("tumour", "normal")) {
      return(" Incorrect sample type assigned. The sample type should be one of ")
    }
    if(! object@variantfiletype %in% c("varfile", "mastervarfile", "vcf")) {
      return(" Incorrect sample type assigned. The sample type should be one of varfile, mastervarfile, vcf")
    }
    
  }
)


initsSamples<-function(samples_df)
{
  sampleList<-list()
  
  # Control of the entries, check the match between the list of properties and the 
  #column 
  # of the samples data frame
  slotList = slotNames("Sample")
  colList =  colnames(samples_df)
  
  presentSlots = intersect(slotList, colList)
  absentSlots = setdiff(slotList, colList)
  toignoreFields = setdiff(colList, slotList)
  
  if (length(presentSlots)>0){
    cat("\n\n\t\t The following properties of samples wil be retrieved : \n\t\t")
    print(presentSlots)
  }
  
  if (length(absentSlots)>0){
    cat("\n\n\t\t The following properties of samples are not present in the fields 
        provided : \n\t\t")
    print(absentSlots)
  }
  
  if (length(toignoreFields)>0){
    cat("\n\n\t\t The following fields provided ar not part of a sample property : \n\t\t")
    print(toignoreFields)
  }
  
  if (length(presentSlots)==0){
    stop("\n\n\t\t Error, No sample property  provided in the data frame, check your input")
  }
  
  
  
  
  
  
  
  for (isample in 1:nrow(samples_df)){
    newSample=Sample(
      name = as.character(samples_df[isample,"name"]),
      ID =  as.character(samples_df[isample,"ID"]),
      type =  as.character(samples_df[isample,"type"]),
      cnvprofile =  as.character(samples_df[isample,"cnvprofile"]),
      group =  as.character(samples_df[isample,"group"]),
      alleletype =  as.character(samples_df[isample,"alleletype"]),
      variantfiletype =  as.character(samples_df[isample,"variantfiletype"]),
      variantfilename =  as.character(samples_df[isample,"variantfilename"]),
      phasinganalysis =  as.logical(samples_df[isample,"phasinganalysis"]),
      prevalenceanalysis =  as.logical(samples_df[isample,"prevalenceanalysis"])
    )
    
    sampleList<-append(sampleList, newSample)
  }
  
  return(sampleList)
}


setGeneric("setProperty",
           def=function(theObject, xname, xvalue){
             standardGeneric("setProperty")
           }
)

setMethod(f="setProperty",
          signature = "Sample",
          definition = function(theObject, xname, xvalue){
            slot(theObject,xname) = xvalue
            return(theObject)
          }
)

#' @export
numeric_column<-function(df,tumoursamples)
{
  for (col in tumoursamples)
  {
    df[col] =  as.numeric(unlist(df[col]))
  }
  df
}



#' Somatic mutations cellular prevalence using haplotype phasing.
#' 
#' This is a generic function to compute the cellular prevalence of somatic mutations in
#'  cancer using haplotype phasing.  The function applies the model to a range of mutations located at a given genomic region or at the whole genome scale. The model computes the prevalence of a somatic
#'   mutation relatively to close and eventually phased germline mutations. It uses three sources
#'    of information as input : The allelic counts, the phasing information and the 
#'    copy number alteration.  Multiple tumor samples can be provided for the prevalence computation.
#' 
#' @param snp_allelecount_df A data frame containing for each mutation the  allelic 
#' counts of the variant at each tumor samples. The data frame should contains at least the following three columns among its firsts columns: Chrom (The mutation
#'  chromosome) , End (The mutation position) and IsGermline (is the mutation a germline
#'   or somatic mutation).
#' @param ref_allelecount_df A data frame containing for each mutation the allelic count
#'  of the reference at each tumor sample. The data frame should contains at least the following three columns among its firsts columns:  Chrom (The mutation
#’ chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic mutation)  
#' @param phasing_association_df A data frame containing for each somatic mutation, 
#' a colon separated list of germline SNP phased to it.
#' @param major_copynumber_df A data frame containing for each mutation, its  major
#’ chromosomal copy number at each tumor samples.
#' @param minor_copynumber_df A data frame containing for each mutation the minor
#'  chromosomal copy number at each tumor samples.
#' @param CNVfraction_df, If provided, represents a data frame containing for each mutation,  the fraction of
#'  cells affected by the copy number alteration. Used only if the mode is "PhasedSNP" and formula is "General".
#' @param nbFirstColumns Number of first columns in snp_allelecount_df to reproduce in
#'  the output dataframe e.g: Chrom, Pos, Vartype. Columns from  nbFirstColumns +1 to the last column should contains the information needed for the prevalence computation at each tumour sample
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end 
#' e.g "chr22:179800-98767 
#' @param tumoursamples : The list of tumor samples to consider for the prevalence
#’ computation.  This samples should be present as column header in the data frame
#'   snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df
#'   and  CNVfraction_df. If not provided, the headers from nbFirstColumns + 1 to 
#'   the last column of snp_allelecount_df is retrieved and its intersection with the
#' other inputted data frames headers is considered.
#' @param mode The mode under which the prevalence is computed  (default : PhasedSNP , alternatives methods  are FlankingSNP, OptimalSNP,and SNVOnly).  Can also be provided as a numeric 0=SNVOnly, 1= PhasedSNP, 2=FlankingSNP and 3 = OptimalSNP
#' #@param formula The formula used to compute the prevalence. can be either "matrix" for the linear equations or "General" for the exact allele count cases. Default : Matrix
#' @param min_cells Minimum number of cells (default 2). In case the estimated number of cells sequenced at the locus of the mutation is less than min_cells, NA is returned.
#' @param min_alleles Minimum number of alleles. (default 4). In case the estimated number of alleles sequenced at the locus of the mutation is less than min_alleles, NA is returned.
#' @param detail when set to TRUE, a detailed output is generated containing, the context and the detailed prevalence for each group of cells (germline cells, cells affected by one of the two genomic alterations (SNV or CNV) but not both, cells affected by  by both copynumber alteration and SNV ). Default : TRUE.
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence dataframe (e.g REF, ALL, ...) }
#'         \item{}{One column per tumour sample reporting the prevalence of the mutation 
#'         at each samples}
#'      }
#'      
#' @examples
#' 
#' #Example 1: Loading a simple example data set with two somatic mutations, 5 germlines SNP, and 3 tumor samples
#' data(simpleExample2)
#' se=simpleExample2
#' prevalence_df=getPrevalenceMultiSamples(se$snp_allelecount_df, se$ref_allelecount_df, se$major_copynumber_df,se$minor_copynumber_df,phasing_association_df=se$phasing_association_df, )
#' print(prevalence_df)
#' 
#' #Chrom     End IsGermline  Tumour1        Tumour2        Tumour3
#' #mutation2  chr2 3003000          0 C2:0|0|1 C2:0.15|0|0.85 C2:0.12|0|0.88
#' #mutation6  chr2 4008000          0 C1:1|0|0       C1:1|0|0 C2:0|0.24|0.76
#' 
#' #Example 2 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient 11152 (data retrieved from a parallel study)
#' 
#' data("chr15_OP1019")
#' ds=chr15_OP1019
#' masterprevalence_df=getPrevalenceMultiSamples(ds$snp_allelecount_df, ds$ref_allelecount_df,  ds$major_copynumber_df,ds$minor_copynumber_df,phasing_association_df = ds$phasing_association_df, cnv_fraction=ds$CNVFraction_df,nbFirstColumns=6,detail=FALSE)
#' print(head(masterprevalence_df))
#' 
#' data("chr10_OP1019")
#' df=chr10_OP1019
#' masterprevalence_df=getPrevalenceMultiSamples(df$snp_allelecount_df, df$ref_allelecount_df, df$major_copynumber_df,df$minor_copynumber_df,phasing_association_df=df$phasing_association_df, cnv_fraction=df$CNVFraction_df,nbFirstColumns=6, region="chr10:50000000-180000000")
#' print(head(masterprevalence_df))
#' 
#' 
#' 
#'@seealso \code{\link{getPrevalence}}
#' @export
getPrevalenceMultiSamples<-function(snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode="PhasedSNP",cnv_fraction=NULL, phasing_association_df=NULL,NormalcellContamination_df=NULL,tumoursamples=NULL,  nbFirstColumns=3, region=NULL,detail=TRUE,  LocusRadius = 10000,NoPrevalence.action="Skip",SameTumour=TRUE)
{
  
  # Extract the somatic mutations 
  
  
  compulsory_columns=c("Chrom","End","IsGermline")
  
  if (length(setdiff(compulsory_columns,colnames(snp_allelecount_df)))>0){
    stop(" The allele count master matrices should have at least the following headers
         columns : ")
    print(compulsory_columns)
  }
  
  
  if (length(setdiff(compulsory_columns,colnames(ref_allelecount_df)))>0){
    stop(" The allele count master matrices should have at least the following 
         headers columns : ")
    print(compulsory_columns)
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
  
  if(!is.null(cnv_fraction))
    tumoursamples =intersect(tumoursamples,colnames(cnv_fraction) )
  
  
  if(length(tumoursamples) ==0)
  {
    stop(" None of the tumour samples provided is present in the five  matrices :
         snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df, cnv_fraction")
  }
  
  snp_allelecount_df=numeric_column(snp_allelecount_df,tumoursamples)
  ref_allelecount_df=numeric_column(ref_allelecount_df,tumoursamples)  
  major_copynumber_df=numeric_column(major_copynumber_df,tumoursamples)
  minor_copynumber_df=numeric_column(minor_copynumber_df,tumoursamples)
  if(!is.null(cnv_fraction)) cnv_fraction=numeric_column(cnv_fraction,tumoursamples)
  
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
  
  
    masterprevalence = getPrevalence_Matrice(snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,mode,cnv_fraction, phasing_association_df,NormalcellContamination_df,tumoursamples,  nbFirstColumns, region,detail,  LocusRadius,NoPrevalence.action,SameTumour)

  masterprevalence
  
  }

