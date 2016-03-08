############################################################
#
# OncoPhase : files Organization 
# data.R - file to generate the Roxygen help documents for the data
# OncoPhase.R -  this file is the main file of the package




#Roxigen help for the package (top help)

#' OncoPhase package for somatic mutations cellular prevalence estimation using haplotype phasing
#' 
#' The main function for somatic mutation cellular prevalence computation is  \code{\link{getPrevalence}}. 
#' See the manual and examples at \code{\link{getPrevalence}} for more details.
#' 
#' Input data for simple case studies can be generated  with the function \code{\link{build_casestudy}}.
#' 
#' The package include experimental data for chromosome 10, 15, 18  and 22 for two patients retrieved from a parallel study.
#'  (see for example  \code{\link{chr15_OP1019}} and \code{\link{chr22_11152}} )
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
#' @name OncoPhase
#' @aliases OncoPhase-package
#' @keywords package
NULL








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
#'  cancer using haplotype phasing. The model computes the prevalence of a somatic
#'   mutation relatively to close and phased germline mutation. It uses three sources
#'    of information as input : The allelic counts, the phasing information and the 
#'    copy number alteration. 
#' 
#' @param snp_allelecount_df A data frame containing for each mutation the  allelic 
#' count of the variant at each tumor samples. The data frame should contains at least the following three columns among its firsts columns : Chrom (The mutation
#'  chromosome) , End (The mutation position) and IsGermline (is the mutation a germline
#'   or somatic).
#' @param ref_allelecount_df A data frame containing for each mutation the allelic count
#'  of the reference at each tumor samples. The data frame should contains at least the following three columns among its firsts columns:  Chrom (The mutation
#’ chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic)  
#' @param phasing_association_df A data frame containing for each somatic mutation, 
#' a colon separated list of germline SNP phased to it.
#' @param major_copynumber_df A data frame containing for each mutation, its  major
#’ chromosomal copy number at each tumor samples.
#' @param minor_copynumber_df A data frame containing for each mutation the minor
#'  chromosomal copy number at each tumor samples.
#' @param normalfraction_df A data frame containing for each mutation the fraction of
#'  normal cell contamination at each tumor samples of the study. Used only if the method is "PhasedSNPgeneral"
#' @param nbFirstColumns Number of first columns in snp_allelecount_df to reproduce in
#'  the output dataframe e.g: Chrom, Pos, Vartype.
#' @param region The region of the genome to consider for the prevalence computation  in the format chrom:start-end 
#' e.g "chr22:179800-98767 
#' @param tumoursamples : The list of tumor samples to consider for the prevalence
#’ computation.  This samples should be present as column header in the data frames 
#'   snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df
#'   and  normalfraction_df. If not provided, the headers from nbFirstColumns + 1 to 
#'   the last column of snp_allelecount_df is retrieved and its intersection with the
#' other inputted data frames headers is considered.
#' @param method The method to be used for prevalence computation  (default : PhasedSNP , alternatives methods  are PhasedSNPGeneral, FlankingSNP,FlankingSNPGeneral)
#' @param min_cells Minimum number of cells (default 2). In case the estimated number of cells sequenced at the locus of the mutation is less than min_cells, NA is returned.
#' @param min_alleles Minimum number of alleles. (default 4). In case the estimated number of alleles sequenced at the locus of the mutation is less than min_alleles, NA is returned.
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other columns to report in the prevalence dataframe (e.g REF, ALL, ...) }
#'         \item{}{One column per tumour sample reporting the prevalences of the mutation 
#'         at each samples}
#'      }
#'      
#' @examples
#' 
#' # Example 1: Loading a simple example data set with two somatic mutations, 5 germlines SNP, and 3 tumor samples
#' data(simpleExample2)
#' attach(simpleExample2)
#' prevalence_df=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df)
#' 
#' 
#' # Example 2: Running a case study as illustrated in the accompanying paper. Available case studies: A, B, C, 1, 2, . . ., 9
#' data(CaseStudy_6)
#' cs=CaseStudy_6
#' prevalence_CaseStudy6=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df, cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#' 
#' data(CaseStudy_A)
#' attach(CaseStudy_A)
#' prevalence_CaseStudy_A=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df)
#' print(prevalence_CaseStudy_A)
#' # 	Chrom  End IsGermline   Tumour1
#' # somaticM  chr3 1000          0 0.6666667
#' 
#' 
#' #Example 3 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient 11152 (data retrieved from a parallel study)
#' 
#' data("chr15_OP1019")
#' ds=chr15_OP1019
#' masterprevalence_df=getPrevalence(ds$snp_allelecount_df, ds$ref_allelecount_df, ds$phasing_association_df, ds$major_copynumber_df,ds$minor_copynumber_df,ds$normalfraction_df,nbFirstColumns=6)
#' 
#' data("chr10_11152")
#' attach(chr10_11152)
#' masterprevalence_df=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df,nbFirstColumns=6, region="chr15:50000000-180000000")
#' 
#' 
#' # Example 4 : Creating a simple example with one somatic mutation and one germline mutation on a single tumor sample
#' 
#' #Empty dataframe
#' snpcount_df=as.data.frame(matrix(ncol=4,nrow=2))
#' names(snpcount_df) = c("Chrom","End","IsGermline","Tumour1")
#' rownames(snpcount_df) = c("mutation1","mutation2")
#' refcount_df = snpcount_df
#' major_cn_df= as.data.frame(matrix(ncol=1,nrow=2))
#' names(major_cn_df) = "Tumour1"
#' rownames(major_cn_df) = c("mutation1","mutation2")
#' minor_cn_df = major_cn_df
#' normalfraction_df = major_cn_df
#' 
#' #Filling the dataframes
#' snpcount_df["mutation1",] = c("chr1", 200100,0,40)
#' snpcount_df["mutation2",] = c("chr1",  200900,1,60)
#' refcount_df["mutation1",] = c("chr1",  200100,0,20)
#' refcount_df["mutation2",] = c("chr1",  200900,1,40)
#' major_cn_df["Tumour1"] = c(1,1)
#' minor_cn_df["Tumour1"] = c(1,1)
#' normalfraction_df["Tumour1"] = c(0.2,0.2)
#' 
#' #Phasing association
#' phasing_association_df = as.data.frame(matrix(ncol=1,nrow=1))
#' colnames(phasing_association_df) = c("PhasedGermline")
#' rownames(phasing_association_df) = c("mutation1")
#' phasing_association_df["mutation1","PhasedMutations"] = "mutation2"
#' 
#' #Computing the prevalence
#' prevalence_df=getPrevalence(snpcount_df, refcount_df, phasing_association_df, major_cn_df, minor_cn_df,normalfraction_df)
#' 
#' print(prevalence_df)
#' 
#' #          Chrom    End IsGermline   Tumour1
#' # mutation1  chr1 200100          0 0.6666667
#' 
#' 
#' @export
getPrevalence<-function(snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,CNV_fraction_df=NULL, nbFirstColumns=3,method="PhasedSNP",tumoursamples=NULL, region=NULL,min_cells=2, min_alleles=4)
{
  
  hg19_dfsize<-list(chr1=249250621,chr2=243199373,chr3=198022430, chr4=191154276, 
                    chr5=180915260, chr6=171115067, chr7=159138663, chr8=146364022, 
                    chr9=141213431, chr10=135534747, chr11=135006516, chr12=133851895, 
                    chr13=115169878, chr14=107349540,chr15=102531392, chr16=90354753,
                    chr17=81195210, chr18=78077248, chr19=59128983, chr20=63025520,
                    chr21=48129895,chr22=51304566, chrX=155270560, chrY=59373566,
                    chrM=16571)
  
  
  
  
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
  
  tumoursamples = Reduce(intersect,list(tumoursamples,colnames(snp_allelecount_df),
                                        colnames(ref_allelecount_df),
                                        colnames(major_copynumber_df), 
                                        colnames(minor_copynumber_df)
                                        #,colnames(normalfraction_df)
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
  
  somatic_snp_allelecount_df = snp_allelecount_df[snp_allelecount_df$IsGermline==0, ]
  
  
  
  # Region to compute the prevalence
  
  if(!is.null(region)){
    region_parts= unlist(strsplit(region,":"))
    
    chrom = region_parts[1]
    startPosition = 1
    endPosition = hg19_dfsize[chrom]
    
    somatic_snp_allelecount_df = somatic_snp_allelecount_df[somatic_snp_allelecount_df$Chrom == chrom , ]
    
    
    if(length(region_parts)>1){
      coordinates = unlist(strsplit(region_parts[2],"-"))
      startPosition = as.numeric(coordinates[1])
      endPosition = as.numeric(coordinates[2])
      somatic_snp_allelecount_df = somatic_snp_allelecount_df[somatic_snp_allelecount_df$Chrom == chrom & somatic_snp_allelecount_df$Start >= startPosition & somatic_snp_allelecount_df$End <= endPosition,]
      
    }
    
  }
  
  masterprevalence<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=nbFirstColumns + length(tumoursamples))
  masterprevalence<-as.data.frame(masterprevalence)
  colnames(masterprevalence) <- c(colnames(somatic_snp_allelecount_df[1:nbFirstColumns]),tumoursamples)
  rownames(masterprevalence) <- rownames(somatic_snp_allelecount_df)
  masterprevalence[1:nbFirstColumns] = somatic_snp_allelecount_df[1:nbFirstColumns]
  


  if(method=="PhasedSNPGeneral")
  masterprevalence = getPrevalence_PhasedSNP( masterprevalence,snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,CNV_fraction_df=CNV_fraction_df,tumoursamples=tumoursamples,min_cells=min_cells, min_alleles=min_alleles,formula="general")

  if(method=="PhasedSNP")
    masterprevalence = getPrevalence_PhasedSNP( masterprevalence,snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,CNV_fraction_df=CNV_fraction_df,tumoursamples=tumoursamples,min_cells=min_cells, min_alleles=min_alleles,formula="matrix")
  
  if(method=="FlankingSNPGeneral")
    masterprevalence = getPrevalence_FlankingSNPGeneral( masterprevalence,snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,CNV_fraction_df=CNV_fraction_df,tumoursamples=tumoursamples,min_cells=min_cells, min_alleles=min_alleles,formula="general")
  
  
  if(method=="FlankingSNP")
    masterprevalence = getPrevalence_FlankingSNPGeneral( masterprevalence,snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,CNV_fraction_df=CNV_fraction_df,tumoursamples=tumoursamples,min_cells=min_cells, min_alleles=min_alleles,formula="matrix")
  
  
  
  masterprevalence
  
  }







#' Build the input data matrices for a case study
#' 
#' This is a generic function to automatically build the five input data frame (snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df if method is PhasedSNPGeneral) for a case study with one somatic mutation, one germline mutation and one or more tumor sample.
#' 
#' @param lambda_G : A count  or a vector of counts (In the case of  multiple tumor samples) of 
#' alleles supporting the variant sequence of  the Germline SNP
#' @param mu_G : A count or a vector of counts (In the case of  multiple tumor samples)  of 
#' alleles  supporting the reference sequence  of  the Germline SNP
#' @param lambda_S : A count or a vector of counts (In the case of  multiple tumor samples)  of 
#' allele supporting the variant sequence of  the somatic mutation 
#' @param mu_S : A count or a vector of counts (In the case of  multiple tumor samples)  of 
#' allele supporting the reference sequence  of  the somatic mutation
#' @param phi_G : Estimated fraction (or a vector of fractions if multiple tumor samples) 
#' of cells affected by the CNV (1- normal genotype cell fraction). Used only in the case of using the PhasedSNPgeneral method
#' @param major_cn: Major copy number (or a vector of copy number if multiple tumor samples)
#' at the locus of the mutation
#' @param minor_cn : Minor copy number (or a vector of copy number if multiple tumor samples)
#'  at the locus of the mutation 
#' @param depthOfCoverage : Coverage depth (or a vector of depth coverage  if multiple tumor samples) 
#' at the locus of the mutation. If not provided the exact value of the counts passed as parameters 
#' are considered. If provided then a binomial sampling with replacement is performed to generate the 
#' counts. For the germline, the sampling is done with the parameters p=lambda_G / (lambda_G + mu_G) 
#' and N= depthOfCoverage and will yield  the count of allele supporting the variant sequence of the 
#' germline and the count of allele supporting the reference.  The same sampling is  apply to the 
#' somatic mutation  with the parameters : p=lambda_S/(lambda_S + mu_S) and N = depthOfCoverage.
#'  
#' @return A list containing the following data frames:
#'  \describe{
#'        \item{snp_allelecount_df}{A data frame containing the  count of allele at each tumour samples suporting the variant sequence at the somatic and germline  mutations. Chrom is set to chr3, Position of the germline and somatic mutations are respectively set to 100 and 1000}
#'        \item{ref_allelecount_df}{A data frame containing the  count of allele  at each tumour samples  suporting the reference sequence at the somatic and germline  mutation. Chrom is set to chr3, Position of the germline and somatic mutations are respectively set to 100 and 1000}
#'        \item{phasing_association_df}{A data frame containing the phasing association between the somatic and the germline mutation}
#'         \item{major_copynumber_df}{A data frame containing the major copy number  at each tumour samples  at the mutation locus}
#'         \item{minor_copynumber_df}{A data frame containing the minor copy number  at each tumour samples  at the mutation locus}
#'        \item{normalfraction_df}{A data frame containing the proportion of cells with a normal genotype  at each tumour samples. Present only if the method is "PhasedSNPgeneral" }
#'      }
#'      
#' @examples
#' 
#' #Example 1
#' # We reproduce here the case study No 6 of the paper
#'  #Build the input data 
#'  cs = build_casestudy(lambda_G=16, mu_G=8,lambda_S=14,mu_S=10,phi_G=4/8,major_cn=3,minor_cn=1  )
#'  #Run the case
#' prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
#'                          cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#' #print the result
#' print(prevalence)
#' #Chrom  End IsGermline Tumour1
#' #somaticM  chr3 1000          0    0.75
#' 
#' 
#' #Example 2
#' #Multiple tumours and stochastic generation of the counts.
#' CaseStudy_10 = build_casestudy(lambda_G=c(8,12,10), mu_G=c(5,4,8),lambda_S=c(3,6,10),mu_S=c(10,8,12),
#' major_cn=c(2,2,3),minor_cn=c(2,1,2) , depthOfCoverage = c(60,100,200))
#' 
#' cs = CaseStudy_10
#' prevalence=getPrevalence(cs$snp_allelecount_df, cs$ref_allelecount_df, cs$phasing_association_df,
#'                          cs$major_copynumber_df,cs$minor_copynumber_df,cs$normalfraction_df)
#' #print the result
#' print(prevalence)
#' 
#' @export
build_casestudy<-function(lambda_G, mu_G,lambda_S, mu_S, major_cn, minor_cn,  cnv_fraction= NULL, depthOfCoverage=NULL)
{
  
  chrom="chr3"
  N=length(lambda_G) # Number of samples
  if((length(mu_G)!=N) || (length(lambda_S) !=N) || (length(mu_S) !=N) ||
      (length(major_cn) !=N) || (length(minor_cn)!=N))
    stop("\n\nThe vectors passed as input should have the same size\n\n")
  
  if(!is.null(depthOfCoverage))
    if(length(depthOfCoverage)!=N)
      stop("\n\n Vector of coverage should have same size as input vector")
  
  
  
  nbTumour = N
  snp_allelecount_df=as.data.frame(matrix(ncol=3+nbTumour,nrow=2))
  names(snp_allelecount_df) = c("Chrom","End","IsGermline",paste("Tumour",1:nbTumour,sep=""))
  rownames(snp_allelecount_df) = c("germlineM","somaticM")
  ref_allelecount_df = snp_allelecount_df
  major_copynumber_df= snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  minor_copynumber_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  
  lambda_G_prime=vector("numeric",nbTumour)
  lambda_S_prime=vector("numeric",nbTumour)
  mu_G_prime =vector("numeric",nbTumour)
  mu_S_prime =vector("numeric",nbTumour)
  
  
  
  if(!is.null(depthOfCoverage))
  {
    
    p_G= lambda_G/(lambda_G + mu_G)
    p_S= lambda_S/(lambda_S + mu_S)
    
    
    for(isample in 1:nbTumour){
      lambda_G_prime[isample] = sum(sample(c(1,0), depthOfCoverage[isample], replace = TRUE,prob=c(p_G[isample],1-p_G[isample])))
      mu_G_prime[isample] = depthOfCoverage[isample] - lambda_G_prime[isample] 
      lambda_S_prime[isample] = sum(sample(c(1,0), depthOfCoverage[isample], replace = TRUE,prob=c(p_S[isample],1-p_S[isample])))
      mu_S_prime[isample] = depthOfCoverage[isample] - lambda_S_prime[isample]
      
    }
    
  }else{
    lambda_G_prime=lambda_G
    lambda_S_prime=lambda_S
    mu_G_prime =mu_G
    mu_S_prime =mu_S
  }
  
  
  snp_allelecount_df["germlineM",] = c(chrom, 100,1,lambda_G_prime)
  snp_allelecount_df["somaticM",] = c(chrom, 1000,0,lambda_S_prime)
  ref_allelecount_df["germlineM",] = c(chrom, 100,1,mu_G_prime)
  ref_allelecount_df["somaticM",] = c(chrom, 1000,0,mu_S_prime)
  
  
  major_copynumber_df["germlineM",paste("Tumour",1:nbTumour,sep="")] = major_cn
  minor_copynumber_df["germlineM",paste("Tumour",1:nbTumour,sep="")] = minor_cn
 # if(!is.null(phi_G)) normalfraction_df["germlineM",paste("Tumour",1:nbTumour,sep="")] = 1-phi_G
  
  major_copynumber_df["somaticM",paste("Tumour",1:nbTumour,sep="")] = major_cn
  minor_copynumber_df["somaticM",paste("Tumour",1:nbTumour,sep="")] = minor_cn
#  if(!is.null(phi_G))normalfraction_df["somaticM",paste("Tumour",1:nbTumour,sep="")] = 1-phi_G
  
  
  
  phasing_association_df = as.data.frame(matrix(ncol=1,nrow=1))
  colnames(phasing_association_df) = c("PhasedMutations")
  rownames(phasing_association_df) = c("somaticM")
  phasing_association_df["somaticM","PhasedMutations"] = "germlineM"
  
 #   stop(1)
  
 
    cs=list(snp_allelecount_df=snp_allelecount_df, ref_allelecount_df=ref_allelecount_df,
         phasing_association_df=phasing_association_df, major_copynumber_df=major_copynumber_df,
         minor_copynumber_df=minor_copynumber_df)
  
  cs
}






































complete_likelihood<-function(lambda_S, lambda_G,theta, C, alpha,beta)
{
  
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=0
  for (i in 1:N)
  {
    
    
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i])&& hatlambda_G[i]>0 )
    {
      # cat (" here PC ")
      #Constrained part
      if(C=="C1")
      {445
        if(is.na(alpha[i])) next
        if(hatlambda_S[i] <= alpha[i]  * hatlambda_G[i])
        {
          if (prob==0)
          {
            prob = dpois(lambda_S[i], hatlambda_S[i])
          }else
          {
            prob = prob * dpois(lambda_S[i], hatlambda_S[i])
          }
          
          for (j in 1:m)
          {
            if(is.na(lambda_G[j,i]) ) next
            prob = prob * dpois(lambda_G[j,i],hatlambda_G[i])
          }
          
          
        }else
        {
          prob = prob * 0
        }
      }else#(C=="C2")
      {
        if (is.na(beta[i] )) next
        if((hatlambda_S[i] >= beta[i]  * hatlambda_G[i]) && (hatlambda_S[i] <=   hatlambda_G[i]))
        {if (prob==0)
        {
          prob = dpois(lambda_S[i], hatlambda_S[i])
        }else
        {
          prob = prob * dpois(lambda_S[i], hatlambda_S[i])
        } 
          for (j in 1:m)
          {
            if(is.na(lambda_G[j,i]) ) next
            prob = prob * dpois(lambda_G[j,i],hatlambda_G[i])
          }
          
          
        }else
        {
          prob = prob * 0
        }
      }
      
    }
    #  cat(prob)
  }
  
  prob
  
  
}



complete_loglikelihood<-function(lambda_S, lambda_G,theta, C, alpha,beta)
{
  
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=-Inf
  for (i in 1:N)
  {
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i] ) && hatlambda_G[i]>0 )
    {
      #Constrained part
      if(C=="C1")
      {
        if(is.na(alpha[i] )) next
        
        if(hatlambda_S[i] <= alpha[i]  * hatlambda_G[i])
        {
          
          if (prob==-Inf)
          {
            prob =  dpois(lambda_S[i], hatlambda_S[i], log=T)
          }else{
            prob = prob + dpois(lambda_S[i], hatlambda_S[i], log=T)            
          }
          
          
          for (j in 1:m)
          {
            if(is.na(lambda_G[j,i]) ) next
            prob = prob + dpois(lambda_G[j,i],hatlambda_G[i], log=T)
          }
          
        }else
        {
          prob = prob  + -Inf 
        }
      }else
      {
        if(is.na(beta[i] )) next
        if((hatlambda_S[i] >= beta[i]  * hatlambda_G[i]) && (hatlambda_S[i] <=   hatlambda_G[i]))
        {
          if (prob==-Inf)
          {
            prob =  dpois(lambda_S[i], hatlambda_S[i], log=T)
          }else{
            prob = prob + dpois(lambda_S[i], hatlambda_S[i], log=T)
          }
          for (j in 1:m)
          {
            if(is.na(lambda_G[j,i]) ) next
            prob = prob + dpois(lambda_G[j,i],hatlambda_G[i], log=T)
          }
          
          
        }else
        {
          prob = prob + -Inf
        }
        
      }
      
      
    }
  }
  
  prob
  
  
}

incomplete_likelihood<-function(lambda_S, lambda_G, theta)
{
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=1
  for (i in 1:N)
  {
    if(is.na(lambda_S[i]) ) next
    if(sum(!is.na(lambda_G[,i]) )==0) next
    
    prob = prob * dpois(lambda_S[i], hatlambda_S[i])
    for (j in 1:m)
    {
      if(is.na(lambda_G[j,i]) ) next
      prob = prob * dpois(lambda_G[j,i],hatlambda_G[i])
    }
    
  }
  
  prob
}

incomplete_loglikelihood<-function(lambda_S, lambda_G, theta)
{
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  prob=0
  for (i in 1:N)
  {
    if(is.na(lambda_S[i]) ) next
    if(sum(!is.na(lambda_G[,i]) )==0) next
    
    prob = prob + dpois(lambda_S[i], hatlambda_S[i])
    for (j in 1:m)
    {
      if(is.na(lambda_G[j,i]) ) next
      prob = prob + dpois(lambda_G[j,i],hatlambda_G[i])
    }
    
  }
  
  prob
  
  
}




Qfunction<-function(lambda_S, lambda_G, theta, theta_m, alpha,beta)
{
  hatlambda_S = theta$hatlambda_S 
  hatlambda_G = theta$hatlambda_G
  
  
  N= length(lambda_S)
  m=nrow(lambda_G)
  logprob=0
  PC=1/2 # Prior
  PC1= complete_loglikelihood(lambda_S, lambda_G, theta_m,"C1", alpha,beta)/ PC
  PC2= complete_loglikelihood(lambda_S, lambda_G, theta_m,"C2", alpha,beta)/ PC
  gamma=1/(PC1 + PC2)
  
  
  #C1_term
  prob =0
  for (i in 1:N)
  {
    
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i]) && !is.na(alpha[i])  && hatlambda_G[i] > 0 )
    {
      
      if(hatlambda_S[i] <= alpha[i]  * hatlambda_G[i])
      {
        prob = prob + dpois(lambda_S[i], hatlambda_S[i], log=T)
        for (j in 1:m)
        {
          if(is.na(lambda_G[j,i]) ) next
          prob = prob + dpois(lambda_G[j,i],hatlambda_G[i], log=T)
        }
        
        if (dpois(lambda_S[i], hatlambda_S[i], log=T)==-Inf) stop(paste("lambda_S[i] ", lambda_S[i], "hatlambda_S[i]", hatlambda_S[i]), "hatlambda_G[i]", hatlambda_G[i])
        
      }else
      {
        prob = prob + -Inf
      }
    }
  }
  
  C1_logprob=prob
  if (PC1==-Inf)
    C1_logprob=0
  
  #C2_term
  prob =0
  for (i in 1:N)
  {
    
    if(!is.na(hatlambda_S[i] ) && !is.na(hatlambda_G[i]) && !is.na(beta[i])  && hatlambda_G[i] > 0)
    {
      if((hatlambda_S[i] >= beta[i]  * hatlambda_G[i]) && (hatlambda_S[i] <=   hatlambda_G[i]))
      {
        prob = prob + dpois(lambda_S[i], hatlambda_S[i], log=T)
        
        for (j in 1:m)
        {
          if(is.na(lambda_G[j,i]) ) next
          prob = prob + dpois(lambda_G[j,i],hatlambda_G[i], log=T)
        }
        
      }else
      {
        prob = prob + -Inf 
      }
    }
  }
  
  C2_logprob=prob
  if (PC2==-Inf)
    C2_logprob=0
  
  
  Qvalue =C1_logprob +    C2_logprob
  Qvalue
  
}



bestAllele<-function(lambda_S,lambda_G, alpha,beta){
  
  
  traceoutput=0
  
  N=length(lambda_S)
  M=nrow(lambda_G)
  
  #if alpha[i] = beta[i], we add a slight increase on beta[i]
  beta[beta==alpha] = beta[beta==alpha]* 1.01
  
  if((N>0) && (M>0))
  {
    
    # Expectation Maximization algorithm
    
    theta<-list()
    theta[[1]] = list(hatlambda_S=lambda_S,hatlambda_G= colMeans(lambda_G,na.rm=T))
    lik0=1
    lik1= incomplete_loglikelihood(lambda_S,lambda_G,theta[[1]])
    
    lik<-list()
    i=0
    while(lik1 -lik0 !=0 )
    {
      i=i+1
      lik[i]=lik1
      lik0=lik1
      #We find theta_C1
      hatlambda_S=c()
      hatlambda_G=c()
      for (isample in 1:N)
      {
        hatlambda_S[isample] = lambda_S[isample]
        hatlambda_G[isample] = as.numeric(colMeans(lambda_G,na.rm=T)[isample])
        if(is.na(hatlambda_S[isample] ) || is.na(hatlambda_G[isample])|| is.na(alpha[isample] )) next
        #Out of boundary cases
        if (hatlambda_S[isample]>alpha[isample] *hatlambda_G[isample])
          hatlambda_S[isample]=alpha[isample] * hatlambda_G[isample]
      }
      theta_C1=list(hatlambda_S=hatlambda_S,hatlambda_G=hatlambda_G )
      
      
      #We find theta_C2
      hatlambda_S=c()
      hatlambda_G=c()
      for (isample in 1:N)
      {
        hatlambda_S[isample] = lambda_S[isample]
        hatlambda_G[isample] = as.numeric(colMeans(lambda_G,na.rm=T)[isample])
        if(is.na(hatlambda_S[isample] ) || is.na(hatlambda_G[isample]) || is.na(beta[isample] )) next
        #Out of boundary cases
        if (hatlambda_S[isample]<beta[isample] *hatlambda_G[isample]) {
          hatlambda_S[isample]=beta[isample] *hatlambda_G[isample]
        }else  if(hatlambda_S[isample]>hatlambda_G[isample])
        {
          hatlambda_S[isample]= hatlambda_G[isample]
        }
        
        
        
      }
      
      theta_C2=list(hatlambda_S=hatlambda_S,hatlambda_G=hatlambda_G )
      
      
      QC2= Qfunction(lambda_S,lambda_G,theta_C2, theta[[i]], alpha,beta )
      QC1= Qfunction(lambda_S,lambda_G,theta_C1, theta[[i]],alpha,beta ) 
      
      if (QC2>=QC1)
        theta[[i+1]] =theta_C2
      else
        theta[[i+1]] =theta_C1
      
      
      lik1= incomplete_loglikelihood(lambda_S,lambda_G,theta[[i+1]])
      
      if(traceoutput)
      {
        cat("\n\n\n\n \t\t One step done \n\t\t############")
        
        cat("\n\t theta[i] :\n"); print(theta[[i]]  )
        cat("\n\t lik[i] : \n"); print(lik[[i]]  )
        cat("\n\t theta_C1 : \n"); print(theta_C1  )
        cat("\n\t theta_C2 :  \n"); print(theta_C2  )
        cat("\n\t QC1 : \n"); print( QC1  )
        cat("\n\t QC2 :  \n"); print(QC2  )
        cat("\n\t theta[i+1] : \n"); print(theta[[i+1]]  )
        cat("\n\t Lik[i+1]:  \n"); print(lik1  )
        
        
        
      }
      
      
      #We find theta_C2 
    }
    
    
    
    LikC2=complete_loglikelihood(lambda_S,lambda_G, theta[[i+1]],"C2",alpha, beta)
    LikC1=complete_loglikelihood(lambda_S,lambda_G, theta[[i+1]],"C1",alpha, beta)
    if(LikC2>=LikC1)
    {
      bestC="C2"
    }else
    {
      bestC="C1"
    }
    
    
    if(traceoutput)
    { 
      cat("\n\t LikC2 :",LikC2  )
      cat("\n\t LikC1 : ", LikC1  )
      
    }
    
    
    list(hatlambda_S=theta[[i+1]]$hatlambda_S, hatlambda_G=theta[[i+1]]$hatlambda_G,bestC=bestC)
    
  }else
  {
    if(traceoutput)
      cat("\n\n No input, check your inputs")
    
    list(hatlambda_S=c(), hatlambda_G=c(),bestC="")
    
  }
  
  
}



























