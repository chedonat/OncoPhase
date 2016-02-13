############################################################
#
# OncoPhase Organisation file
# data.R file to generate the Roxygen help on data
# OncoPhase.R this file, main file of the package



#Roxigen help for the package (top help)

#' OncoPhase package for somatic mutation cellular prevalence estimation using haplotype phasing
#' 
#' The main functions for somatic mutation cellular prevalence computation is  \code{\link{getPrevalence}}. 
#' See the examples at \code{\link{getPrevalence}} for basic steps.
#' 
#' Input data for simple case studies as the ones illustrated in the accompagning paper can be generated 
#' with the function \code{\link{build_casestudy}}.
#' 
#' The package include experimental data for chromosome 15 and 22 for two patients of a parallele study.
#'  (see \code{\link{chr15_OP1019}} and \code{\link{chr22_11152}} )
#' 
#' For more detailed information on usage, see the package vignette, by typing
#' \code{vignette("OncoPhase")}. All support questions should be emailed to  the authors.
#'
#' @references
#'
#' OncoPhase reference:
#' 
#' OncoPhase: A package for computing Somatic Mutation cellular Prevalence in cancer using haplotype phasing. Bioinformatics 2016. Submitted
#'
#' OncoPhase reference:
#' 
#' SOX2 paper. Detail of the reference to be added
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


numeric_column<-function(df,tumoursamples)
{
    for (col in tumoursamples)
    {
      df[col] =  as.numeric(unlist(df[col]))
    }
  df
}



#' Somatic mutation cellular prevalence using haplotype phasing.
#' 
#' This is a generic function computing the cellular prevalence of somatic mutations in
#'  cancer using haplotype phasing. The model compute the prevalence of a somatic
#'   mutation relatively to close and phased germline mutation. It uses three sources
#'    of information as input : The allelic counts, the phasing information and the 
#'    copy number alteration. 
#' 
#' @param snp_allelecount_df A data frame containing for each mutation the  allelic 
#' count of the variant at each tumour. Should contains at least Chrom (The mutation
#'  chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'   or Somatic)  among the first columnns.
#' @param ref_allelecount_df A data frame containing for each mutation the allelic count
#'  of the reference at each tumour.Should contains at least Chrom (The mutation
#'   chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic)  among the first columnns.
#' @param phasing_association_df A data frame containing for each somatic mutation, 
#' a colon separated list of germline SNP phased to it.
#' @param major_copynumber_df A data frame containing for each mutation the major
#'  chromosonal copy number at each tumour sample.
#' @param minor_copynumber_df A data frame containing for each mutation the minor
#'  chromosonal copy number at each tumour sample.
#' @param normalfraction_df A data frame containing for each mutation the fraction of
#'  normal cell contamination at each tumour/sample of the study.
#' @param nbFirstColumns Number of first columns in snp_allelecount_df to reproduce in
#'  the output dataframe e.g: Chrom, Pos, Vartype.
#' @param region The region of the genome to consider in the format chrom:start-end 
#' e.g "chr22:179800-98767 
#' @param tumoursamples The list of tumour samples to consider for the prevalence
#'  computation. This samples should be present as column header in the data frames 
#'   snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df
#'   and  normalfraction_df. If not provided, the headers from nbFirstColumns + 1 to 
#'   the last column of snp_allelecount_df is retrieved and its intersection with the
#'    headers of the four othe rmatrices is considered as the tumour samples
#'   of the allelecount or copy number matrices)
#' @param method The method to be used for prevalence computation  (default : General , alternatives are FlankingGeneral 
#' , Linear, FlankingLinear.)
#' @param min_cells Minimum number of cells (default 2)  In case the estimated number of cells sequenced at the locus of the mutation is less than min_cells, NA is returned.
#' @param min_alleles Minimum number of alleles (default 4). In case the estimated number of alleles sequenced at the locus of the mutation is less than min_alleles, NA is returned.
#' @return A data frame containing :
#'  \describe{
#'        \item{}{Column 1 to NbFirstcolumn of the input data frame snp_allelecount_df. 
#'        This will generally include the chromosome and the position of the mutation plus
#'        any other column you will like to report in the prevalence (e.g REF, ALL, ...) }
#'         \item{}{One column per tumour sample reporting the prevalences of the mutation 
#'         at each samples}
#'      }
#'      
#' @examples
#' 
#' # Example 1: Loading a simple  example data set  with two somatic mutations, 5 germlines SNP, and 3 tumour samples
#' data(simpleExample2)
#' attach(simpleExample2)
#' prevalence_df=getPrevalence(snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df)
#' 
#' 
#' # Example 2: Running a case study as illustrated in the accompagning paper. Available case studies : A, B, C, 1 ,2, . . ., 9
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
#' #Example 3 : Computing somatic mutation cellular prevalence on chromosome 15 of  patient 11152 (data retrieved from a parallele study)
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
#' # Example 4 : Creating a simple example with one somatic mutation and one germline mutation on a single tumour sample
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
getPrevalence<-function(snp_allelecount_df, ref_allelecount_df,phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df,nbFirstColumns=3,method="General",tumoursamples=NULL, region=NULL,min_cells=2, min_alleles=4)
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
                            colnames(minor_copynumber_df),
                            colnames(normalfraction_df)))
  if(length(tumoursamples) ==0)
  {
    stop(" None of the tumour samples provided is present in the five  matrices :
snp_allelecount_df, ref_allelecount_df, major_copynumber_df,minor_copynumber_df,
         normalfraction_df")
  }
  
  snp_allelecount_df=numeric_column(snp_allelecount_df,tumoursamples)
  ref_allelecount_df=numeric_column(ref_allelecount_df,tumoursamples)  
  major_copynumber_df=numeric_column(major_copynumber_df,tumoursamples)
  minor_copynumber_df=numeric_column(minor_copynumber_df,tumoursamples)
  normalfraction_df=numeric_column(normalfraction_df,tumoursamples)
  
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
  
  
  
  
  
  
  imutation=0
  
  
  masterprevalence<-matrix(nrow=nrow(somatic_snp_allelecount_df),ncol=nbFirstColumns + length(tumoursamples))
  masterprevalence<-as.data.frame(masterprevalence)
  colnames(masterprevalence) <- c(colnames(somatic_snp_allelecount_df[1:nbFirstColumns]),tumoursamples)
  rownames(masterprevalence) <- rownames(somatic_snp_allelecount_df)
  masterprevalence[1:nbFirstColumns] = somatic_snp_allelecount_df[1:nbFirstColumns]
  
  for (imut in 1:nrow(masterprevalence))
  {
    
    
    imutation =imutation+1; 
    mut <- rownames(masterprevalence[imut,]); mut_pos=as.numeric(masterprevalence[imut,"End"])
    
    
    #  output=(((imutation-1) %%100)==0)
    
    
    # if(output)
    #   cat (" \n\n\n\n  mut ",imutation,"/ ",nrow(masterprevalence)," (", mut,") ")
    
    
    
    ############################
    #####Source of information 1: The list of phased germline
    ##############
    
    # Retrieving the phased germline, if not phased germline skip
    
    phased_germline=""
    phased_germline=  as.character(phasing_association_df[mut,"PhasedMutations"])
    
    if(is.null(phased_germline))
      next
    phased_list<-unlist(strsplit(phased_germline,":"));  
    if (length(unlist(phased_list))==0)
      next
    
    #Keep only the germline having allele counts
    phased_list=intersect(phased_list, rownames(snp_allelecount_df[snp_allelecount_df$IsGermline==1,]))
    
    if (length(unlist(phased_list))==0)
      next
    
    #  if(output) 
    #   cat ( " :  ", length(phased_list), "Phased present in the Alelle matrice")  
    
    
    
    
    
    ############################
    #####Source of information 2: The Copy Number 
    ##############
    
    
    #For the somatic mutation
    phi_cn_sample_somatic_list = 1-normalfraction_df[mut,tumoursamples,drop=F]
    Major_cn_sample_somatic_list = major_copynumber_df[mut,tumoursamples, drop=F]
    Minor_cn_sample_somatic_list = minor_copynumber_df[mut,tumoursamples,drop=F]
    Total_cn_sample_somatic_list = Major_cn_sample_somatic_list + Minor_cn_sample_somatic_list 
    
    
    #For the phased germline mutations
    phi_cn_sample_PhasedGermline_df = 1-normalfraction_df[phased_list,tumoursamples,drop=F]
    Major_cn_sample_PhasedGermline_df = major_copynumber_df[phased_list,tumoursamples,drop=F ]
    Minor_cn_sample_PhasedGermline_df = minor_copynumber_df[phased_list,tumoursamples,drop=F]
    Total_cn_sample_PhasedGermline_df = Major_cn_sample_PhasedGermline_df + Minor_cn_sample_PhasedGermline_df
    
    
    
    ############################
    #####Source of information 3: The Allele Count 
    ##############
    ############################
    
    
    #Somatic
    #wellfraction_somatic =snp_allelecount_df[mut,cifs:nbcolumns_wellfraction]
    snpwellcount_somatic =  data.matrix(snp_allelecount_df[mut,tumoursamples,drop=F])
    refwellcount_somatic = data.matrix(ref_allelecount_df[mut,tumoursamples,drop=F])
    
    #Germline
    # wellfraction_germlines=snp_allelecount_df[phased_list,cifs:nbcolumns_wellfraction]
    snpwellcount_germlines=  data.matrix(snp_allelecount_df[phased_list,tumoursamples, drop=F])
    refwellcount_germlines = data.matrix(ref_allelecount_df[phased_list,tumoursamples, drop=F])
    
    
    
    ############################
    ##### Mutations locus restriction
    ##############
    ############################
    
    # We need to pass trough each samples,  and assign NA to the entries of any phased germline appearing to not be on the same locus with the somatic mutation ( not having same phi, major_cn and minor_cn as the somatic and with no  cn breakpoint between the germline and the somatic.)
    
    submatrix_phased = snp_allelecount_df[phased_list, 1:3 , drop=F] #retrieve a submatrix of the phased with the chromosome and the position.
    
    submatrix_phased_leftside = submatrix_phased[submatrix_phased$End <mut_pos,1:3 , drop=F ] #germlines at the left of the mutation 
    submatrix_phased_rightside = submatrix_phased[submatrix_phased$End >mut_pos,1:3 , drop=F ] # germline at the right
    
    #order the leftside from highest to smallest position, leave the rightside from smallest to highest
    orders_phase=order(submatrix_phased_leftside["End"],decreasing=T)
    submatrix_phased_leftside=submatrix_phased_leftside[orders_phase,]
    
    
    #Germline
    #wellfraction_samelocus_germlines=wellfraction_germlines
    snpwellcount_samelocus_germlines= snpwellcount_germlines
    refwellcount_samelocus_germlines = refwellcount_germlines 
    
    
    at_least_one_good_germline=F # is there atleast one good germline kept?
    
    for(sample in tumoursamples)
    {
      samelocus_germline=c()
      phi_som=phi_cn_sample_somatic_list[sample]
      major_som=Major_cn_sample_somatic_list[sample]
      minor_som=Minor_cn_sample_somatic_list[sample]
      allelecount_som=snpwellcount_somatic[mut,sample]
      if (is.na(allelecount_som)) #Nothing wont be done on this sample anyway withut the snp count of the somatic.
        next
      
      absence_copynumberprofile=c()
      count_lower_than_somatic<-c() # For more accuracy in case of abundance of germline, someone can choose to exclude germline having an allele count less than the somatic allele count.
      
      for (germ in rownames(submatrix_phased_leftside))
      {
        phi_germ=phi_cn_sample_PhasedGermline_df[germ , sample]
        major_germ=Major_cn_sample_PhasedGermline_df[germ , sample]
        minor_germ=Minor_cn_sample_PhasedGermline_df[germ, sample]
        # allelecount_germ= snpwellcount_samelocus_germlines[germ,sample]
        
        if ( (phi_germ==phi_som || is.na(phi_germ)|| is.na(phi_som))&& 
             (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) &&
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          next
        }
        
        if(is.na(phi_germ) || is.na(major_germ))# We dont want NA to be considered as a breakpoint, but if the copy number information is absent we need to not consider this germline at that particular sample
          absence_copynumberprofile=c(absence_copynumberprofile,germ)
        
        # if(allelecount_germ< 1 * allelecount_som & !is.na(allelecount_germ))
        #   count_lower_than_somatic<-c( count_lower_than_somatic, germ)
        
        
      }
      
      for (germ in rownames(submatrix_phased_rightside))
      {
        phi_germ=phi_cn_sample_PhasedGermline_df[germ , sample]
        major_germ=Major_cn_sample_PhasedGermline_df[germ , sample]
        minor_germ=Minor_cn_sample_PhasedGermline_df[germ, sample]
        # allelecount_germ= snpwellcount_samelocus_germlines[germ,sample]
        
        if ( (phi_germ==phi_som || is.na(phi_germ)|| is.na(phi_som))&& 
             (major_germ==major_som || is.na(major_germ)|| is.na(major_som)) && 
             (minor_germ == minor_som || is.na(minor_germ)|| is.na(minor_som)))
        {
          samelocus_germline<-c(samelocus_germline,germ)
        }else
        {
          next
        }
        
        if(is.na(phi_germ) || is.na(major_germ))
          absence_copynumberprofile=c(absence_copynumberprofile,germ)
        
        # if(allelecount_germ< 1 * allelecount_som & !is.na(allelecount_germ))
        #   count_lower_than_somatic<-c( count_lower_than_somatic, germ)
        
        
      }
      
      #if we assign NA to the alelle count information of germline mutation not on the same somatic mutation, we are removing them from the list at the considered sample
      notsamelocus_germline = setdiff(phased_list, samelocus_germline)
      germline_to_exclude_at_this_sample=c(notsamelocus_germline,absence_copynumberprofile)
      
      if(!at_least_one_good_germline)
        at_least_one_good_germline=(length(germline_to_exclude_at_this_sample)< length(phased_list))
      
      
      if(length(germline_to_exclude_at_this_sample) ==0) next
      
      
      
      #  wellfraction_samelocus_germlines[germline_to_exclude_at_this_sample,sample]= rep(NA, length(germline_to_exclude_at_this_sample))
      snpwellcount_samelocus_germlines[germline_to_exclude_at_this_sample,sample]= rep(NA, length(germline_to_exclude_at_this_sample))
      refwellcount_samelocus_germlines [germline_to_exclude_at_this_sample,sample]= rep(NA, length(germline_to_exclude_at_this_sample))
      
    }
    
    if(!  at_least_one_good_germline) next # there is no good gerline
    
    ############################
    ##### We prepare the input for the formulw
    ##############
    ############################
    
    #Allele Count
    ############
    # The somatic mutation
    lambda_somatic_list=snpwellcount_somatic[mut, tumoursamples,drop=F] # \lambda(S) across the samples (see paper)
    
    # The germline  mutation, approximation according to poisson distribution, see paper
    lambda_PhasedGermline_list<-colMeans(data.matrix(snpwellcount_samelocus_germlines), na.rm=T)  #  for \lambda(G)
    mu_PhasedGermline_list<-colMeans(data.matrix(refwellcount_samelocus_germlines), na.rm=T) # for \mu(G)
    
    
    
    omega_PhasedGermline_mean_list<-lambda_PhasedGermline_list/(lambda_PhasedGermline_list+mu_PhasedGermline_list) # for \omega(G)
    
    #Copy number profile (simply the one of the somatic mutation locus)
    ############
    phi_cn_list= unlist(phi_cn_sample_somatic_list)
    Major_cn_list= unlist(Major_cn_sample_somatic_list)
    Minor_cn_list = unlist(Minor_cn_sample_somatic_list)
    
    
    
    
    
    ##########
    #THE Allele count MODEL
    #########
    ########
    
    
    
    #First, we separate the samples according to their copy number profile. 
    CNV_groups<-list()
    for (phival in unique(phi_cn_list)) 
      for(majorval in unique(Major_cn_list)) 
        for (minorval in unique(Minor_cn_list))
        {
          
          if(!is.na(phival) && !is.na(majorval) && !is.na(minorval))
          {
            samplecnvgroups=c()
            
            for (sample in tumoursamples)
              if(!is.na(phi_cn_list[sample]) && !is.na(Major_cn_list[sample]) && !is.na(Minor_cn_list[sample]))
                if ((phi_cn_list[sample] ==phival) &&(Major_cn_list[sample]==majorval) && (Minor_cn_list[sample]== minorval ))
                  samplecnvgroups=c(samplecnvgroups,sample)
                if(length(samplecnvgroups)!=0)
                  CNV_groups[length(CNV_groups) +1] = list(samplecnvgroups)
          }
        }
    
    
    
    
    
    sigma_PhasedGermline_list=Major_cn_list-Major_cn_list
    tau_PhasedGermline_list= phi_cn_list * sigma_PhasedGermline_list + (1-phi_cn_list) # for \tau(G)
    hatlambda_somatic_list=lambda_somatic_list-lambda_somatic_list
    hatlambda_PhasedGermline_list=lambda_PhasedGermline_list-lambda_PhasedGermline_list
    condition_list=rep(NA,length(lambda_somatic_list))
    names(condition_list)= names(lambda_somatic_list)
    
    
    
    
    #We run the model on each group
    for (icnvgroup in 1 : length(CNV_groups))
    {
      
      sample_cnvgroup=CNV_groups[[icnvgroup]]
      
      
      # Estimate of sigma. See the paper
      #################################
      A=0
      B=0
      Asnp=0
      Bref=0
      for (sample_i in sample_cnvgroup)
      {
        
        
        for ( germlinemut in rownames(snpwellcount_samelocus_germlines) )
        {
          if (!is.na(snpwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (snpwellcount_samelocus_germlines[germlinemut,sample_i]>0)
          )
            Asnp= Asnp + snpwellcount_samelocus_germlines[germlinemut,sample_i]
          if (!is.na(refwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (refwellcount_samelocus_germlines[germlinemut,sample_i]>0)
          )
            Bref=Bref+refwellcount_samelocus_germlines[germlinemut,sample_i]
          
          if(!is.na(snpwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (snpwellcount_samelocus_germlines[germlinemut,sample_i]>0) 
             &&
             !is.na(refwellcount_samelocus_germlines[germlinemut,sample_i]) #&& (refwellcount_samelocus_germlines[germlinemut,sample_i]>0)
          )
          {
            A= A + snpwellcount_samelocus_germlines[germlinemut,sample_i]
            B=B+refwellcount_samelocus_germlines[germlinemut,sample_i]
          }
        }
      }
      
      if((is.na(A)) || (is.na(B)))
      {
        A=Asnp
        B=Bref
      }
      if(A>=B)
        sigma_PhasedGermline_list[sample_cnvgroup] = Major_cn_list[sample_cnvgroup] 
      if(A<B)
        sigma_PhasedGermline_list[sample_cnvgroup] = Minor_cn_list[sample_cnvgroup] 
      
      
      tau_PhasedGermline_list[sample_cnvgroup]= phi_cn_list[sample_cnvgroup] * sigma_PhasedGermline_list[sample_cnvgroup] + (1-phi_cn_list[sample_cnvgroup]) # for \tau(G)
      
      
      #Allele_Count
      alpha_list=1/tau_PhasedGermline_list[sample_cnvgroup]
      beta_list=(phi_cn_list[sample_cnvgroup]  * sigma_PhasedGermline_list[sample_cnvgroup] )/ tau_PhasedGermline_list[sample_cnvgroup]
      
      #condition_list2=lambda_somatic_list - alpha_condition * lambda_PhasedGermline_list# The condition, if >=0 then case 1 else case 2
      lambda_S=as.numeric(lambda_somatic_list[mut,sample_cnvgroup])
      lambda_G=data.matrix(snpwellcount_samelocus_germlines[phased_list,sample_cnvgroup,drop=F])
      alpha=as.numeric(alpha_list[sample_cnvgroup])
      beta=as.numeric(beta_list[sample_cnvgroup])
      
      EM_parameters=bestAllele(lambda_S, lambda_G, alpha, beta)
      
      
      
      hatlambda_somatic_list[mut, sample_cnvgroup] =EM_parameters$hatlambda_S
      hatlambda_PhasedGermline_list[sample_cnvgroup] =EM_parameters$hatlambda_G
      condition=EM_parameters$bestC
      condition_list[sample_cnvgroup]=rep(condition, length(lambda_S))
      
      
    }
    
    
    ############################
    #####   ############################
    ##### Intermediary computation
    ##############
    ############################
    
    u_list =phi_cn_list * hatlambda_PhasedGermline_list/tau_PhasedGermline_list # for u
    v_list=(1-phi_cn_list) * hatlambda_PhasedGermline_list/tau_PhasedGermline_list # for v
    
    sum_cells_list<-u_list+v_list # estimation of the number of cells if well counts and the K * the number of cells of reads count with k = amplification factor
    sum_allele_list= hatlambda_PhasedGermline_list +  mu_PhasedGermline_list # estimation of the number of allele
    
    
    
    ############################
    ##### Prevalence computation
    ##############
    ############################
    
    #`The user should specify a minimum count of cells and allele for the prevalence to be computed. By default we take 6 cells and 8 alleles.
#     min_cells=4
#     min_alleles=6
#     
    
    
    prev_somatic_list =vector("numeric", length=length(hatlambda_somatic_list))
    names(prev_somatic_list) = colnames(hatlambda_somatic_list)
    prev_somatic_list[prev_somatic_list==0]<-NA
    
    
    for(sample in tumoursamples)
    {
      
      if (is.na(condition_list[sample]))
        next
      
      if(is.na(sum_cells_list[sample]) || is.na(sum_allele_list[sample]))
        next
      if (sum_cells_list[sample]< min_cells)
        next
      if(sum_allele_list[sample]< min_alleles)
        next
      
      if (condition_list[sample]=="C2")#CONTEXTE 2
        prev_somatic_list[sample] = as.numeric(unlist(   phi_cn_list[sample] + (1-phi_cn_list[sample]) * ((hatlambda_somatic_list[mut,sample] - u_list[sample] * sigma_PhasedGermline_list[sample])/v_list[sample])   ))
      else if (condition_list[sample]=="C1")
        prev_somatic_list[sample] = as.numeric(unlist( (hatlambda_somatic_list[mut,sample] / hatlambda_PhasedGermline_list[sample])* tau_PhasedGermline_list[sample]    ))
      
      #if(sample=="O13_A_ABpre")
      #   exit()
      #if(!is.na(prev_somatic_list[sample]) && prev_somatic_list[sample] > 1.00000001 ) exit()
      #if(!is.na(prev_somatic_list[sample]) && prev_somatic_list[sample] < 0 ) exit()
      
    }
    
    
    masterprevalence[mut,names(prev_somatic_list)] = prev_somatic_list
    
    
    
    masterprevalence[mut,names(prev_somatic_list)] = prev_somatic_list
    
    list_prev=unlist(prev_somatic_list)
    list_prev=list_prev[!is.na(list_prev)]
    
    
    
    
  }
  
  
  masterprevalence
  
}







#' Build the input data matrix for a case study
#' 
#' This is a generic function to automatically build the five input data frame (snp_allelecount_df, ref_allelecount_df, phasing_association_df, major_copynumber_df,minor_copynumber_df,normalfraction_df) for a case study with one somatic mutation, one germline mutation and one tumour sample.
#' 
#' @param lambda_G : count of allele suporting the variant sequence of  the Germline SNP
#' @param mu_G : count of allele suporting the reference sequence  of  the Germline SNP
#' @param lambda_S : count of allele suporting the variant sequence of  the somatic mutation 
#' @param mu_S : count of allele suporting the reference sequence  of  the somatic mutation
#' @param phi_G : Estimated fraction of cells affected by the CNV (1- normal genotype cell fraction)
#' @param major_cn: Major copy number at the locus of the mutation
#' @param minor_cn : Minor copy number at the locus of the mutation 
#' @return A list containing the following data frames:
#' 
#' @param snp_allelecount_df A data frame containing for each mutation the  allelic 
#' count of the variant at each tumour. Should contains at least Chrom (The mutation
#'  chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'   or Somatic)  among the first columnns.
#' @param ref_allelecount_df A data frame containing for each mutation the allelic count
#'  of the reference at each tumour.Should contains at least Chrom (The mutation
#'   chromosome) , End (The mutation position) and IsGermline (is the mutation a Germline
#'    or Somatic)  among the first columnns.
#' @param phasing_association_df A data frame containing for each somatic mutation, 
#' a colon separated list of germline SNP phased to it.
#' @param major_copynumber_df A data frame containing for each mutation the major
#'  chromosonal copy number at each tumour sample.
#' @param minor_copynumber_df A data frame containing for each mutation the minor
#'  chromosonal copy number at each tumour sample.
#' @param normalfraction_df A data frame containing for each mutation the fraction of
#'  normal cell contamination at each tumour/sample of the study.
#' @return A list containing the following dataframes :
#'  \describe{
#'        \item{snp_allelecount_df}{A data frame containing the  count of allele suporting the variant sequence at the somatic and germline  mutations. Chrom is set to chr3, Position of the germline and somatic mutations are respectively set to 100 and 1000}
#'        \item{ref_allelecount_df}{A data frame containing the  count of allele suporting the reference sequence at the somatic and germline  mutation. Chrom is set to chr3, Position of the germline and somatic mutations are respectively set to 100 and 1000}
#'        \item{phasing_association_df}{A data frame containing phasing association}
#'         \item{major_copynumber_df}{A data frame containing the major copy number at the mutation locus}
#'         \item{minor_copynumber_df}{A data frame containing the minor copy number at the mutation locus}
#'        \item{normalfraction_df}{A data frame containing the proportion of cells with a normal genotype}
#'      }
#'      
#' @examples
#' 
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
#' @export
build_casestudy<-function(lambda_G, mu_G,lambda_S, mu_S, phi_G, major_cn, minor_cn)
{
  nbTumour = 1
  chrom="chr3"
  
  snp_allelecount_df=as.data.frame(matrix(ncol=3+nbTumour,nrow=2))
  names(snp_allelecount_df) = c("Chrom","End","IsGermline",paste("Tumour",1:nbTumour,sep=""))
  rownames(snp_allelecount_df) = c("germlineM","somaticM")
  ref_allelecount_df = snp_allelecount_df
  major_copynumber_df= snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  minor_copynumber_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  normalfraction_df = snp_allelecount_df[paste("Tumour",1:nbTumour,sep="")]
  
  snp_allelecount_df["germlineM",] = c(chrom, 100,1,lambda_G)
  snp_allelecount_df["somaticM",] = c(chrom, 1000,0,lambda_S)
  ref_allelecount_df["germlineM",] = c(chrom, 100,1,mu_G)
  ref_allelecount_df["somaticM",] = c(chrom, 1000,0,mu_S)
  
  
  major_copynumber_df["Tumour1"] = c(major_cn,major_cn)
  minor_copynumber_df["Tumour1"] = c(minor_cn,minor_cn)
  normalfraction_df["Tumour1"] = c(1-phi_G,1-phi_G)
  
  phasing_association_df = as.data.frame(matrix(ncol=1,nrow=1))
  colnames(phasing_association_df) = c("PhasedMutations")
  rownames(phasing_association_df) = c("somaticM")
  phasing_association_df["somaticM","PhasedMutations"] = "germlineM"
  
  
  # nbFirstColumns = 3
  # tumoursamples = "Tumour1"
  # region = chrom
  list(snp_allelecount_df=snp_allelecount_df, ref_allelecount_df=ref_allelecount_df,
       phasing_association_df=phasing_association_df, major_copynumber_df=major_copynumber_df,
       minor_copynumber_df=minor_copynumber_df,normalfraction_df=normalfraction_df)
  
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


























